# Sample sex, imputed from expression, and used as a model covariate.
#
# Sex is a large, genome-wide source of variance that sample sheets often
# omit or get wrong. run_qc() therefore always imputes it from the
# sex-chromosome marker genes (via biowrenchR::impute_sex()), after the
# low-expression filter and before any X/Y genes are removed, and stores it
# as `colData(dds)$imputed_sex`. The pipeline then adds it to every model
# whenever doing so is estimable -- i.e. sex is not already confounded with
# `var`, `pair` or `batch` -- so a sex imbalance between groups cannot
# masquerade as a treatment effect.
#
# biowrenchR::impute_sex() always splits the samples into two k-means
# clusters, so a single-sex cohort would still be divided in two. The
# clusters are therefore only trusted when they differ clearly on the marker
# genes (`min_gap`, in log2 units); otherwise sex is "undetermined" and is
# not used.

#' Column of `colData(dds)` holding the imputed sex
#' @noRd
.sex_col <- "imputed_sex"

#' Impute sample sex from expression
#'
#' Clusters samples on sex-chromosome marker genes with
#' [biowrenchR::impute_sex()] (female: X-inactivation escape genes such as
#' *XIST*; male: Y-chromosome genes) on log2(CPM + 1). [run_qc()] calls this
#' on every run, after low-expression filtering and before X/Y gene removal.
#'
#' The two clusters are accepted only if the female cluster is at least
#' `min_gap` (log2 units) higher on the female markers *and* the male cluster
#' at least `min_gap` higher on the Y markers (each set that is present).
#' Otherwise the cohort is most likely single-sex (k-means always returns two
#' clusters), so every sample is reported as `NA` with `status =
#' "undetermined"`, and sex is not used as a covariate.
#'
#' @param dds A `DESeqDataSet` (raw counts, or a limma carrier from
#'   [expr_to_dds()]).
#' @param org `"human"`, `"mouse"`, or `NULL` (default) to detect from gene
#'   symbols.
#' @param min_gap Minimum log2 difference, in the expected direction, between
#'   the two clusters' mean marker expression on each marker set for the split
#'   to be trusted. Default `1` (two-fold).
#' @param seed Random seed for the k-means starts. Default `1`.
#' @return A list with `table` (data frame: `sample`, `imputed_sex`,
#'   `female_score`, `male_score`), `status` (`"imputed"`, `"undetermined"`
#'   or `"failed"`), `message`, `org`, and the marker genes used
#'   (`female_genes`, `male_genes`).
#' @examples
#' \dontrun{
#' sex <- impute_sample_sex(dds)
#' sex$table
#' }
#' @export
impute_sample_sex <- function(dds, org = NULL, min_gap = 1, seed = 1L) {
    dds <- validate_dds(dds)
    if (!is.numeric(min_gap) || length(min_gap) != 1L || is.na(min_gap) || min_gap < 0) {
        stop("`min_gap` must be a single non-negative number.", call. = FALSE)
    }
    blank <- data.frame(sample = colnames(dds), imputed_sex = NA_character_,
                        female_score = NA_real_, male_score = NA_real_,
                        stringsAsFactors = FALSE)
    if (is.null(org)) {
        org <- tryCatch(detect_organism(rownames(dds)), error = function(e) NA_character_)
        if (is.na(org)) {
            # Never stop QC for this: sex is an aid, and non-symbol gene IDs
            # (Ensembl, Entrez) simply cannot be matched to the markers.
            return(list(table = blank, status = "failed", org = NA_character_,
                        female_genes = character(0), male_genes = character(0),
                        message = paste0("Sex not imputed: the organism could not be detected ",
                                         "from rownames(dds) (gene symbols are needed).")))
        }
    }
    org <- validate_choice(org, c("human", "mouse"), "org")

    raw <- .raw_count_matrix(dds)
    lib <- colSums(raw)
    lib[lib == 0] <- 1
    logcpm <- log2(sweep(raw, 2, lib, "/") * 1e6 + 1)

    markers <- tryCatch(biowrenchR::sex_marker_genes(org),
                        error = function(e) list(female = character(0), male = character(0)))
    fem <- intersect(markers$female, rownames(logcpm))
    mal <- intersect(markers$male, rownames(logcpm))
    score <- function(genes) {
        if (!length(genes)) return(rep(NA_real_, ncol(logcpm)))
        colMeans(logcpm[genes, , drop = FALSE])
    }
    tbl <- data.frame(sample = colnames(dds), imputed_sex = NA_character_,
                      female_score = score(fem), male_score = score(mal),
                      stringsAsFactors = FALSE)
    out <- list(table = tbl, status = "failed", message = "", org = org,
                female_genes = fem, male_genes = mal)

    # Keep biowrenchR's own warning (e.g. "markers disagree") for the report.
    warn <- character(0)
    sex <- tryCatch(
        withCallingHandlers(
            withr::with_seed(seed, biowrenchR::impute_sex(
                logcpm, female_genes = fem, male_genes = mal, org = org)),
            warning = function(w) {
                warn <<- c(warn, conditionMessage(w))
                invokeRestart("muffleWarning")
            }),
        error = function(e) {
            out$message <<- paste0("Sex imputation failed: ", conditionMessage(e))
            NULL
        }
    )
    if (is.null(sex) || all(is.na(sex))) {
        if (!nzchar(out$message)) {
            disagree <- any(grepl("disagree", warn))
            out$status <- if (disagree) "undetermined" else "failed"
            out$message <- paste0(
                if (disagree) {
                    "Female and male markers do not split the samples the same way; the cohort is probably single-sex. "
                } else "Sex could not be imputed. ",
                "(", length(fem), " female and ", length(mal), " male marker genes after filtering",
                if (length(warn)) paste0("; biowrenchR: ", paste(sub("[.]+$", "", unique(warn)), collapse = "; ")),
                ".) Sex is not used as a covariate.")
        }
        return(out)
    }

    sex <- as.character(sex[colnames(dds)])
    # Directional gaps: how much higher the female cluster is on the female
    # markers, and the male cluster on the Y markers. k-means picks the split
    # that maximises separation, so on a single-sex cohort noise alone can
    # open a gap on one marker set; requiring every available set to agree
    # in the expected direction is what tells real sexes from noise.
    dgap <- function(s, high) {
        if (all(is.na(s))) return(NA_real_)
        low <- setdiff(c("female", "male"), high)
        mean(s[sex == high], na.rm = TRUE) - mean(s[sex == low], na.rm = TRUE)
    }
    gaps <- c(female = dgap(tbl$female_score, "female"),
              male = dgap(tbl$male_score, "male"))
    gaps <- gaps[!is.na(gaps)]
    if (length(unique(sex)) < 2L || !length(gaps) || min(gaps) < min_gap) {
        out$status <- "undetermined"
        out$message <- paste0(
            "Sex markers do not separate the samples (cluster difference ",
            paste(names(gaps), round(gaps, 2), collapse = ", "), " log2; need >= ", min_gap,
            " on each); the cohort is probably single-sex. Sex is not used as a covariate.")
        return(out)
    }
    out$table$imputed_sex <- sex
    out$status <- "imputed"
    out$message <- paste0("Imputed ", sum(sex == "female"), " female and ",
                          sum(sex == "male"), " male sample(s) (marker difference: ",
                          paste(names(gaps), round(gaps, 2), collapse = ", "), " log2).")
    out
}

#' Model covariates for a run: `batch` / `pair`, plus imputed sex when estimable
#'
#' Sex is added when `colData(dds)$imputed_sex` has two levels with no
#' missing values and adding it raises the rank of every design the run fits
#' (additive, and the product design when interactions are on) by one. If
#' sex is confounded with something already in the design -- `var` itself, a
#' `pair` whose units each have one sex, or a `batch` -- it is left out.
#'
#' @param dds QC-filtered `DESeqDataSet`.
#' @param var,time Tested factor columns (`time` for LRT).
#' @param covars Blocking covariates already in the model (`c(batch, pair)`).
#' @param include_interactions Whether the run fits a product design.
#' @return A list: `covars` (character vector of covariate columns, or
#'   `NULL`), `sex_used` (logical) and `sex_reason` (why sex was or was not
#'   used).
#' @noRd
.model_covariates <- function(dds, var, covars = NULL, include_interactions = FALSE,
                              time = NULL) {
    keep <- covars
    reason <- NULL
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    sx <- cd[[.sex_col]]
    tested <- c(var, time)
    if (is.null(sx) || anyNA(sx)) {
        reason <- "not used: sex was not imputed"
    } else if (.sex_col %in% c(tested, covars)) {
        reason <- "not used: already in the design"
    } else if (length(unique(sx)) < 2L) {
        reason <- "not used: all samples have the same sex"
    } else {
        terms <- c(list(paste(tested, collapse = " + ")),
                   if (isTRUE(include_interactions) && length(tested) > 1L) {
                       list(paste(tested, collapse = " * "))
                   })
        estimable <- vapply(terms, function(rhs) {
            .full_rank(cd, c(covars, .sex_col, rhs))
        }, logical(1))
        if (all(estimable)) {
            keep <- c(covars, .sex_col)
            reason <- "used: added to every model"
        } else {
            reason <- "not used: confounded with variables already in the design"
        }
    }
    message("Imputed sex as a model covariate -- ", reason, ".")
    list(covars = keep, sex_used = .sex_col %in% keep, sex_reason = reason)
}

#' Drop covariates that cannot be estimated within a subset of samples
#'
#' Within-stratum fits use a subset, where a covariate that is fine on all
#' samples can be constant or aliased (e.g. every sample of the stratum has
#' the same sex). Those are dropped, with a message, rather than failing the
#' fit.
#' @noRd
.estimable_covariates <- function(dds, covars, tested) {
    if (!length(covars)) return(covars)
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    keep <- character(0)
    for (b in covars) {
        if (length(unique(as.character(cd[[b]]))) < 2L) next
        if (.full_rank(cd, c(keep, b, tested))) keep <- c(keep, b)
    }
    dropped <- setdiff(covars, keep)
    if (length(dropped)) {
        message("  Covariate(s) not estimable in this subset, dropped: ",
                paste(dropped, collapse = ", "))
    }
    if (length(keep)) keep else NULL
}

#' Is the additive model `~ terms` full rank on these samples?
#' @noRd
.full_rank <- function(cd, terms) {
    cd[] <- lapply(cd, function(x) if (is.character(x)) factor(x) else x)
    for (nm in names(cd)) if (is.factor(cd[[nm]])) cd[[nm]] <- droplevels(cd[[nm]])
    mm <- stats::model.matrix(stats::as.formula(paste("~", paste(terms, collapse = " + "))), cd)
    qr(mm)$rank == ncol(mm)
}

#' Plot female- vs male-marker expression per sample
#'
#' One point per sample: mean log2(CPM + 1) of the female markers (x) against
#' the male markers (y), coloured by imputed sex. Two well-separated groups
#' (top-left and bottom-right) mean the imputation is clear-cut; a single
#' cloud means a single-sex cohort.
#'
#' @param sex_table The `table` element of [impute_sample_sex()].
#' @return A `ggplot`, or `NULL` when no marker scores are available.
#' @export
plot_sex_markers <- function(sex_table) {
    df <- as.data.frame(sex_table)
    if (!nrow(df) || all(is.na(df$female_score)) || all(is.na(df$male_score))) return(NULL)
    df$sex <- ifelse(is.na(df$imputed_sex), "undetermined", df$imputed_sex)
    df$sex <- factor(df$sex, levels = c("female", "male", "undetermined"))
    ggplot2::ggplot(df, ggplot2::aes(x = .data$female_score, y = .data$male_score,
                                     colour = .data$sex,
                                     text = paste0("Sample: ", .data$sample))) +
        ggplot2::geom_point(size = 3, alpha = 0.85) +
        ggplot2::scale_colour_manual(
            values = c(female = "#D6604D", male = "#4393C3", undetermined = "grey50"),
            drop = TRUE, name = "Imputed sex") +
        ggplot2::labs(x = "Female markers, mean log2(CPM + 1)",
                      y = "Male (Y) markers, mean log2(CPM + 1)") +
        ggplot2::theme_bw(base_size = 12) +
        square_panel_theme()
}

#' `imputed_sex` if it holds a usable (non-missing, two-level) sex, else `NULL`
#' @noRd
.sex_assoc_col <- function(dds) {
    x <- SummarizedExperiment::colData(dds)[[.sex_col]]
    if (is.null(x) || anyNA(x) || length(unique(as.character(x))) < 2L) return(NULL)
    .sex_col
}

#' Variables to colour the PCA by, one report tab each
#'
#' Group (the visualisation factor), then pair, imputed sex and each batch
#' column, when present.
#' @noRd
.pca_colour_vars <- function(dds, viz_var, pair = NULL, batch = NULL) {
    out <- c(Group = viz_var)
    if (!is.null(pair)) out <- c(out, Pair = pair)
    sx <- SummarizedExperiment::colData(dds)[[.sex_col]]
    if (!is.null(sx) && !all(is.na(sx))) out <- c(out, Sex = .sex_col)
    if (length(batch) == 1L) {
        out <- c(out, Batch = batch)
    } else if (length(batch) > 1L) {
        out <- c(out, stats::setNames(batch, paste0("Batch (", batch, ")")))
    }
    out
}
