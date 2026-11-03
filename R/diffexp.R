#' Run differential expression for a single comparison
#'
#' Fits (or reuses a pre-fit) DESeq2 model and extracts results for one
#' comparison string (`"A_vs_B"` or a pooled `"A_vs_B+C"` one-vs-rest
#' comparison). For pooled comparisons, group levels are merged by exact
#' match (never regex/substring matching) on a private copy of `dds`, so the
#' caller's original grouping factor and any other in-flight comparisons are
#' never affected. Before fitting, the design's model matrix is checked for
#' full column rank via [assert_full_rank_model_matrix()] and a clear error
#' is raised early if it is not.
#'
#' @param dds A `DESeqDataSet` object, QC-filtered and with `design(dds)`
#'   already set (typically the output of [run_qc()]).
#' @param var Column name in `colData(dds)` used for grouping.
#' @param comparison A single comparison string, e.g. `"Treated_vs_Control"`
#'   or `"A_vs_B+C"`.
#' @param shrinkage LFC shrinkage method: one of `"ashr"` (default; the only
#'   method that supports the pooled contrasts this package generates),
#'   `"apeglm"` (only valid for simple two-group contrasts, not pooled ones),
#'   `"normal"`, or `"none"`.
#' @param alpha Significance threshold (adjusted p-value) used for the
#'   independent filtering / `results()` call. Default `0.05`.
#' @param lfc_threshold Log2 fold-change threshold used only to flag
#'   "significant" genes in the returned table (`sig = padj < alpha &
#'   abs(log2FoldChange) > lfc_threshold`); does not affect the DESeq2 test
#'   itself. Default `0`.
#' @param ranking Ranking metric added as a `rank_score` column, via
#'   [compute_ranking_score()]. Default `"log2FoldChange"`.
#' @param gene_col Name to give the gene-identifier column in the output
#'   (row names of `dds` are used as the values). Default `"gene"`.
#' @return A data frame of DESeq2 results (one row per gene) with additional
#'   columns `comparison`, `mean_cpm` (mean counts-per-million across all
#'   samples in `dds`, for use as the MA plot's x-axis in place of
#'   `baseMean`), `rank_score`, and `sig`.
#' @examples
#' \dontrun{
#' res <- run_diffexp(dds, var = "condition", comparison = "Treated_vs_Control")
#' }
#' @export
run_diffexp <- function(
    dds,
    var,
    comparison,
    shrinkage = "ashr",
    alpha = 0.05,
    lfc_threshold = 0.5,
    ranking = "log2FoldChange",
    gene_col = "gene") {

    dds <- validate_dds(dds)
    dds <- validate_var(dds, var)
    shrinkage <- validate_shrinkage(shrinkage)
    alpha <- validate_alpha(alpha)
    lfc_threshold <- validate_lfc_threshold(lfc_threshold)
    ranking <- validate_ranking(ranking)
    levels_var <- levels(SummarizedExperiment::colData(dds)[[var]])
    validate_comparison(comparison, levels_var)

    pooled <- grepl("+", comparison, fixed = TRUE)
    if (pooled && shrinkage == "apeglm") {
        stop("`shrinkage = \"apeglm\"` does not support pooled one-vs-rest comparisons ",
             "(\"", comparison, "\"); use \"ashr\", \"normal\", or \"none\" instead.", call. = FALSE)
    }

    if (pooled) {
        merged <- merge_comparison_levels(dds, var, comparison)
        fit_dds <- merged$dds
        fit_var <- merged$merged_var
        test_level <- merged$test
        ref_level <- merged$ref
        # rebuild the design formula, substituting the original grouping
        # term for the merged one -- done on the parsed formula terms
        # (not string substring replacement) to avoid mismatches when `var`
        # is itself a substring of another design term.
        orig_design <- DESeq2::design(fit_dds)
        design_terms <- all.vars(orig_design)
        new_terms <- ifelse(design_terms == var, fit_var, design_terms)
        rhs <- paste(new_terms, collapse = " + ")
        new_design <- stats::as.formula(paste("~", rhs))
        DESeq2::design(fit_dds) <- new_design
    } else {
        parts <- strsplit(comparison, "_vs_", fixed = TRUE)[[1]]
        test_level <- parts[1]
        ref_level <- parts[2]
        fit_dds <- dds
        fit_var <- var
    }

    assert_full_rank_model_matrix(fit_dds, comparison, var)

    fit_dds <- DESeq2::DESeq(fit_dds, quiet = TRUE)
    contrast <- c(fit_var, test_level, ref_level)

    res <- DESeq2::results(fit_dds, contrast = contrast, alpha = alpha)
    # the Wald `stat` column is preserved from the unshrunken results and
    # re-attached below -- `lfcShrink()` (all types except "normal") drops
    # it from its return value, but it's needed downstream for the "stat"
    # ranking metric and for user inspection. Keep the unshrunken LFC too
    # for the volcano x-axis (shrunk values stay in `log2FoldChange`).
    stat_col <- res$stat
    unshrunken_lfc <- stats::setNames(as.numeric(res$log2FoldChange), rownames(res))
    if (shrinkage != "none") {
        res <- tryCatch({
            if (shrinkage == "apeglm") {
                coef_name <- paste0(fit_var, "_", test_level, "_vs_", ref_level)
                DESeq2::lfcShrink(fit_dds, coef = coef_name, type = "apeglm", quiet = TRUE)
            } else {
                DESeq2::lfcShrink(fit_dds, contrast = contrast, type = shrinkage, res = res, quiet = TRUE)
            }
        }, error = function(e) {
            warning("LFC shrinkage (\"", shrinkage, "\") failed for comparison \"", comparison,
                    "\": ", conditionMessage(e), ". Returning unshrunken results.", call. = FALSE)
            res
        })
    }

    res_df <- as.data.frame(res)
    if (!"stat" %in% colnames(res_df)) {
        res_df$stat <- stat_col
    }
    res_df$log2FoldChange_unshrunken <- unshrunken_lfc[rownames(res_df)]
    res_df[[gene_col]] <- rownames(res_df)
    rownames(res_df) <- NULL
    res_df$comparison <- comparison
    res_df$rank_score <- compute_ranking_score(res_df, ranking)
    res_df$sig <- !is.na(res_df$padj) & res_df$padj < alpha & abs(res_df$log2FoldChange) > lfc_threshold

    # Mean CPM across *all* samples in the original (pre-merge) `dds`, for
    # the MA plot's x-axis -- library-size-normalised and on a scale
    # comparable across comparisons/datasets, unlike DESeq2's own
    # `baseMean` (the mean of size-factor-normalised counts, which is
    # comparison-specific for pooled contrasts since `fit_dds`'s size
    # factors are re-estimated on the merged grouping).
    raw_counts <- DESeq2::counts(dds)
    lib_sizes <- colSums(raw_counts)
    cpm_mat <- sweep(raw_counts, 2, lib_sizes, "/") * 1e6
    mean_cpm <- rowMeans(cpm_mat)
    res_df$mean_cpm <- mean_cpm[res_df[[gene_col]]]

    res_df[, c(gene_col, "comparison", "baseMean", "mean_cpm", "log2FoldChange",
               "log2FoldChange_unshrunken", "lfcSE", "stat",
               "pvalue", "padj", "rank_score", "sig")]
}

#' Stop early when a cell of the `var` factorial has no samples
#'
#' An interaction contrast needs the full product design, which is
#' rank-deficient as soon as one level combination is unobserved. Naming the
#' empty cells is far more actionable than the model-matrix rank message.
#' @noRd
.assert_full_factorial <- function(dds, var) {
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    fac <- lapply(var, function(v) droplevels(factor(cd[[v]])))
    names(fac) <- var
    cells <- table(fac)
    if (all(cells > 0)) return(invisible(TRUE))
    grid <- expand.grid(dimnames(cells), stringsAsFactors = FALSE)
    empty <- grid[as.vector(cells) == 0, , drop = FALSE]
    named <- apply(empty, 1, function(row) {
        paste(paste0(var, "=", row), collapse = ", ")
    })
    stop("Interaction contrasts need every combination of ",
         paste(var, collapse = " \u00d7 "), " to have samples, but ",
         nrow(empty), " combination(s) are empty: ",
         paste(named, collapse = "; "),
         ". Drop a factor from `var`, or run without `include_interactions`.",
         call. = FALSE)
}

# `.clean_coef_names()` used to live here; it is now in R/contrast_spec.R with
# the other spec-field codecs.

#' Run differential expression on a combination of named model coefficients
#'
#' Used for interaction contrasts from `~ A * B [* C]`. A single `coef` with
#' nothing subtracted is extracted by name; several are compared as a list
#' contrast, `results(contrast = list(coef, coef_subtract))`, summing the
#' first vector's coefficients against the second's. `coef_weights` covers
#' the contrasts a list cannot express -- those that *average* coefficients,
#' as a pooled one-vs-rest comparison does -- by passing DESeq2 a numeric
#' weight vector over `resultsNames()` instead. Same return schema as
#' [run_diffexp()].
#'
#' @inheritParams run_diffexp
#' @param coef Coefficient name(s) as in `resultsNames(dds)`; several names
#'   are summed. `NULL` when `coef_weights` is given.
#' @param coef_subtract Optional coefficient name(s) whose sum is subtracted
#'   from `coef`. Default `NULL`.
#' @param coef_weights Optional named numeric vector (or its
#'   `"name=weight;..."` encoding) of coefficient weights, used instead of
#'   `coef` / `coef_subtract`. Default `NULL`.
#' @param prefit Optional already-`DESeq()`-ed object for `dds`, so a grid of
#'   contrasts sharing one design is fitted once. Default `NULL` (fit here).
#' @param comparison Name recorded in the `comparison` column of the result.
#' @return A results data frame, with the same columns as [run_diffexp()].
#' @examples
#' \dontrun{
#' dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, ~ tissue * genotype)
#' DESeq2::resultsNames(DESeq2::DESeq(dds))
#' # Does the genotype effect differ between tissues?
#' res <- run_diffexp_coef(dds, coef = "tissueSP.genotypeKO")
#' # The average of two interaction terms against a third:
#' res2 <- run_diffexp_coef(
#'   dds,
#'   coef_weights = c(tissueSP.genotypeKO = 1,
#'                    tissueLN.genotypeKO = -0.5,
#'                    tissueBM.genotypeKO = -0.5)
#' )
#' }
#' @seealso [contrast_spec_from_coefs()] to put the same contrast into a grid
#'   that [run_diffexp_pipeline()] runs.
#' @export
run_diffexp_coef <- function(
    dds,
    coef = NULL,
    coef_subtract = NULL,
    coef_weights = NULL,
    shrinkage = "ashr",
    alpha = 0.05,
    lfc_threshold = 0.5,
    ranking = "log2FoldChange",
    gene_col = "gene",
    comparison = paste(coef, collapse = "+"),
    prefit = NULL) {

    dds <- validate_dds(dds)
    shrinkage <- validate_shrinkage(shrinkage)
    alpha <- validate_alpha(alpha)
    lfc_threshold <- validate_lfc_threshold(lfc_threshold)
    ranking <- validate_ranking(ranking)

    coef <- .clean_coef_names(coef)
    coef_subtract <- .clean_coef_names(coef_subtract)
    if (is.character(coef_weights)) coef_weights <- .parse_coef_weights(coef_weights)
    weighted <- length(coef_weights) > 0L
    if (weighted && (length(coef) || length(coef_subtract))) {
        stop("Contrast \"", comparison, "\": give either `coef` / ",
             "`coef_subtract` or `coef_weights`, not both.", call. = FALSE)
    }
    if (weighted && (is.null(names(coef_weights)) || anyNA(coef_weights))) {
        stop("Contrast \"", comparison,
             "\": `coef_weights` must be a fully named, non-`NA` numeric vector.",
             call. = FALSE)
    }
    if (!weighted && !length(coef)) {
        stop("Contrast \"", comparison, "\": no numerator coefficient given.",
             call. = FALSE)
    }

    wanted <- if (weighted) names(coef_weights) else c(coef, coef_subtract)
    if (is.null(prefit)) {
        assert_full_rank_model_matrix(dds, comparison, wanted[1])
        fit_dds <- DESeq2::DESeq(dds, quiet = TRUE)
    } else {
        fit_dds <- prefit
    }
    rn <- DESeq2::resultsNames(fit_dds)
    unknown <- setdiff(wanted, rn)
    if (length(unknown)) {
        stop("Coefficient(s) ", paste0("\"", unknown, "\"", collapse = ", "),
             " not in resultsNames() for contrast \"", comparison, "\". Available: ",
             paste(rn, collapse = ", "), call. = FALSE)
    }

    # A weight vector spans every coefficient in resultsNames() order; a list
    # contrast names only the coefficients it combines.
    contrast_arg <- if (weighted) {
        w <- stats::setNames(numeric(length(rn)), rn)
        w[names(coef_weights)] <- unname(coef_weights)
        unname(w)
    } else if (length(coef) > 1L || length(coef_subtract)) {
        if (length(coef_subtract)) list(coef, coef_subtract) else list(coef)
    } else {
        NULL
    }

    if (!is.null(contrast_arg)) {
        res <- DESeq2::results(fit_dds, contrast = contrast_arg, alpha = alpha)
        shrink_type <- shrinkage
        if (shrinkage == "apeglm") {
            warning("apeglm cannot shrink a combined contrast; using ashr for \"",
                    comparison, "\".", call. = FALSE)
            shrink_type <- "ashr"
        }
        stat_col <- res$stat
        unshrunken_lfc <- stats::setNames(as.numeric(res$log2FoldChange), rownames(res))
        if (shrink_type != "none") {
            res <- tryCatch({
                DESeq2::lfcShrink(fit_dds, contrast = contrast_arg, type = shrink_type,
                                  res = res, quiet = TRUE)
            }, error = function(e) {
                warning("LFC shrinkage (\"", shrink_type, "\") failed for \"", comparison,
                        "\": ", conditionMessage(e), ". Returning unshrunken results.",
                        call. = FALSE)
                res
            })
        }
    } else {
        res <- DESeq2::results(fit_dds, name = coef, alpha = alpha)
        stat_col <- res$stat
        unshrunken_lfc <- stats::setNames(as.numeric(res$log2FoldChange), rownames(res))
        if (shrinkage != "none") {
            res <- tryCatch({
                DESeq2::lfcShrink(fit_dds, coef = coef, type = shrinkage, res = res, quiet = TRUE)
            }, error = function(e) {
                warning("LFC shrinkage (\"", shrinkage, "\") failed for \"", comparison,
                        "\": ", conditionMessage(e), ". Returning unshrunken results.",
                        call. = FALSE)
                res
            })
        }
    }

    res_df <- as.data.frame(res)
    if (!"stat" %in% colnames(res_df)) res_df$stat <- stat_col
    res_df$log2FoldChange_unshrunken <- unshrunken_lfc[rownames(res_df)]
    res_df[[gene_col]] <- rownames(res_df)
    rownames(res_df) <- NULL
    res_df$comparison <- comparison
    res_df$rank_score <- compute_ranking_score(res_df, ranking)
    res_df$sig <- !is.na(res_df$padj) & res_df$padj < alpha & abs(res_df$log2FoldChange) > lfc_threshold

    raw_counts <- DESeq2::counts(dds)
    lib_sizes <- colSums(raw_counts)
    cpm_mat <- sweep(raw_counts, 2, lib_sizes, "/") * 1e6
    res_df$mean_cpm <- rowMeans(cpm_mat)[res_df[[gene_col]]]

    res_df[, c(gene_col, "comparison", "baseMean", "mean_cpm", "log2FoldChange",
               "log2FoldChange_unshrunken", "lfcSE", "stat",
               "pvalue", "padj", "rank_score", "sig")]
}

#' Run differential expression for all requested comparisons
#'
#' Loops [run_diffexp()] over every comparison in `comparisons` (or every
#' auto-generated comparison if `comparisons` is `NULL`), writes one results
#' CSV per comparison plus MA and volcano plots, and returns all results
#' combined. With `var` of length 2 or 3 the grid comes from
#' [expand_contrast_specs()], unless `comparisons` subsets those labels.
#'
#' How a spec row is fitted depends on its `class`:
#' \describe{
#'   \item{`"single"`}{The comparison on `dds` as handed in.}
#'   \item{`"main"`}{All samples, design rebuilt additively as
#'     `~ [batch +] other factors + tested factor`, so a factor's main effect
#'     is adjusted for the other factors rather than measured at their
#'     reference levels (which is what an interaction design's main-effect
#'     coefficient would give).}
#'   \item{`"interaction"`}{A difference of differences on all samples: the
#'     full product design (`~ [batch +] A * B [* C]`) is fitted once for
#'     the whole interaction grid, and each row compares the sum of its
#'     `comparison` coefficients against the sum of its `coef_subtract`
#'     ones — or, where a pooled one-vs-rest comparison needs coefficients
#'     averaged, weights them by its `coef_weights`. Fails naming the empty
#'     cells if the factorial is incomplete, since the product design would
#'     then be rank-deficient.}
#'   \item{`"within_*"`}{The samples in that stratum only, design
#'     `~ [batch +] tested factor`.}
#' }
#'
#' @inheritParams run_diffexp
#' @param comparisons Character vector of comparisons / contrast labels to
#'   run, a contrast-specs data frame from [expand_contrast_specs()], or
#'   `NULL` (default) to auto-generate the full grid.
#' @param batch Optional covariate column name(s) kept in the design when it
#'   is rebuilt for a main effect or a within-stratum subset. Default `NULL`.
#' @param include_interactions Logical. When `var` has length ≥ 2, add the
#'   difference-of-differences interaction contrasts to the main effects.
#'   Default `FALSE`. Passed to [expand_contrast_specs()].
#' @param stratified Logical. Add the within-stratum contrasts, fitted on
#'   each subset alone. Default `FALSE`; not exposed by [run_deseq2_wald()].
#'   Passed to [expand_contrast_specs()].
#' @param top_label_n Number of top genes (by `rank_score`) to label on the
#'   volcano plot. Default `10`.
#' @param save_dir Directory to write per-comparison subfolders to. Default
#'   `getwd()`.
#' @param save_plot Logical. If `TRUE`, saves MA and volcano plots per
#'   comparison. Default `TRUE`.
#' @return A list with elements `results` (combined data frame across all
#'   comparisons, with a `comparison` column), `by_comparison` (a named
#'   list of per-comparison result data frames), `contrast_specs` (the
#'   specs data frame used), and `manifest` ([contrast_manifest()]: what each
#'   contrast asks, the design fitted, the samples behind it and the
#'   coefficients combined, also written to `contrasts.csv`).
#' @examples
#' \dontrun{
#' res <- run_diffexp_pipeline(dds, var = "condition")
#' }
#' @export
run_diffexp_pipeline <- function(
    dds,
    var,
    comparisons = NULL,
    batch = NULL,
    include_interactions = FALSE,
    stratified = FALSE,
    shrinkage = "ashr",
    alpha = 0.05,
    lfc_threshold = 0.5,
    ranking = "log2FoldChange",
    gene_col = "gene",
    top_label_n = 10,
    save_dir = getwd(),
    save_plot = TRUE) {

    dds <- validate_dds(dds)
    dds <- validate_var(dds, var)
    comparisons <- validate_comparisons(
        dds, var, comparisons,
        include_interactions = include_interactions,
        stratified = stratified
    )
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)

    # If `batch` was requested but is missing from design(dds) (e.g. the
    # caller built the DESeqDataSet with `~ var` only), rebuild so every
    # contrast — not just within-stratum subsets — is `~ batch + var`.
    if (!is.null(batch)) {
        des_vars <- all.vars(DESeq2::design(dds))
        if (!all(as.character(batch) %in% des_vars)) {
            DESeq2::design(dds) <- resolve_design(
                var, batch = batch, design = NULL,
                include_interactions = include_interactions
            )
        }
    }

    # One spec table, whatever the caller passed: a character vector of
    # comparisons, an auto-generated grid, or a hand-built table.
    specs <- as_contrast_specs(comparisons, var = var[1])
    if (!"coef_subtract" %in% names(specs)) specs$coef_subtract <- NA_character_
    specs <- validate_contrast_specs(specs, dds = dds, var = var)

    # Every interaction row is a different contrast of the *same* product
    # model, so fit it once rather than per row.
    ix_fit <- NULL
    if (any(specs$fit_method %in% "coefs")) {
        .assert_full_factorial(dds, var)
        dds_ix <- dds
        DESeq2::design(dds_ix) <- resolve_design(
            var, batch = batch, design = NULL, include_interactions = TRUE
        )
        assert_full_rank_model_matrix(dds_ix, "interaction contrasts",
                                      paste(var, collapse = " * "))
        message("Fitting the interaction model (", deparse1(DESeq2::design(dds_ix)),
                ") ...")
        ix_fit <- DESeq2::DESeq(dds_ix, quiet = TRUE)
    }

    # Say what is about to be compared, and leave the same statement on disk
    # next to the results.
    manifest <- contrast_manifest(specs, dds = dds, var = var, batch = batch)
    .announce_contrasts(manifest)
    write_output(manifest, file.path(save_dir, "contrasts.csv"), format = "csv")

    by_comparison <- list()
    for (i in seq_len(nrow(specs))) {
        spec <- specs[i, , drop = FALSE]
        label <- spec$label[[1]]
        message("Running differential expression: ", label, " ...")

        if (identical(spec$fit_method[[1]], "coefs")) {
            # A pooled one-vs-rest comparison averages coefficients, which
            # only the weight vector can express; the rest sum them.
            coefs <- contrast_coefs(spec)
            weighted <- length(coefs$weights) > 0L
            res_df <- run_diffexp_coef(
                dds = dds,
                coef = if (weighted) NULL else coefs$pos,
                coef_subtract = if (weighted) NULL else coefs$neg,
                coef_weights = if (weighted) coefs$weights else NULL,
                shrinkage = shrinkage, alpha = alpha, lfc_threshold = lfc_threshold,
                ranking = ranking, gene_col = gene_col, comparison = label,
                prefit = ix_fit
            )
        } else if (identical(spec$fit_method[[1]], "levels_additive")) {
            # Main effect: every sample, additive design, so the factors not
            # under test are adjusted for instead of inflating the residual.
            tv <- spec$test_var[[1]]
            dds_main <- dds
            DESeq2::design(dds_main) <- resolve_design(
                c(setdiff(var, tv), tv), batch = batch, design = NULL,
                include_interactions = FALSE
            )
            res_df <- run_diffexp(
                dds = dds_main, var = tv, comparison = spec$comparison[[1]],
                shrinkage = shrinkage, alpha = alpha, lfc_threshold = lfc_threshold,
                ranking = ranking, gene_col = gene_col
            )
            res_df$comparison <- label
        } else if (identical(spec$fit_method[[1]], "levels_subset")) {
            stratum <- contrast_stratum(spec)
            if (!length(stratum)) {
                stop("Contrast \"", label, "\": fit_method \"levels_subset\" needs ",
                     "`stratum_var` / `stratum_level`.", call. = FALSE)
            }
            tv <- contrast_test_vars(spec)[[1]]
            cd <- SummarizedExperiment::colData(dds)
            keep <- rep(TRUE, ncol(dds))
            for (sv in names(stratum)) {
                keep <- keep & as.character(cd[[sv]]) == stratum[[sv]]
            }
            dds_sub <- dds[, keep]
            DESeq2::design(dds_sub) <- resolve_design(tv, batch = batch, design = NULL)
            res_df <- run_diffexp(
                dds = dds_sub, var = tv, comparison = spec$comparison[[1]],
                shrinkage = shrinkage, alpha = alpha, lfc_threshold = lfc_threshold,
                ranking = ranking, gene_col = gene_col
            )
            res_df$comparison <- label
        } else {
            tv <- contrast_test_vars(spec)
            tv <- if (length(tv)) tv[[1]] else var[1]
            res_df <- run_diffexp(
                dds = dds, var = tv, comparison = spec$comparison[[1]],
                shrinkage = shrinkage, alpha = alpha, lfc_threshold = lfc_threshold,
                ranking = ranking, gene_col = gene_col
            )
            res_df$comparison <- label
        }

        # `label` is the contrast's identity (report tab, by_comparison key);
        # `path` is where it is written. They coincide unless a hand-written
        # spec overrides `path`.
        comp_dir <- file.path(save_dir, spec$path[[1]])
        stem <- contrast_file_stem(spec$path[[1]])
        write_output(res_df, file.path(comp_dir, paste0(stem, "_diffexp.csv")), format = "csv")

        if (save_plot) {
            ma_plot <- plot_ma(res_df, alpha = alpha, lfc_threshold = lfc_threshold)
            save_ggplot(ma_plot, file.path(comp_dir, paste0(stem, "_ma_plot.pdf")), width = 6, height = 5)

            volcano_plot <- plot_volcano(res_df, alpha = alpha, lfc_threshold = lfc_threshold,
                                          gene_col = gene_col, top_label_n = top_label_n)
            save_ggplot(volcano_plot, file.path(comp_dir, paste0(stem, "_volcano_plot.pdf")), width = 7, height = 6)
        }

        by_comparison[[label]] <- res_df
    }

    combined <- do.call(rbind, by_comparison)
    rownames(combined) <- NULL
    list(results = combined, by_comparison = by_comparison, contrast_specs = specs,
         manifest = manifest)
}

#' Classify genes as up/down/not-significant for direction-coloured plots
#'
#' @param res_df A data frame as returned by [run_diffexp()] (must contain
#'   `padj` and `log2FoldChange`; uses the pre-computed `sig` column when
#'   present, else recomputes it from `alpha`/`lfc_threshold`).
#' @param alpha Significance threshold, used only if `res_df$sig` is absent.
#' @param lfc_threshold LFC threshold, used only if `res_df$sig` is absent.
#' @return Factor with levels `"up"`, `"down"`, `"ns"` (in that order).
#' @noRd
classify_de_direction <- function(res_df, alpha = 0.05, lfc_threshold = 0) {
    # Always recompute from padj + |LFC| so plot-level `lfc_threshold`
    # (e.g. volcano colouring at 0.5) is honoured even when `res_df$sig`
    # was produced under a different threshold.
    sig <- !is.na(res_df$padj) & res_df$padj < alpha &
        abs(res_df$log2FoldChange) > lfc_threshold
    sig[is.na(sig)] <- FALSE
    direction <- ifelse(sig & res_df$log2FoldChange > 0, "up",
                  ifelse(sig & res_df$log2FoldChange < 0, "down", "ns"))
    factor(direction, levels = c("up", "down", "ns"))
}

#' Robust symmetric |log2FC| axis limit that only trims *super*-extreme
#' outliers (e.g. one gene at |LFC| ≫ the rest of the cloud).
#'
#' Uses a high quantile of `|lfc|`; if the absolute maximum is at least
#' `outlier_factor` times that bulk quantile, the quantile is returned so
#' the axis is not blown out by a handful of extremes. Mild outliers are
#' left alone. Squished points still draw at the panel edge via
#' `scales::oob_squish`.
#' @noRd
.robust_lfc_limit <- function(lfc, bulk_quantile = 0.999, outlier_factor = 5) {
    abs_lfc <- abs(lfc[!is.na(lfc)])
    if (!length(abs_lfc)) {
        return(1)
    }
    max_lfc <- max(abs_lfc)
    bulk <- as.numeric(stats::quantile(abs_lfc, probs = bulk_quantile, names = FALSE, type = 7))
    # Only crop truly extreme outliers (e.g. |LFC|=40 while 99.9% sit below ~6).
    if (is.finite(bulk) && bulk > 0 && max_lfc >= outlier_factor * bulk) {
        return(max(bulk, 1))
    }
    max_lfc
}

#' Parse test/reference level labels from a contrast comparison string
#'
#' Handles plain `"A_vs_B"`, stratum-prefixed `"BM__IgG_vs_IgM"`, and
#' interaction labels (falls back to positive/negative). An interaction label
#' such as `"C_vs_noC_x_K_vs_noK_in_CB"` does contain `"_vs_"`, but neither
#' side of it is a group of samples, so it gets the neutral wording too.
#' @noRd
.parse_contrast_levels <- function(comparison) {
    comparison <- as.character(comparison)[1]
    if (is.na(comparison) || !nzchar(comparison)) {
        return(list(test = "up", ref = "down"))
    }
    if (!grepl("_vs_", comparison, fixed = TRUE) ||
        grepl("_x_", comparison, fixed = TRUE)) {
        return(list(test = "positive LFC", ref = "negative LFC"))
    }
    core <- if (grepl("__", comparison, fixed = TRUE)) {
        sub("^.*?__", "", comparison)
    } else {
        comparison
    }
    parts <- strsplit(core, "_vs_", fixed = TRUE)[[1]]
    list(test = parts[1], ref = if (length(parts) >= 2) parts[2] else "ref")
}

# Shared direction -> colour/size mapping for MA and volcano plots: DE-up
# genes red, DE-down genes blue, non-significant genes grey, with
# significant points drawn slightly larger so they stand out.
.de_direction_colours <- c(up = "#B2182B", down = "#2166AC", ns = "grey70")
.de_direction_sizes <- c(up = 1.6, down = 1.6, ns = 0.8)

#' MA plot (mean expression vs. log2 fold-change)
#'
#' Points are coloured by direction: red = significantly up, blue =
#' significantly down, grey = not significant (using the `sig` column
#' already computed by [run_diffexp()]); significant points are drawn
#' slightly larger. The number of up- and down-regulated genes is
#' annotated directly on the plot. A `text` aesthetic (gene name, log2
#' fold-change, adjusted p-value) is included on every point so the plot
#' can be converted to an interactive, hoverable widget with
#' [make_interactive()].
#'
#' The x-axis is `log2(mean CPM + 1)` (mean counts-per-million across all
#' samples, from the `mean_cpm` column added by [run_diffexp()]) rather
#' than DESeq2's own `baseMean` -- CPM is library-size-normalised on a
#' scale that stays comparable across comparisons and datasets, whereas
#' `baseMean` is the mean of size-factor-normalised counts and so is
#' re-scaled per comparison for pooled one-vs-rest contrasts.
#'
#' @param res_df A data frame as returned by [run_diffexp()] (must contain
#'   a `mean_cpm` column).
#' @param alpha Significance threshold used to colour points (only used if
#'   `res_df` lacks a `sig` column). Default `0.05`.
#' @param gene_col Name of the gene-identifier column, used for the hover
#'   tooltip. Default `"gene"`.
#' @return A `ggplot` object.
#' @export
plot_ma <- function(res_df, alpha = 0.05, lfc_threshold = 0.5, gene_col = "gene") {
    if (!"mean_cpm" %in% colnames(res_df)) {
        stop("`res_df` must contain a `mean_cpm` column (as produced by run_diffexp()).", call. = FALSE)
    }
    res_df$direction <- classify_de_direction(res_df, alpha = alpha, lfc_threshold = lfc_threshold)
    n_up <- sum(res_df$direction == "up")
    n_down <- sum(res_df$direction == "down")

    res_df$log2_mean_cpm <- log2(res_df$mean_cpm + 1)

    res_df$.tooltip <- paste0(
        "Gene: ", res_df[[gene_col]],
        "<br>log2FC: ", formatC(res_df$log2FoldChange, format = "f", digits = 2),
        "<br>padj: ", formatC(res_df$padj, format = "e", digits = 2)
    )

    comparison <- unique(res_df$comparison)
    levels <- .parse_contrast_levels(comparison)
    test_level <- levels$test
    ref_level <- levels$ref

    # Draw ns first so grey points sit behind red/blue.
    res_df <- res_df[order(match(as.character(res_df$direction), c("ns", "up", "down"))), ]

    p <- ggplot2::ggplot(res_df, ggplot2::aes(x = .data$log2_mean_cpm, y = .data$log2FoldChange,
                                          colour = .data$direction, size = .data$direction,
                                          text = .data$.tooltip)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_x_continuous(expand = zero_expansion()) +
        # A little headroom (10% of the data range on each side) rather
        # than zero data-padding, so points near the log2FC extremes --
        # and the up/down count and direction-arrow annotations anchored
        # near them -- aren't crowded against the panel border.
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
        ggplot2::scale_colour_manual(values = .de_direction_colours, drop = FALSE,
                                      breaks = c("up", "down", "ns")) +
        ggplot2::scale_size_manual(values = .de_direction_sizes, guide = "none", drop = FALSE) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
        # Counts (inside panel) + up/down arrows on the *right* naming which
        # group is higher in each direction. (Also mirrored to plotly via
        # attr below -- ggplotly drops annotate().)
        ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.15, vjust = 1.5,
                           label = paste0("\u2191 ", test_level, "  (n=", n_up, ")"),
                           colour = .de_direction_colours[["up"]], size = 3.5, fontface = "bold") +
        ggplot2::annotate("text", x = Inf, y = -Inf, hjust = 1.15, vjust = -0.5,
                           label = paste0("\u2193 ", ref_level, "  (n=", n_down, ")"),
                           colour = .de_direction_colours[["down"]], size = 3.5, fontface = "bold") +
        unclipped_coord() +
        ggplot2::theme_bw(base_size = 12) +
        square_panel_theme() +
        ggplot2::theme(
            legend.position = "right",
            legend.justification = "center",
            plot.margin = ggplot2::margin(t = 12, r = 14, b = 5.5, l = 5.5)
        ) +
        ggplot2::xlab("log2 mean CPM") +
        ggplot2::ylab("log2 fold change") +
        ggplot2::labs(colour = NULL) +
        ggplot2::ggtitle(comparison)

    attr(p, "plotly_annotations") <- list(
        list(x = 0.98, y = 0.98, xref = "paper", yref = "paper",
             text = paste0("<b>\u2191 ", test_level, "  (n=", n_up, ")</b>"),
             showarrow = FALSE, xanchor = "right", yanchor = "top",
             font = list(color = .de_direction_colours[["up"]], size = 12)),
        list(x = 0.98, y = 0.02, xref = "paper", yref = "paper",
             text = paste0("<b>\u2193 ", ref_level, "  (n=", n_down, ")</b>"),
             showarrow = FALSE, xanchor = "right", yanchor = "bottom",
             font = list(color = .de_direction_colours[["down"]], size = 12))
    )
    p
}

#' Volcano plot (log2 fold-change vs. -log10 p-value), with top genes labelled
#'
#' Styled after the conventions popularised by the `EnhancedVolcano`
#' package (dashed significance/LFC threshold lines, points coloured by
#' regulation direction, repelled gene-symbol labels), implemented as a
#' plain `ggplot2` build so the up/down/not-significant colour legend
#' requested for this pipeline is preserved (rather than EnhancedVolcano's
#' own four-way NS/LFC/p/both legend). Points are coloured red
#' (significantly up), blue (significantly down), or grey (not
#' significant), using the `sig` column already computed by
#' [run_diffexp()]; significant points are drawn slightly larger. The
#' top `top_label_n` up- and top `top_label_n` down-regulated genes (by
#' `rank_score`) are labelled, and the number of up-/down-regulated genes
#' is annotated in the top corners of the plot -- the y-axis is expanded
#' with extra headroom above the highest point (via
#' `scale_y_continuous(expand = )`) so these corner labels don't overlap
#' the data, and the x-axis is set symmetric around 0
#' (`max(abs(xmin), abs(xmax))` on each side) so up- and down-regulation
#' are visually comparable. A `text` aesthetic (gene name, log2
#' fold-change, adjusted p-value) is included on every point so the plot
#' can be converted to an interactive, hoverable widget with
#' [make_interactive()].
#'
#' @param res_df A data frame as returned by [run_diffexp()].
#' @param alpha Significance threshold (only used if `res_df` lacks a `sig`
#'   column, and for the horizontal threshold line). Default `0.05`.
#' @param lfc_threshold Log2 fold-change threshold used for the vertical
#'   threshold lines and the significance overlay. Default `0`.
#' @param gene_col Name of the gene-identifier column. Default `"gene"`.
#' @param top_label_n Number of top **up** genes and top **down** genes
#'   (by `rank_score`) to label. Default `5`.
#' @return A `ggplot` object.
#' @export
plot_volcano <- function(res_df, alpha = 0.05, lfc_threshold = 0.5, gene_col = "gene", top_label_n = 5) {
    # Volcano uses the primary (shrunken) log2FoldChange column — same
    # values as the MA plot and the DEG table. Unshrunken MLE values stay
    # available as `log2FoldChange_unshrunken` for inspection only.
    res_df$direction <- classify_de_direction(res_df, alpha = alpha, lfc_threshold = lfc_threshold)
    # Floor padj away from 0 so -log10 stays finite; NA padj -> NA y (dropped).
    res_df$neglog10padj <- -log10(pmax(res_df$padj, .Machine$double.xmin))
    n_up <- sum(res_df$direction == "up")
    n_down <- sum(res_df$direction == "down")

    comparison <- unique(res_df$comparison)
    levels <- .parse_contrast_levels(comparison)
    test_level <- levels$test
    ref_level <- levels$ref

    res_df$.tooltip <- paste0(
        "Gene: ", res_df[[gene_col]],
        "<br>log2FC: ", formatC(res_df$log2FoldChange, format = "f", digits = 2),
        "<br>padj: ", formatC(res_df$padj, format = "e", digits = 2)
    )

    up_df <- res_df[res_df$direction == "up", ]
    up_df <- up_df[order(-up_df$rank_score), ][seq_len(min(top_label_n, nrow(up_df))), ]
    down_df <- res_df[res_df$direction == "down", ]
    down_df <- down_df[order(down_df$rank_score), ][seq_len(min(top_label_n, nrow(down_df))), ]
    label_df <- rbind(up_df, down_df)

    # Draw ns first so grey points sit behind red/blue (geom_point uses
    # data-row order). Legend level order stays up/down/ns via the factor.
    res_df <- res_df[order(match(as.character(res_df$direction), c("ns", "up", "down"))), ]

    # Symmetric x-axis: only trim super-extreme |log2FC| outliers; milder
    # extremes are kept so the cloud is not over-cropped.
    x_bound <- .robust_lfc_limit(res_df$log2FoldChange)

    p <- ggplot2::ggplot(res_df, ggplot2::aes(x = .data$log2FoldChange, y = .data$neglog10padj,
                                               colour = .data$direction, size = .data$direction,
                                               text = .data$.tooltip)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_colour_manual(values = .de_direction_colours, drop = FALSE,
                                      breaks = c("up", "down", "ns")) +
        ggplot2::scale_size_manual(values = .de_direction_sizes, guide = "none", drop = FALSE) +
        # A little headroom (10% of the data range on each side) rather
        # than zero data-padding, on both axes, so points near the
        # extremes -- and the count/direction-arrow annotations anchored
        # near the top/sides -- aren't crowded against the panel border.
        ggplot2::scale_x_continuous(limits = c(-x_bound, x_bound),
                                     expand = ggplot2::expansion(mult = c(0.1, 0.1)),
                                     oob = scales::oob_squish) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
        ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", colour = "grey40") +
        ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "grey40") +
        # Left/right arrows on top naming the groups, with DEG counts.
        # (Mirrored to plotly via attr below -- ggplotly drops annotate().)
        ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,
                           label = paste0(test_level, " \u2192  (up n=", n_up, ")"),
                           colour = .de_direction_colours[["up"]], size = 3.5, fontface = "bold") +
        ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
                           label = paste0("(down n=", n_down, ")  \u2190 ", ref_level),
                           colour = .de_direction_colours[["down"]], size = 3.5, fontface = "bold") +
        unclipped_coord() +
        ggplot2::theme_bw(base_size = 12) +
        square_panel_theme() +
        ggplot2::theme(
            legend.position = "right",
            legend.justification = "center",
            plot.margin = ggplot2::margin(t = 12, r = 5.5, b = 5.5, l = 5.5)
        ) +
        ggplot2::xlab("log2 fold change") +
        ggplot2::ylab("-log10(padj)") +
        ggplot2::labs(colour = NULL) +
        ggplot2::ggtitle(comparison)

    if (requireNamespace("ggrepel", quietly = TRUE) && nrow(label_df) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = label_df,
            ggplot2::aes(label = .data[[gene_col]]),
            colour = "black", size = 3, max.overlaps = 30
        )
    }
    attr(p, "plotly_annotations") <- list(
        list(x = 0.98, y = 0.98, xref = "paper", yref = "paper",
             text = paste0("<b>", test_level, " \u2192  (up n=", n_up, ")</b>"),
             showarrow = FALSE, xanchor = "right", yanchor = "top",
             font = list(color = .de_direction_colours[["up"]], size = 12)),
        list(x = 0.02, y = 0.98, xref = "paper", yref = "paper",
             text = paste0("<b>(down n=", n_down, ")  \u2190 ", ref_level, "</b>"),
             showarrow = FALSE, xanchor = "left", yanchor = "top",
             font = list(color = .de_direction_colours[["down"]], size = 12))
    )
    p
}

#' Top and bottom differentially expressed genes by ranking score
#'
#' Selects the `top_n` genes with the highest `rank_score` ("top", most
#' up-regulated in the test direction) and the `top_n` with the lowest
#' `rank_score` ("bottom", most down-regulated), for reporting alongside
#' the MA/volcano plots for a comparison.
#'
#' @param res_df A data frame as returned by [run_diffexp()] (must contain
#'   `rank_score`).
#' @param top_n Number of genes *per direction* to return. Default `15`.
#' @return A data frame: up to `2 * top_n` rows, ordered by descending
#'   `rank_score`, with all columns of `res_df`.
#' @export
diffexp_top_bottom_table <- function(res_df, top_n = 15) {
    df <- res_df[order(-res_df$rank_score), ]
    n <- nrow(df)
    top_n <- min(top_n, ceiling(n / 2))
    top_idx <- seq_len(min(top_n, n))
    bottom_idx <- seq(to = n, length.out = min(top_n, n))
    idx <- union(top_idx, bottom_idx)
    df[idx, ]
}

#' Interactive, sortable/searchable DEG table
#'
#' Wraps a differential expression result table (typically a top/bottom
#' subset for HTML reports, or the full table for interactive exploration)
#' in a `DT::datatable()`. Numeric columns are rounded for display only.
#'
#' @param res_df A data frame as returned by [run_diffexp()].
#' @param digits Number of decimal places to round numeric columns to for
#'   display. Default `4`.
#' @param page_length Number of rows shown per page. Default `10`.
#' @return A `DT::datatable` htmlwidget.
#' @export
make_deg_datatable <- function(res_df, digits = 4, page_length = 10) {
    if (!requireNamespace("DT", quietly = TRUE)) {
        stop("Package \"DT\" is required for interactive tables. Install it with ",
             "`install.packages(\"DT\")`.", call. = FALSE)
    }
    display_df <- res_df
    num_cols <- vapply(display_df, is.numeric, logical(1))
    display_df[num_cols] <- lapply(display_df[num_cols], round, digits = digits)

    DT::datatable(
        display_df,
        filter = "top",
        extensions = "Buttons",
        rownames = FALSE,
        options = list(
            pageLength = page_length,
            lengthMenu = c(10, 25, 50, 100, nrow(display_df)),
            scrollX = TRUE,
            dom = "Blfrtip",
            buttons = c("csv")
        )
    )
}

#' Summarise significant up/down DEG counts per contrast
#'
#' For each contrast in `diffexp_out$by_comparison`, counts significant
#' genes (`padj < alpha` and `|log2FoldChange| > lfc_threshold` when an
#' LFC column is present; LRT results use `padj < alpha` only) and reports
#' n / % up and down relative to all genes with non-`NA` `padj`.
#'
#' @param diffexp_out Output of [run_diffexp_pipeline()],
#'   [run_diffexp_lrt_pipeline()], or [prepare_diffexp_from_deg()].
#' @param alpha Significance threshold. Default `0.05`.
#' @param lfc_threshold Absolute LFC threshold for Wald / DEG-table
#'   contrasts. Ignored for LRT (no fold-change direction). Default `0.5`.
#' @return A data frame with one row per contrast.
#' @export
summarize_deg_contrasts <- function(diffexp_out, alpha = 0.05, lfc_threshold = 0.5) {
    by_comp <- diffexp_out$by_comparison
    if (is.null(by_comp) || !length(by_comp)) {
        return(data.frame(
            contrast = character(0),
            n_tested = integer(0),
            n_sig = integer(0),
            n_up = integer(0),
            pct_up = numeric(0),
            n_down = integer(0),
            pct_down = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    rows <- lapply(names(by_comp), function(comp) {
        res_df <- by_comp[[comp]]
        tested <- !is.na(res_df$padj)
        n_tested <- sum(tested)
        has_lfc <- "log2FoldChange" %in% colnames(res_df) &&
            any(is.finite(as.numeric(res_df$log2FoldChange)))
        if (has_lfc) {
            dir <- classify_de_direction(res_df, alpha = alpha, lfc_threshold = lfc_threshold)
            n_up <- sum(dir == "up", na.rm = TRUE)
            n_down <- sum(dir == "down", na.rm = TRUE)
            n_sig <- n_up + n_down
        } else {
            n_sig <- sum(tested & res_df$padj < alpha, na.rm = TRUE)
            n_up <- NA_integer_
            n_down <- NA_integer_
        }
        pct <- function(n) {
            if (is.na(n) || n_tested == 0L) NA_real_ else round(100 * n / n_tested, 2)
        }
        data.frame(
            contrast = comp,
            n_tested = as.integer(n_tested),
            n_sig = as.integer(n_sig),
            n_up = as.integer(n_up),
            pct_up = pct(n_up),
            n_down = as.integer(n_down),
            pct_down = pct(n_down),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}
