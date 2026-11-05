#' Build a design formula string from model terms plus an optional batch prefix
#'
#' Batch covariates are always prefixed (never tested), so every model in an
#' LRT pair carries the same batch terms and the likelihood-ratio test is
#' driven purely by the biological terms that differ.
#'
#' @param terms Character vector of biological model terms (may be empty).
#' @param batch Optional character vector of covariate column names.
#' @return A single character string, e.g. `"~ batch + genotype + time"`,
#'   or `"~ 1"` when nothing is left.
#' @noRd
.lrt_formula_string <- function(terms, batch = NULL) {
    rhs <- c(batch, terms)
    rhs <- rhs[!is.na(rhs) & nzchar(rhs)]
    if (!length(rhs)) return("~ 1")
    paste("~", paste(rhs, collapse = " + "))
}

#' Sanitize model terms into a filename/label-safe fragment
#'
#' `a * b` / `a:b` become `a_x_b`, `a + b` becomes `a_plus_b`, and an empty
#' term set becomes `"intercept"`.
#' @noRd
.lrt_sanitize_terms <- function(terms) {
    terms <- terms[!is.na(terms) & nzchar(terms)]
    if (!length(terms)) return("intercept")
    s <- paste(terms, collapse = " + ")
    s <- gsub("[[:space:]]*[*:][[:space:]]*", "_x_", s)
    s <- gsub("[[:space:]]*\\+[[:space:]]*", "_plus_", s)
    s <- gsub("[^A-Za-z0-9_]+", "_", s)
    s <- gsub("_+", "_", s)
    sub("_$", "", sub("^_", "", s))
}

#' Resolve abstract model components to concrete design terms
#'
#' Components are `"G"` (the `var` block, additive across all of `var`),
#' `"T"` (the time column), and `"GxT"` (their interaction, which already
#' implies both main effects and is therefore emitted as a single term).
#' @noRd
.lrt_component_terms <- function(comps, var, time) {
    if ("GxT" %in% comps) {
        return(if (length(var) == 1L) {
            paste(var, time, sep = " * ")
        } else {
            paste0("(", paste(var, collapse = " + "), ") * ", time)
        })
    }
    terms <- character(0)
    if ("G" %in% comps && length(var)) terms <- c(terms, paste(var, collapse = " + "))
    if ("T" %in% comps) terms <- c(terms, time)
    terms
}

#' Every strictly nested, hierarchy-respecting reduced model for one full model
#'
#' Returns component sets ordered from largest (closest to the full model) to
#' the empty set. Interaction components are never kept in a reduced model on
#' their own, since a model with `G:T` but without both main effects is not a
#' meaningful nested alternative.
#' @noRd
.lrt_reduced_components <- function(full_comps) {
    base <- intersect(c("G", "T"), full_comps)
    subs <- list()
    if (length(base)) {
        for (k in rev(seq_along(base))) {
            subs <- c(subs, utils::combn(base, k, simplify = FALSE))
        }
    }
    subs <- c(subs, list(character(0)))
    Filter(function(s) !setequal(s, full_comps), subs)
}

#' Human-readable description of the terms an LRT drops
#' @noRd
.lrt_tested_terms <- function(full_comps, reduced_comps, var, time) {
    dropped <- setdiff(full_comps, reduced_comps)
    labels <- vapply(dropped, function(cmp) {
        switch(cmp,
               G = paste(var, collapse = " + "),
               T = time,
               GxT = if (length(var) == 1L) {
                   paste0(var, ":", time)
               } else {
                   paste0("(", paste(var, collapse = " + "), "):", time)
               },
               cmp)
    }, character(1))
    paste(labels, collapse = ", ")
}

#' Expand the automatic full-vs-reduced model grid for likelihood-ratio tests
#'
#' Enumerates every candidate DESeq2 likelihood-ratio test for a time-course
#' design: each *full* model is paired with every strictly nested *reduced*
#' model, so the test isolates the biological terms that differ between them.
#' Batch covariates are never tested -- they are prefixed onto both models of
#' every pair (and therefore omitted from the label, where they would carry no
#' information).
#'
#' Full models are built from additive blocks: `G` (all of `var`, additive)
#' and `T` (`time`). With `var` non-empty the candidates are `~ B + G + T`,
#' `~ B + G` and `~ B + T`; with `var` empty (time-only designs) the single
#' candidate is `~ B + T`. When `include_interactions = TRUE` and `var` is
#' non-empty, `~ B + G * T` is added as well (`var * time` for a single
#' factor, `(var1 + var2) * time` for several).
#'
#' @param var Character vector of non-time factor column names in
#'   `colData(dds)`, or `NULL` / `character(0)` for a time-only design.
#' @param time Column name (character scalar) holding the time variable.
#' @param batch Optional character vector of covariate column names, always
#'   kept in both models of every pair. Default `NULL`.
#' @param include_interactions Logical; add the `var * time` full model.
#'   Default `FALSE` (ignored when `var` is empty).
#' @return A data frame with one row per test and columns `label`, `class`
#'   (always `"lrt"`), `type` (always `"lrt"`), `full` (full design formula
#'   as a string), `reduced` (reduced design formula as a string),
#'   `test_terms` (the terms the test drops) and `comparison` (a copy of
#'   `label`, so the table can be consumed like a contrast-spec table).
#' @examples
#' expand_lrt_specs(var = "genotype", time = "time", include_interactions = TRUE)
#' @export
expand_lrt_specs <- function(var, time, batch = NULL, include_interactions = FALSE) {
    if (!is.character(time) || length(time) != 1L || is.na(time) || !nzchar(time)) {
        stop("`time` must be a single non-empty character string naming a colData column.",
             call. = FALSE)
    }
    if (is.null(var)) var <- character(0)
    if (!is.character(var) || anyNA(var)) {
        stop("`var` must be `NULL`, `character(0)`, or a character vector of colData column names.",
             call. = FALSE)
    }
    var <- var[nzchar(var)]
    if (anyDuplicated(var)) {
        stop("`var` must not contain duplicated column names.", call. = FALSE)
    }
    if (time %in% var) {
        stop("`time` (\"", time, "\") must not also appear in `var`.", call. = FALSE)
    }
    if (!is.null(batch)) {
        if (!is.character(batch) || !length(batch) || anyNA(batch)) {
            stop("`batch` must be `NULL` or a non-empty character vector of column names.",
                 call. = FALSE)
        }
        if (length(intersect(batch, c(var, time)))) {
            stop("`batch` columns must not overlap with `var` or `time`.", call. = FALSE)
        }
    }
    validate_logical(include_interactions, "include_interactions")

    if (length(var)) {
        fulls <- list(c("G", "T"), "G", "T")
        if (isTRUE(include_interactions)) fulls <- c(list(c("G", "T", "GxT")), fulls)
    } else {
        if (isTRUE(include_interactions)) {
            message("`include_interactions = TRUE` ignored: `var` is empty, so there is ",
                    "nothing to interact with `", time, "`.")
        }
        fulls <- list("T")
    }

    rows <- list()
    for (full_comps in fulls) {
        full_terms <- .lrt_component_terms(full_comps, var, time)
        for (red_comps in .lrt_reduced_components(full_comps)) {
            red_terms <- .lrt_component_terms(red_comps, var, time)
            label <- paste0("lrt__full_", .lrt_sanitize_terms(full_terms),
                            "_vs_reduced_", .lrt_sanitize_terms(red_terms))
            rows[[length(rows) + 1L]] <- data.frame(
                label = label,
                class = "lrt",
                type = "lrt",
                full = .lrt_formula_string(full_terms, batch),
                reduced = .lrt_formula_string(red_terms, batch),
                test_terms = .lrt_tested_terms(full_comps, red_comps, var, time),
                comparison = label,
                stringsAsFactors = FALSE
            )
        }
    }
    specs <- do.call(rbind, rows)
    specs <- specs[!duplicated(specs$label), , drop = FALSE]
    rownames(specs) <- NULL
    new_contrast_specs(specs)
}

#' Validate that column names exist in `colData(dds)` and coerce them to factors
#' @noRd
.lrt_factorize_cols <- function(dds, cols, arg) {
    if (!length(cols)) return(dds)
    cd_names <- colnames(SummarizedExperiment::colData(dds))
    missing <- setdiff(cols, cd_names)
    if (length(missing)) {
        stop("`", arg, "` column(s) not in colData(dds): ", paste(missing, collapse = ", "),
             ". Available columns: ", paste(cd_names, collapse = ", "), call. = FALSE)
    }
    for (v in cols) {
        col <- SummarizedExperiment::colData(dds)[[v]]
        if (!is.factor(col)) col <- factor(col)
        if (nlevels(col) < 2L) {
            stop("`", arg, "` (\"", v, "\") must have at least 2 levels; found ",
                 nlevels(col), ".", call. = FALSE)
        }
        SummarizedExperiment::colData(dds)[[v]] <- col
    }
    dds
}

#' Sanitize a label into a directory-safe name
#'
#' Underscores are preserved so LRT labels (already `[A-Za-z0-9_]` only) map
#' onto the same folder names the LRT results were written to.
#' @noRd
.lrt_safe_dir <- function(x) gsub("[^A-Za-z0-9_]+", "_", as.character(x))

#' Run DESeq2 likelihood-ratio tests across the full/reduced model grid
#'
#' Fits every full-vs-reduced model pair from [expand_lrt_specs()] with
#' `DESeq2::DESeq(test = "LRT")` and returns results tables shaped exactly
#' like [run_diffexp_pipeline()]'s, so downstream consumers (TF activity,
#' PPI, CSPA, reporting) can treat them interchangeably. An LRT asks whether
#' *any* of the dropped terms explains variation, so it produces no single
#' fold change: `log2FoldChange`, `log2FoldChange_unshrunken` and `lfcSE` are
#' `NA`, and `rank_score` is the (always non-negative) LRT chi-square `stat`.
#'
#' @param dds A QC-filtered `DESeqDataSet`.
#' @param var Character vector of non-time factor columns, or `NULL` /
#'   `character(0)` for a time-only design.
#' @param time Column name holding the time variable.
#' @param batch Optional covariate column name(s), kept in every model and
#'   never tested. Default `NULL`.
#' @param include_interactions Logical; also test the `var * time` full
#'   model. Default `FALSE`.
#' @param alpha Adjusted p-value threshold used for `results()` independent
#'   filtering and the `sig` column. Default `0.05`.
#' @param gene_col Name of the gene-identifier column. Default `"gene"`.
#' @param save_dir Directory to write per-test subfolders to. Default
#'   `getwd()`.
#' @return A list with `results` (all tests row-bound), `by_comparison`
#'   (named list of per-test data frames), `contrast_specs` / `lrt_specs`
#'   (the spec table actually fitted) and `alpha`.
#' @examples
#' \dontrun{
#' lrt_out <- run_diffexp_lrt_pipeline(dds_qc, var = "genotype", time = "day",
#'                                     include_interactions = TRUE,
#'                                     save_dir = "results/lrt")
#' }
#' @export
run_diffexp_lrt_pipeline <- function(
    dds,
    var,
    time,
    batch = NULL,
    include_interactions = FALSE,
    alpha = 0.05,
    gene_col = "gene",
    save_dir = getwd()) {

    dds <- validate_dds(dds)
    alpha <- validate_alpha(alpha)
    save_dir <- validate_save_dir(save_dir)
    if (is.null(var)) var <- character(0)

    dds <- .lrt_factorize_cols(dds, time, "time")
    dds <- .lrt_factorize_cols(dds, var, "var")
    dds <- .lrt_factorize_cols(dds, batch, "batch")

    specs <- expand_lrt_specs(var = var, time = time, batch = batch,
                              include_interactions = include_interactions)

    raw_counts <- DESeq2::counts(dds)
    lib_sizes <- colSums(raw_counts)
    cpm_mat <- sweep(raw_counts, 2, lib_sizes, "/") * 1e6
    mean_cpm <- rowMeans(cpm_mat)

    by_comparison <- list()
    fitted <- rep(FALSE, nrow(specs))
    message("Running ", nrow(specs), " likelihood-ratio test(s) ...")
    for (i in seq_len(nrow(specs))) {
        label <- specs$label[[i]]
        full_fml <- stats::as.formula(specs$full[[i]], env = baseenv())
        reduced_fml <- stats::as.formula(specs$reduced[[i]], env = baseenv())
        message("  ", label, ": ", specs$full[[i]], " vs ", specs$reduced[[i]])

        res <- tryCatch({
            fit_dds <- dds
            DESeq2::design(fit_dds) <- full_fml
            fit_dds <- DESeq2::DESeq(fit_dds, test = "LRT", reduced = reduced_fml,
                                     quiet = TRUE)
            DESeq2::results(fit_dds, alpha = alpha)
        }, error = function(e) {
            warning("Skipping LRT \"", label, "\": ", conditionMessage(e),
                    call. = FALSE)
            NULL
        })
        if (is.null(res)) next

        res_df <- as.data.frame(res)
        genes <- rownames(res_df)
        out <- data.frame(
            gene = genes,
            comparison = label,
            baseMean = as.numeric(res_df$baseMean),
            mean_cpm = unname(mean_cpm[genes]),
            log2FoldChange = NA_real_,
            log2FoldChange_unshrunken = NA_real_,
            lfcSE = NA_real_,
            stat = as.numeric(res_df$stat),
            pvalue = as.numeric(res_df$pvalue),
            padj = as.numeric(res_df$padj),
            stringsAsFactors = FALSE
        )
        colnames(out)[1] <- gene_col
        out$rank_score <- out$stat
        out$sig <- !is.na(out$padj) & out$padj < alpha
        rownames(out) <- NULL

        write_output(out, file.path(save_dir, label, paste0(contrast_file_stem(label), "_lrt.csv")))
        message("    ", sum(out$sig), " significant gene(s) at padj < ", alpha, ".")
        by_comparison[[label]] <- out
        fitted[i] <- TRUE
    }

    specs <- specs[fitted, , drop = FALSE]
    rownames(specs) <- NULL
    combined <- if (length(by_comparison)) do.call(rbind, by_comparison) else NULL
    if (!is.null(combined)) rownames(combined) <- NULL

    list(
        results = combined,
        by_comparison = by_comparison,
        contrast_specs = specs,
        lrt_specs = specs,
        alpha = alpha
    )
}

#' Cluster LRT-significant genes into temporal expression patterns
#'
#' Runs `DEGreport::degPatterns()` on the variance-stabilized profiles of the
#' genes each likelihood-ratio test called significant, grouping them into
#' clusters of genes that change over `time` in the same way. Tests with
#' fewer than `minc` significant genes are skipped (a cluster cannot be
#' smaller than `minc`, so there is nothing to find).
#'
#' Requires Imports package `DEGreport`.
#'
#' @param dds A QC-filtered `DESeqDataSet` (the same object passed to
#'   [run_diffexp_lrt_pipeline()]).
#' @param lrt_out Output of [run_diffexp_lrt_pipeline()].
#' @param time Column name holding the time variable.
#' @param var Optional non-time factor column(s); the first is passed to
#'   `degPatterns(col = )` so trajectories are split by group. Default
#'   `NULL`.
#' @param alpha Adjusted p-value threshold, used when a results table has no
#'   `sig` column. Default `0.05`.
#' @param minc Minimum cluster size, passed to `degPatterns()`. Also the
#'   minimum number of significant genes required to attempt clustering.
#'   Default `15`.
#' @param gene_col Name of the gene-identifier column. Default `"gene"`.
#' @param save_dir Directory to write per-test subfolders to. Default
#'   `getwd()`.
#' @param save_plot Logical; write the cluster-trajectory PDF. Default
#'   `TRUE`.
#' @return A named list (one element per LRT label) of lists with `df` (the
#'   gene/cluster assignment table, `NULL` when skipped), `plot` (the
#'   `degPatterns` trajectory `ggplot`, or `NULL`), `pass` (logical: did
#'   clustering run and return clusters), `n_sig` (significant genes) and
#'   `n_clusters`.
#' @examples
#' \dontrun{
#' dp <- run_degpatterns_pipeline(dds_qc, lrt_out, time = "day", var = "genotype")
#' }
#' @export
run_degpatterns_pipeline <- function(
    dds,
    lrt_out,
    time,
    var = NULL,
    alpha = 0.05,
    minc = 15,
    gene_col = "gene",
    save_dir = getwd(),
    save_plot = TRUE) {

    if (!requireNamespace("DEGreport", quietly = TRUE)) {
        stop("Package \"DEGreport\" is required for run_degpatterns_pipeline() ",
             "but is not installed. Install with BiocManager::install(\"DEGreport\").",
             call. = FALSE)
    }
    dds <- validate_dds(dds)
    alpha <- validate_alpha(alpha)
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)
    minc <- as.integer(minc)[1]
    if (is.na(minc) || minc < 2L) {
        stop("`minc` must be an integer >= 2.", call. = FALSE)
    }
    if (is.null(var)) var <- character(0)
    dds <- .lrt_factorize_cols(dds, time, "time")
    dds <- .lrt_factorize_cols(dds, var, "var")

    by_comp <- lrt_out$by_comparison
    if (is.null(by_comp) || !length(by_comp)) {
        stop("`lrt_out$by_comparison` is empty; nothing to cluster.", call. = FALSE)
    }

    message("Computing blind VST for temporal clustering ...")
    vsd <- quiet(safe_vst(dds, blind = TRUE))
    mat <- as.matrix(SummarizedExperiment::assay(vsd))
    storage.mode(mat) <- "double"

    metadata <- as.data.frame(SummarizedExperiment::colData(dds))
    rownames(metadata) <- colnames(dds)
    col_var <- if (length(var)) var[[1L]] else NULL

    out <- list()
    for (comp in names(by_comp)) {
        res <- by_comp[[comp]]
        sig <- if ("sig" %in% colnames(res)) {
            !is.na(res$sig) & as.logical(res$sig)
        } else {
            !is.na(res$padj) & res$padj < alpha
        }
        genes <- unique(as.character(res[[gene_col]][sig]))
        genes <- genes[genes %in% rownames(mat)]
        n_sig <- length(genes)

        if (n_sig < minc) {
            message("  Skipping \"", comp, "\": ", n_sig,
                    " significant gene(s) < minc = ", minc, ".")
            out[[comp]] <- list(df = NULL, plot = NULL, pass = FALSE,
                                n_sig = n_sig, n_clusters = 0L)
            next
        }

        message("  ", comp, ": clustering ", n_sig, " significant gene(s) ...")
        ma <- mat[genes, , drop = FALSE]
        dp <- tryCatch(
            quiet(DEGreport::degPatterns(
                ma, metadata, time = time, col = col_var,
                minc = minc, plot = FALSE
            )),
            error = function(e) {
                warning("degPatterns() failed for \"", comp, "\": ",
                        conditionMessage(e), call. = FALSE)
                NULL
            }
        )
        if (is.null(dp) || is.null(dp$df) || !nrow(dp$df)) {
            message("    No clusters returned for \"", comp, "\".")
            out[[comp]] <- list(df = NULL, plot = NULL, pass = FALSE,
                                n_sig = n_sig, n_clusters = 0L)
            next
        }

        df <- as.data.frame(dp$df)
        if (!gene_col %in% colnames(df)) {
            src <- intersect(c("genes", "gene", "rownames"), colnames(df))
            df[[gene_col]] <- if (length(src)) as.character(df[[src[[1L]]]]) else rownames(df)
        }
        df[[gene_col]] <- as.character(df[[gene_col]])
        df$cluster <- as.character(df$cluster)
        rownames(df) <- NULL

        comp_dir <- file.path(save_dir, .lrt_safe_dir(comp))
        write_output(df, file.path(comp_dir, "degpatterns_clusters.csv"))

        # `plot = FALSE` above keeps degPatterns() from drawing to the
        # current device; depending on the DEGreport version it may then
        # return no plot object at all, so rebuild one from its normalized
        # long table when that happens.
        p <- dp$plot
        if (!inherits(p, "ggplot") && !is.null(dp$normalized) &&
            "degPlotCluster" %in% getNamespaceExports("DEGreport")) {
            p <- tryCatch(
                quiet(DEGreport::degPlotCluster(dp$normalized, time = time,
                                                color = col_var)),
                error = function(e) NULL
            )
        }
        if (save_plot && inherits(p, "ggplot")) {
            save_ggplot(p, file.path(comp_dir, "degpatterns_clusters.pdf"),
                        width = 8, height = 6)
        }
        n_clusters <- length(unique(df$cluster))
        message("    ", n_clusters, " cluster(s) covering ", nrow(df), " gene(s).")
        out[[comp]] <- list(df = df, plot = p, pass = TRUE,
                            n_sig = n_sig, n_clusters = n_clusters)
    }
    out
}

#' Over-representation analysis for every temporal expression cluster
#'
#' Runs `fgsea::fora()` for each `degPatterns` cluster against each requested
#' MSigDB collection, using the same universe/plot conventions as
#' [run_wgcna_pipeline()]'s per-module ORA.
#'
#' @param degpatterns_out Output of [run_degpatterns_pipeline()].
#' @param universe Character vector of all tested genes (typically
#'   `rownames(dds)`), used as the ORA background.
#' @param collections MSigDB collection alias(es), see
#'   `validate_msigdb_collections()`. Default matches the GSEA/WGCNA default.
#' @param top_n Number of top pathways (by adjusted p-value) shown per dot
#'   plot. Default `10`.
#' @param alpha Adjusted p-value threshold for flagging a pathway
#'   significant (and for deciding whether a dot plot is worth drawing).
#'   Default `0.05`.
#' @param gene_col Name of the gene-identifier column in the cluster tables.
#'   Default `"gene"`.
#' @param save_dir Directory to write `<lrt>/<cluster>/ora_<collection>.csv`
#'   under. Default `getwd()`.
#' @param save_plot Logical; write one dot-plot PDF per cluster/collection.
#'   Default `TRUE`.
#' @return A list with `org`, `collections`, `ora_by_lrt`
#'   (`[[lrt]][[cluster]][[collection]]` data frames) and `plots_by_lrt`
#'   (same keying; `NULL` where nothing was significant).
#' @examples
#' \dontrun{
#' ora <- run_cluster_ora_pipeline(dp, universe = rownames(dds_qc),
#'                                 collections = c("HALLMARK", "GOBP"))
#' }
#' @export
run_cluster_ora_pipeline <- function(
    degpatterns_out,
    universe,
    collections = c("GOBP", "GOCC", "GOMF", "KEGG",
                    "REACTOME", "PID"),
    top_n = 10,
    alpha = 0.05,
    gene_col = "gene",
    save_dir = getwd(),
    save_plot = TRUE) {

    collections <- validate_msigdb_collections(collections)
    alpha <- validate_alpha(alpha)
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)
    universe <- unique(as.character(universe))
    universe <- universe[!is.na(universe) & nzchar(universe)]
    if (length(universe) < 2L) {
        stop("`universe` must contain at least 2 gene identifiers.", call. = FALSE)
    }

    org <- detect_organism(universe)
    message("Fetching ", length(collections), " MSigDB collection(s) for organism=\"",
            org, "\" ...")
    pathway_sets <- stats::setNames(
        lapply(collections, function(coll) import_msigdbr(org, coll)),
        collections
    )
    can_plot <- exists("plot_module_ora_dotplot", mode = "function")

    ora_by_lrt <- list()
    plots_by_lrt <- list()
    for (lrt in names(degpatterns_out)) {
        entry <- degpatterns_out[[lrt]]
        df <- entry$df
        if (is.null(df) || !nrow(df)) next
        if (!gene_col %in% colnames(df)) {
            warning("Skipping \"", lrt, "\": no `", gene_col, "` column in the cluster table.",
                    call. = FALSE)
            next
        }
        clusters <- unique(as.character(df$cluster))
        message("Running cluster ORA for \"", lrt, "\" (", length(clusters),
                " cluster(s)) ...")

        ora_by_cluster <- list()
        plots_by_cluster <- list()
        for (cl in clusters) {
            cl_genes <- unique(as.character(df[[gene_col]][as.character(df$cluster) == cl]))
            cl_genes <- cl_genes[!is.na(cl_genes) & nzchar(cl_genes)]
            cl_dir <- file.path(save_dir, .lrt_safe_dir(lrt),
                                paste0("cluster_", .lrt_safe_dir(cl)))
            ora_by_coll <- list()
            plots_by_coll <- list()
            for (coll in collections) {
                res <- as.data.frame(fgsea::fora(
                    pathway_sets[[coll]], genes = cl_genes, universe = universe
                ))
                if (nrow(res)) {
                    res <- res[order(res$padj), , drop = FALSE]
                    res$overlapGenes <- vapply(res$overlapGenes, paste, character(1),
                                               collapse = ",")
                    res$sig <- !is.na(res$padj) & res$padj < alpha
                    rownames(res) <- NULL
                }
                write_output(res, file.path(cl_dir, paste0("ora_", coll, ".csv")))
                ora_by_coll[[coll]] <- res

                p <- NULL
                sig_res <- if (nrow(res)) res[!is.na(res$padj) & res$padj < alpha, , drop = FALSE] else res
                if (can_plot && nrow(sig_res)) {
                    p <- plot_module_ora_dotplot(
                        sig_res, top_n = top_n,
                        title = paste0("Cluster ", cl, ": ", lrt)
                    )
                    if (save_plot && !is.null(p)) {
                        save_ggplot(p, file.path(cl_dir, paste0("ora_dotplot_", coll, ".pdf")),
                                    width = 9 * 0.65, height = 4.8 * 0.65)
                    }
                }
                plots_by_coll[[coll]] <- p
            }
            ora_by_cluster[[as.character(cl)]] <- ora_by_coll
            plots_by_cluster[[as.character(cl)]] <- plots_by_coll
        }
        ora_by_lrt[[lrt]] <- ora_by_cluster
        plots_by_lrt[[lrt]] <- plots_by_cluster
    }

    list(
        org = org,
        collections = collections,
        ora_by_lrt = ora_by_lrt,
        plots_by_lrt = plots_by_lrt
    )
}

#' Reshape LRT results into a Wald-style `diffexp_out` for TF-activity analysis
#'
#' [run_decoupler_tf_pipeline()] and [run_viper_pipeline()] consume
#' `diffexp_out$by_comparison` and rank genes by the `stat` column. LRT
#' tables carry a chi-square `stat` for every gene, including the vast
#' majority that show no temporal signal, so this helper subsets each table
#' to its significant genes (keeping their LRT `stat`) before handing it
#' over.
#'
#' @param lrt_out Output of [run_diffexp_lrt_pipeline()].
#' @param alpha Adjusted p-value threshold, used when a table has no `sig`
#'   column. Default `0.05`.
#' @param gene_col Name of the gene-identifier column. Default `"gene"`.
#' @return A list shaped like [run_diffexp_pipeline()]'s output: `results`,
#'   `by_comparison` and `contrast_specs`.
#' @examples
#' \dontrun{
#' tf_in <- lrt_sig_diffexp_for_tf(lrt_out)
#' run_decoupler_tf_pipeline(dds_qc, tf_in, save_dir = "results/lrt/tf")
#' }
#' @export
lrt_sig_diffexp_for_tf <- function(lrt_out, alpha = 0.05, gene_col = "gene") {
    alpha <- validate_alpha(alpha)
    by_comp <- lrt_out$by_comparison
    if (is.null(by_comp) || !length(by_comp)) {
        stop("`lrt_out$by_comparison` is empty; nothing to hand to TF-activity analysis.",
             call. = FALSE)
    }

    by_comparison <- list()
    for (comp in names(by_comp)) {
        res <- by_comp[[comp]]
        sig <- if ("sig" %in% colnames(res)) {
            !is.na(res$sig) & as.logical(res$sig)
        } else {
            !is.na(res$padj) & res$padj < alpha
        }
        sub <- res[sig & !is.na(res$stat), , drop = FALSE]
        if (!nrow(sub)) {
            message("  Dropping \"", comp, "\": no significant genes with an LRT statistic.")
            next
        }
        sub <- sub[order(-sub$stat), , drop = FALSE]
        rownames(sub) <- NULL
        by_comparison[[comp]] <- sub
    }
    if (!length(by_comparison)) {
        stop("No LRT comparison had significant genes; cannot build a TF-activity input.",
             call. = FALSE)
    }

    combined <- do.call(rbind, by_comparison)
    rownames(combined) <- NULL
    specs <- lrt_out$contrast_specs
    if (!is.null(specs) && "label" %in% colnames(specs)) {
        specs <- specs[specs$label %in% names(by_comparison), , drop = FALSE]
        rownames(specs) <- NULL
    }
    list(results = combined, by_comparison = by_comparison,
         contrast_specs = specs, gene_col = gene_col)
}

#' Reshape degPatterns clusters into a Wald-style `diffexp_out` for TF analysis
#'
#' Builds one "contrast" per temporal cluster (`"<lrt> cluster <n>"`) whose
#' result table holds the LRT `stat` for genes in that cluster. Pass the
#' result to [run_decoupler_tf_pipeline()] / [run_viper_pipeline()] so ULM /
#' msVIPER score TFs associated with each cluster.
#'
#' @param degpatterns_out Output of [run_degpatterns_pipeline()].
#' @param lrt_out Output of [run_diffexp_lrt_pipeline()] (supplies LRT
#'   statistics). If `NULL`, cluster tables must already contain a `stat`
#'   column.
#' @param gene_col Gene-id column. Default `"gene"`.
#' @param min_genes Minimum genes per cluster to keep. Default `5`.
#' @return A list shaped like [run_diffexp_pipeline()] output, with
#'   `source = "lrt_clusters"`.
#' @export
lrt_clusters_diffexp_for_tf <- function(degpatterns_out,
                                        lrt_out = NULL,
                                        gene_col = "gene",
                                        min_genes = 5L) {
    min_genes <- max(1L, as.integer(min_genes)[1])
    by_comparison <- list()
    for (lrt in names(degpatterns_out)) {
        dp <- degpatterns_out[[lrt]]
        df <- dp$df
        if (is.null(df) || !nrow(df) || !gene_col %in% colnames(df)) next
        if (!"cluster" %in% colnames(df)) next

        stats_map <- NULL
        if (!is.null(lrt_out) && !is.null(lrt_out$by_comparison[[lrt]])) {
            res <- lrt_out$by_comparison[[lrt]]
            gcol <- if (gene_col %in% colnames(res)) gene_col else {
                if ("gene" %in% colnames(res)) "gene" else colnames(res)[[1]]
            }
            if ("stat" %in% colnames(res)) {
                stats_map <- stats::setNames(
                    as.numeric(res$stat), as.character(res[[gcol]])
                )
            }
        }

        for (cl in sort(unique(as.character(df$cluster)))) {
            sub <- df[as.character(df$cluster) == cl, , drop = FALSE]
            genes <- as.character(sub[[gene_col]])
            ok <- !is.na(genes) & nzchar(genes)
            genes <- genes[ok]
            if (length(genes) < min_genes) next

            if (!is.null(stats_map)) {
                st <- unname(stats_map[genes])
            } else if ("stat" %in% colnames(sub)) {
                st <- as.numeric(sub$stat[ok])
            } else {
                st <- rep(1, length(genes))
            }
            keep <- !is.na(st)
            if (sum(keep) < min_genes) next
            label <- paste0(lrt, " cluster ", cl)
            tab <- data.frame(
                gene = genes[keep],
                stat = st[keep],
                padj = 0,
                sig = TRUE,
                comparison = label,
                stringsAsFactors = FALSE
            )
            tab <- tab[order(-tab$stat), , drop = FALSE]
            rownames(tab) <- NULL
            by_comparison[[label]] <- tab
        }
    }
    if (!length(by_comparison)) {
        stop("No degPatterns clusters had enough genes for TF-activity input.",
             call. = FALSE)
    }
    combined <- do.call(rbind, by_comparison)
    rownames(combined) <- NULL
    specs <- data.frame(
        label = names(by_comparison),
        class = "lrt_cluster",
        type = "lrt_cluster",
        stratum_var = NA_character_,
        stratum_level = NA_character_,
        stringsAsFactors = FALSE
    )
    list(
        results = combined,
        by_comparison = by_comparison,
        contrast_specs = specs,
        gene_col = gene_col,
        source = "lrt_clusters"
    )
}

#' Present temporal clusters as a WGCNA-like module object
#'
#' [run_ppi_pipeline()] reads gene/module assignments from
#' `wgcna_out$module_df` (columns `gene` and `module`, optionally `kME`).
#' This helper flattens [run_degpatterns_pipeline()]'s per-test cluster
#' tables into exactly that shape, so temporal clusters can be fed straight
#' into the STRING PPI pipeline as if they were co-expression modules.
#' Cluster "modules" are named `"<lrt label> cluster <n>"`; there is no
#' membership statistic, so no `kME` column is emitted (PPI seeds then fall
#' back to unranked cluster membership).
#'
#' @param degpatterns_out Output of [run_degpatterns_pipeline()].
#' @param gene_col Name of the gene-identifier column. Default `"gene"`.
#' @return A list with `module_df` (all clusters row-bound, columns `gene`,
#'   `module`, `cluster`, `lrt`), `by_lrt` (the same table split by LRT
#'   label) and `source = "degPatterns"`.
#' @examples
#' \dontrun{
#' clusters <- lrt_clusters_as_wgcna_like(dp)
#' run_ppi_pipeline(diffexp_out, wgcna_out = clusters, save_dir = "results/lrt/ppi")
#' }
#' @export
lrt_clusters_as_wgcna_like <- function(degpatterns_out, gene_col = "gene") {
    by_lrt <- list()
    for (lrt in names(degpatterns_out)) {
        df <- degpatterns_out[[lrt]]$df
        if (is.null(df) || !nrow(df) || !gene_col %in% colnames(df)) next
        mod <- data.frame(
            gene = as.character(df[[gene_col]]),
            module = paste0(lrt, " cluster ", as.character(df$cluster)),
            cluster = as.character(df$cluster),
            lrt = lrt,
            stringsAsFactors = FALSE
        )
        colnames(mod)[1] <- gene_col
        rownames(mod) <- NULL
        by_lrt[[lrt]] <- mod
    }
    module_df <- if (length(by_lrt)) do.call(rbind, by_lrt) else data.frame(
        gene = character(0), module = character(0),
        cluster = character(0), lrt = character(0),
        stringsAsFactors = FALSE
    )
    if (!length(by_lrt)) colnames(module_df)[1] <- gene_col
    rownames(module_df) <- NULL
    list(module_df = module_df, by_lrt = by_lrt, source = "degPatterns")
}

#' Row-scaled VST heatmap of CSPA surface genes within temporal clusters
#'
#' For every LRT with temporal clusters, intersects the clustered genes with
#' the Cell Surface Protein Atlas symbol set (see [load_cspa_genes()]) and
#' draws one heatmap per test. Up to `n_per_cluster` CSPA genes are kept
#' **per degPatterns cluster** (ranked by LRT `stat` when `lrt_out` is
#' supplied, otherwise by VST row variance). Genes are ordered by cluster
#' (row annotation) and scaled by row across all samples.
#'
#' @param dds A QC-filtered `DESeqDataSet`.
#' @param degpatterns_out Output of [run_degpatterns_pipeline()].
#' @param lrt_out Optional output of [run_diffexp_lrt_pipeline()]; used to
#'   rank CSPA genes within each cluster by LRT `stat` / `rank_score`.
#'   Default `NULL` (rank by VST variance).
#' @param cspa_genes Character vector of CSPA symbols. Default `NULL` loads
#'   them via [load_cspa_genes()] for the organism detected from
#'   `rownames(dds)`.
#' @param n_per_cluster Maximum CSPA genes shown per temporal cluster.
#'   Default `10`.
#' @param annotation_var Optional `colData(dds)` column name(s) used as
#'   column annotation (e.g. the time variable). Default `NULL`.
#' @param batch Optional batch column(s), including the `pair` column of a
#'   paired design. When set, heatmaps use
#'   batch-corrected expression (protecting `annotation_var`).
#' @param gene_col Name of the gene-identifier column in the cluster tables.
#'   Default `"gene"`.
#' @param save_dir Directory for per-LRT outputs. Files go in
#'   `<save_dir>/<lrt>/cspa/` (pipeline uses `deseq2/`). Default `getwd()`.
#' @param save_plot Logical; write the heatmap PDF. Default `TRUE`.
#' @return A named list (one element per LRT with at least one CSPA gene) of
#'   lists with `heatmap` (a `pheatmap` object or `NULL`), `genes` (the
#'   gene/cluster table plotted) and `n` (gene count).
#' @examples
#' \dontrun{
#' run_lrt_cspa_heatmap(dds_qc, dp, lrt_out = lrt_out, n_per_cluster = 10)
#' }
#' @export
run_lrt_cspa_heatmap <- function(
    dds,
    degpatterns_out,
    lrt_out = NULL,
    cspa_genes = NULL,
    n_per_cluster = 10,
    annotation_var = NULL,
    gene_col = "gene",
    batch = NULL,
    save_dir = getwd(),
    save_plot = TRUE) {

    dds <- validate_dds(dds)
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)
    n_per_cluster <- as.integer(n_per_cluster)[1]
    if (is.na(n_per_cluster) || n_per_cluster < 1L) {
        stop("`n_per_cluster` must be a positive integer.", call. = FALSE)
    }

    if (is.null(cspa_genes)) {
        org <- detect_organism(rownames(dds))
        message("Loading CSPA surface markers for organism=\"", org, "\" ...")
        cspa_genes <- load_cspa_genes(organism = org)
    }
    cspa_genes <- unique(as.character(cspa_genes))
    cspa_genes <- cspa_genes[!is.na(cspa_genes) & nzchar(cspa_genes)]
    if (!length(cspa_genes)) {
        stop("`cspa_genes` is empty; nothing to intersect with the temporal clusters.",
             call. = FALSE)
    }

    anno_col <- NULL
    if (!is.null(annotation_var)) {
        cd <- as.data.frame(SummarizedExperiment::colData(dds))
        missing <- setdiff(annotation_var, colnames(cd))
        if (length(missing)) {
            stop("`annotation_var` column(s) not in colData(dds): ",
                 paste(missing, collapse = ", "), call. = FALSE)
        }
        anno_col <- cd[, annotation_var, drop = FALSE]
        rownames(anno_col) <- colnames(dds)
    }

    if (!is.null(batch)) dds <- validate_batch(dds, batch)

    message("Computing sample-level expression for CSPA cluster heatmaps ...")
    protect <- if (!is.null(annotation_var) && length(annotation_var)) {
        annotation_var
    } else if (!is.null(batch)) {
        stop("`annotation_var` must be supplied when `batch` is set.", call. = FALSE)
    } else {
        colnames(SummarizedExperiment::colData(dds))[[1L]]
    }
    mat <- sample_level_expr_matrix(dds, var = protect, batch = batch)
    expr_correction <- attr(mat, "expr_correction")
    storage.mode(mat) <- "double"

    palette_length <- 100L
    colours <- grDevices::colorRampPalette(c("Darkblue", "white", "red"))(palette_length)
    breaks <- c(
        seq(-3, 0, length.out = ceiling(palette_length / 2) + 1),
        seq(0.05, 3, length.out = floor(palette_length / 2))
    )

    out <- list()
    for (lrt in names(degpatterns_out)) {
        df <- degpatterns_out[[lrt]]$df
        if (is.null(df) || !nrow(df) || !gene_col %in% colnames(df)) next
        hits <- .match_cspa_genes(df[[gene_col]], cspa_genes)
        sub <- df[as.character(df[[gene_col]]) %in% hits, , drop = FALSE]
        sub <- sub[as.character(sub[[gene_col]]) %in% rownames(mat), , drop = FALSE]
        sub <- sub[!duplicated(as.character(sub[[gene_col]])), , drop = FALSE]
        if (!nrow(sub)) {
            message("  ", lrt, ": no CSPA surface genes among the clustered genes.")
            next
        }

        # Rank score: LRT statistic when available, else VST row variance.
        score <- rep(0, nrow(sub))
        lrt_tbl <- if (!is.null(lrt_out) && !is.null(lrt_out$by_comparison[[lrt]])) {
            lrt_out$by_comparison[[lrt]]
        } else {
            NULL
        }
        if (!is.null(lrt_tbl) && gene_col %in% colnames(lrt_tbl)) {
            sc_col <- if ("rank_score" %in% colnames(lrt_tbl)) {
                "rank_score"
            } else if ("stat" %in% colnames(lrt_tbl)) {
                "stat"
            } else {
                NA_character_
            }
            if (!is.na(sc_col)) {
                idx <- match(as.character(sub[[gene_col]]),
                             as.character(lrt_tbl[[gene_col]]))
                score <- abs(as.numeric(lrt_tbl[[sc_col]][idx]))
                score[!is.finite(score)] <- 0
            }
        }
        if (!any(score > 0)) {
            score <- apply(mat[as.character(sub[[gene_col]]), , drop = FALSE], 1, stats::var,
                           na.rm = TRUE)
            score[!is.finite(score)] <- 0
        }
        sub$cspa_rank_score <- as.numeric(score)

        # Top n_per_cluster CSPA genes within each temporal cluster.
        parts <- split(sub, as.character(sub$cluster))
        parts <- lapply(parts, function(cl_df) {
            cl_df <- cl_df[order(-cl_df$cspa_rank_score, cl_df[[gene_col]]), ,
                           drop = FALSE]
            utils::head(cl_df, n_per_cluster)
        })
        sub <- do.call(rbind, parts)
        rownames(sub) <- NULL
        sub <- sub[order(as.character(sub$cluster), -sub$cspa_rank_score,
                         as.character(sub[[gene_col]])), , drop = FALSE]
        genes <- as.character(sub[[gene_col]])

        lrt_dir <- file.path(save_dir, lrt, "cspa")
        write_output(sub, file.path(lrt_dir, "cspa_cluster_genes.csv"))

        hm <- NULL
        if (length(genes) >= 2L && ncol(mat) >= 2L) {
            scaled <- t(scale(t(mat[genes, , drop = FALSE])))
            scaled[!is.finite(scaled)] <- 0
            anno_row <- data.frame(cluster = factor(sub$cluster),
                                   stringsAsFactors = FALSE)
            rownames(anno_row) <- genes
            hm <- pheatmap::pheatmap(
                scaled,
                color = colours,
                breaks = breaks,
                border_color = NA,
                cluster_rows = FALSE,
                cluster_cols = TRUE,
                annotation_row = anno_row,
                annotation_col = anno_col,
                main = paste0("Top ", n_per_cluster,
                              " CSPA genes / cluster: ", lrt),
                fontsize_row = 8,
                silent = TRUE
            )
            if (save_plot) {
                grDevices::pdf(
                    file.path(lrt_dir, "cspa_cluster_heatmap.pdf"),
                    width = max(6, 0.25 * ncol(mat) + 2) * 0.65,
                    height = max(4, 0.22 * length(genes) + 2) * 0.65
                )
                grid::grid.newpage()
                grid::grid.draw(hm$gtable)
                grDevices::dev.off()
            }
        }
        message("  ", lrt, ": ", length(genes), " CSPA surface gene(s) ",
                "(top ", n_per_cluster, " per cluster; ",
                length(unique(sub$cluster)), " cluster(s)).")
        out[[lrt]] <- list(heatmap = hm, genes = sub, n = length(genes),
                           n_per_cluster = n_per_cluster)
    }
    attr(out, "expr_correction") <- expr_correction
    out
}
