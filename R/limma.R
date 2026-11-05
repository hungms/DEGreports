# limma-trend engine for TPM / FPKM (non-count) expression.
#
# Every module reads a `DESeqDataSet` and asks it for a VST matrix. Instead of
# teaching each module a second object type, a limma run wraps the data in a
# *carrier* `DESeqDataSet`:
#
# * `"counts"` holds the values rounded to integers, only because a
#   `DESeqDataSet` accepts nothing else -- nothing in a limma run reads it;
# * `"expr_values"` holds the values as supplied; QC, the count export and
#   limma read this ([.raw_count_matrix()]);
# * where a module would take a VST / rlog it gets log2(x + prior_count)
#   ([safe_vst()], [safe_rlog()]), and batch correction is ComBat on that
#   matrix, since ComBat-seq needs integer counts.
#
# DE is limma-trend, returned in the DESeq2 column schema so GSEA, decoupleR,
# VIPER, PPI, CSPA, multihit and the report run unchanged.

#' Assay holding the supplied values on a carrier
#' @noRd
.expr_assay <- "expr_values"

#' Is this a carrier built by [expr_to_dds()]?
#' @noRd
.is_expr_dds <- function(dds) {
    methods::is(dds, "SummarizedExperiment") &&
        .expr_assay %in% SummarizedExperiment::assayNames(dds)
}

#' The values to treat as "raw counts": supplied values on a carrier, else
#' the DESeq2 counts
#' @noRd
.raw_count_matrix <- function(dds) {
    if (.is_expr_dds(dds)) {
        as.matrix(SummarizedExperiment::assay(dds, .expr_assay))
    } else {
        DESeq2::counts(dds, normalized = FALSE)
    }
}

#' Mean CPM per gene across every sample (MA-plot x-axis)
#' @noRd
.mean_cpm <- function(dds) {
    raw <- .raw_count_matrix(dds)
    rowMeans(sweep(raw, 2, colSums(raw), "/") * 1e6)
}

#' log2(x + prior_count) of a carrier: limma's input and the VST stand-in
#' @noRd
.log_expr <- function(dds) {
    prior <- S4Vectors::metadata(dds)$degreports_prior_count
    log2(.raw_count_matrix(dds) + if (is.null(prior)) 1 else prior)
}

#' [.log_expr()] wrapped like a `DESeqTransform`, so `assay()` callers work
#' @noRd
.log_expr_transform <- function(dds) {
    SummarizedExperiment::SummarizedExperiment(
        assays = list(log2expr = .log_expr(dds)),
        colData = SummarizedExperiment::colData(dds)
    )
}

#' rlog, or log2 expression for a carrier
#' @noRd
safe_rlog <- function(dds, blind = TRUE) {
    if (.is_expr_dds(dds)) return(.log_expr_transform(dds))
    DESeq2::rlog(dds, blind = blind)
}

#' Wrap TPM / FPKM data as a carrier `DESeqDataSet`
#'
#' @param object A `SummarizedExperiment` (values in `assay`, sample metadata
#'   in `colData`) or an `edgeR::DGEList` (values in `$counts`, metadata in
#'   `$samples`). Values are on the linear scale (not logged), finite and
#'   non-negative, with gene symbols as row names.
#' @param assay Assay name or index to use from a `SummarizedExperiment`.
#' @param prior_count Offset added before log2. Default `1`.
#' @return A `DESeqDataSet` with an extra `"expr_values"` assay.
#' @export
expr_to_dds <- function(object, assay = 1L, prior_count = 1) {
    if (inherits(object, "DGEList")) {
        x <- object$counts
        cd <- object$samples[, setdiff(colnames(object$samples),
                                       c("lib.size", "norm.factors")), drop = FALSE]
    } else if (methods::is(object, "SummarizedExperiment")) {
        x <- SummarizedExperiment::assay(object, assay)
        cd <- SummarizedExperiment::colData(object)
    } else {
        stop("`object` must be a SummarizedExperiment or an edgeR DGEList.",
             call. = FALSE)
    }
    x <- as.matrix(x)
    if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0)) {
        stop("Expression values must be finite and non-negative (TPM / FPKM ",
             "on the linear scale, not log-transformed).", call. = FALSE)
    }
    if (is.null(rownames(x)) || anyDuplicated(rownames(x))) {
        stop("Expression values need unique row names (gene symbols).", call. = FALSE)
    }
    if (!is.numeric(prior_count) || length(prior_count) != 1L ||
            is.na(prior_count) || prior_count <= 0) {
        stop("`prior_count` must be a single positive number.", call. = FALSE)
    }
    stand_in <- round(x)
    storage.mode(stand_in) <- "integer"
    dds <- quiet(DESeq2::DESeqDataSetFromMatrix(
        stand_in, colData = as.data.frame(cd), design = ~1
    ))
    SummarizedExperiment::assay(dds, .expr_assay, withDimnames = FALSE) <- x
    S4Vectors::metadata(dds)$degreports_prior_count <- prior_count
    dds
}

#' DESeq2-style coefficient names for a model matrix
#'
#' The contrast grid names coefficients as `DESeq2::resultsNames()` does:
#' `Intercept`, factor main effects as `<var>_<level>_vs_<ref>`, interaction
#' terms through [make.names()] (`a:b` -> `a.b`).
#' @noRd
.deseq2_coef_names <- function(mm, design, coldata) {
    nm <- make.names(sub("^\\(Intercept\\)$", "Intercept", colnames(mm)))
    for (v in all.vars(design)) {
        lv <- levels(coldata[[v]])
        if (length(lv) < 2L) next
        hit <- match(nm, make.names(paste0(v, lv[-1L])))
        nm[!is.na(hit)] <- make.names(paste0(v, "_", lv[-1L], "_vs_", lv[1L]))[hit[!is.na(hit)]]
    }
    nm
}

#' Fit the limma linear model to a carrier under its design
#' @return List with `fit` (`MArrayLM`), `coef_names` (DESeq2-style).
#' @noRd
.limma_fit <- function(dds) {
    des <- DESeq2::design(dds)
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    mm <- stats::model.matrix(des, data = cd)
    list(fit = limma::lmFit(.log_expr(dds), mm),
         coef_names = .deseq2_coef_names(mm, des, cd))
}

#' limma-trend test of one coefficient contrast, as a DESeq2-schema table
#'
#' `stat` is the moderated t, `lfcSE` its standard error, `padj` BH.
#' `log2FoldChange` is ashr-shrunk unless `shrinkage = "none"`; the limma
#' estimate stays in `log2FoldChange_unshrunken`.
#' @param weights Contrast weights over `prefit$coef_names`.
#' @param dds The full object, for `baseMean` / `mean_cpm`.
#' @noRd
.limma_results <- function(prefit, weights, dds, shrinkage, alpha, lfc_threshold,
                           ranking, gene_col, comparison) {
    fit <- limma::contrasts.fit(prefit$fit, contrasts = matrix(weights))
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    tab <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
    lfc <- tab$logFC
    se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, 1]

    shrunk <- lfc
    if (shrinkage != "none") {
        if (shrinkage != "ashr") {
            warning("limma uses ashr shrinkage, not \"", shrinkage, "\".", call. = FALSE)
        }
        shrunk <- tryCatch(
            ashr::ash(lfc, se, mixcompdist = "normal", method = "shrink")$result$PosteriorMean,
            error = function(e) {
                warning("LFC shrinkage (ashr) failed for \"", comparison, "\": ",
                        conditionMessage(e), ". Returning unshrunken results.",
                        call. = FALSE)
                lfc
            })
    }

    genes <- rownames(tab)
    res_df <- data.frame(
        gene = genes, comparison = comparison,
        baseMean = rowMeans(.raw_count_matrix(dds))[genes],
        mean_cpm = .mean_cpm(dds)[genes],
        log2FoldChange = shrunk, log2FoldChange_unshrunken = lfc,
        lfcSE = unname(se), stat = tab$t, pvalue = tab$P.Value, padj = tab$adj.P.Val,
        row.names = NULL, stringsAsFactors = FALSE
    )
    names(res_df)[1L] <- gene_col
    res_df$rank_score <- compute_ranking_score(res_df, ranking)
    res_df$sig <- !is.na(res_df$padj) & res_df$padj < alpha &
        abs(res_df$log2FoldChange) > lfc_threshold
    res_df
}

#' limma-trend for one level comparison -- [run_diffexp()]'s counterpart
#' @noRd
run_limma_diffexp <- function(dds, var, comparison, shrinkage = "ashr",
                              alpha = 0.05, lfc_threshold = 0.5,
                              ranking = "log2FoldChange", gene_col = "gene") {
    setup <- .comparison_fit_dds(dds, var, comparison)
    fit_dds <- setup$dds
    # Reference level first: the test level is then one treatment coefficient.
    cd <- SummarizedExperiment::colData(fit_dds)
    cd[[setup$var]] <- stats::relevel(droplevels(factor(cd[[setup$var]])), ref = setup$ref)
    SummarizedExperiment::colData(fit_dds) <- cd
    assert_full_rank_model_matrix(fit_dds, comparison, var)

    prefit <- .limma_fit(fit_dds)
    coef <- make.names(paste0(setup$var, "_", setup$test, "_vs_", setup$ref))
    weights <- as.numeric(prefit$coef_names == coef)
    if (!any(weights)) {
        stop("No \"", setup$test, "\" coefficient for comparison \"", comparison,
             "\" in the limma model.", call. = FALSE)
    }
    .limma_results(prefit, weights, dds, shrinkage, alpha, lfc_threshold,
                   ranking, gene_col, comparison)
}

#' limma-trend on named coefficients -- [run_diffexp_coef()]'s counterpart
#' @param prefit List from [.limma_fit()]; fitted here when `NULL`.
#' @noRd
run_limma_diffexp_coef <- function(dds, coef = NULL, coef_subtract = NULL,
                                   coef_weights = NULL, shrinkage = "ashr",
                                   alpha = 0.05, lfc_threshold = 0.5,
                                   ranking = "log2FoldChange", gene_col = "gene",
                                   comparison = paste(coef, collapse = "+"),
                                   prefit = NULL) {
    if (is.null(prefit)) prefit <- .limma_fit(dds)
    rn <- prefit$coef_names
    w <- stats::setNames(numeric(length(rn)), rn)
    if (length(coef_weights)) {
        if (is.character(coef_weights)) coef_weights <- .parse_coef_weights(coef_weights)
        wanted <- names(coef_weights)
    } else {
        wanted <- c(.clean_coef_names(coef), .clean_coef_names(coef_subtract))
    }
    unknown <- setdiff(wanted, rn)
    if (length(unknown)) {
        stop("Coefficient(s) ", paste0("\"", unknown, "\"", collapse = ", "),
             " not in the limma model for contrast \"", comparison, "\". Available: ",
             paste(rn, collapse = ", "), call. = FALSE)
    }
    if (length(coef_weights)) {
        w[names(coef_weights)] <- coef_weights
    } else {
        w[.clean_coef_names(coef)] <- 1
        sub <- .clean_coef_names(coef_subtract)
        w[sub] <- w[sub] - 1
    }
    .limma_results(prefit, unname(w), dds, shrinkage, alpha, lfc_threshold,
                   ranking, gene_col, comparison)
}

#' Run the complete pipeline on TPM / FPKM data with limma-trend
#'
#' The limma-trend counterpart of [run_deseq2()] for expression that is
#' not raw counts (TPM, FPKM, RPKM). Same stages, files and HTML report as
#' [run_deseq2()]; only the statistics differ:
#'
#' \describe{
#'   \item{Differential expression}{limma-trend on log2(x + `prior_count`):
#'     `lmFit()`, `contrasts.fit()`, `eBayes(trend = TRUE, robust = TRUE)`,
#'     for every contrast in the grid (pairwise, one-vs-rest, main effects,
#'     interactions). Results keep the DESeq2 columns: `stat` is the
#'     moderated t, `padj` the BH FDR, and `log2FoldChange` is ashr-shrunk
#'     unless `opts$shrinkage = "none"` (`log2FoldChange_unshrunken` keeps
#'     limma's `logFC`).}
#'   \item{Sample-level expression}{The same log2(x + `prior_count`) matrix
#'     wherever the DESeq2 pipeline uses a blind VST or rlog: QC, PCA, GSVA,
#'     WGCNA, decoupleR / VIPER sample scores, CSPA heatmaps, count exports.}
#'   \item{Batch correction}{With `batch`, ComBat on that matrix (ComBat-seq
#'     needs integer counts); `batch` also enters every limma design.}
#' }
#'
#' There is no LRT / time-course mode; use [`run_deseq2(mode = "lrt")`][run_deseq2] on counts.
#'
#' @inheritParams run_deseq2
#' @inheritParams expr_to_dds
#' @return Invisibly, the same list as [run_deseq2()]; `dds` is the
#'   QC-filtered carrier, whose `"expr_values"` assay holds the input values.
#' @examples
#' \dontrun{
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(tpm = tpm_matrix), colData = coldata)
#' out <- run_limma(se, var = "condition", save_dir = "results_limma")
#' }
#' @export
run_limma <- function(
    object,
    var,
    assay = 1L,
    prior_count = 1,
    design = NULL,
    comparisons = NULL,
    include_interactions = FALSE,
    stratified = FALSE,
    batch = NULL,
    export_counts = TRUE,
    filter_method = "cpm",
    cpm_percentile = 0.05,
    min_count = 10,
    min_samples_per_group = NULL,
    remove_mt = FALSE,
    remove_xy = FALSE,
    mt_pattern = "^MT-|^mt-",
    gene_chr = NULL,
    remove_vdj = NULL,
    ranking = "log2FoldChange",
    collections = c("GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "PID"),
    pals = NULL,
    shape_var = NULL,
    n_scree_pcs = 10,
    n_loading_pcs = 6,
    n_top_genes = 10,
    top_label_n = 5,
    top_n_plot = 10,
    run_clusterprofiler = FALSE,
    run_gsva = FALSE,
    run_wgcna = FALSE,
    run_decoupler_tf = FALSE,
    run_decoupler_ccc = FALSE,
    run_viper = FALSE,
    run_cspa = TRUE,
    render_html = TRUE,
    parallel_render = FALSE,
    render_workers = NULL,
    use_targets = FALSE,
    opts = degreports_opts(),
    save_dir = getwd()) {

    call <- match.call()
    args <- mget(names(formals()), envir = environment())
    .write_repro_script(
        "run_limma", data_arg = "object", data_var = "object",
        data_comment = "# The TPM / FPKM object passed to run_limma(), before QC",
        args = args, call = call, save_dir = save_dir,
        run_comment = "Run QC, PCA, limma-trend, GSEA and the enabled modules, then render the HTML report")
    dds <- expr_to_dds(object, assay = assay, prior_count = prior_count)
    engine_args <- args[setdiff(names(args), c("object", "assay", "prior_count"))]
    do.call(.run_wald_engine, c(
        list(dds = dds, platform = "rnaseq"), engine_args,
        list(engine = "limma", call = call,
             engine_params = list(prior_count = prior_count))
    ), quote = TRUE)
}
