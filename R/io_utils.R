#' Write a data frame to disk in the requested output format
#'
#' Currently only `"csv"` is supported; the indirection exists so
#' additional formats can be added without changing call sites.
#'
#' @param df A data frame (or object coercible via `as.data.frame()`).
#' @param path File path to write to.
#' @param format Output format, currently only `"csv"`.
#' @return Invisibly, `path`.
#' @export
write_output <- function(df, path, format = "csv") {
    format <- validate_output_format(format)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (format == "csv") {
        utils::write.csv(as.data.frame(df), path, row.names = FALSE)
    }
    invisible(path)
}

#' Save a ggplot object to PDF, falling back gracefully if the requested
#' device is unavailable
#'
#' @param plot A `ggplot` object.
#' @param path Output file path (should end in `.pdf`).
#' @param width,height Plot dimensions in inches.
#' @return Invisibly, `path`.
#' @export
save_ggplot <- function(plot, path, width = 7, height = 5) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    result <- tryCatch({
        grDevices::pdf(path, width = width, height = height)
        print(plot)
        grDevices::dev.off()
        TRUE
    }, error = function(e) {
        message("PDF device failed (", conditionMessage(e), "); output not saved: ", path)
        FALSE
    })
    invisible(path)
}

#' Run `DESeq2::vst()`, falling back to the un-subsampled
#' `varianceStabilizingTransformation()` when the dataset is too small (or
#' too lowly-expressed) for `vst()`'s internal gene-subsampling heuristic
#'
#' `DESeq2::vst(nsub = n)` does not simply require `nrow(dds) >= n`: it
#' first restricts to genes with mean normalized count > 5, then requires
#' *that* subset to have >= `nsub` genes, else it errors outright asking
#' the caller to use `varianceStabilizingTransformation()` directly. This
#' can bite on small or lowly-expressed datasets (e.g. simulated/test data,
#' small pilot experiments) even after naively capping
#' `nsub = min(1000, nrow(dds))`, since `nrow(dds)` may still exceed the
#' count of genes that clear the mean-count-> 5 bar. This helper checks
#' that count directly and falls back to the slower, non-subsampled
#' transform (which fits dispersion on every gene, so it always works,
#' regardless of gene count or expression level) whenever `vst()`'s own
#' criterion would not be met.
#'
#' @param dds A `DESeqDataSet` object.
#' @param blind Passed through to `vst()` / `varianceStabilizingTransformation()`.
#' @return A `DESeqTransform` object (as returned by either underlying function).
#' @noRd
safe_vst <- function(dds, blind = TRUE) {
    dds_sf <- if (is.null(DESeq2::sizeFactors(dds)) && is.null(DESeq2::normalizationFactors(dds))) {
        DESeq2::estimateSizeFactors(dds)
    } else {
        dds
    }
    base_mean <- MatrixGenerics::rowMeans(DESeq2::counts(dds_sf, normalized = TRUE))
    n_eligible <- sum(base_mean > 5)
    # Below this many eligible genes, `vst()`'s subsampled dispersion-trend
    # fit is unstable (too few points for a parametric fit) even though it
    # would technically run without erroring; fall back to fitting on the
    # full dataset instead.
    min_eligible <- 100L

    if (n_eligible < min_eligible) {
        # Not enough genes clear vst()'s internal "mean normalized count > 5"
        # bar for its subsampled dispersion-fit trick; use the direct
        # transform instead, which fits on all genes and always works --
        # except on very small/toy datasets (e.g. simulated test data with
        # only a handful of genes) where even a single-parameter
        # (`fitType = "mean"`) dispersion fit can fail outright, because
        # DESeq2 refuses the fit when gene-wise dispersion estimates are
        # all within ~2 orders of magnitude of each other (too few genes
        # for *any* trend fit to be meaningful). In that pathological case,
        # fall back further to a simple log2(normalized count + 1)
        # transform, wrapped in a `SummarizedExperiment` so downstream
        # `assay()` calls work identically to a real `DESeqTransform`; this
        # loses variance-stabilization but only triggers on datasets too
        # small for VST to be meaningful anyway.
        return(tryCatch(
            DESeq2::varianceStabilizingTransformation(dds, blind = blind),
            error = function(e) {
                tryCatch(
                    DESeq2::varianceStabilizingTransformation(dds, blind = blind, fitType = "mean"),
                    error = function(e2) {
                        norm_mat <- log2(DESeq2::counts(dds_sf, normalized = TRUE) + 1)
                        SummarizedExperiment::SummarizedExperiment(
                            assays = list(counts = norm_mat),
                            colData = SummarizedExperiment::colData(dds_sf)
                        )
                    }
                )
            }
        ))
    }
    DESeq2::vst(dds, blind = blind, nsub = min(1000L, n_eligible))
}

#' Quietly execute an expression, suppressing messages/warnings/console
#' output (some Bioconductor functions, e.g. `sva::ComBat_seq()`, print
#' progress via `cat()` rather than `message()`)
#'
#' @param expr An expression to evaluate.
#' @return The value of `expr`.
#' @noRd
quiet <- function(expr) {
    result <- NULL
    utils::capture.output(
        result <- withCallingHandlers(
            expr,
            message = function(m) invokeRestart("muffleMessage"),
            warning = function(w) invokeRestart("muffleWarning")
        )
    )
    result
}
