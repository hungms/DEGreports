#' Export count tables in multiple normalizations
#'
#' Writes raw counts, variance-stabilized (VST) counts, regularized-log
#' (rlog) counts, and counts-per-million (CPM) as CSV files to `save_dir`.
#' Optionally also writes a batch-corrected table: the VST matrix (log2
#' expression for a limma carrier) after `limma::removeBatchEffect()`,
#' protecting `var`. It is on the log scale, for clustering and plotting --
#' not a corrected count matrix to model. TPM is not produced --
#' it requires a gene-length table that a `DESeqDataSet` does not carry, and
#' is out of scope for this function; compute it separately if you have gene
#' lengths available.
#'
#' All exports other than the batch-corrected table are derived directly
#' from `dds` as given (typically the QC-filtered object returned by
#' [run_qc()]) and are independent of the differential-expression model.
#'
#' @param dds A `DESeqDataSet` object (typically post-[run_qc()]).
#' @param batch Column name(s) in `colData(dds)` to correct for (batch
#'   and/or the `pair` column of a paired design), or `NULL` (default) to
#'   skip the batch-corrected export entirely.
#' @param var Column name in `colData(dds)` identifying the biological
#'   group of interest. Required if `batch` is supplied: the correction
#'   protects this signal rather than removing it along with the batch
#'   effect. Ignored if `batch` is `NULL`.
#' @param transform One of `"vst"` (default) or `"rlog"` -- which
#'   variance-stabilizing transform to use for the VST/rlog export. Both
#'   are always computed at `blind = TRUE` (i.e. ignoring the design), since
#'   these tables are for general-purpose downstream use (clustering,
#'   visualization, exporting), not the DE test itself.
#' @param save_dir Directory to write outputs to. Default `getwd()`.
#' @return Invisibly, a named list of the data frames written (`raw`,
#'   `vst`, `rlog`, `cpm`, and `batch_corrected` if `batch` was supplied),
#'   each with genes as rows (first column `gene`) and samples as columns.
#' @examples
#' \dontrun{
#' out <- run_counts_export(dds, save_dir = "results")
#' out_bc <- run_counts_export(dds, batch = "batch", var = "condition",
#'                              save_dir = "results")
#' }
#' @export
run_counts_export <- function(
    dds,
    batch = NULL,
    var = NULL,
    transform = "vst",
    save_dir = getwd()) {

    dds <- validate_dds(dds)
    save_dir <- validate_save_dir(save_dir)
    transform <- validate_choice(transform, c("vst", "rlog"), "transform")
    if (!is.null(batch)) {
        dds <- validate_batch(dds, batch)
        if (is.null(var) || !is.character(var) || length(var) != 1) {
            stop("`var` must be supplied (as a single column name in colData(dds)) when `batch` is set, ",
                 "so the batch correction can protect this biological signal while removing the batch effect.",
                 call. = FALSE)
        }
        dds <- validate_var(dds, var)
    }

    raw_counts <- .raw_count_matrix(dds)

    message("Computing VST counts...")
    vst_mat <- SummarizedExperiment::assay(quiet(safe_vst(dds, blind = TRUE)))

    message("Computing rlog counts...")
    rlog_mat <- SummarizedExperiment::assay(quiet(safe_rlog(dds, blind = TRUE)))

    message("Computing CPM...")
    lib_sizes <- colSums(raw_counts)
    cpm_mat <- sweep(raw_counts, 2, lib_sizes, "/") * 1e6

    to_df <- function(mat) {
        out <- as.data.frame(mat, check.names = FALSE)
        out <- cbind(gene = rownames(mat), out)
        rownames(out) <- NULL
        out
    }

    result <- list(
        raw  = to_df(raw_counts),
        vst  = to_df(vst_mat),
        rlog = to_df(rlog_mat),
        cpm  = to_df(cpm_mat)
    )

    write_output(result$raw,  file.path(save_dir, "counts_raw.csv"))
    write_output(result$vst,  file.path(save_dir, "counts_vst.csv"))
    write_output(result$rlog, file.path(save_dir, "counts_rlog.csv"))
    write_output(result$cpm,  file.path(save_dir, "counts_cpm.csv"))

    if (!is.null(batch)) {
        message("Running removeBatchEffect (removing \"", paste(batch, collapse = "\", \""),
                "\", protecting \"", var, "\")...")
        corrected <- batch_corrected_matrix(dds, batch = batch, var = var,
                                            transform = transform)
        outfile <- if (.is_expr_dds(dds)) {
            "counts_batch_corrected_log2expr.csv"
        } else {
            paste0("counts_batch_corrected_", transform, ".csv")
        }
        result$batch_corrected <- to_df(corrected)
        write_output(result$batch_corrected, file.path(save_dir, outfile))
    }

    invisible(result)
}
