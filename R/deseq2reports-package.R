#' deseq2reports: Modular DESeq2 Pipeline with Automated HTML Reports
#'
#' A self-contained, modular pipeline for RNA-seq differential expression
#' analysis built on DESeq2. See [run_deseq2_pipeline()] for the one-command
#' entry point, or use the standalone stage functions [run_qc()],
#' [run_distance_pca()], [run_diffexp_pipeline()], [run_gsea_pipeline()],
#' [run_wgcna_pipeline()], [run_decoupler_pipeline()], [run_viper_pipeline()],
#' [run_ppi_pipeline()], and [render_report()] independently.
#'
#' Package vignettes (HTML):
#' \itemize{
#'   \item \code{vignette("deseq2reports")} -- one-command and modular usage
#'   \item \code{vignette("deseq2reports-modules")} -- detailed module reference
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats model.matrix as.formula relevel prcomp dist hclust setNames
#' @importFrom utils packageVersion write.csv sessionInfo capture.output
#' @importFrom methods is
#' @importFrom rlang .data
## usethis namespace: end
NULL
