#' Fixed analysis settings for the pipeline entry points
#'
#' These used to be arguments of [run_deseq2()], [`run_deseq2(mode = "lrt")`][run_deseq2], and
#' [run_deg()]. They are now fixed so that every run is comparable and
#' the "Pipeline arguments" table stays about study design rather than
#' tuning knobs. The underlying module functions
#' ([run_diffexp_pipeline()], [run_wgcna_pipeline()], ...) still expose them,
#' so a direct caller can deviate deliberately.
#'
#' `cspa_n_markers` differs by mode: the Wald / DEG-table pipelines pick
#' markers per contrast, the LRT pipeline picks them per degPatterns cluster
#' (of which there are many more), hence the smaller per-cluster budget.
#' @noRd
.pipeline_const <- list(
    shrinkage = "ashr",
    alpha = 0.05,
    lfc_threshold = 0.5,
    gsva_n_heatmap = 50L,
    wgcna_min_module_size = 100L,
    decoupler_n_tfs = 25L,
    decoupler_n_lr = 10L,
    viper_n_tfs = 10L,
    ppi_top_n_hubs = 10L,
    cspa_n_markers = 25L,
    cspa_n_markers_lrt = 10L
)
