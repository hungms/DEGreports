#' Run the complete RNA-seq DESeq2 pipeline in one command
#'
#' Chains quality control ([run_qc()]), sample distance/PCA
#' ([run_distance_pca()]), differential expression, gene set enrichment, the
#' optional modules and an HTML report ([render_report()]) into a single call,
#' in one of two modes:
#'
#' \describe{
#'   \item{`mode = "wald"` (default)}{Wald tests for every contrast in the grid
#'     built from `var` ([run_diffexp_pipeline()]), ranked GSEA
#'     ([run_gsea_pipeline()]).}
#'   \item{`mode = "lrt"`}{Likelihood-ratio tests of `time` and/or `var`
#'     over the nested full/reduced grid from [expand_lrt_specs()]
#'     ([run_diffexp_lrt_pipeline()]). Significant genes (`padj < alpha`) are
#'     clustered with `DEGreport::degPatterns()` (x = `time`, colour = first
#'     of `var`); the clusters drive per-cluster ORA
#'     ([run_cluster_ora_pipeline()]), STRING PPI, CSPA heatmaps and the
#'     per-cluster decoupleR / VIPER signatures, scored by the LRT statistic.}
#' }
#'
#' For TPM / FPKM data use [run_limma()]; for a pre-computed DEG table use
#' [run_deg()].
#'
#' Optional modules (default off): `run_gsva`, `run_wgcna`, `run_decoupler_tf`,
#' `run_decoupler_ccc`, `run_viper`. Other stages always run: GSEA / cluster
#' ORA, CSPA (`run_cspa = TRUE`), and STRING PPI via [run_ppi_pipeline()]
#' for every DE contrast and (when enabled) every WGCNA module. PPI needs
#' network access on first run to download STRING files; failures are caught
#' with a warning so the rest of the pipeline still finishes.
#'
#' Statistical and plotting thresholds (LFC shrinkage, `alpha`, the LFC
#' cutoff, heatmap/barplot sizes, WGCNA module size, hub/TF/marker counts)
#' are fixed -- see `.pipeline_const`. Call the module functions directly
#' ([run_diffexp_pipeline()], [run_wgcna_pipeline()], ...) to change them.
#'
#' @param dds A `DESeqDataSet` object containing raw counts.
#' @param mode `"wald"` (default) or `"lrt"`; see Description.
#' @param var Column name(s) in `colData(dds)` used for grouping. In LRT mode
#'   these are the non-time factors, and `NULL` gives a time-only design.
#'   Length 1
#'   = single-factor mode; length 2 = two-factor mode; length 3 = three-factor
#'   mode (e.g. `celltype * C_mutation * K_mutation`). Every factor's main
#'   effect is always tested across all levels of the others;
#'   `include_interactions = TRUE` makes the design the full product of the
#'   factors and adds the difference-of-differences interaction contrasts.
#'   Comparisons / the design formula are auto-built from these columns unless
#'   overridden.
#' @param time LRT mode only (required there): `colData` column holding the
#'   time / ordered trajectory factor (discrete levels).
#' @param degpatterns_minc LRT mode only: minimum genes per degPatterns
#'   cluster (passed to [DEGreport::degPatterns()]). Default `15`.
#' @param design Wald mode only. Design formula (or character string) for the DESeq2 model,
#'   or `NULL` (default) to auto-build via [resolve_design()] from `var`,
#'   `batch`, and `include_interactions` (`~ var`, `~ batch + var`,
#'   additive `~ A + B` / `~ A + B + C`, or `~ A * B` / `~ A * B * C` when
#'   interactions are requested).
#' @param comparisons Wald mode only. Character vector of comparisons / contrast labels to
#'   run (`"A_vs_B"` / `"A_vs_B+C"`, or multi-factor labels like
#'   `"C_vs_noC_x_K_vs_noK_in_CB"`), or `NULL` (default) to auto-generate
#'   the full grid via [expand_contrast_specs()].
#' @param include_interactions Logical. When `var` has length ≥ 2 and
#'   `design = NULL`, use `~ A * B` / `~ A * B * C` and add the
#'   difference-of-differences interaction contrasts — "does the K effect
#'   differ between C and noC, within CB"
#'   (`C_vs_noC_x_K_vs_noK_in_CB`), plus the three-way row for three
#'   factors — to the main effects that are always expanded. Each combines
#'   coefficients from one fit of the product design on all samples. A factor
#'   with three or more levels is compared one-vs-rest
#'   (`CA_vs_CB+CC_x_K_vs_noK_in_C`), its pooled side an average with equal
#'   weight per level. Default `FALSE`: additive design, main effects only.
#'   Ignored for single-factor `var`. When `design` is supplied explicitly,
#'   only the contrast grid is affected.
#' @param stratified Wald mode only. Logical. Also expand the within-stratum contrasts --
#'   each factor compared inside one fixed level combination of the others,
#'   fitted on that subset alone (`CB__K_vs_noK` for two factors,
#'   `K_vs_noK/CB__C` for three). Default `FALSE`: an interaction contrast
#'   asks the same question with every sample behind it, so these are opt-in,
#'   but they are the right tool when the strata are not expected to share a
#'   dispersion. Ignored for single-factor `var`.
#' @param batch Optional covariate column name(s) in `colData(dds)`. Added
#'   to every model (`~ batch + var`), and removed from the expression matrix
#'   behind the PCA, sample similarity, count export and the sample-level
#'   modules (WGCNA / GSVA / sample-level ULM and VIPER / CSPA heatmaps) with
#'   `limma::removeBatchEffect()` on the VST matrix, protecting `var`.
#'   DESeq2 itself always uses raw counts with `batch` in the design.
#'   Default `NULL`.
#' @param pair Optional `colData(dds)` column naming the unit measured
#'   under each condition in a paired design -- a donor, patient or animal.
#'   It is added to every model as a blocking term (`~ [batch +] pair +
#'   var`), so each comparison is made within units, and it is removed from
#'   the sample-level matrix together with `batch`. Every tested factor must
#'   vary within at least one unit (a factor that is constant within every
#'   unit cannot be estimated in a paired design, and is an error); units
#'   seen in only one condition trigger a warning. Default `NULL`.
#' @param export_counts Logical. If `TRUE` (default), also writes raw,
#'   VST, rlog, and CPM count tables via [run_counts_export()]. If `batch`
#'   is set, a batch-corrected table is written as well (protecting `var`
#'   as the biological signal of interest). These tables are independent of
#'   the DE model and provided for general downstream use.
#' @param filter_method Low-expression gene-filtering strategy: `"cpm"`
#'   (default, percentile-based), `"count"`, or `NULL` (no expression
#'   filtering). See [run_qc()].
#' @param cpm_percentile Fraction of the *expressed* genes dropped on top of
#'   the genes with no reads at all, when `filter_method = "cpm"`. See
#'   [run_qc()].
#' @param min_count,min_samples_per_group Used when `filter_method =
#'   "count"`. See [run_qc()].
#' @param remove_mt,remove_xy,mt_pattern,gene_chr,remove_vdj See [run_qc()].
#' @param ranking DEG/GSEA ranking metric: `"log2FoldChange"`, `"padj"`,
#'   `"pxfc"`, or `"stat"`. Default `NULL`: `"log2FoldChange"` in Wald mode,
#'   `"stat"` in LRT mode.
#' @param collections MSigDB collections for GSEA: one or more of
#'   `"HALLMARK"`, `"GOBP"`, `"GOCC"`, `"GOMF"`, `"KEGG"`,
#'   `"REACTOME"`, `"PID"`. Default `c("GOBP", "GOCC",
#'   "GOMF", "KEGG", "REACTOME", "PID")` --
#'   `"HALLMARK"` remains available but is not forced on by default. GSEA
#'   runs independently against every requested collection, for every
#'   comparison. Organism (human/mouse) is auto-detected from gene symbols
#'   -- there is no `org` argument.
#' @param pals Optional colour palette (character vector), one colour per
#'   level of the visualization grouping factor. Default `NULL`.
#' @param shape_var Optional column name for a second grouping variable
#'   mapped to point shape in the PCA plot. Default `NULL`: the `pair`
#'   column when `pair` is set, else the first `var` in multi-factor mode.
#' @param n_scree_pcs,n_loading_pcs,n_top_genes See [run_distance_pca()].
#' @param top_label_n Number of top **up** and top **down** genes (by
#'   `rank_score`) labelled per volcano plot. Default `5`.
#' @param top_n_plot Number of top pathways per direction shown per GSEA
#'   barplot (regardless of significance). Default `10`.
#' @param run_clusterprofiler Logical. If `TRUE` (default `FALSE`), also
#'   runs a second, independent GSEA pass via `clusterProfiler::GSEA()`
#'   alongside the primary `fgsea`-based GSEA, exporting the raw
#'   `clusterProfiler` `gseaResult` objects. See [run_gsea_pipeline()].
#' @param run_gsva Logical. If `TRUE` (default `FALSE`), runs sample-level GSVA
#'   via [run_gsva_pipeline()] and adds a **GSVA** HTML page (Imports:
#'   `GSVA`). Independent of clusterProfiler / fgsea GSEA. HALLMARK is
#'   always scored and heatmapped in full. Large collections (GO / KEGG /
#'   REACTOME / PID) are scored in full (CSV) but the heatmap shows the
#'   most-associated pathways for the grouping `var` (Kruskal–Wallis on
#'   GSVA scores; variance fallback). GSEA `padj` is not used — every set
#'   is scored, so a significance cutoff would still leave hundreds of rows.
#' @param run_wgcna Logical. If `TRUE` (default `FALSE`), also runs
#'   weighted gene co-expression network analysis via
#'   [run_wgcna_pipeline()] -- soft-threshold power selection, module
#'   detection, condition- and sample-level module-trait correlation, hub
#'   genes, and per-module over-representation analysis against
#'   `"GOBP"`, `"GOCC"`, `"GOMF"`, `"KEGG"`, and
#'   `"REACTOME"`.
#'   Requires the `WGCNA` package. This is a heavier, exploratory
#'   companion analysis independent of the `var`-based differential
#'   expression comparisons above; it uses every sample's full expression
#'   profile rather than testing one gene at a time against a design.
#' @param run_decoupler_tf Logical. If `TRUE` (default `FALSE`), runs
#'   CollecTRI + ULM TF activity analysis via [run_decoupler_tf_pipeline()]
#'   and adds a **decoupleR-TF** HTML page (Imports: `decoupleR`,
#'   `OmnipathR` for network download when no local prior is present).
#' @param run_decoupler_ccc Logical. If `TRUE` (default `FALSE`), runs
#'   ligand–receptor (cell–cell communication) ULM via
#'   [run_decoupler_ccc_pipeline()] and adds a **decoupleR-CCC** hub plus
#'   one HTML page per comparison. Independent of `run_decoupler_tf`.
#' @param run_viper Logical. If `TRUE` (default `FALSE`), runs
#'   [VIPER](https://califano.c2b2.columbia.edu/viper) / msVIPER with the
#'   human B-cell regulon from [bcellViper](https://bioconductor.org/packages/bcellViper/)
#'   via [run_viper_pipeline()] and adds a **VIPER** HTML page (Imports:
#'   `viper`, `bcellViper`, and for mouse `babelgene`).
#' @param run_cspa Logical. If `TRUE` (default), selects top Cell Surface
#'   Protein Atlas (CSPA) markers among DEGs per contrast and draws
#'   all-sample VST heatmaps (plus a combined heatmap) via
#'   [run_cspa_marker_pipeline()] (Imports: `readxl`).
#' @param render_html Logical. If `TRUE` (default), renders an HTML summary
#'   report via [render_report()] after all stages complete.
#' @param parallel_render If `TRUE`, knit per-contrast HTML pages in
#'   parallel on Unix (`parallel::mclapply()`). Hub pages stay sequential
#'   in the current R session. Passed to [render_report()]. Default `FALSE`.
#' @param render_workers Parallel worker count for `parallel_render`.
#'   Default `NULL` (auto).
#' @param use_targets If `TRUE`, render HTML via [render_report_targets()]
#'   (requires Suggested package `targets`) instead of [render_report()].
#'   Default `FALSE`.
#' @param opts Analysis settings from [degreports_opts()] -- LFC
#'   shrinkage, `alpha`, the LFC cutoff, and the per-module top-N budgets.
#'   The defaults reproduce the values these entry points used to hard-code,
#'   and only settings that differ from them are listed in the report's
#'   "Pipeline arguments" table. Default `degreports_opts()`.
#' @param save_dir Directory where all outputs will be written. Default
#'   `getwd()`.
#' @return Invisibly, a list with elements `dds` (the QC-filtered object),
#'   `dds_path` (path to `dds.rds` under `save_dir`),
#'   `distance` (output of [run_distance_pca()]), `diffexp` (output of
#'   [run_diffexp_pipeline()]), `gsea` (output of [run_gsea_pipeline()]),
#'   `wgcna` (output of [run_wgcna_pipeline()], or `NULL` if `run_wgcna =
#'   FALSE` / the stage failed), `ppi` (output of [run_ppi_pipeline()], or
#'   `NULL` if the stage failed), `cspa` (output of
#'   [run_cspa_marker_pipeline()], or `NULL` if `run_cspa = FALSE`), `viper`
#'   (output of [run_viper_pipeline()], or `NULL` if `run_viper = FALSE`), and
#'   `report_path` (a named character vector of the paths to the rendered
#'   HTML report pages -- `"index"`, `"qc"`, `"deseq2"`,
#'   `"clusterProfiler"`, and (if enabled) `"wgcna"` / `"decoupleR"` /
#'   `"viper"`; see [render_report()] -- or `NULL` if `render_html =
#'   FALSE`). Module artifacts and their HTML pages are written under
#'   per-module folders (`qc/` incl. `counts/` and `distance/`, `deseq2/`,
#'   `clusterprofiler/`, `wgcna/` incl. per-module `ppi/`, `decoupler/`,
#'   `viper/`; per-contrast CSPA/PPI under `deseq2/<contrast>/cspa/` and
#'   `deseq2/<contrast>/ppi/`; LRT-cluster PPI under `ppi/<cluster>/ppi/`;
#'   GSEA overlap under `clusterprofiler/<contrast>/gsea_overlap/`) with
#'   `index.html` and the QC-filtered `dds.rds` at the root of `save_dir`.
#'   LRT mode also returns `degpatterns` and `cluster_ora`, with `gsea` `NULL`
#'   (cluster ORA replaces ranked GSEA).
#' @examples
#' \dontrun{
#' out <- run_deseq2(dds, var = "condition", save_dir = "results")
#'
#' out2 <- run_deseq2(dds, var = c("Tissue", "Isotype"),
#'                    save_dir = "results_2factor")
#'
#' out3 <- run_deseq2(dds, mode = "lrt", var = "genotype", time = "day",
#'                    include_interactions = TRUE, save_dir = "results_lrt")
#' }
#' @export
run_deseq2 <- function(
    dds,
    var = NULL,
    mode = "wald",
    time = NULL,
    design = NULL,
    comparisons = NULL,
    include_interactions = FALSE,
    stratified = FALSE,
    batch = NULL,
    pair = NULL,
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
    ranking = NULL,
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
    degpatterns_minc = 15,
    render_html = TRUE,
    parallel_render = FALSE,
    render_workers = NULL,
    use_targets = FALSE,
    opts = degreports_opts(),
    save_dir = getwd()) {

    call <- match.call()
    mode <- validate_choice(mode, c("wald", "lrt"), "mode")
    # Arguments that belong to the other mode are an error, not silently dropped.
    only <- list(wald = c("design", "comparisons", "stratified"),
                 lrt = c("time", "degpatterns_minc"))
    other <- setdiff(names(only), mode)
    wrong <- intersect(names(as.list(call)), only[[other]])
    if (length(wrong)) {
        stop("`", paste(wrong, collapse = "`, `"), "` only appl",
             if (length(wrong) == 1L) "ies" else "y", " to mode = \"", other, "\".",
             call. = FALSE)
    }
    if (mode == "lrt" && is.null(time)) {
        stop("`time` is required when mode = \"lrt\".", call. = FALSE)
    }
    if (is.null(ranking)) ranking <- if (mode == "lrt") "stat" else "log2FoldChange"

    args <- mget(names(formals()), envir = environment())
    .write_repro_script(
        "run_deseq2", data_arg = "dds", data_var = "dds",
        data_comment = "# The DESeqDataSet passed to run_deseq2(), before QC",
        args = args, call = call, save_dir = save_dir,
        run_comment = if (mode == "lrt") {
            "Run QC, PCA, DESeq2 LRT, degPatterns clustering and the enabled modules, then render the HTML report"
        } else {
            "Run QC, PCA, DESeq2, GSEA and the enabled modules, then render the HTML report"
        })
    engine <- if (mode == "lrt") .run_lrt_engine else .run_wald_engine
    engine_args <- args[setdiff(names(args), c("mode", only[[other]]))]
    do.call(engine, c(engine_args, list(call = call)), quote = TRUE)
}

#' The Wald-style pipeline shared by [run_deseq2()] and [run_limma()]
#'
#' `engine = "limma"` expects the carrier from [expr_to_dds()]: the DE stage
#' becomes limma-trend, and every VST / rlog step uses log2 expression
#' instead (see R/limma.R). Everything else -- stages,
#' files, report -- is identical.
#' @param engine `"deseq2"` or `"limma"`.
#' @param call The user's call, for the run log.
#' @param engine_params Extra settings to list in the report's argument table.
#' @noRd
.run_wald_engine <- function(
    dds,
    var,
    design = NULL,
    comparisons = NULL,
    include_interactions = FALSE,
    stratified = FALSE,
    batch = NULL,
    pair = NULL,
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
    save_dir = getwd(),
    engine = "deseq2",
    call = NULL,
    engine_params = NULL) {

    save_dir <- validate_save_dir(save_dir)
    opts <- .as_opts(opts)
    alpha <- opts$alpha
    lfc_threshold <- opts$lfc_threshold
    # Blocking covariates: in every model, and removed from the sample-level
    # matrix. `pair` sits next to `var` so the design reads `~ batch + pair + var`.
    covars <- c(batch, pair)

    design_resolved <- resolve_design(
        var, batch = covars, design = design,
        include_interactions = include_interactions
    )
    params <- list(
        var = var,
        design = paste(deparse(design_resolved), collapse = ""),
        comparisons = comparisons,
        include_interactions = include_interactions,
        stratified = stratified,
        batch = batch, pair = pair,
        export_counts = export_counts, filter_method = filter_method,
        cpm_percentile = cpm_percentile, min_count = min_count,
        min_samples_per_group = min_samples_per_group,
        remove_mt = remove_mt, remove_xy = remove_xy, remove_vdj = remove_vdj,
        ranking = ranking, collections = collections, run_clusterprofiler = run_clusterprofiler,
        run_gsva = run_gsva,
        run_wgcna = run_wgcna,
        run_decoupler_tf = run_decoupler_tf,
        run_decoupler_ccc = run_decoupler_ccc,
        run_viper = run_viper,
        run_cspa = run_cspa,
        render_html = render_html,
        parallel_render = parallel_render,
        render_workers = render_workers,
        use_targets = use_targets
    )
    params <- c(params, .opts_diff(opts))
    is_limma <- identical(engine, "limma")
    if (is_limma) {
        params <- c(list(analysis_mode = "limma"), params, engine_params)
    }

    log_output(call, save_dir, params = params, expr = quote({

        message("Running DEGreports ", if (is_limma) "limma-trend" else "Wald",
                " pipeline v", utils::packageVersion("DEGreports"), "...")

        .warn_if_no_pandoc(render_html)

        dds <- validate_dds(dds)
        dds <- validate_var(dds, var)
        dds <- validate_batch(dds, batch)
        dds <- validate_pair(dds, pair, var = var, batch = batch)
        dds <- validate_design(dds, design_resolved)
        contrast_specs <- validate_comparisons(
            dds, var, comparisons,
            include_interactions = include_interactions,
            stratified = stratified
        )
        n_contrasts <- nrow(contrast_specs)

        # QC / PCA / WGCNA colour by a single factor: the first `var`, or
        # an interaction paste when two+ factors are supplied.
        viz_var <- var[1]
        if (length(var) >= 2L) {
            viz_var <- paste(var, collapse = "_x_")
            SummarizedExperiment::colData(dds)[[viz_var]] <- do.call(
                interaction,
                c(
                    lapply(var, function(v) SummarizedExperiment::colData(dds)[[v]]),
                    list(sep = "_", drop = TRUE)
                )
            )
            if (is.null(shape_var) && is.null(pair)) shape_var <- var[1]
        }
        # Paired design: shape the PCA points by unit (donor / patient).
        if (is.null(shape_var) && !is.null(pair)) shape_var <- pair

        # The run as data: one descriptor per stage, in order. The driver
        # owns the progress counter and the failure policy, so a module can be
        # added, reordered or made optional without touching either.
        stages <- list(
            .stage("dds", "Quality control", optional = FALSE, run = function(st) {
                run_qc(
                    dds = st$dds, var = viz_var, filter_method = filter_method,
                    cpm_percentile = cpm_percentile, min_count = min_count,
                    min_samples_per_group = min_samples_per_group,
                    remove_mt = remove_mt, mt_pattern = mt_pattern,
                    remove_xy = remove_xy, gene_chr = gene_chr,
                    remove_vdj = remove_vdj,
                    save_dir = file.path(save_dir, "qc"))
            }),
            # Covariates for every later model and the batch-corrected matrix:
            # batch / pair, plus imputed sex when it is not confounded.
            .stage("model", "Model covariates", optional = FALSE, silent = TRUE,
                   run = function(st) {
                .model_covariates(st$dds, var = var, covars = covars,
                                  include_interactions = include_interactions)
            }),
            .stage(NULL, "Exporting count tables (raw / VST / rlog / CPM / batch-corrected)",
                enabled = isTRUE(export_counts), run = function(st) {
                run_counts_export(
                    dds = st$dds, batch = st$model$covars, var = viz_var,
                    save_dir = file.path(save_dir, "qc", "counts"))
            }),
            .stage("distance", "Sample similarity / PCA", optional = FALSE,
                   run = function(st) {
                run_distance_pca(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    pals = pals, shape_var = shape_var,
                    n_scree_pcs = n_scree_pcs, n_loading_pcs = n_loading_pcs,
                    n_top_genes = n_top_genes,
                    assoc_vars = unique(c(var, st$model$covars, .sex_assoc_col(st$dds))),
                    colour_vars = .pca_colour_vars(st$dds, viz_var, pair, batch),
                    save_dir = file.path(save_dir, "qc", "distance"))
            }),
            .stage("diffexp", paste0(if (is_limma) "limma-trend" else "DESeq2",
                                     " analysis (", n_contrasts, " contrast(s))"),
                   optional = FALSE, run = function(st) {
                run_diffexp_pipeline(
                    dds = st$dds, var = var, comparisons = contrast_specs, batch = st$model$covars,
                    include_interactions = include_interactions, stratified = stratified,
                    shrinkage = opts$shrinkage,
                    alpha = alpha, lfc_threshold = lfc_threshold,
                    ranking = ranking, top_label_n = top_label_n,
                    save_dir = file.path(save_dir, "deseq2"),
                    engine = engine)
            }),
            .stage("cspa", "CSPA surface-marker heatmaps", enabled = isTRUE(run_cspa),
                   run = function(st) {
                run_cspa_marker_pipeline(
                    dds = st$dds, diffexp_out = st$diffexp, var = var, batch = st$model$covars,
                    n_markers = opts$cspa_n_markers,
                    save_dir = file.path(save_dir, "deseq2"))
            }),
            .stage("gsea", "clusterProfiler analysis", run = function(st) {
                run_gsea_pipeline(
                    diffexp_results = st$diffexp$results, collections = collections,
                    ranking = ranking, top_n_plot = top_n_plot,
                    run_clusterprofiler = run_clusterprofiler,
                    save_dir = file.path(save_dir, "clusterprofiler"))
            }),
            .stage("gsea_overlap", "GSEA leading-edge overlap", run = function(st) {
                run_gsea_overlap_pipeline(
                    gsea_out = st$gsea, alpha = alpha,
                    save_dir = file.path(save_dir, "clusterprofiler"))
            }),
            .stage("gsva", "GSVA (sample-level gene-set scores)",
                   enabled = isTRUE(run_gsva), run = function(st) {
                run_gsva_pipeline(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    collections = unique(c("HALLMARK", collections)),
                    n_heatmap = opts$gsva_n_heatmap,
                    save_dir = file.path(save_dir, "gsva"))
            }),
            .stage("wgcna", "Weighted gene co-expression network analysis (WGCNA)",
                   enabled = isTRUE(run_wgcna), run = function(st) {
                run_wgcna_pipeline(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    min_module_size = opts$wgcna_min_module_size,
                    ora_collection = collections,
                    save_dir = file.path(save_dir, "wgcna"))
            }),
            .stage("decoupler_tf", "DecoupleR (CollecTRI + ULM)",
                   enabled = isTRUE(run_decoupler_tf), run = function(st) {
                run_decoupler_tf_pipeline(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    diffexp_out = st$diffexp, n_tfs = opts$decoupler_n_tfs,
                    save_dir = file.path(save_dir, "decoupler-tf"))
            }),
            .stage("decoupler_ccc", "DecoupleR-CCC (ligand\u2013receptor ULM)",
                   enabled = isTRUE(run_decoupler_ccc), run = function(st) {
                run_decoupler_ccc_pipeline(
                    diffexp_out = st$diffexp, n_lr = opts$decoupler_n_lr,
                    save_dir = file.path(save_dir, "decoupler-ccc"))
            }),
            .stage("viper", "VIPER (bcellViper)", enabled = isTRUE(run_viper),
                   run = function(st) {
                run_viper_pipeline(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    diffexp_out = st$diffexp, n_tfs = opts$viper_n_tfs,
                    save_dir = file.path(save_dir, "viper"))
            }),
            .stage("ppi", "PPI connectivity (STRINGdb)", run = function(st) {
                run_ppi_pipeline(
                    diffexp_out = st$diffexp, wgcna_out = st$wgcna,
                    top_n_hubs = opts$ppi_top_n_hubs,
                    alpha = alpha, lfc_threshold = lfc_threshold,
                    save_dir = file.path(save_dir, "deseq2"),
                    module_save_dir = file.path(save_dir, "ppi"),
                    wgcna_module_save_dir = file.path(save_dir, "wgcna"))
            }),
            .stage("multihit", "Cross-analysis gene hit ranking", run = function(st) {
                run_multihit_pipeline(
                    diffexp_out = st$diffexp, gsea_out = st$gsea,
                    viper_out = st$viper, decoupler_tf_out = st$decoupler_tf,
                    decoupler_ccc_out = st$decoupler_ccc, cspa_out = st$cspa,
                    alpha = alpha, lfc_threshold = lfc_threshold, top_n = 10,
                    save_dir = file.path(save_dir, "multihit"))
            })
        )

        st <- .run_stages(stages, list(dds = dds))
        dds <- st$dds
        params$imputed_sex <- st$model$sex_reason
        distance_out <- st$distance
        diffexp_out <- st$diffexp
        cspa_out <- st$cspa
        gsea_out <- st$gsea
        gsea_overlap_out <- st$gsea_overlap
        gsva_out <- st$gsva
        wgcna_out <- st$wgcna
        decoupler_tf_out <- st$decoupler_tf
        decoupler_ccc_out <- st$decoupler_ccc
        viper_out <- st$viper
        ppi_out <- st$ppi
        multihit_out <- st$multihit

        report_path <- NULL
        if (render_html) {
            message("Rendering HTML report...")
            ppi_for_report <- ppi_out
            render_fn <- if (isTRUE(use_targets)) render_report_targets else render_report
            report_path <- render_fn(
                dds = dds, var = var, distance_out = distance_out,
                diffexp_out = diffexp_out, gsea_out = gsea_out,
                gsea_overlap_out = gsea_overlap_out,
                multihit_out = multihit_out,
                wgcna_out = wgcna_out, decoupler_tf_out = decoupler_tf_out,
                decoupler_ccc_out = decoupler_ccc_out,
                ppi_out = ppi_for_report, viper_out = viper_out,
                cspa_out = cspa_out,
                gsva_out = gsva_out,
                alpha = alpha, lfc_threshold = lfc_threshold,
                run_args = params,
                save_dir = save_dir,
                parallel_render = parallel_render,
                render_workers = render_workers)
        }

        cleanup_ppi_dir(save_dir)

        dds_rds <- file.path(save_dir, "dds.rds")
        message("Saving QC-filtered DESeqDataSet to: ", dds_rds)
        saveRDS(dds, dds_rds)

        message("Pipeline complete. Outputs written to: ", save_dir)
        invisible(list(
            dds = dds, distance = distance_out, diffexp = diffexp_out,
            gsea = gsea_out, gsea_overlap = gsea_overlap_out,
            multihit = multihit_out,
            wgcna = wgcna_out, decoupler_tf = decoupler_tf_out,
            decoupler_ccc = decoupler_ccc_out, viper = viper_out,
            ppi = ppi_out, cspa = cspa_out, gsva = gsva_out,
            report_path = report_path,
            dds_path = dds_rds
        ))
    }))
}

#' The LRT (time-course) pipeline behind `run_deseq2(mode = "lrt")`
#'
#' Arguments as documented for [run_deseq2()].
#' @param call The user's call, for the run log.
#' @noRd
.run_lrt_engine <- function(
    dds,
    var = NULL,
    time,
    include_interactions = FALSE,
    batch = NULL,
    pair = NULL,
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
    ranking = "stat",
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
    degpatterns_minc = 15,
    render_html = TRUE,
    parallel_render = FALSE,
    render_workers = NULL,
    use_targets = FALSE,
    opts = degreports_opts(),
    save_dir = getwd(),
    call = NULL) {

    save_dir <- validate_save_dir(save_dir)
    opts <- .as_opts(opts)
    covars <- c(batch, pair)
    alpha <- opts$alpha
    lfc_threshold <- opts$lfc_threshold

    if (is.null(var)) var <- character(0)
    if (!is.character(time) || length(time) != 1L || is.na(time) || !nzchar(time)) {
        stop("`time` must be a single non-empty character string.", call. = FALSE)
    }

    params <- list(
        analysis_mode = "lrt",
        var = var,
        time = time,
        include_interactions = include_interactions,
        batch = covars,
        export_counts = export_counts, filter_method = filter_method,
        cpm_percentile = cpm_percentile, min_count = min_count,
        min_samples_per_group = min_samples_per_group,
        remove_mt = remove_mt, remove_xy = remove_xy, remove_vdj = remove_vdj,
        ranking = ranking, collections = collections,
        run_clusterprofiler = run_clusterprofiler,
        run_gsva = run_gsva,
        run_wgcna = run_wgcna,
        run_decoupler_tf = run_decoupler_tf,
        run_decoupler_ccc = run_decoupler_ccc,
        run_viper = run_viper,
        run_cspa = run_cspa,
        degpatterns_minc = degpatterns_minc,
        render_html = render_html,
        parallel_render = parallel_render,
        render_workers = render_workers,
        use_targets = use_targets
    )
    params <- c(params, .opts_diff(opts))

    log_output(call, save_dir, params = params, expr = quote({

        message("Running DEGreports LRT pipeline v",
                utils::packageVersion("DEGreports"), "...")

        .warn_if_no_pandoc(render_html)

        dds <- validate_dds(dds)
        dds <- .lrt_factorize_cols(dds, time, "time")
        dds <- .lrt_factorize_cols(dds, var, "var")
        dds <- .lrt_factorize_cols(dds, batch, "batch")
        dds <- validate_pair(dds, pair, var = var, batch = batch, time = time)

        lrt_specs <- expand_lrt_specs(
            var = var, time = time, batch = covars,
            include_interactions = include_interactions
        )
        n_tests <- nrow(lrt_specs)
        # Stamp the largest full model (always includes `batch` / `pair` when
        # set) onto dds so QC / report / any design-aware downstream step see
        # `~ batch + pair + time` (or `... + var + time`), not the input
        # formula the caller built the DESeqDataSet with.
        dds <- validate_design(
            dds, stats::as.formula(lrt_specs$full[[1L]], env = baseenv())
        )
        params$design <- lrt_specs$full[[1L]]

        # Visualization factor for QC / PCA / WGCNA: time alone, or
        # var × time paste when grouping factors are present.
        if (length(var)) {
            viz_cols <- c(var, time)
            viz_var <- paste(viz_cols, collapse = "_x_")
            SummarizedExperiment::colData(dds)[[viz_var]] <- interaction(
                as.data.frame(SummarizedExperiment::colData(dds))[viz_cols],
                sep = "_", drop = TRUE
            )
            if (is.null(shape_var) && is.null(pair)) shape_var <- var[[1L]]
        } else {
            viz_var <- time
        }
        if (is.null(shape_var) && !is.null(pair)) shape_var <- pair

        # The LRT run as data; see .stage() / .run_stages(). Same driver,
        # same failure policy and progress counter as the Wald pipeline --
        # only the stages differ.
        stages <- list(
            .stage("dds", "Quality control", optional = FALSE, run = function(st) {
                run_qc(
                    dds = st$dds, var = viz_var, filter_method = filter_method,
                    cpm_percentile = cpm_percentile, min_count = min_count,
                    min_samples_per_group = min_samples_per_group,
                    remove_mt = remove_mt, mt_pattern = mt_pattern,
                    remove_xy = remove_xy, gene_chr = gene_chr,
                    remove_vdj = remove_vdj,
                    save_dir = file.path(save_dir, "qc"))
            }),
            # Covariates for every later model and the batch-corrected matrix:
            # batch / pair, plus imputed sex when it is not confounded.
            .stage("model", "Model covariates", optional = FALSE, silent = TRUE,
                   run = function(st) {
                .model_covariates(st$dds, var = var, covars = covars,
                                  include_interactions = include_interactions, time = time)
            }),
            .stage(NULL, "Exporting count tables", enabled = isTRUE(export_counts),
                   run = function(st) {
                run_counts_export(
                    dds = st$dds, batch = st$model$covars, var = viz_var,
                    save_dir = file.path(save_dir, "qc", "counts"))
            }),
            .stage("distance", "Sample similarity / PCA", optional = FALSE,
                   run = function(st) {
                run_distance_pca(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    pals = pals, shape_var = shape_var,
                    n_scree_pcs = n_scree_pcs, n_loading_pcs = n_loading_pcs,
                    n_top_genes = n_top_genes,
                    assoc_vars = unique(c(var, time, st$model$covars, .sex_assoc_col(st$dds))),
                    colour_vars = .pca_colour_vars(st$dds, viz_var, pair, batch),
                    save_dir = file.path(save_dir, "qc", "distance"))
            }),
            .stage("diffexp", paste0("DESeq2 LRT (", n_tests, " test(s))"),
                   optional = FALSE, run = function(st) {
                run_diffexp_lrt_pipeline(
                    dds = st$dds, var = var, time = time, batch = st$model$covars,
                    include_interactions = include_interactions,
                    alpha = alpha,
                    save_dir = file.path(save_dir, "deseq2"))
            }),
            .stage("degpatterns", "degPatterns temporal clustering", optional = FALSE,
                   run = function(st) {
                run_degpatterns_pipeline(
                    dds = st$dds, lrt_out = st$diffexp, time = time, var = var,
                    alpha = alpha, minc = degpatterns_minc,
                    save_dir = file.path(save_dir, "deseq2"))
            }),
            .stage("cspa", "CSPA heatmaps by temporal cluster",
                   enabled = isTRUE(run_cspa), run = function(st) {
                hm <- run_lrt_cspa_heatmap(
                    dds = st$dds, degpatterns_out = st$degpatterns,
                    lrt_out = st$diffexp, n_per_cluster = opts$cspa_n_markers_lrt,
                    annotation_var = unique(c(time, var)), batch = st$model$covars,
                    save_dir = file.path(save_dir, "deseq2"))
                list(
                    source = "lrt_clusters",
                    org = detect_organism(rownames(st$dds)),
                    by_lrt = hm,
                    heatmaps = lapply(hm, function(x) x$heatmap),
                    n_markers = opts$cspa_n_markers_lrt,
                    expr_correction = attr(hm, "expr_correction")
                )
            }),
            .stage("cluster_ora", "Cluster ORA (clusterProfiler page)",
                   run = function(st) {
                run_cluster_ora_pipeline(
                    degpatterns_out = st$degpatterns, universe = rownames(st$dds),
                    collections = collections, top_n = top_n_plot, alpha = alpha,
                    save_dir = file.path(save_dir, "clusterprofiler"))
            }),
            .stage("gsva", "GSVA (sample-level gene-set scores)",
                   enabled = isTRUE(run_gsva), run = function(st) {
                run_gsva_pipeline(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    collections = unique(c("HALLMARK", collections)),
                    n_heatmap = opts$gsva_n_heatmap,
                    save_dir = file.path(save_dir, "gsva"))
            }),
            .stage("wgcna", "WGCNA", enabled = isTRUE(run_wgcna), run = function(st) {
                run_wgcna_pipeline(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    min_module_size = opts$wgcna_min_module_size,
                    ora_collection = collections,
                    save_dir = file.path(save_dir, "wgcna"))
            }),
            # One TF signature per degPatterns cluster, shared by the three
            # regulator modules below, rather than the whole-LRT gene set.
            .stage("tf_diffexp", "Per-cluster TF signatures", silent = TRUE,
                   enabled = any(c(run_decoupler_tf, run_decoupler_ccc, run_viper)),
                   run = function(st) {
                lrt_clusters_diffexp_for_tf(
                    degpatterns_out = st$degpatterns, lrt_out = st$diffexp)
            }),
            .stage("decoupler_tf", "DecoupleR (sample + per-cluster ULM)",
                   enabled = isTRUE(run_decoupler_tf), run = function(st) {
                if (is.null(st$tf_diffexp)) {
                    warning("Skipping decoupleR: no temporal clusters for TF input.",
                            call. = FALSE)
                    return(NULL)
                }
                out <- run_decoupler_tf_pipeline(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    diffexp_out = st$tf_diffexp, n_tfs = opts$decoupler_n_tfs,
                    save_dir = file.path(save_dir, "decoupler-tf"))
                if (!is.null(out)) out$tf_level <- "lrt_cluster"
                out
            }),
            .stage("decoupler_ccc", "DecoupleR-CCC (ligand\u2013receptor ULM)",
                   enabled = isTRUE(run_decoupler_ccc), run = function(st) {
                run_decoupler_ccc_pipeline(
                    diffexp_out = st$diffexp, cluster_diffexp = st$tf_diffexp,
                    n_lr = opts$decoupler_n_lr,
                    save_dir = file.path(save_dir, "decoupler-ccc"))
            }),
            .stage("viper", "VIPER (sample + per-cluster msVIPER)",
                   enabled = isTRUE(run_viper), run = function(st) {
                if (is.null(st$tf_diffexp)) {
                    warning("Skipping VIPER: no temporal clusters for TF input.",
                            call. = FALSE)
                    return(NULL)
                }
                out <- run_viper_pipeline(
                    dds = st$dds, var = viz_var, batch = st$model$covars,
                    diffexp_out = st$tf_diffexp, n_tfs = opts$viper_n_tfs,
                    save_dir = file.path(save_dir, "viper"))
                if (!is.null(out)) out$tf_level <- "lrt_cluster"
                out
            }),
            .stage("ppi", "PPI connectivity (clusters + WGCNA)", run = function(st) {
                clusters_wgcna <- tryCatch(
                    lrt_clusters_as_wgcna_like(st$degpatterns),
                    error = function(e) NULL
                )
                # degPatterns clusters and WGCNA modules are merged into one
                # module_df so PPI covers both; the WGCNA page only displays
                # real WGCNA keys. Keep the columns PPI seeds on.
                .ppi_module_cols <- function(df) {
                    gene_col <- if ("gene" %in% colnames(df)) {
                        "gene"
                    } else {
                        setdiff(colnames(df), c("module", "module_label",
                                                "module_color", "kME",
                                                "cluster", "lrt"))[1]
                    }
                    keep <- c(gene_col, "module", intersect("kME", colnames(df)))
                    out <- df[, keep, drop = FALSE]
                    colnames(out)[1] <- "gene"
                    if (!"kME" %in% colnames(out)) out$kME <- NA_real_
                    out
                }
                mod_parts <- list()
                if (!is.null(clusters_wgcna) && nrow(clusters_wgcna$module_df)) {
                    mod_parts[[length(mod_parts) + 1L]] <-
                        .ppi_module_cols(clusters_wgcna$module_df)
                }
                if (!is.null(st$wgcna) && !is.null(st$wgcna$module_df) &&
                    nrow(st$wgcna$module_df)) {
                    mod_parts[[length(mod_parts) + 1L]] <-
                        .ppi_module_cols(st$wgcna$module_df)
                }
                ppi_modules <- NULL
                if (length(mod_parts)) {
                    merged <- do.call(rbind, mod_parts)
                    rownames(merged) <- NULL
                    ppi_modules <- list(module_df = merged)
                }
                run_ppi_pipeline(
                    diffexp_out = st$diffexp, wgcna_out = ppi_modules,
                    top_n_hubs = opts$ppi_top_n_hubs, max_plot_nodes = 50L,
                    alpha = alpha, lfc_threshold = 0,
                    save_dir = file.path(save_dir, "deseq2"),
                    module_save_dir = file.path(save_dir, "ppi"),
                    wgcna_module_save_dir = file.path(save_dir, "wgcna"))
            }),
            .stage("multihit", "Cross-analysis gene hit ranking (per cluster)",
                   run = function(st) {
                run_multihit_lrt_pipeline(
                    degpatterns_out = st$degpatterns, diffexp_out = st$diffexp,
                    wgcna_out = st$wgcna, viper_out = st$viper,
                    decoupler_tf_out = st$decoupler_tf,
                    decoupler_ccc_out = st$decoupler_ccc,
                    alpha = alpha, top_n = 10,
                    save_dir = file.path(save_dir, "multihit"))
            })
        )

        st <- .run_stages(stages, list(dds = dds))
        dds <- st$dds
        params$imputed_sex <- st$model$sex_reason
        distance_out <- st$distance
        diffexp_out <- st$diffexp
        degpatterns_out <- st$degpatterns
        cspa_out <- st$cspa
        cluster_ora_out <- st$cluster_ora
        # Ranked GSEA is not used in LRT mode (no signed LFC ranking).
        gsea_out <- NULL
        gsea_overlap_out <- NULL
        gsva_out <- st$gsva
        wgcna_out <- st$wgcna
        decoupler_tf_out <- st$decoupler_tf
        decoupler_ccc_out <- st$decoupler_ccc
        viper_out <- st$viper
        ppi_out <- st$ppi
        multihit_out <- st$multihit

        report_path <- NULL
        if (render_html) {
            message("Rendering HTML report...")
            render_fn <- if (isTRUE(use_targets)) render_report_targets else render_report
            report_path <- render_fn(
                dds = dds, var = if (length(var)) var else time,
                distance_out = distance_out,
                diffexp_out = diffexp_out, gsea_out = gsea_out,
                gsea_overlap_out = gsea_overlap_out,
                multihit_out = multihit_out,
                wgcna_out = wgcna_out, decoupler_tf_out = decoupler_tf_out,
                decoupler_ccc_out = decoupler_ccc_out,
                ppi_out = ppi_out, viper_out = viper_out,
                cspa_out = cspa_out,
                gsva_out = gsva_out,
                degpatterns_out = degpatterns_out,
                cluster_ora_out = cluster_ora_out,
                alpha = alpha, lfc_threshold = lfc_threshold,
                run_args = params,
                save_dir = save_dir,
                parallel_render = parallel_render,
                render_workers = render_workers)
        }

        cleanup_ppi_dir(save_dir)

        dds_rds <- file.path(save_dir, "dds.rds")
        message("Saving QC-filtered DESeqDataSet to: ", dds_rds)
        saveRDS(dds, dds_rds)

        message("Pipeline complete. Outputs written to: ", save_dir)
        invisible(list(
            dds = dds, distance = distance_out, diffexp = diffexp_out,
            degpatterns = degpatterns_out, cluster_ora = cluster_ora_out,
            gsea = gsea_out, gsea_overlap = gsea_overlap_out,
            multihit = multihit_out,
            wgcna = wgcna_out, decoupler_tf = decoupler_tf_out,
            decoupler_ccc = decoupler_ccc_out, viper = viper_out,
            ppi = ppi_out, cspa = cspa_out, gsva = gsva_out,
            report_path = report_path,
            dds_path = dds_rds
        ))
    }))
}
