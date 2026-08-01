#' Run the complete RNA-seq DESeq2 pipeline in one command
#'
#' Chains quality control ([run_qc()]), sample distance/PCA
#' ([run_distance_pca()]), differential expression across all requested
#' comparisons ([run_diffexp_pipeline()]), gene set enrichment analysis
#' ([run_gsea_pipeline()]), and an HTML report ([render_report()]) into a
#' single call. All arguments that commonly differ between analyses --
#' filtering thresholds, batch correction, LFC shrinkage method,
#' significance thresholds, the DEG ranking metric, and MSigDB collections
#' -- are exposed directly; every other stage can also be run standalone
#' with the same arguments if you need to inspect or rerun a single step.
#'
#' @param dds A `DESeqDataSet` object containing raw counts.
#' @param var Column name(s) in `colData(dds)` used for grouping. Length 1
#'   = single-factor mode; length 2 = two-factor mode with interaction
#'   contrasts. Comparisons / the design formula are auto-built from these
#'   columns unless overridden.
#' @param design Design formula (or character string) for the DESeq2 model,
#'   or `NULL` (default) to auto-build via [resolve_design()] from `var`
#'   and `batch` (`~ var`, `~ batch + var`, `~ A * B`, or
#'   `~ batch + A * B`).
#' @param comparisons Character vector of comparisons / contrast labels to
#'   run (`"A_vs_B"` / `"A_vs_B+C"`, or two-factor labels like
#'   `"BM__IgG_vs_IgM"`), or `NULL` (default) to auto-generate the full
#'   grid via [expand_contrast_specs()].
#' @param batch Optional covariate column name(s) in `colData(dds)`. Added
#'   to the auto design when `design = NULL`, passed through to within-
#'   stratum subset fits in two-factor mode, and used for PCA visualization
#'   correction. Default `NULL`.
#' @param batch_method Batch-correction method for PCA visualization only:
#'   `"limma"` (default), `"combat"`, or `"none"`.
#' @param export_counts Logical. If `TRUE` (default), also writes raw,
#'   VST, rlog, and CPM count tables via [run_counts_export()]. If `batch`
#'   is set, a ComBat_seq-corrected raw-count table is written as well
#'   (protecting `var` as the biological signal of interest). These tables
#'   are independent of the DE model and provided for general downstream
#'   use.
#' @param filter_method Low-expression gene-filtering strategy: `"cpm"`
#'   (default, percentile-based), `"count"`, or `NULL` (no expression
#'   filtering). See [run_qc()].
#' @param cpm_percentile Used when `filter_method = "cpm"`. See [run_qc()].
#' @param min_count,min_samples_per_group Used when `filter_method =
#'   "count"`. See [run_qc()].
#' @param remove_mt,remove_xy,mt_pattern,gene_chr,remove_vdj See [run_qc()].
#' @param shrinkage LFC shrinkage method, one of `"ashr"` (default),
#'   `"apeglm"`, `"normal"`, `"none"`. See [run_diffexp()].
#' @param alpha Significance threshold (adjusted p-value). Default `0.05`.
#' @param lfc_threshold Log2 fold-change threshold for significance
#'   flagging / up-down colouring. Default `0.5`.
#' @param ranking DEG/GSEA ranking metric: `"pxfc"` (default), `"padj"`,
#'   `"log2FoldChange"`, or `"stat"`.
#' @param collections MSigDB collections for GSEA: one or more of
#'   `"HALLMARK"`, `"GOBP"`, `"GOCC"`, `"GOMF"`, `"KEGG_MEDICUS"`,
#'   `"KEGG_LEGACY"`, `"REACTOME"`, `"PID"`. Default `c("GOBP", "GOCC",
#'   "GOMF", "KEGG_MEDICUS", "KEGG_LEGACY", "REACTOME", "PID")` --
#'   `"HALLMARK"` remains available but is not forced on by default. GSEA
#'   runs independently against every requested collection, for every
#'   comparison. Organism (human/mouse) is auto-detected from gene symbols
#'   -- there is no `org` argument.
#' @param pals Optional colour palette (character vector), one colour per
#'   level of the visualization grouping factor. Default `NULL`.
#' @param shape_var Optional column name for a second grouping variable
#'   mapped to point shape in the PCA plot. Default `NULL`.
#' @param n_scree_pcs,n_loading_pcs,n_top_genes See [run_distance_pca()].
#' @param top_label_n Number of top **up** and top **down** genes (by
#'   `rank_score`) labelled per volcano plot. Default `5`.
#' @param top_n_plot Number of top pathways per direction shown per GSEA
#'   barplot (regardless of significance). Default `10`.
#' @param run_clusterprofiler Logical. If `TRUE` (default `FALSE`), also
#'   runs a second, independent GSEA pass via `clusterProfiler::GSEA()`
#'   alongside the primary `fgsea`-based GSEA, exporting the raw
#'   `clusterProfiler` `gseaResult` objects. See [run_gsea_pipeline()].
#' @param run_wgcna Logical. If `TRUE` (default `FALSE`), also runs
#'   weighted gene co-expression network analysis via
#'   [run_wgcna_pipeline()] -- soft-threshold power selection, module
#'   detection, condition- and sample-level module-trait correlation, hub
#'   genes, and per-module over-representation analysis against
#'   `"GOBP"`, `"GOCC"`, `"GOMF"`, `"KEGG_MEDICUS"`, `"KEGG_LEGACY"`, and
#'   `"REACTOME"`.
#'   Requires the `WGCNA` package. This is a heavier, exploratory
#'   companion analysis independent of the `var`-based differential
#'   expression comparisons above; it uses every sample's full expression
#'   profile rather than testing one gene at a time against a design.
#' @param wgcna_power Soft-thresholding power for WGCNA. Default `NULL`
#'   auto-selects. Ignored if `run_wgcna = FALSE`. See
#'   [run_wgcna_pipeline()].
#' @param wgcna_min_module_size Minimum WGCNA module size. Default `30`.
#'   Ignored if `run_wgcna = FALSE`.
#' @param run_decoupler Logical. If `TRUE` (default `FALSE`), runs
#'   CollecTRI + ULM TF activity analysis via [run_decoupler_pipeline()]
#'   and adds a **TF activity** HTML page. Requires Suggests packages
#'   `decoupleR` and `OmnipathR`.
#' @param decoupler_n_tfs Number of top TFs shown in the sample heatmap
#'   and contrast barplots. Default `25`. Ignored if `run_decoupler = FALSE`.
#' @param render_html Logical. If `TRUE` (default), renders an HTML summary
#'   report via [render_report()] after all stages complete.
#' @param save_dir Directory where all outputs will be written. Default
#'   `getwd()`.
#' @return Invisibly, a list with elements `dds` (the QC-filtered object),
#'   `distance` (output of [run_distance_pca()]), `diffexp` (output of
#'   [run_diffexp_pipeline()]), `gsea` (output of [run_gsea_pipeline()]),
#'   `wgcna` (output of [run_wgcna_pipeline()], or `NULL` if `run_wgcna =
#'   FALSE`), and `report_path` (a named character vector of the paths to
#'   the rendered HTML report pages -- `"index"`, `"qc"`, `"diffexp"`,
#'   `"pathway_enrichment"`, and (if enabled) `"wgcna"` / `"tf_activity"`;
#'   see [render_report()] -- or `NULL` if `render_html = FALSE`). Module
#'   artifacts and their HTML pages are written under per-module folders
#'   (`qc/`, `distance/`, `diffexp/`, `pathway/`, `wgcna/`, `tf_activity/`)
#'   with `index.html` at the root of `save_dir`.
#' @examples
#' \dontrun{
#' out <- run_deseq2_pipeline(
#'   dds, var = "condition",
#'   save_dir = "results")
#'
#' out2 <- run_deseq2_pipeline(
#'   dds, var = c("Tissue", "Isotype"),
#'   save_dir = "results_2factor")
#' }
#' @export
run_deseq2_pipeline <- function(
    dds,
    var,
    design = NULL,
    comparisons = NULL,
    batch = NULL,
    batch_method = "limma",
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
    shrinkage = "ashr",
    alpha = 0.05,
    lfc_threshold = 0.5,
    ranking = "pxfc",
    collections = c("GOBP", "GOCC", "GOMF", "KEGG_MEDICUS", "KEGG_LEGACY", "REACTOME", "PID"),
    pals = NULL,
    shape_var = NULL,
    n_scree_pcs = 10,
    n_loading_pcs = 6,
    n_top_genes = 10,
    top_label_n = 5,
    top_n_plot = 10,
    run_clusterprofiler = FALSE,
    run_wgcna = FALSE,
    wgcna_power = NULL,
    wgcna_min_module_size = 30,
    run_decoupler = FALSE,
    decoupler_n_tfs = 25,
    render_html = TRUE,
    save_dir = getwd()) {

    call <- match.call()
    save_dir <- validate_save_dir(save_dir)

    design_resolved <- resolve_design(var, batch = batch, design = design)
    params <- list(
        var = var,
        design = paste(deparse(design_resolved), collapse = ""),
        comparisons = comparisons, batch = batch, batch_method = batch_method,
        export_counts = export_counts, filter_method = filter_method,
        cpm_percentile = cpm_percentile, min_count = min_count,
        min_samples_per_group = min_samples_per_group,
        remove_mt = remove_mt, remove_xy = remove_xy, remove_vdj = remove_vdj,
        shrinkage = shrinkage, alpha = alpha, lfc_threshold = lfc_threshold,
        ranking = ranking, collections = collections, run_clusterprofiler = run_clusterprofiler,
        run_wgcna = run_wgcna, wgcna_power = wgcna_power,
        wgcna_min_module_size = wgcna_min_module_size,
        run_decoupler = run_decoupler, decoupler_n_tfs = decoupler_n_tfs,
        render_html = render_html
    )

    log_output(call, save_dir, params = params, expr = quote({

        message("Running deseq2reports pipeline v", utils::packageVersion("deseq2reports"), "...")

        dds <- validate_dds(dds)
        dds <- validate_var(dds, var)
        dds <- validate_design(dds, design_resolved)
        contrast_specs <- validate_comparisons(dds, var, comparisons)
        n_contrasts <- if (is.data.frame(contrast_specs)) nrow(contrast_specs) else length(contrast_specs)

        # QC / PCA / WGCNA colour by a single factor: the first `var`, or
        # an interaction paste when two factors are supplied.
        viz_var <- var[1]
        if (length(var) == 2) {
            viz_var <- paste(var, collapse = "_x_")
            SummarizedExperiment::colData(dds)[[viz_var]] <- interaction(
                SummarizedExperiment::colData(dds)[[var[1]]],
                SummarizedExperiment::colData(dds)[[var[2]]],
                sep = "_", drop = TRUE
            )
            if (is.null(shape_var)) shape_var <- var[1]
        }

        n_steps <- 4L + as.integer(isTRUE(export_counts)) +
            as.integer(isTRUE(run_wgcna)) + as.integer(isTRUE(run_decoupler))
        step <- 0L
        step <- step + 1L
        message("[", step, "/", n_steps, "] Quality control...")
        dds <- run_qc(
            dds = dds, var = viz_var, filter_method = filter_method,
            cpm_percentile = cpm_percentile, min_count = min_count,
            min_samples_per_group = min_samples_per_group,
            remove_mt = remove_mt, mt_pattern = mt_pattern,
            remove_xy = remove_xy, gene_chr = gene_chr,
            remove_vdj = remove_vdj,
            save_dir = file.path(save_dir, "qc"))

        if (export_counts) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] Exporting count tables (raw / VST / rlog / CPM", if (!is.null(batch)) " / ComBat_seq" else "", ")...")
            run_counts_export(
                dds = dds, batch = batch, var = viz_var,
                save_dir = file.path(save_dir, "counts"))
        }

        step <- step + 1L
        message("[", step, "/", n_steps, "] Sample similarity / PCA...")
        distance_out <- run_distance_pca(
            dds = dds, var = viz_var, batch = batch, batch_method = batch_method,
            pals = pals, shape_var = shape_var,
            n_scree_pcs = n_scree_pcs, n_loading_pcs = n_loading_pcs, n_top_genes = n_top_genes,
            save_dir = file.path(save_dir, "distance"))

        step <- step + 1L
        message("[", step, "/", n_steps, "] Differential expression (", n_contrasts, " contrast(s))...")
        diffexp_out <- run_diffexp_pipeline(
            dds = dds, var = var, comparisons = comparisons, batch = batch,
            shrinkage = shrinkage, alpha = alpha, lfc_threshold = lfc_threshold,
            ranking = ranking, top_label_n = top_label_n,
            save_dir = file.path(save_dir, "diffexp"))

        step <- step + 1L
        message("[", step, "/", n_steps, "] Gene set enrichment analysis...")
        gsea_out <- run_gsea_pipeline(
            diffexp_results = diffexp_out$results, collections = collections,
            ranking = ranking, top_n_plot = top_n_plot,
            run_clusterprofiler = run_clusterprofiler,
            save_dir = file.path(save_dir, "pathway"))

        wgcna_out <- NULL
        if (run_wgcna) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] Weighted gene co-expression network analysis (WGCNA)...")
            wgcna_out <- run_wgcna_pipeline(
                dds = dds, var = viz_var, power = wgcna_power,
                min_module_size = wgcna_min_module_size,
                ora_collection = collections,
                save_dir = file.path(save_dir, "wgcna"))
        }

        decoupler_out <- NULL
        if (run_decoupler) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] TF activity (CollecTRI + ULM)...")
            decoupler_out <- run_decoupler_pipeline(
                dds = dds, var = viz_var, diffexp_out = diffexp_out,
                n_tfs = decoupler_n_tfs,
                save_dir = file.path(save_dir, "tf_activity"))
        }

        report_path <- NULL
        if (render_html) {
            message("Rendering HTML report...")
            report_path <- render_report(
                dds = dds, var = var, distance_out = distance_out,
                diffexp_out = diffexp_out, gsea_out = gsea_out,
                wgcna_out = wgcna_out, decoupler_out = decoupler_out,
                alpha = alpha, lfc_threshold = lfc_threshold,
                run_args = params,
                save_dir = save_dir)
        }

        message("Pipeline complete. Outputs written to: ", save_dir)
        invisible(list(dds = dds, distance = distance_out, diffexp = diffexp_out,
                       gsea = gsea_out, wgcna = wgcna_out,
                       decoupler = decoupler_out, report_path = report_path))
    }))
}
