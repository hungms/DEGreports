#' Run the complete RNA-seq DESeq2 Wald pipeline in one command
#'
#' Chains quality control ([run_qc()]), sample distance/PCA
#' ([run_distance_pca()]), differential expression across all requested
#' comparisons ([run_diffexp_pipeline()]), gene set enrichment analysis
#' ([run_gsea_pipeline()]), and an HTML report ([render_report()]) into a
#' single call. For time-course likelihood-ratio testing see
#' [run_deseq2_lrt()]. The previous name [run_deseq2_pipeline()] is retained
#' as a deprecated alias.
#'
#' @param dds A `DESeqDataSet` object containing raw counts.
#' @param var Column name(s) in `colData(dds)` used for grouping. Length 1
#'   = single-factor mode; length 2 = two-factor mode (within-factor contrasts;
#'   interaction design + contrasts when `include_interactions = TRUE`).
#'   Comparisons / the design formula are auto-built from these columns
#'   unless overridden.
#' @param design Design formula (or character string) for the DESeq2 model,
#'   or `NULL` (default) to auto-build via [resolve_design()] from `var`,
#'   `batch`, and `include_interactions` (`~ var`, `~ batch + var`,
#'   additive `~ A + B`, or `~ A * B` when interactions are requested).
#' @param comparisons Character vector of comparisons / contrast labels to
#'   run (`"A_vs_B"` / `"A_vs_B+C"`, or two-factor labels like
#'   `"BM__IgG_vs_IgM"`), or `NULL` (default) to auto-generate the full
#'   grid via [expand_contrast_specs()].
#' @param include_interactions Logical. When `var` has length 2 and
#'   `design = NULL`, use `~ A * B` and expand interaction contrasts
#'   (reference level × non-ref and pairwise non-ref differences, e.g.
#'   mutant vs mutant × cell type). Default `FALSE` (additive design; no
#'   interaction section). Ignored for single-factor `var`. When `design`
#'   is supplied explicitly, only the contrast grid is affected.
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
#' @param ranking DEG/GSEA ranking metric: `"log2FoldChange"` (default), `"padj"`,
#'   `"pxfc"`, or `"stat"`.
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
#'   and adds a **decoupleR** HTML page. Requires Suggests packages
#'   `decoupleR` and `OmnipathR`.
#' @param decoupler_n_tfs Number of top TFs shown in the sample heatmap
#'   and contrast barplots. Default `25`. Ignored if `run_decoupler = FALSE`.
#' @param run_viper Logical. If `TRUE` (default `FALSE`), runs
#'   [VIPER](https://califano.c2b2.columbia.edu/viper) / msVIPER with the
#'   human B-cell regulon from [bcellViper](https://bioconductor.org/packages/bcellViper/)
#'   via [run_viper_pipeline()] and adds a **VIPER** HTML page. Requires
#'   Suggests packages `viper`, `bcellViper`, and (for mouse) `babelgene`.
#' @param viper_n_tfs Number of top master regulators shown in each
#'   contrast barcode plot. Default `10`. Ignored if `run_viper = FALSE`.
#' @param run_ppi Logical. If `TRUE` (default), runs STRING PPI
#'   connectivity analysis via [run_ppi_pipeline()] for every DE contrast
#'   and (when enabled) every WGCNA module, embedding BioNERO-style
#'   community/hub network plots (plus native STRINGdb PNGs on disk) in
#'   the DE and WGCNA HTML pages. Requires Suggests packages `STRINGdb`
#'   and `igraph` and network access on first run to download STRING
#'   files. Failures are caught with a warning so the rest of the
#'   pipeline still finishes.
#' @param ppi_top_n_hubs Number of highest-degree genes highlighted per
#'   PPI plot. Default `10`. Ignored if `run_ppi = FALSE`.
#' @param run_cspa Logical. If `TRUE` (default), selects top Cell Surface
#'   Protein Atlas (CSPA) markers among DEGs per contrast and draws
#'   all-sample VST heatmaps (plus a combined heatmap) via
#'   [run_cspa_marker_pipeline()]. Requires Suggests package `readxl`.
#' @param cspa_n_markers Top CSPA genes per contrast. Default `25`.
#' @param render_html Logical. If `TRUE` (default), renders an HTML summary
#'   report via [render_report()] after all stages complete.
#' @param save_dir Directory where all outputs will be written. Default
#'   `getwd()`.
#' @return Invisibly, a list with elements `dds` (the QC-filtered object),
#'   `dds_path` (path to `dds.rds` under `save_dir`),
#'   `distance` (output of [run_distance_pca()]), `diffexp` (output of
#'   [run_diffexp_pipeline()]), `gsea` (output of [run_gsea_pipeline()]),
#'   `wgcna` (output of [run_wgcna_pipeline()], or `NULL` if `run_wgcna =
#'   `run_ppi = FALSE` / the stage failed), `cspa` (output of
#'   [run_cspa_marker_pipeline()], or `NULL` if `run_cspa = FALSE`), `viper`
#'   (output of [run_viper_pipeline()], or `NULL` if `run_viper = FALSE`), and
#'   `report_path` (a named character vector of the paths to the rendered
#'   HTML report pages -- `"index"`, `"qc"`, `"deseq2"`,
#'   `"clusterProfiler"`, and (if enabled) `"wgcna"` / `"decoupleR"` /
#'   `"viper"`; see [render_report()] -- or `NULL` if `render_html =
#'   FALSE`). Module artifacts and their HTML pages are written under
#'   per-module folders (`qc/`, `distance/`, `deseq2/`, `clusterprofiler/`,
#'   `wgcna/`, `decoupler/`, `viper/`, `ppi/`; CSPA under `deseq2/cspa/`,
#'   GSEA overlap under `clusterprofiler/<contrast>/gsea_overlap/`) with
#'   `index.html` and the QC-filtered `dds.rds` at the root of `save_dir`.
#' @examples
#' \dontrun{
#' out <- run_deseq2_wald(
#'   dds, var = "condition",
#'   save_dir = "results")
#'
#' out2 <- run_deseq2_wald(
#'   dds, var = c("Tissue", "Isotype"),
#'   save_dir = "results_2factor")
#' }
#' @export
run_deseq2_wald <- function(
    dds,
    var,
    design = NULL,
    comparisons = NULL,
    include_interactions = FALSE,
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
    ranking = "log2FoldChange",
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
    run_viper = FALSE,
    viper_n_tfs = 10,
    run_ppi = TRUE,
    ppi_top_n_hubs = 10,
    run_cspa = TRUE,
    cspa_n_markers = 25,
    render_html = TRUE,
    save_dir = getwd()) {

    call <- match.call()
    save_dir <- validate_save_dir(save_dir)

    design_resolved <- resolve_design(
        var, batch = batch, design = design,
        include_interactions = include_interactions
    )
    params <- list(
        var = var,
        design = paste(deparse(design_resolved), collapse = ""),
        comparisons = comparisons,
        include_interactions = include_interactions,
        batch = batch, batch_method = batch_method,
        export_counts = export_counts, filter_method = filter_method,
        cpm_percentile = cpm_percentile, min_count = min_count,
        min_samples_per_group = min_samples_per_group,
        remove_mt = remove_mt, remove_xy = remove_xy, remove_vdj = remove_vdj,
        shrinkage = shrinkage, alpha = alpha, lfc_threshold = lfc_threshold,
        ranking = ranking, collections = collections, run_clusterprofiler = run_clusterprofiler,
        run_wgcna = run_wgcna, wgcna_power = wgcna_power,
        wgcna_min_module_size = wgcna_min_module_size,
        run_decoupler = run_decoupler, decoupler_n_tfs = decoupler_n_tfs,
        run_viper = run_viper, viper_n_tfs = viper_n_tfs,
        run_ppi = run_ppi, ppi_top_n_hubs = ppi_top_n_hubs,
        run_cspa = run_cspa, cspa_n_markers = cspa_n_markers,
        render_html = render_html
    )

    log_output(call, save_dir, params = params, expr = quote({

        message("Running deseq2reports Wald pipeline v", utils::packageVersion("deseq2reports"), "...")

        dds <- validate_dds(dds)
        dds <- validate_var(dds, var)
        dds <- validate_design(dds, design_resolved)
        contrast_specs <- validate_comparisons(
            dds, var, comparisons,
            include_interactions = include_interactions
        )
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

        n_steps <- 6L + as.integer(isTRUE(export_counts)) +
            as.integer(isTRUE(run_wgcna)) + as.integer(isTRUE(run_decoupler)) +
            as.integer(isTRUE(run_viper)) + as.integer(isTRUE(run_ppi)) +
            as.integer(isTRUE(run_cspa))
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
        message("[", step, "/", n_steps, "] DESeq2 analysis (", n_contrasts, " contrast(s))...")
        diffexp_out <- run_diffexp_pipeline(
            dds = dds, var = var, comparisons = comparisons, batch = batch,
            include_interactions = include_interactions,
            shrinkage = shrinkage, alpha = alpha, lfc_threshold = lfc_threshold,
            ranking = ranking, top_label_n = top_label_n,
            save_dir = file.path(save_dir, "deseq2"))

        cspa_out <- NULL
        if (run_cspa) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] CSPA surface-marker heatmaps...")
            cspa_out <- tryCatch(
                run_cspa_marker_pipeline(
                    dds = dds,
                    diffexp_out = diffexp_out,
                    var = var,
                    n_markers = cspa_n_markers,
                    save_dir = file.path(save_dir, "deseq2", "cspa")
                ),
                error = function(e) {
                    warning("CSPA marker analysis failed (",
                            conditionMessage(e), "); continuing without CSPA heatmaps.",
                            call. = FALSE)
                    NULL
                }
            )
        }

        step <- step + 1L
        message("[", step, "/", n_steps, "] clusterProfiler analysis...")
        gsea_out <- run_gsea_pipeline(
            diffexp_results = diffexp_out$results, collections = collections,
            ranking = ranking, top_n_plot = top_n_plot,
            run_clusterprofiler = run_clusterprofiler,
            save_dir = file.path(save_dir, "clusterprofiler"))

        step <- step + 1L
        message("[", step, "/", n_steps, "] GSEA leading-edge overlap...")
        gsea_overlap_out <- tryCatch(
            run_gsea_overlap_pipeline(
                gsea_out = gsea_out,
                alpha = alpha,
                # Per comparison: clusterprofiler/<comp>/gsea_overlap/
                save_dir = file.path(save_dir, "clusterprofiler")
            ),
            error = function(e) {
                warning("GSEA overlap analysis failed (",
                        conditionMessage(e),
                        "); continuing without leading-edge overlap section.",
                        call. = FALSE)
                NULL
            }
        )

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
            message("[", step, "/", n_steps, "] DecoupleR (CollecTRI + ULM)...")
            decoupler_out <- run_decoupler_pipeline(
                dds = dds, var = viz_var, diffexp_out = diffexp_out,
                n_tfs = decoupler_n_tfs,
                save_dir = file.path(save_dir, "decoupler"))
        }

        viper_out <- NULL
        if (run_viper) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] VIPER (bcellViper)...")
            viper_out <- tryCatch(
                run_viper_pipeline(
                    dds = dds, var = viz_var, diffexp_out = diffexp_out,
                    n_tfs = viper_n_tfs,
                    save_dir = file.path(save_dir, "viper")
                ),
                error = function(e) {
                    warning("VIPER analysis failed (",
                            conditionMessage(e), "); continuing without VIPER page.",
                            call. = FALSE)
                    NULL
                }
            )
        }

        ppi_out <- NULL
        if (run_ppi) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] PPI connectivity (STRINGdb)...")
            ppi_out <- tryCatch(
                run_ppi_pipeline(
                    diffexp_out = diffexp_out,
                    wgcna_out = wgcna_out,
                    top_n_hubs = ppi_top_n_hubs,
                    alpha = alpha,
                    lfc_threshold = lfc_threshold,
                    save_dir = file.path(save_dir, "ppi")
                ),
                error = function(e) {
                    warning("PPI connectivity analysis failed (",
                            conditionMessage(e), "); continuing without PPI plots.",
                            call. = FALSE)
                    NULL
                }
            )
        }

        step <- step + 1L
        message("[", step, "/", n_steps, "] Cross-analysis gene hit ranking...")
        multihit_out <- tryCatch(
            run_multihit_pipeline(
                diffexp_out = diffexp_out,
                gsea_out = gsea_out,
                viper_out = viper_out,
                decoupler_out = decoupler_out,
                cspa_out = cspa_out,
                alpha = alpha,
                top_n = 10,
                save_dir = file.path(save_dir, "multihit")
            ),
            error = function(e) {
                warning("Multi-hit gene ranking failed (",
                        conditionMessage(e), "); continuing without index barplots.",
                        call. = FALSE)
                NULL
            }
        )

        report_path <- NULL
        if (render_html) {
            message("Rendering HTML report...")
            ppi_for_report <- ppi_out
            report_path <- render_report(
                dds = dds, var = var, distance_out = distance_out,
                diffexp_out = diffexp_out, gsea_out = gsea_out,
                gsea_overlap_out = gsea_overlap_out,
                multihit_out = multihit_out,
                wgcna_out = wgcna_out, decoupler_out = decoupler_out,
                ppi_out = ppi_for_report, viper_out = viper_out,
                cspa_out = cspa_out,
                alpha = alpha, lfc_threshold = lfc_threshold,
                run_args = params,
                save_dir = save_dir)
        }

        dds_rds <- file.path(save_dir, "dds.rds")
        message("Saving QC-filtered DESeqDataSet to: ", dds_rds)
        saveRDS(dds, dds_rds)

        message("Pipeline complete. Outputs written to: ", save_dir)
        invisible(list(
            dds = dds, distance = distance_out, diffexp = diffexp_out,
            gsea = gsea_out, gsea_overlap = gsea_overlap_out,
            multihit = multihit_out,
            wgcna = wgcna_out, decoupler = decoupler_out, viper = viper_out,
            ppi = ppi_out, cspa = cspa_out,
            report_path = report_path,
            dds_path = dds_rds
        ))
    }))
}

#' @rdname run_deseq2_wald
#' @export
run_deseq2_pipeline <- function(...) {
    .Deprecated("run_deseq2_wald")
    run_deseq2_wald(...)
}

#' Run the complete RNA-seq DESeq2 LRT (time-course) pipeline in one command
#'
#' Parallel entry point to [run_deseq2_wald()] for likelihood-ratio testing of
#' time and/or grouping factors. Shares QC, counts export, PCA, WGCNA, and
#' sample-level TF modules with the Wald pipeline, but replaces pairwise
#' Wald contrasts with the nested full/reduced grid from
#' [expand_lrt_specs()] / [run_diffexp_lrt_pipeline()].
#'
#' Significant LRT genes (`padj < alpha`) are clustered with
#' `DEGreport::degPatterns()` (x = `time`, colour = first of `var` when
#' present). Those clusters drive STRING PPI (top 50 nodes), CSPA heatmaps
#' with cluster row annotation, and per-cluster ORA
#' ([run_cluster_ora_pipeline()]). Contrast-level decoupleR / VIPER use the
#' **union of all LRT-significant genes** (not clusters), scored by the
#' unsigned LRT statistic; directional patterns are read from degPatterns.
#'
#' @param dds A `DESeqDataSet` of raw counts.
#' @param var Non-time biological factor column(s) in `colData(dds)` (e.g.
#'   `"genotype"`), or `NULL` / `character(0)` for a time-only design.
#' @param time Required `colData` column for the time / ordered trajectory
#'   factor (discrete levels).
#' @param design,comparisons Unused for LRT (grid is auto-built). Passing
#'   non-`NULL` values raises an error.
#' @param include_interactions If `TRUE` and `var` is non-empty, also test
#'   the `var * time` full model against nested reduced models. Default
#'   `FALSE`.
#' @param batch Optional covariate column name(s). Prefixed onto every LRT
#'   model and never tested.
#' @inheritParams run_deseq2_wald
#' @param cspa_n_markers Top CSPA genes **per degPatterns cluster** on the
#'   LRT CSPA heatmap. Default `10`.
#' @param degpatterns_minc Minimum genes per degPatterns cluster (passed to
#'   [DEGreport::degPatterns()]). Default `15`.
#' @return Invisibly, a list with `dds`, `dds_path`, `distance`, `diffexp`
#'   (LRT results), `degpatterns`, `cluster_ora`, `gsea` (`NULL`; ORA
#'   replaces ranked GSEA), `wgcna`, `decoupler`, `viper`, `ppi`, `cspa`,
#'   and `report_path`.
#' @examples
#' \dontrun{
#' out <- run_deseq2_lrt(
#'   dds, var = "genotype", time = "day",
#'   include_interactions = TRUE,
#'   save_dir = "results_lrt")
#' }
#' @export
run_deseq2_lrt <- function(
    dds,
    var = NULL,
    time,
    design = NULL,
    comparisons = NULL,
    include_interactions = FALSE,
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
    ranking = "stat",
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
    run_viper = FALSE,
    viper_n_tfs = 10,
    run_ppi = TRUE,
    ppi_top_n_hubs = 10,
    run_cspa = TRUE,
    cspa_n_markers = 10,
    degpatterns_minc = 15,
    render_html = TRUE,
    save_dir = getwd()) {

    call <- match.call()
    save_dir <- validate_save_dir(save_dir)

    if (!is.null(design)) {
        stop("`design` is unused in run_deseq2_lrt(); LRT models are built from ",
             "`var`, `time`, `batch`, and `include_interactions`.", call. = FALSE)
    }
    if (!is.null(comparisons)) {
        stop("`comparisons` is unused in run_deseq2_lrt(); the LRT grid is ",
             "auto-built by expand_lrt_specs().", call. = FALSE)
    }
    if (is.null(var)) var <- character(0)
    if (!is.character(time) || length(time) != 1L || is.na(time) || !nzchar(time)) {
        stop("`time` must be a single non-empty character string.", call. = FALSE)
    }

    params <- list(
        analysis_mode = "lrt",
        var = var,
        time = time,
        include_interactions = include_interactions,
        batch = batch, batch_method = batch_method,
        export_counts = export_counts, filter_method = filter_method,
        cpm_percentile = cpm_percentile, min_count = min_count,
        min_samples_per_group = min_samples_per_group,
        remove_mt = remove_mt, remove_xy = remove_xy, remove_vdj = remove_vdj,
        shrinkage = shrinkage, alpha = alpha, lfc_threshold = lfc_threshold,
        ranking = ranking, collections = collections,
        run_clusterprofiler = run_clusterprofiler,
        run_wgcna = run_wgcna, wgcna_power = wgcna_power,
        wgcna_min_module_size = wgcna_min_module_size,
        run_decoupler = run_decoupler, decoupler_n_tfs = decoupler_n_tfs,
        run_viper = run_viper, viper_n_tfs = viper_n_tfs,
        run_ppi = run_ppi, ppi_top_n_hubs = ppi_top_n_hubs,
        run_cspa = run_cspa, cspa_n_markers = cspa_n_markers,
        degpatterns_minc = degpatterns_minc,
        render_html = render_html
    )

    log_output(call, save_dir, params = params, expr = quote({

        message("Running deseq2reports LRT pipeline v",
                utils::packageVersion("deseq2reports"), "...")

        dds <- validate_dds(dds)
        dds <- .lrt_factorize_cols(dds, time, "time")
        dds <- .lrt_factorize_cols(dds, var, "var")
        dds <- .lrt_factorize_cols(dds, batch, "batch")

        lrt_specs <- expand_lrt_specs(
            var = var, time = time, batch = batch,
            include_interactions = include_interactions
        )
        n_tests <- nrow(lrt_specs)

        # Visualization factor for QC / PCA / WGCNA: time alone, or
        # var × time paste when grouping factors are present.
        if (length(var)) {
            viz_cols <- c(var, time)
            viz_var <- paste(viz_cols, collapse = "_x_")
            SummarizedExperiment::colData(dds)[[viz_var]] <- interaction(
                as.data.frame(SummarizedExperiment::colData(dds))[viz_cols],
                sep = "_", drop = TRUE
            )
            if (is.null(shape_var)) shape_var <- var[[1L]]
        } else {
            viz_var <- time
        }

        n_steps <- 8L + as.integer(isTRUE(export_counts)) +
            as.integer(isTRUE(run_wgcna)) + as.integer(isTRUE(run_decoupler)) +
            as.integer(isTRUE(run_viper)) + as.integer(isTRUE(run_ppi)) +
            as.integer(isTRUE(run_cspa))
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
            message("[", step, "/", n_steps, "] Exporting count tables...")
            run_counts_export(
                dds = dds, batch = batch, var = viz_var,
                save_dir = file.path(save_dir, "counts"))
        }

        step <- step + 1L
        message("[", step, "/", n_steps, "] Sample similarity / PCA...")
        distance_out <- run_distance_pca(
            dds = dds, var = viz_var, batch = batch, batch_method = batch_method,
            pals = pals, shape_var = shape_var,
            n_scree_pcs = n_scree_pcs, n_loading_pcs = n_loading_pcs,
            n_top_genes = n_top_genes,
            save_dir = file.path(save_dir, "distance"))

        step <- step + 1L
        message("[", step, "/", n_steps, "] DESeq2 LRT (", n_tests, " test(s))...")
        diffexp_out <- run_diffexp_lrt_pipeline(
            dds = dds, var = var, time = time, batch = batch,
            include_interactions = include_interactions,
            alpha = alpha,
            save_dir = file.path(save_dir, "deseq2"))

        step <- step + 1L
        message("[", step, "/", n_steps, "] degPatterns temporal clustering...")
        degpatterns_out <- run_degpatterns_pipeline(
            dds = dds,
            lrt_out = diffexp_out,
            time = time,
            var = var,
            alpha = alpha,
            minc = degpatterns_minc,
            save_dir = file.path(save_dir, "deseq2")
        )

        cspa_out <- NULL
        if (run_cspa) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] CSPA heatmaps by temporal cluster...")
            hm <- run_lrt_cspa_heatmap(
                dds = dds,
                degpatterns_out = degpatterns_out,
                lrt_out = diffexp_out,
                n_per_cluster = cspa_n_markers,
                annotation_var = unique(c(time, var)),
                save_dir = file.path(save_dir, "deseq2", "cspa")
            )
            cspa_out <- list(
                source = "lrt_clusters",
                org = detect_organism(rownames(dds)),
                by_lrt = hm,
                heatmaps = lapply(hm, function(x) x$heatmap),
                n_markers = cspa_n_markers
            )
        }

        step <- step + 1L
        message("[", step, "/", n_steps, "] Cluster ORA (clusterProfiler page)...")
        cluster_ora_out <- run_cluster_ora_pipeline(
            degpatterns_out = degpatterns_out,
            universe = rownames(dds),
            collections = collections,
            top_n = top_n_plot,
            alpha = alpha,
            save_dir = file.path(save_dir, "clusterprofiler")
        )
        # Ranked GSEA is not used in LRT mode (no signed LFC ranking).
        gsea_out <- NULL
        gsea_overlap_out <- NULL

        wgcna_out <- NULL
        if (run_wgcna) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] WGCNA...")
            wgcna_out <- run_wgcna_pipeline(
                dds = dds, var = viz_var, power = wgcna_power,
                min_module_size = wgcna_min_module_size,
                ora_collection = collections,
                save_dir = file.path(save_dir, "wgcna"))
        }

        # One TF signature per degPatterns cluster (not whole-LRT sig set).
        tf_diffexp <- tryCatch(
            lrt_clusters_diffexp_for_tf(
                degpatterns_out = degpatterns_out,
                lrt_out = diffexp_out
            ),
            error = function(e) {
                warning("No degPatterns clusters for TF activity (",
                        conditionMessage(e), ").", call. = FALSE)
                NULL
            }
        )

        decoupler_out <- NULL
        if (run_decoupler) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] DecoupleR (sample + per-cluster ULM)...")
            if (is.null(tf_diffexp)) {
                warning("Skipping decoupleR: no temporal clusters for TF input.",
                        call. = FALSE)
            } else {
                decoupler_out <- tryCatch(
                    run_decoupler_pipeline(
                        dds = dds, var = viz_var, diffexp_out = tf_diffexp,
                        n_tfs = decoupler_n_tfs,
                        save_dir = file.path(save_dir, "decoupler")
                    ),
                    error = function(e) {
                        warning("decoupleR failed (", conditionMessage(e), ").",
                                call. = FALSE)
                        NULL
                    }
                )
                if (!is.null(decoupler_out)) {
                    decoupler_out$tf_level <- "lrt_cluster"
                }
            }
        }

        viper_out <- NULL
        if (run_viper) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] VIPER (sample + per-cluster msVIPER)...")
            if (is.null(tf_diffexp)) {
                warning("Skipping VIPER: no temporal clusters for TF input.",
                        call. = FALSE)
            } else {
                viper_out <- tryCatch(
                    run_viper_pipeline(
                        dds = dds, var = viz_var, diffexp_out = tf_diffexp,
                        n_tfs = viper_n_tfs,
                        save_dir = file.path(save_dir, "viper")
                    ),
                    error = function(e) {
                        warning("VIPER failed (", conditionMessage(e), ").",
                                call. = FALSE)
                        NULL
                    }
                )
                if (!is.null(viper_out)) {
                    viper_out$tf_level <- "lrt_cluster"
                }
            }
        }

        ppi_out <- NULL
        if (run_ppi) {
            step <- step + 1L
            message("[", step, "/", n_steps, "] PPI connectivity (clusters + WGCNA)...")
            clusters_wgcna <- tryCatch(
                lrt_clusters_as_wgcna_like(degpatterns_out),
                error = function(e) NULL
            )
            # Merge degPatterns clusters and WGCNA modules into one module_df
            # so PPI covers both; WGCNA HTML only displays real WGCNA keys.
            # Keep only columns PPI seeds on (gene + module [+ optional kME]).
            .ppi_module_cols <- function(df) {
                gene_col <- if ("gene" %in% colnames(df)) {
                    "gene"
                } else {
                    setdiff(colnames(df), c("module", "module_label",
                                            "module_color", "kME",
                                            "cluster", "lrt"))[1]
                }
                keep <- c(gene_col, "module",
                          intersect("kME", colnames(df)))
                out <- df[, keep, drop = FALSE]
                colnames(out)[1] <- "gene"
                if (!"kME" %in% colnames(out)) out$kME <- NA_real_
                out
            }
            ppi_modules <- NULL
            mod_parts <- list()
            if (!is.null(clusters_wgcna) && nrow(clusters_wgcna$module_df)) {
                mod_parts[[length(mod_parts) + 1L]] <-
                    .ppi_module_cols(clusters_wgcna$module_df)
            }
            if (!is.null(wgcna_out) && !is.null(wgcna_out$module_df) &&
                nrow(wgcna_out$module_df)) {
                mod_parts[[length(mod_parts) + 1L]] <-
                    .ppi_module_cols(wgcna_out$module_df)
            }
            if (length(mod_parts)) {
                merged <- do.call(rbind, mod_parts)
                rownames(merged) <- NULL
                ppi_modules <- list(module_df = merged)
            }
            ppi_out <- tryCatch(
                run_ppi_pipeline(
                    diffexp_out = diffexp_out,
                    wgcna_out = ppi_modules,
                    top_n_hubs = ppi_top_n_hubs,
                    max_plot_nodes = 50L,
                    alpha = alpha,
                    lfc_threshold = 0,
                    save_dir = file.path(save_dir, "ppi")
                ),
                error = function(e) {
                    warning("PPI connectivity analysis failed (",
                            conditionMessage(e), "); continuing without PPI plots.",
                            call. = FALSE)
                    NULL
                }
            )
        }

        multihit_out <- NULL
        step <- step + 1L
        message("[", step, "/", n_steps, "] Cross-analysis gene hit ranking (per cluster)...")
        multihit_out <- tryCatch(
            run_multihit_lrt_pipeline(
                degpatterns_out = degpatterns_out,
                diffexp_out = diffexp_out,
                wgcna_out = wgcna_out,
                viper_out = viper_out,
                decoupler_out = decoupler_out,
                alpha = alpha,
                top_n = 10,
                save_dir = file.path(save_dir, "multihit")
            ),
            error = function(e) {
                warning("Multi-hit gene ranking failed (",
                        conditionMessage(e), "); continuing without index tables.",
                        call. = FALSE)
                NULL
            }
        )

        message("Rendering HTML report...")
        report_path <- render_report(
            dds = dds, var = if (length(var)) var else time,
            distance_out = distance_out,
            diffexp_out = diffexp_out, gsea_out = gsea_out,
            gsea_overlap_out = gsea_overlap_out,
            multihit_out = multihit_out,
            wgcna_out = wgcna_out, decoupler_out = decoupler_out,
            ppi_out = ppi_out, viper_out = viper_out,
            cspa_out = cspa_out,
            degpatterns_out = degpatterns_out,
            cluster_ora_out = cluster_ora_out,
            alpha = alpha, lfc_threshold = lfc_threshold,
            run_args = params,
            save_dir = save_dir)

        dds_rds <- file.path(save_dir, "dds.rds")
        message("Saving QC-filtered DESeqDataSet to: ", dds_rds)
        saveRDS(dds, dds_rds)

        message("Pipeline complete. Outputs written to: ", save_dir)
        invisible(list(
            dds = dds, distance = distance_out, diffexp = diffexp_out,
            degpatterns = degpatterns_out, cluster_ora = cluster_ora_out,
            gsea = gsea_out, gsea_overlap = gsea_overlap_out,
            multihit = multihit_out,
            wgcna = wgcna_out, decoupler = decoupler_out, viper = viper_out,
            ppi = ppi_out, cspa = cspa_out,
            report_path = report_path,
            dds_path = dds_rds
        ))
    }))
}
