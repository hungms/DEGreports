#' Render the HTML report with the **targets** package (parallel, resumable)
#'
#' Drop-in alternative to [render_report()] that uses **targets** to render
#' independent hub pages in parallel (resumable via `tar_make()`). Per-contrast
#' pages are rendered afterward via [render_report()] with `contrasts_only =
#' TRUE`.
#'
#' Requires Suggested packages **targets** and **future**.
#'
#' @inheritParams render_report
#' @param script_path Path for the generated `_targets.R` script. Default
#'   `file.path(save_dir, "_targets.R")`.
#' @return Invisibly, the same named path vector as [render_report()], or
#'   `NULL` when pandoc is unavailable.
#' @seealso [render_report()], [run_deseq2()]
#' @export
render_report_targets <- function(
    dds,
    var,
    distance_out,
    diffexp_out,
    gsea_out,
    gsea_overlap_out = NULL,
    multihit_out = NULL,
    wgcna_out = NULL,
    decoupler_tf_out = NULL,
    ppi_out = NULL,
    viper_out = NULL,
    cspa_out = NULL,
    gsva_out = NULL,
    decoupler_ccc_out = NULL,
    degpatterns_out = NULL,
    cluster_ora_out = NULL,
    alpha = 0.05,
    lfc_threshold = 0.5,
    run_args = NULL,
    save_dir = getwd(),
    parallel_render = TRUE,
    render_workers = NULL,
    script_path = NULL) {

    if (!.pandoc_ok()) {
        warning(.pandoc_missing_msg(), " Skipping the report; the analysis outputs ",
                "already written to disk are unaffected.", call. = FALSE)
        return(invisible(NULL))
    }

    if (!requireNamespace("targets", quietly = TRUE)) {
        stop("Package \"targets\" is required for render_report_targets(). ",
             "Install with install.packages(\"targets\").", call. = FALSE)
    }

    save_dir <- validate_save_dir(save_dir)
    if (is.null(script_path)) {
        script_path <- file.path(save_dir, "_targets.R")
    }

    meta_rds <- file.path(save_dir, ".degreports_render_meta.rds")
    meta <- list(
        dds = dds,
        var = var,
        distance_out = distance_out,
        diffexp_out = diffexp_out,
        gsea_out = gsea_out,
        gsea_overlap_out = gsea_overlap_out,
        multihit_out = multihit_out,
        wgcna_out = wgcna_out,
        decoupler_tf_out = decoupler_tf_out,
        ppi_out = ppi_out,
        viper_out = viper_out,
        cspa_out = cspa_out,
        gsva_out = gsva_out,
        decoupler_ccc_out = decoupler_ccc_out,
        degpatterns_out = degpatterns_out,
        cluster_ora_out = cluster_ora_out,
        alpha = alpha,
        lfc_threshold = lfc_threshold,
        run_args = run_args,
        save_dir = save_dir
    )
    saveRDS(meta, meta_rds)

    meta_q <- dQuote(meta_rds, q = FALSE)
    repo_q <- dQuote(normalizePath(file.path(save_dir, "_targets"), mustWork = FALSE), q = FALSE)
    script_lines <- c(
        "library(targets)",
        "library(DEGreports)",
        sprintf("meta <- readRDS(%s)", meta_q),
        "render_order <- local({",
        "  ord <- character(0)",
        "  if (!is.null(meta$dds)) ord <- c(ord, \"qc\")",
        "  ord <- c(ord, \"deseq2\", \"clusterProfiler\")",
        "  if (!is.null(meta$gsva_out)) ord <- c(ord, \"gsva\")",
        "  if (!is.null(meta$wgcna_out)) ord <- c(ord, \"wgcna\")",
        "  if (!is.null(meta$decoupler_tf_out)) ord <- c(ord, \"decoupleR-TF\")",
        "  if (!is.null(meta$decoupler_ccc_out)) ord <- c(ord, \"decoupleR-CCC\")",
        "  if (!is.null(meta$viper_out)) ord <- c(ord, \"viper\")",
        "  if (!is.null(meta$multihit_out)) ord <- c(ord, \"multihit\")",
        "  c(ord, \"index\")",
        "})",
        "hub_keys <- render_order",
        "idx <- which(hub_keys == \"index\")",
        "if (length(idx)) hub_keys <- c(hub_keys[-idx], \"index\")",
        sprintf("tar_option_set(packages = \"DEGreports\", format = \"rds\", repository = %s)", repo_q),
        ".render_hub_page <- function(page_key) {",
        "  disabled <- character(0)",
        "  if (is.null(meta$dds)) disabled <- c(disabled, \"QC\")",
        "  if (is.null(meta$gsva_out)) disabled <- c(disabled, \"GSVA\")",
        "  if (is.null(meta$wgcna_out)) disabled <- c(disabled, \"WGCNA\")",
        "  if (is.null(meta$decoupler_tf_out)) disabled <- c(disabled, \"decoupleR-TF\")",
        "  if (is.null(meta$decoupler_ccc_out)) disabled <- c(disabled, \"decoupleR-CCC\")",
        "  if (is.null(meta$viper_out)) disabled <- c(disabled, \"VIPER\")",
        "  if (is.null(meta$multihit_out)) disabled <- c(disabled, \"Multihit\")",
        "  shared <- DEGreports:::.report_shared_params(",
        "    meta$dds, meta$var, meta$distance_out, meta$diffexp_out, meta$gsea_out,",
        "    meta$wgcna_out, meta$decoupler_tf_out, meta$ppi_out, meta$viper_out, meta$cspa_out,",
        "    meta$alpha, meta$lfc_threshold, meta$run_args, disabled,",
        "    gsea_overlap_out = meta$gsea_overlap_out, multihit_out = meta$multihit_out,",
        "    degpatterns_out = meta$degpatterns_out, cluster_ora_out = meta$cluster_ora_out,",
        "    gsva_out = meta$gsva_out, decoupler_ccc_out = meta$decoupler_ccc_out",
        "  )",
        "  params <- DEGreports:::.report_params_for_page(page_key, shared)",
        "  DEGreports:::.render_report_page(page_key, params, meta$save_dir)",
        "}",
        "list(",
        "  tar_target(render_meta, meta, format = \"object\"),",
        "  tar_map(",
        "    values = tibble::tibble(page_key = hub_keys),",
        "    tar_target(",
        "      hub_html,",
        "      .render_hub_page(page_key),",
        "      format = \"file\"",
        "    )",
        "  )",
        ")"
    )
    writeLines(script_lines, script_path)

    n_hubs <- length({
        ord <- character(0)
        if (!is.null(dds)) ord <- c(ord, "qc")
        ord <- c(ord, "deseq2", "clusterProfiler")
        if (!is.null(gsva_out)) ord <- c(ord, "gsva")
        if (!is.null(wgcna_out)) ord <- c(ord, "wgcna")
        if (!is.null(decoupler_tf_out)) ord <- c(ord, "decoupleR-TF")
        if (!is.null(decoupler_ccc_out)) ord <- c(ord, "decoupleR-CCC")
        if (!is.null(viper_out)) ord <- c(ord, "viper")
        if (!is.null(multihit_out)) ord <- c(ord, "multihit")
        c(ord, "index")
    })
    workers <- render_workers
    if (is.null(workers)) {
        workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
    }
    if (requireNamespace("future", quietly = TRUE)) {
        future::plan(future::multisession, workers = workers)
    }

    message("Running targets render pipeline (", n_hubs, " hub page(s))...")
    targets::tar_make(script = script_path, callr_function = NULL)

    render_report(
        dds = dds,
        var = var,
        distance_out = distance_out,
        diffexp_out = diffexp_out,
        gsea_out = gsea_out,
        gsea_overlap_out = gsea_overlap_out,
        multihit_out = multihit_out,
        wgcna_out = wgcna_out,
        decoupler_tf_out = decoupler_tf_out,
        ppi_out = ppi_out,
        viper_out = viper_out,
        cspa_out = cspa_out,
        gsva_out = gsva_out,
        decoupler_ccc_out = decoupler_ccc_out,
        degpatterns_out = degpatterns_out,
        cluster_ora_out = cluster_ora_out,
        alpha = alpha,
        lfc_threshold = lfc_threshold,
        run_args = run_args,
        save_dir = save_dir,
        parallel_render = parallel_render,
        render_workers = render_workers,
        contrasts_only = TRUE
    )
}
