#' Normalize and validate a merged DEG results table
#'
#' Accepts a single data frame that stacks every contrast, distinguished by
#' a `comparison` column (e.g. `"groupA_vs_groupB"`, `"groupA_vs_groupC"`).
#'
#' **Required columns** (aliases accepted; see Details):
#' \itemize{
#'   \item `gene` -- gene symbol (HGNC / MGI)
#'   \item `comparison` -- contrast label
#'   \item `log2FoldChange` -- effect size
#'   \item `padj` -- adjusted p-value (or `pvalue` if `padj` missing)
#' }
#'
#' **Strongly recommended:**
#' \itemize{
#'   \item `stat` -- Wald / signed statistic (used by GSEA `"stat"` ranking,
#'     decoupleR ULM, and msVIPER). If absent, derived as
#'     `sign(log2FoldChange) * abs(qnorm(padj/2))`.
#'   \item `pvalue` -- raw p-value (derived from `padj` if missing)
#'   \item `baseMean` / `mean_cpm` -- for MA plots (`mean_cpm` falls back to
#'     `baseMean` when present)
#' }
#'
#' Optional aliases remapped automatically:
#' `Gene` / `gene_name` / `symbol` → `gene`;
#' `contrast` / `Contrast` → `comparison`;
#' `log2FC` / `logFC` / `LFC` → `log2FoldChange`;
#' `FDR` / `fdr` / `qvalue` / `adj.P.Val` → `padj`;
#' `P.Value` / `pval` / `p_value` → `pvalue`;
#' `AveExpr` → `baseMean`.
#'
#' @param deg_table Data frame (or path to CSV/TSV readable by
#'   [utils::read.delim()] / [utils::read.csv()]).
#' @param gene_col,comparison_col Optional overrides if names are non-standard
#'   after alias remapping. Default `NULL`.
#' @param alpha,lfc_threshold Used to (re)compute `sig` when missing.
#'   Defaults `0.05` / `0.5`.
#' @param ranking Ranking metric for `rank_score`
#'   (see [compute_ranking_score()]). Default `"log2FoldChange"`. When
#'   required score columns are incomplete, a simpler fallback is used.
#' @return A list shaped like [run_diffexp_pipeline()] output:
#'   `results`, `by_comparison`, `contrast_specs`, `gene_col`,
#'   `source = "deg_table"`.
#' @examples
#' \dontrun{
#' deg <- read.csv("all_contrasts_deg.csv")
#' diffexp_out <- prepare_diffexp_from_deg(deg)
#' }
#' @export
prepare_diffexp_from_deg <- function(
    deg_table,
    gene_col = NULL,
    comparison_col = NULL,
    alpha = 0.05,
    lfc_threshold = 0.5,
    ranking = "log2FoldChange") {

    alpha <- validate_alpha(alpha)
    lfc_threshold <- as.numeric(lfc_threshold)[1]
    if (is.na(lfc_threshold) || lfc_threshold < 0) {
        stop("`lfc_threshold` must be a non-negative number.", call. = FALSE)
    }
    ranking <- validate_ranking(ranking)

    df <- .read_deg_table_input(deg_table)
    df <- .remap_deg_colnames(df)

    if (!is.null(gene_col)) {
        if (!gene_col %in% colnames(df)) {
            stop("`gene_col` \"", gene_col, "\" not found in deg_table.",
                 call. = FALSE)
        }
        if (!identical(gene_col, "gene")) {
            df$gene <- df[[gene_col]]
        }
    }
    if (!is.null(comparison_col)) {
        if (!comparison_col %in% colnames(df)) {
            stop("`comparison_col` \"", comparison_col,
                 "\" not found in deg_table.", call. = FALSE)
        }
        if (!identical(comparison_col, "comparison")) {
            df$comparison <- df[[comparison_col]]
        }
    }

    need <- c("gene", "comparison", "log2FoldChange")
    missing <- setdiff(need, colnames(df))
    if (length(missing)) {
        stop(
            "DEG table is missing required column(s): ",
            paste(missing, collapse = ", "), ".\n",
            "Required: gene, comparison, log2FoldChange, and padj (or pvalue).\n",
            "Example comparison labels: \"groupA_vs_groupB\", \"groupA_vs_groupC\".\n",
            "Present columns: ", paste(colnames(df), collapse = ", "),
            call. = FALSE
        )
    }
    if (!"padj" %in% colnames(df) && !"pvalue" %in% colnames(df)) {
        stop("DEG table needs `padj` and/or `pvalue`.", call. = FALSE)
    }

    df$gene <- as.character(df$gene)
    df$comparison <- as.character(df$comparison)
    df$log2FoldChange <- as.numeric(df$log2FoldChange)
    if ("padj" %in% colnames(df)) df$padj <- as.numeric(df$padj)
    if ("pvalue" %in% colnames(df)) df$pvalue <- as.numeric(df$pvalue)

    keep <- !is.na(df$gene) & nzchar(df$gene) &
        !is.na(df$comparison) & nzchar(df$comparison)
    if (!all(keep)) {
        message("Dropping ", sum(!keep),
                " row(s) with missing gene/comparison.")
        df <- df[keep, , drop = FALSE]
    }
    if (!nrow(df)) {
        stop("No usable rows remain in deg_table.", call. = FALSE)
    }

    if (!"padj" %in% colnames(df)) {
        message("No `padj` column; using `pvalue` as padj.")
        df$padj <- df$pvalue
    }
    if (!"pvalue" %in% colnames(df)) {
        df$pvalue <- df$padj
    }
    if (!"stat" %in% colnames(df)) {
        message("No `stat` column; deriving a signed z-like score from padj ",
                "and log2FoldChange for TF / optional GSEA-stat ranking.")
        p <- df$padj
        p[!is.finite(p) | p <= 0] <- .Machine$double.xmin
        p[p > 1] <- 1
        df$stat <- sign(df$log2FoldChange) * abs(stats::qnorm(p / 2))
        df$stat[!is.finite(df$stat)] <- 0
    } else {
        df$stat <- as.numeric(df$stat)
    }

    if (!"mean_cpm" %in% colnames(df)) {
        if ("baseMean" %in% colnames(df)) {
            df$mean_cpm <- as.numeric(df$baseMean)
        } else {
            # MA plot needs a numeric x; volcano does not.
            df$mean_cpm <- 1
            message("No `mean_cpm`/`baseMean`; MA plots will use a flat x-axis.")
        }
    } else {
        df$mean_cpm <- as.numeric(df$mean_cpm)
    }
    if (!"baseMean" %in% colnames(df)) {
        df$baseMean <- df$mean_cpm
    }

    # rank_score: prefer full compute_ranking_score; fall back if columns thin
    df$rank_score <- tryCatch(
        compute_ranking_score(df, ranking = ranking),
        error = function(e) {
            message("Could not compute full rank_score (", conditionMessage(e),
                    "); using log2FoldChange.")
            replace(df$log2FoldChange, is.na(df$log2FoldChange), 0)
        }
    )

    if (!"sig" %in% colnames(df)) {
        df$sig <- !is.na(df$padj) & df$padj < alpha &
            abs(df$log2FoldChange) > lfc_threshold
    } else {
        df$sig <- as.logical(df$sig)
        df$sig[is.na(df$sig)] <- FALSE
    }

    comps <- unique(df$comparison)
    if (length(comps) < 1L) {
        stop("No comparison labels found.", call. = FALSE)
    }
    message("DEG table: ", nrow(df), " row(s), ", length(comps),
            " contrast(s): ", paste(comps, collapse = ", "))

    by_comparison <- split(df, df$comparison)
    by_comparison <- lapply(by_comparison, function(x) {
        rownames(x) <- NULL
        x
    })

    specs <- data.frame(
        label = names(by_comparison),
        class = "single",
        type = classify_comparison_type(names(by_comparison)),
        stratum_var = NA_character_,
        stratum_level = NA_character_,
        comparison = names(by_comparison),
        stringsAsFactors = FALSE
    )

    list(
        results = df,
        by_comparison = by_comparison,
        contrast_specs = specs,
        gene_col = "gene",
        source = "deg_table",
        alpha = alpha,
        lfc_threshold = lfc_threshold,
        ranking = ranking
    )
}

#' @noRd
.read_deg_table_input <- function(deg_table) {
    if (is.character(deg_table) && length(deg_table) == 1L &&
            file.exists(deg_table)) {
        path <- normalizePath(deg_table, mustWork = TRUE)
        if (grepl("\\.csv$", path, ignore.case = TRUE)) {
            return(utils::read.csv(path, stringsAsFactors = FALSE,
                                   check.names = FALSE))
        }
        return(utils::read.delim(path, stringsAsFactors = FALSE,
                                 check.names = FALSE))
    }
    if (is.data.frame(deg_table)) {
        return(as.data.frame(deg_table, stringsAsFactors = FALSE,
                             check.names = FALSE))
    }
    stop("`deg_table` must be a data frame or a path to a CSV/TSV file.",
         call. = FALSE)
}

#' @noRd
.remap_deg_colnames <- function(df) {
    alias <- c(
        Gene = "gene", gene_name = "gene", GeneSymbol = "gene",
        genesymbol = "gene", symbol = "gene", SYMBOL = "gene",
        contrast = "comparison", Contrast = "comparison",
        COMPARISON = "comparison",
        log2FC = "log2FoldChange", logFC = "log2FoldChange",
        Log2FC = "log2FoldChange", LFC = "log2FoldChange",
        log2foldchange = "log2FoldChange",
        FDR = "padj", fdr = "padj", qvalue = "padj",
        `adj.P.Val` = "padj", adj.P.Val = "padj",
        `P.Value` = "pvalue", pval = "pvalue", p_value = "pvalue",
        PValue = "pvalue", AveExpr = "baseMean", aveExpr = "baseMean"
    )
    cn <- colnames(df)
    for (i in seq_along(cn)) {
        if (cn[i] %in% names(alias) && !alias[[cn[i]]] %in% cn) {
            cn[i] <- alias[[cn[i]]]
        }
    }
    colnames(df) <- cn
    # If both gene and Gene existed, keep first gene.
    if (anyDuplicated(colnames(df))) {
        df <- df[, !duplicated(colnames(df)), drop = FALSE]
    }
    df
}

#' CSPA marker tables from DEG results (no expression heatmap)
#' @noRd
.run_cspa_markers_deg_only <- function(
    diffexp_out,
    n_markers = 25,
    cspa_file = NULL,
    gene_col = "gene",
    significant_only = TRUE,
    save_dir = getwd()) {

    save_dir <- validate_save_dir(save_dir)
    by_comp <- diffexp_out$by_comparison
    sample_genes <- unlist(lapply(by_comp, function(df) {
        as.character(df[[gene_col]])
    }), use.names = FALSE)
    org <- detect_organism(sample_genes)
    message("Loading CSPA surface markers for organism=\"", org, "\" ...")
    cspa_symbols <- load_cspa_genes(organism = org, file = cspa_file)
    message("  ", length(cspa_symbols), " CSPA symbols (marker selection only; ",
            "no VST heatmap without counts).")

    markers_by_comparison <- list()
    for (comp in names(by_comp)) {
        sel <- select_cspa_markers(
            by_comp[[comp]], cspa_symbols,
            n_markers = n_markers,
            gene_col = gene_col,
            significant_only = significant_only
        )
        markers_by_comparison[[comp]] <- sel
        if (nrow(sel)) {
            write_output(sel, file.path(save_dir, comp, "cspa", "cspa_markers.csv"))
            message("  ", comp, ": ", nrow(sel), " CSPA marker(s).")
        }
    }
    combined <- unique(do.call(rbind, Filter(function(x) nrow(x) > 0, markers_by_comparison)))
    if (!is.null(combined) && nrow(combined)) {
        write_output(combined, file.path(save_dir, "cspa_markers_combined.csv"))
    }
    list(
        source = "deg_table",
        org = org,
        cspa_n = length(cspa_symbols),
        n_markers = as.integer(n_markers)[1],
        markers_by_comparison = markers_by_comparison,
        combined_markers = combined,
        heatmaps = list()
    )
}

#' Downstream analysis from a merged DEG table (no counts)
#'
#' Starts from a **merged** differential-expression table with a
#' `comparison` column and runs every analysis that does not need raw
#' counts / a `DESeqDataSet`:
#' \enumerate{
#'   \item GSEA ([run_gsea_pipeline()])
#'   \item GSEA leading-edge overlap ([run_gsea_overlap_pipeline()])
#'   \item CSPA marker selection (no expression heatmap)
#'   \item Contrast-level CollecTRI / ULM ([run_decoupler_tf_pipeline()]
#'     without sample heatmap)
#'   \item Contrast-level msVIPER ([run_viper_pipeline()] without sample
#'     heatmap)
#'   \item STRING PPI among significant DEGs ([run_ppi_pipeline()])
#'   \item Cross-analysis multihit ranking ([run_multihit_pipeline()])
#'   \item HTML report ([render_report()]; QC / WGCNA pages disabled)
#' }
#'
#' See [prepare_diffexp_from_deg()] for the required table format.
#'
#' @param deg_table Merged DEG data frame or path to CSV/TSV. Must contain
#'   at least `gene`, `comparison`, `log2FoldChange`, and `padj` (or
#'   `pvalue`). Contrast labels should look like `"groupA_vs_groupB"`.
#' @param gene_col,comparison_col Optional column-name overrides.
#' @param ranking DEG/GSEA ranking metric (see [compute_ranking_score()]).
#'   Default `"log2FoldChange"`.
#' @param collections MSigDB collections for GSEA. Default includes
#'   HALLMARK, GO, KEGG, REACTOME, PID.
#' @param top_n_plot Pathways per NES direction on GSEA barplots. Default `10`.
#' @param run_gsea Logical. Default `TRUE`.
#' @param run_decoupler_tf Logical; contrast-level ULM. Default `FALSE`
#'   (optional module).
#' @param run_decoupler_ccc Logical; ligand–receptor ULM. Default `FALSE`
#'   (optional module).
#' @param run_viper Logical; contrast-level msVIPER. Default `FALSE`
#'   (optional module).
#' @param run_cspa Logical; CSPA **marker tables** only. Default `TRUE`.
#' @param render_html Logical; knit the multi-page HTML report. Default `TRUE`.
#' @param opts Analysis settings from [deseq2reports_opts()] -- `alpha`, the
#'   LFC cutoff, and the per-module top-N budgets. The defaults reproduce the
#'   values this entry point used to hard-code, and only settings that differ
#'   from them are listed in the report's "Pipeline arguments" table.
#' @param save_dir Output root directory. Default `getwd()`.
#' @return Invisibly, a list with `diffexp`, `gsea`, `gsea_overlap`,
#'   `cspa`, `decoupler`, `viper`, `ppi`, `multihit`, and `report_path`.
#' @examples
#' \dontrun{
#' out <- run_deseq2_deg(
#'   deg_table = "all_contrasts.csv",
#'   collections = c("HALLMARK", "GOBP"),
#'   save_dir = "results_from_deg"
#' )
#' }
#' @export
run_deseq2_deg <- function(
    deg_table,
    gene_col = NULL,
    comparison_col = NULL,
    ranking = "log2FoldChange",
    collections = c("HALLMARK", "GOBP", "GOCC", "GOMF",
                    "KEGG", "REACTOME", "PID"),
    top_n_plot = 10,
    run_gsea = TRUE,
    run_decoupler_tf = FALSE,
    run_decoupler_ccc = FALSE,
    run_viper = FALSE,
    run_cspa = TRUE,
    render_html = TRUE,
    opts = deseq2reports_opts(),
    save_dir = getwd()) {

    call <- match.call()
    save_dir <- validate_save_dir(save_dir)
    opts <- .as_opts(opts)
    alpha <- opts$alpha
    lfc_threshold <- opts$lfc_threshold
    params <- list(
        analysis_mode = "deg_table",
        ranking = ranking,
        collections = collections, top_n_plot = top_n_plot,
        run_gsea = run_gsea, run_decoupler_tf = run_decoupler_tf,
        run_decoupler_ccc = run_decoupler_ccc,
        run_viper = run_viper,
        run_cspa = run_cspa,
        render_html = render_html,
        save_dir = save_dir
    )
    params <- c(params, .opts_diff(opts))

    log_output(call, save_dir, expr = {
        .warn_if_no_pandoc(render_html)

        # The DEG-table run as data; see .stage() / .run_stages(). The same
        # driver as the Wald and LRT entry points, minus every stage that
        # needs counts.
        stages <- list(
            .stage("diffexp", "Preparing DEG table", optional = FALSE,
                   run = function(st) {
                out <- prepare_diffexp_from_deg(
                    deg_table = deg_table, gene_col = gene_col,
                    comparison_col = comparison_col, alpha = alpha,
                    lfc_threshold = lfc_threshold, ranking = ranking)
                # Per-contrast CSVs under deseq2/, for the report downloads.
                deseq2_dir <- file.path(save_dir, "deseq2")
                dir.create(deseq2_dir, recursive = TRUE, showWarnings = FALSE)
                for (comp in names(out$by_comparison)) {
                    write_output(
                        out$by_comparison[[comp]],
                        file.path(deseq2_dir, comp,
                                  paste0(contrast_file_stem(comp), "_diffexp.csv")))
                }
                write_output(out$results,
                             file.path(deseq2_dir, "all_comparisons_diffexp.csv"))
                out
            }),
            .stage("gsea", "GSEA", enabled = isTRUE(run_gsea), run = function(st) {
                run_gsea_pipeline(
                    diffexp_results = st$diffexp$results, collections = collections,
                    ranking = ranking, top_n_plot = top_n_plot,
                    save_dir = file.path(save_dir, "clusterprofiler"))
            }),
            .stage("gsea_overlap", "GSEA leading-edge overlap",
                   enabled = isTRUE(run_gsea), run = function(st) {
                run_gsea_overlap_pipeline(
                    gsea_out = st$gsea, alpha = alpha,
                    save_dir = file.path(save_dir, "clusterprofiler"))
            }),
            .stage("cspa", "CSPA marker selection", enabled = isTRUE(run_cspa),
                   run = function(st) {
                .run_cspa_markers_deg_only(
                    diffexp_out = st$diffexp, n_markers = opts$cspa_n_markers,
                    save_dir = file.path(save_dir, "deseq2"))
            }),
            .stage("decoupler_tf", "DecoupleR (contrast-level ULM)",
                   enabled = isTRUE(run_decoupler_tf), run = function(st) {
                run_decoupler_tf_pipeline(
                    dds = NULL, var = NULL, diffexp_out = st$diffexp,
                    n_tfs = opts$decoupler_n_tfs,
                    save_dir = file.path(save_dir, "decoupler-tf"))
            }),
            .stage("decoupler_ccc", "DecoupleR-CCC (ligand\u2013receptor ULM)",
                   enabled = isTRUE(run_decoupler_ccc), run = function(st) {
                run_decoupler_ccc_pipeline(
                    diffexp_out = st$diffexp, n_lr = opts$decoupler_n_lr,
                    save_dir = file.path(save_dir, "decoupler-ccc"))
            }),
            .stage("viper", "VIPER (contrast-level msVIPER)",
                   enabled = isTRUE(run_viper), run = function(st) {
                run_viper_pipeline(
                    dds = NULL, var = NULL, diffexp_out = st$diffexp,
                    n_tfs = opts$viper_n_tfs,
                    save_dir = file.path(save_dir, "viper"))
            }),
            .stage("ppi", "PPI connectivity", run = function(st) {
                run_ppi_pipeline(
                    diffexp_out = st$diffexp, wgcna_out = NULL,
                    top_n_hubs = opts$ppi_top_n_hubs, alpha = alpha,
                    lfc_threshold = lfc_threshold,
                    save_dir = file.path(save_dir, "deseq2"),
                    module_save_dir = file.path(save_dir, "ppi"))
            }),
            .stage("multihit", "Cross-analysis gene hit ranking", run = function(st) {
                run_multihit_pipeline(
                    diffexp_out = st$diffexp, gsea_out = st$gsea,
                    viper_out = st$viper, decoupler_tf_out = st$decoupler_tf,
                    decoupler_ccc_out = st$decoupler_ccc, cspa_out = st$cspa,
                    alpha = alpha, top_n = 10,
                    save_dir = file.path(save_dir, "multihit"))
            })
        )

        st <- .run_stages(stages)
        diffexp_out <- st$diffexp
        gsea_out <- st$gsea
        gsea_overlap_out <- st$gsea_overlap
        cspa_out <- st$cspa
        decoupler_tf_out <- st$decoupler_tf
        decoupler_ccc_out <- st$decoupler_ccc
        viper_out <- st$viper
        ppi_out <- st$ppi
        multihit_out <- st$multihit

        report_path <- NULL
        if (render_html) {
            message("Rendering HTML report...")
            report_path <- render_report(
                dds = NULL,
                var = "comparison",
                distance_out = NULL,
                diffexp_out = diffexp_out,
                gsea_out = gsea_out,
                gsea_overlap_out = gsea_overlap_out,
                multihit_out = multihit_out,
                wgcna_out = NULL,
                decoupler_tf_out = decoupler_tf_out,
                decoupler_ccc_out = decoupler_ccc_out,
                ppi_out = ppi_out,
                viper_out = viper_out,
                cspa_out = cspa_out,
                alpha = alpha,
                lfc_threshold = lfc_threshold,
                run_args = params,
                save_dir = save_dir
            )
        }

        message("DEG-table pipeline complete. Outputs written to: ", save_dir)
        invisible(list(
            diffexp = diffexp_out,
            gsea = gsea_out,
            gsea_overlap = gsea_overlap_out,
            cspa = cspa_out,
            decoupler_tf = decoupler_tf_out,
            decoupler_ccc = decoupler_ccc_out,
            viper = viper_out,
            ppi = ppi_out,
            multihit = multihit_out,
            report_path = report_path
        ))
    })
}
