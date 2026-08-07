#' Load CSPA cell-surface gene symbols (Human / Mouse)
#'
#' Reads the Cell Surface Protein Atlas supplementary tables from
#' [Bausch-Fluck et al.](https://doi.org/10.1371/journal.pone.0121314)
#' (PLOS ONE S3 File). Sheet 1 (`Table A`) = human; sheet 2 (`Table B`) =
#' mouse. By default **all** CSPA categories are kept (high confidence,
#' putative, and unspecific). Symbols come from the `ENTREZ gene symbol`
#' column.
#'
#' @param organism `"human"` or `"mouse"`.
#' @param file Path to the CSPA Excel workbook. Default `NULL` uses the
#'   file shipped under `inst/extdata/` (or downloads it from the GitHub
#'   repo into the user cache).
#' @param categories Optional character vector of `CSPA category` values to
#'   keep. Default `NULL` keeps every category.
#' @return A character vector of unique gene symbols (as stored in the
#'   table -- Title Case for mouse, UPPER for human).
#' @export
load_cspa_genes <- function(organism = c("human", "mouse"),
                            file = NULL,
                            categories = NULL) {
    organism <- match.arg(organism)
    if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package \"readxl\" is required to read the CSPA Excel table. ",
             "Install with install.packages(\"readxl\").", call. = FALSE)
    }
    file <- if (is.null(file)) {
        .cspa_default_path(download = TRUE)
    } else {
        normalizePath(file, mustWork = TRUE)
    }
    sheets <- readxl::excel_sheets(file)
    if (length(sheets) < 2L) {
        stop("CSPA workbook must have at least 2 sheets (human, mouse); found ",
             length(sheets), " in ", file, ".", call. = FALSE)
    }
    sheet <- if (identical(organism, "human")) sheets[[1L]] else sheets[[2L]]
    raw <- as.data.frame(readxl::read_excel(file, sheet = sheet))
    sym_col <- intersect(c("ENTREZ gene symbol", "ENTREZ_gene_symbol", "gene_symbol", "gene"),
                         colnames(raw))
    if (!length(sym_col)) {
        stop("Could not find an ENTREZ gene symbol column in CSPA sheet \"",
             sheet, "\". Columns: ", paste(colnames(raw), collapse = ", "),
             call. = FALSE)
    }
    sym_col <- sym_col[[1L]]
    if (!is.null(categories)) {
        cat_col <- intersect(c("CSPA category", "CSPA_category"), colnames(raw))
        if (!length(cat_col)) {
            stop("`categories` was set but no CSPA category column was found.",
                 call. = FALSE)
        }
        keep <- as.character(raw[[cat_col[[1L]]]]) %in% as.character(categories)
        raw <- raw[keep, , drop = FALSE]
    }
    genes <- unique(as.character(raw[[sym_col]]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!length(genes)) {
        stop("No CSPA gene symbols remaining after filtering for organism=\"",
             organism, "\".", call. = FALSE)
    }
    genes
}

#' Match gene IDs to a CSPA symbol set (case-insensitive)
#' @noRd
.match_cspa_genes <- function(genes, cspa_symbols) {
    genes <- as.character(genes)
    ok <- !is.na(genes) & nzchar(genes)
    genes <- genes[ok]
    if (!length(genes) || !length(cspa_symbols)) return(character(0))
    key <- toupper(genes)
    hit <- key %in% toupper(as.character(cspa_symbols))
    unique(genes[hit])
}

#' Pick top CSPA surface markers from one contrast result table
#'
#' Prefers significant DEGs (`sig == TRUE` when present), ranked by
#' absolute `rank_score` (falling back to `log2FoldChange`). Returns at
#' most `n_markers` gene IDs that intersect the CSPA list.
#'
#' @param res_df A per-contrast DE data frame from [run_diffexp()].
#' @param cspa_symbols Character vector from [load_cspa_genes()].
#' @param n_markers Maximum markers to return. Default `25`.
#' @param gene_col Gene ID column. Default `"gene"`.
#' @param significant_only If `TRUE`, restrict to `sig` genes when that
#'   column exists and is non-empty; otherwise fall back to all ranked
#'   CSPA genes in `res_df`. Default `TRUE`.
#' @return A data frame of selected markers (subset of `res_df` rows),
#'   ordered by decreasing `|rank_score|`.
#' @export
select_cspa_markers <- function(res_df,
                                cspa_symbols,
                                n_markers = 25,
                                gene_col = "gene",
                                significant_only = TRUE) {
    if (is.null(res_df) || !nrow(res_df)) {
        return(res_df[0, , drop = FALSE])
    }
    if (!gene_col %in% colnames(res_df)) {
        stop("`gene_col` \"", gene_col, "\" not in res_df.", call. = FALSE)
    }
    n_markers <- as.integer(n_markers)[1]
    if (is.na(n_markers) || n_markers < 1L) {
        stop("`n_markers` must be a positive integer.", call. = FALSE)
    }
    df <- res_df
    df[[gene_col]] <- as.character(df[[gene_col]])
    hit <- .match_cspa_genes(df[[gene_col]], cspa_symbols)
    df <- df[df[[gene_col]] %in% hit, , drop = FALSE]
    if (!nrow(df)) return(df)

    if (isTRUE(significant_only) && "sig" %in% colnames(df)) {
        sig_df <- df[!is.na(df$sig) & as.logical(df$sig), , drop = FALSE]
        if (nrow(sig_df)) df <- sig_df
    }

    score <- if ("rank_score" %in% colnames(df)) {
        abs(as.numeric(df$rank_score))
    } else if ("log2FoldChange" %in% colnames(df)) {
        abs(as.numeric(df$log2FoldChange))
    } else {
        rep(0, nrow(df))
    }
    score[!is.finite(score)] <- 0
    df <- df[order(-score, df[[gene_col]]), , drop = FALSE]
    # One row per gene symbol
    df <- df[!duplicated(df[[gene_col]]), , drop = FALSE]
    utils::head(df, n_markers)
}

#' Sample×gene VST heatmap for CSPA surface markers (all samples)
#'
#' @param mat Gene × sample numeric matrix (e.g. VST).
#' @param genes Character vector of gene IDs (rows of `mat`).
#' @param annotation_col Optional data.frame of column annotations for
#'   [pheatmap::pheatmap()] (`rownames` = sample IDs).
#' @param title Optional plot title.
#' @return A `pheatmap` object, or `NULL` if fewer than 2 genes/samples.
#' @export
plot_cspa_marker_heatmap <- function(mat,
                                     genes,
                                     annotation_col = NULL,
                                     title = NULL) {
    if (is.null(mat) || !length(genes)) return(NULL)
    genes <- as.character(genes)
    genes <- genes[genes %in% rownames(mat)]
    genes <- unique(genes)
    if (length(genes) < 1L || ncol(mat) < 2L) return(NULL)

    sub <- as.matrix(mat[genes, , drop = FALSE])
    storage.mode(sub) <- "double"
    # Row-scale for display
    sub <- t(scale(t(sub)))
    sub[!is.finite(sub)] <- 0

    if (!is.null(annotation_col)) {
        annotation_col <- as.data.frame(annotation_col)
        annotation_col <- annotation_col[colnames(sub), , drop = FALSE]
    }

    palette_length <- 100L
    my_color <- grDevices::colorRampPalette(c("Darkblue", "white", "red"))(palette_length)
    my_breaks <- c(
        seq(-3, 0, length.out = ceiling(palette_length / 2) + 1),
        seq(0.05, 3, length.out = floor(palette_length / 2))
    )
    pheatmap::pheatmap(
        sub,
        color = my_color,
        breaks = my_breaks,
        border_color = NA,
        cluster_rows = length(genes) > 1L,
        cluster_cols = TRUE,
        annotation_col = annotation_col,
        main = title,
        fontsize_row = 8,
        silent = TRUE
    )
}

#' CSPA surface-marker heatmaps from DEGs (per contrast + combined)
#'
#' For each DE contrast, selects up to `n_markers` CSPA genes among that
#' contrast's DEGs, then plots a **row-scaled VST heatmap across all
#' samples** (not just the contrast groups). Also builds a combined
#' heatmap from the union of per-contrast selections.
#'
#' Requires Suggests package `readxl` to read the CSPA workbook.
#'
#' @param dds A QC-filtered `DESeqDataSet`.
#' @param diffexp_out Output of [run_diffexp_pipeline()].
#' @param var Column name(s) in `colData(dds)` used for heatmap column
#'   annotation. Length 1–2; two factors are pasted for annotation.
#' @param n_markers Top CSPA markers per contrast. Default `25`.
#' @param cspa_file Optional override path to the CSPA Excel file.
#' @param gene_col Gene ID column in DE tables. Default `"gene"`.
#' @param significant_only Passed to [select_cspa_markers()]. Default `TRUE`.
#' @param save_dir Output directory (pipeline uses `deseq2/cspa/` if the caller
#'   already chose that folder). Default `getwd()`.
#' @param save_plot Logical; write PDFs. Default `TRUE`.
#' @return A list with `org`, `cspa_n`, `markers_by_comparison`,
#'   `combined_markers`, `heatmaps` (named list including `"combined"`),
#'   and `annotation_col`.
#' @export
run_cspa_marker_pipeline <- function(
    dds,
    diffexp_out,
    var,
    n_markers = 25,
    cspa_file = NULL,
    gene_col = "gene",
    significant_only = TRUE,
    save_dir = getwd(),
    save_plot = TRUE) {

    dds <- validate_dds(dds)
    save_dir <- validate_save_dir(save_dir)
    n_markers <- as.integer(n_markers)[1]
    if (is.na(n_markers) || n_markers < 1L) {
        stop("`n_markers` must be a positive integer.", call. = FALSE)
    }
    if (!is.character(var) || !length(var) || anyNA(var)) {
        stop("`var` must be a non-empty character vector of colData column names.",
             call. = FALSE)
    }
    missing_var <- setdiff(var, colnames(SummarizedExperiment::colData(dds)))
    if (length(missing_var)) {
        stop("`var` not in colData(dds): ", paste(missing_var, collapse = ", "),
             call. = FALSE)
    }
    by_comp <- diffexp_out$by_comparison
    if (is.null(by_comp) || !length(by_comp)) {
        stop("`diffexp_out$by_comparison` is empty; cannot select CSPA markers.",
             call. = FALSE)
    }

    org <- detect_organism(rownames(dds))
    message("Loading CSPA surface markers for organism=\"", org, "\" ...")
    cspa_symbols <- load_cspa_genes(organism = org, file = cspa_file)
    message("  ", length(cspa_symbols), " CSPA gene symbols loaded (all categories).")

    message("Computing blind VST for marker heatmaps ...")
    vsd <- quiet(safe_vst(dds, blind = TRUE))
    mat <- as.matrix(SummarizedExperiment::assay(vsd))
    storage.mode(mat) <- "double"

    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    anno <- cd[, var, drop = FALSE]
    rownames(anno) <- colnames(dds)

    markers_by_comparison <- list()
    heatmaps <- list()
    combined_genes <- character(0)

    message("Selecting top ", n_markers, " CSPA markers for ",
            length(by_comp), " contrast(s) ...")
    for (comp in names(by_comp)) {
        sel <- select_cspa_markers(
            by_comp[[comp]], cspa_symbols,
            n_markers = n_markers,
            gene_col = gene_col,
            significant_only = significant_only
        )
        markers_by_comparison[[comp]] <- sel
        genes <- as.character(sel[[gene_col]])
        combined_genes <- union(combined_genes, genes)

        comp_dir <- file.path(save_dir, comp)
        if (nrow(sel)) {
            write_output(sel, file.path(comp_dir, "cspa_markers.csv"))
        } else {
            message("  Skipping heatmap for \"", comp, "\": no CSPA DEGs selected.")
            heatmaps[[comp]] <- NULL
            next
        }
        hm <- plot_cspa_marker_heatmap(
            mat, genes, annotation_col = anno,
            title = paste0("CSPA surface markers: ", comp)
        )
        heatmaps[[comp]] <- hm
        if (save_plot && !is.null(hm)) {
            dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
            grDevices::pdf(file.path(comp_dir, "cspa_marker_heatmap.pdf"),
                           width = max(6, 0.25 * ncol(mat) + 2) * 0.65,
                           height = max(4, 0.22 * length(genes) + 2) * 0.65)
            grid::grid.newpage()
            grid::grid.draw(hm$gtable)
            grDevices::dev.off()
        }
        message("  ", comp, ": ", length(genes), " CSPA marker(s).")
    }

    combined_markers <- data.frame(
        gene = combined_genes,
        stringsAsFactors = FALSE
    )
    # Attach best |rank_score| / mean abs LFC across contrasts when available
    if (length(combined_genes)) {
        score_map <- stats::setNames(rep(NA_real_, length(combined_genes)), combined_genes)
        for (comp in names(markers_by_comparison)) {
            sel <- markers_by_comparison[[comp]]
            if (!nrow(sel)) next
            g <- as.character(sel[[gene_col]])
            sc <- if ("rank_score" %in% colnames(sel)) {
                abs(as.numeric(sel$rank_score))
            } else if ("log2FoldChange" %in% colnames(sel)) {
                abs(as.numeric(sel$log2FoldChange))
            } else {
                rep(0, length(g))
            }
            for (i in seq_along(g)) {
                if (!is.finite(score_map[[g[i]]]) || sc[i] > score_map[[g[i]]]) {
                    score_map[[g[i]]] <- sc[i]
                }
            }
        }
        combined_markers$best_abs_score <- unname(score_map[combined_markers$gene])
        combined_markers <- combined_markers[order(-combined_markers$best_abs_score,
                                                    combined_markers$gene), , drop = FALSE]
        write_output(combined_markers, file.path(save_dir, "cspa_markers_combined.csv"))

        hm_c <- plot_cspa_marker_heatmap(
            mat, combined_markers$gene, annotation_col = anno,
            title = "CSPA surface markers (combined contrasts)"
        )
        heatmaps[["combined"]] <- hm_c
        if (save_plot && !is.null(hm_c)) {
            grDevices::pdf(file.path(save_dir, "cspa_marker_heatmap_combined.pdf"),
                           width = max(6, 0.25 * ncol(mat) + 2) * 0.65,
                           height = max(4, 0.22 * nrow(combined_markers) + 2) * 0.65)
            grid::grid.newpage()
            grid::grid.draw(hm_c$gtable)
            grDevices::dev.off()
        }
        message("Combined CSPA heatmap: ", nrow(combined_markers), " unique marker(s).")
    } else {
        heatmaps[["combined"]] <- NULL
        message("No CSPA markers selected in any contrast.")
    }

    list(
        org = org,
        cspa_n = length(cspa_symbols),
        n_markers = n_markers,
        markers_by_comparison = markers_by_comparison,
        combined_markers = combined_markers,
        heatmaps = heatmaps,
        annotation_col = anno
    )
}
