#' Run GSVA on VST expression for requested MSigDB collections
#'
#' Sample-level gene-set variation analysis ([GSVA](https://doi.org/10.1186/1471-2105-14-7))
#' on blind VST counts, batch-corrected (`limma::removeBatchEffect()`) when `batch` is set.
#' Independent of clusterProfiler / fgsea GSEA (those remain ranked-list
#' tests per contrast). **HALLMARK** is always scored
#' when present in `collections`. Large collections (GO / KEGG / REACTOME /
#' PID) are scored in full and written to CSV; heatmaps show only the
#' top `n_heatmap` pathways most associated with `var` (Kruskal-Wallis on
#' GSVA scores), because a GSEA-style `padj` cutoff does not apply: every
#' set is scored, so filtering by significance still leaves hundreds of
#' rows. HALLMARK is always plotted in full.
#'
#' @param dds A QC-filtered `DESeqDataSet`.
#' @param var Column in `colData(dds)` used to annotate samples and to
#'   rank pathways for large-collection heatmaps.
#' @param batch Optional batch column(s), including the `pair` column of a
#'   paired design. When set, GSVA uses batch-corrected
#'   expression (protecting `var`) instead of uncorrected VST.
#' @param collections MSigDB aliases (see [validate_msigdb_collections()]).
#'   Default `c("HALLMARK", "GOBP", "GOCC", "GOMF", "KEGG", "REACTOME",
#'   "PID")`. HALLMARK is prepended by the Wald/LRT runners even when it
#'   is absent from the GSEA `collections` argument.
#' @param n_heatmap Pathways shown on heatmaps for **non-HALLMARK**
#'   collections. Default `50`. HALLMARK always plots every scored set.
#' @param min_size Minimum overlapping genes per set (passed to GSVA).
#'   Default `10`.
#' @param gene_sets Optional named list of collections, each a named list
#'   of gene-symbol vectors. When supplied, `collections` is ignored and
#'   MSigDB is not queried (used in tests).
#' @param save_dir Output directory. Default `getwd()`.
#' @param save_plot Logical; write per-collection heatmap PDFs. Default
#'   `TRUE`.
#' @return A list with `org`, `collections`, `scores` (named list of
#'   pathway × sample matrices), `stats` (per-collection ranking tables),
#'   `heatmaps` (named list of `pheatmap` objects), `n_heatmap`,
#'   `min_size`, and `var`.
#' @export
run_gsva_pipeline <- function(
    dds,
    var,
    collections = c("HALLMARK", "GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "PID"),
    n_heatmap = 50,
    min_size = 10,
    gene_sets = NULL,
    batch = NULL,
    save_dir = getwd(),
    save_plot = TRUE) {

    if (!requireNamespace("GSVA", quietly = TRUE)) {
        stop("Package \"GSVA\" is required for GSVA analysis. ",
             "Install with BiocManager::install(\"GSVA\").", call. = FALSE)
    }
    dds <- validate_dds(dds)
    dds <- validate_var(dds, var)
    if (!is.null(batch)) dds <- validate_batch(dds, batch)
    n_heatmap <- as.integer(n_heatmap)[1]
    min_size <- as.integer(min_size)[1]
    if (is.na(n_heatmap) || n_heatmap < 1L) {
        stop("`n_heatmap` must be a positive integer.", call. = FALSE)
    }
    if (is.na(min_size) || min_size < 1L) {
        stop("`min_size` must be a positive integer.", call. = FALSE)
    }
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)

    if (is.null(gene_sets)) {
        collections <- validate_msigdb_collections(collections)
        collections <- unique(c("HALLMARK", collections))
        org <- detect_organism(rownames(dds))
        message("Fetching ", length(collections),
                " MSigDB collection(s) for GSVA (organism=\"", org, "\") ...")
        gene_sets <- stats::setNames(
            lapply(collections, function(coll) import_msigdbr(org, coll)),
            collections
        )
    } else {
        if (!is.list(gene_sets) || is.null(names(gene_sets)) ||
                any(!nzchar(names(gene_sets)))) {
            stop("`gene_sets` must be a named list of collections ",
                 "(each a named list of gene vectors).", call. = FALSE)
        }
        collections <- names(gene_sets)
        org <- tryCatch(detect_organism(rownames(dds)),
                        error = function(e) NA_character_)
    }

    message("Computing sample-level expression for GSVA ...")
    expr <- sample_level_expr_matrix(dds, var = var, batch = batch)
    expr_correction <- attr(expr, "expr_correction")
    storage.mode(expr) <- "double"
    group <- SummarizedExperiment::colData(dds)[[var]]
    if (!is.factor(group)) group <- factor(group)

    scores <- list()
    stats_by <- list()
    heatmaps <- list()
    for (coll in collections) {
        sets <- gene_sets[[coll]]
        if (is.null(sets) || !length(sets)) {
            warning("Skipping GSVA collection \"", coll, "\": empty gene sets.",
                    call. = FALSE)
            next
        }
        message("  GSVA: ", coll, " (", length(sets), " gene set(s)) ...")
        sc <- tryCatch(
            .gsva_score_matrix(expr, sets, min_size = min_size),
            error = function(e) {
                warning("GSVA failed for \"", coll, "\": ",
                        conditionMessage(e), call. = FALSE)
                NULL
            }
        )
        if (is.null(sc) || !nrow(sc) || !ncol(sc)) next

        scores[[coll]] <- sc
        st <- .gsva_pathway_stats(sc, group)
        stats_by[[coll]] <- st
        write_output(
            data.frame(pathway = rownames(sc), sc, check.names = FALSE,
                       stringsAsFactors = FALSE),
            file.path(save_dir, coll, "gsva_scores.csv")
        )
        write_output(st, file.path(save_dir, coll, "gsva_pathway_stats.csv"))

        plot_all <- identical(coll, "HALLMARK")
        keep <- .gsva_heatmap_pathways(st, n = n_heatmap, plot_all = plot_all)
        hm <- plot_gsva_heatmap(
            sc, group = group, var_label = var,
            keep = keep, collection = coll
        )
        heatmaps[[coll]] <- hm
        if (save_plot && !is.null(hm)) {
            save_pheatmap(hm, file.path(save_dir, coll, "gsva_heatmap.pdf"))
        }
        message("    ", nrow(sc), " scored set(s); heatmap shows ",
                length(keep),
                if (plot_all) " (all HALLMARK)." else
                    paste0(" (top ", n_heatmap, " by group association)."))
    }

    if (!length(scores)) {
        stop("GSVA produced no scores for any collection.", call. = FALSE)
    }

    list(
        org = org,
        collections = names(scores),
        scores = scores,
        stats = stats_by,
        heatmaps = heatmaps,
        n_heatmap = n_heatmap,
        min_size = min_size,
        var = var,
        expr_correction = expr_correction
    )
}

#' GSVA score matrix (pathways × samples) for one collection
#' @noRd
.gsva_score_matrix <- function(expr, gene_sets, min_size = 10L) {
    ns <- asNamespace("GSVA")
    if (exists("gsvaParam", envir = ns, inherits = FALSE)) {
        param <- GSVA::gsvaParam(
            exprData = expr,
            geneSets = gene_sets,
            minSize = min_size,
            kcdf = "Gaussian"
        )
        out <- GSVA::gsva(param, verbose = FALSE)
    } else {
        out <- GSVA::gsva(
            expr, gene_sets,
            method = "gsva",
            kcdf = "Gaussian",
            min.sz = min_size,
            verbose = FALSE
        )
    }
    mat <- as.matrix(out)
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
        stop("GSVA did not return a named score matrix.", call. = FALSE)
    }
    mat
}

#' Per-pathway variance and Kruskal-Wallis p-value vs grouping factor
#' @noRd
.gsva_pathway_stats <- function(scores, group) {
    group <- factor(group)
    pathways <- rownames(scores)
    rows <- lapply(pathways, function(pw) {
        y <- as.numeric(scores[pw, ])
        ok <- is.finite(y) & !is.na(group)
        v <- if (sum(ok) >= 2L) stats::var(y[ok]) else NA_real_
        p <- NA_real_
        g <- droplevels(group[ok])
        if (nlevels(g) >= 2L && sum(ok) >= 3L) {
            p <- tryCatch(
                stats::kruskal.test(y[ok] ~ g)$p.value,
                error = function(e) NA_real_
            )
        }
        data.frame(
            pathway = pw,
            variance = v,
            p_group = p,
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)
    df$padj_group <- stats::p.adjust(df$p_group, method = "BH")
    df <- df[order(df$p_group, -df$variance, na.last = TRUE), , drop = FALSE]
    rownames(df) <- NULL
    df
}

#' Pathway names to show on a GSVA heatmap
#' @noRd
.gsva_heatmap_pathways <- function(stats_df, n = 50L, plot_all = FALSE) {
    if (is.null(stats_df) || !nrow(stats_df)) return(character(0))
    if (isTRUE(plot_all)) return(as.character(stats_df$pathway))
    n <- min(as.integer(n)[1], nrow(stats_df))
    if (any(is.finite(stats_df$p_group))) {
        ord <- order(stats_df$p_group, -stats_df$variance, na.last = TRUE)
    } else {
        ord <- order(-stats_df$variance, na.last = TRUE)
    }
    as.character(stats_df$pathway[ord[seq_len(n)]])
}

#' Strip the MSigDB collection prefix from pathway names
#' @noRd
.gsva_display_names <- function(pathways, collection) {
    spec <- tryCatch(
        resolve_msigdb_alias(collection),
        error = function(e) NULL
    )
    lab <- as.character(pathways)
    if (!is.null(spec) && !is.null(spec$prefix) && nzchar(spec$prefix)) {
        lab <- sub(paste0("^", spec$prefix), "", lab)
    }
    lab <- gsub("_", " ", lab)
    lab
}

#' Wrap pathway labels at `width` characters (word-aware)
#' @noRd
.gsva_wrap_labels <- function(labels, width = 40L) {
    width <- as.integer(width)[1]
    if (is.na(width) || width < 8L) width <- 40L
    vapply(as.character(labels), function(lab) {
        if (is.na(lab) || !nzchar(lab) || nchar(lab) <= width) return(lab)
        words <- unlist(strsplit(lab, "\\s+"), use.names = FALSE)
        words <- words[nzchar(words)]
        if (!length(words)) return(lab)
        lines <- character(0)
        cur <- ""
        flush_long <- function(w) {
            out <- character(0)
            while (nchar(w) > width) {
                out <- c(out, substr(w, 1L, width))
                w <- substr(w, width + 1L, nchar(w))
            }
            list(parts = out, rest = w)
        }
        for (w in words) {
            if (nchar(w) > width) {
                if (nzchar(cur)) {
                    lines <- c(lines, cur)
                    cur <- ""
                }
                split_w <- flush_long(w)
                lines <- c(lines, split_w$parts)
                w <- split_w$rest
                if (!nzchar(w)) next
            }
            trial <- if (nzchar(cur)) paste(cur, w) else w
            if (nchar(trial) <= width) {
                cur <- trial
            } else {
                if (nzchar(cur)) lines <- c(lines, cur)
                cur <- w
            }
        }
        if (nzchar(cur)) lines <- c(lines, cur)
        paste(lines, collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
}

#' Shrink the pheatmap row-name column to the longest wrapped line
#' @noRd
.gsva_fit_heatmap_rownames <- function(hm, fontsize_row) {
    if (is.null(hm) || is.null(hm$gtable)) return(hm)
    gt <- hm$gtable
    hit <- which(gt$layout$name == "row_names")
    if (!length(hit)) return(hm)
    grob <- gt$grobs[[hit[[1]]]]
    labs <- grob$label
    if (is.null(labs) && inherits(grob, "gTree") && length(grob$children)) {
        for (ch in grob$children) {
            if (!is.null(ch$label) && length(ch$label)) {
                labs <- ch$label
                break
            }
        }
    }
    if (is.null(labs) || !length(labs)) return(hm)
    lines <- unlist(strsplit(as.character(labs), "\n", fixed = TRUE),
                    use.names = FALSE)
    lines <- lines[!is.na(lines) & nzchar(lines)]
    if (!length(lines)) return(hm)
    longest <- lines[[which.max(nchar(lines))]]
    new_w <- grid::grobWidth(
        grid::textGrob(longest, gp = grid::gpar(fontsize = fontsize_row))
    ) + grid::unit(6, "bigpts")
    col <- gt$layout$l[[hit[[1]]]]
    gt$widths[col] <- new_w
    hm$gtable <- gt
    hm
}

#' Sample×pathway GSVA heatmap (row-scaled)
#'
#' Labels wrap at 40 characters. If any plotted name wraps, row font size
#' is reduced; [gtable_page_size()] then gives each row enough height for the
#' wrapped label. Note that the returned gtable's matrix row and column are
#' `npc`, so it has no natural size -- never scale its `widths`/`heights`, which
#' would push the layout past whatever page it is drawn on and clip it.
#'
#' @param scores Pathway × sample numeric matrix.
#' @param group Sample grouping factor (same order as columns of `scores`).
#' @param var_label Annotation legend title. Default `"group"`.
#' @param keep Pathway names to plot. Default all rows.
#' @param collection Collection alias (for prefix stripping). Default `NULL`.
#' @return A `pheatmap` object, or `NULL` if fewer than 2 pathways/samples.
#' @export
plot_gsva_heatmap <- function(scores, group, var_label = "group",
                              keep = NULL, collection = NULL) {
    if (is.null(scores) || !length(scores)) return(NULL)
    mat <- as.matrix(scores)
    if (!is.null(keep) && length(keep)) {
        keep <- intersect(as.character(keep), rownames(mat))
        mat <- mat[keep, , drop = FALSE]
    }
    if (nrow(mat) < 1L || ncol(mat) < 2L) return(NULL)

    display <- .gsva_display_names(rownames(mat), collection)
    # Keep labels unique after prefix stripping.
    if (anyDuplicated(display)) {
        display <- make.unique(display, sep = " ")
    }
    display <- .gsva_wrap_labels(display, width = 40L)
    rownames(mat) <- display
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0

    group <- factor(group)
    if (length(group) != ncol(mat)) {
        stop("`group` length must match ncol(scores).", call. = FALSE)
    }
    anno <- data.frame(group = group, row.names = colnames(mat),
                       check.names = FALSE)
    colnames(anno) <- var_label

    palette_length <- 100L
    my_color <- grDevices::colorRampPalette(c("Darkblue", "white", "red"))(palette_length)
    my_breaks <- c(
        seq(-3, 0, length.out = ceiling(palette_length / 2) + 1),
        seq(0.05, 3, length.out = floor(palette_length / 2))
    )
    fs_row <- if (nrow(mat) > 40L) 6 else if (nrow(mat) > 25L) 7 else 8
    if (any(grepl("\n", display, fixed = TRUE))) {
        fs_row <- max(4, fs_row - 2L)
    }
    hm <- pheatmap::pheatmap(
        mat,
        color = my_color,
        breaks = my_breaks,
        border_color = NA,
        cluster_rows = nrow(mat) > 1L,
        cluster_cols = TRUE,
        annotation_col = anno,
        fontsize_row = fs_row,
        fontsize_col = 8,
        silent = TRUE
    )
    .gsva_fit_heatmap_rownames(hm, fontsize_row = fs_row)
}
