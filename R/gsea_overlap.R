#' Parse a comma-separated leading-edge gene string into a character vector
#' @param x Character scalar (possibly NA / empty).
#' @return Character vector of gene symbols (possibly length 0).
#' @noRd
.parse_leading_edge <- function(x) {
    if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(x[[1]])) {
        return(character(0))
    }
    genes <- strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]]
    genes <- trimws(genes)
    genes[nzchar(genes)]
}

#' Clean pathway labels into tokens for word-cloud annotation
#' @noRd
.pathway_terms_for_cloud <- function(labels) {
    labs <- as.character(labels)
    # Strip "COLLECTION | " prefix when present
    labs <- sub("^[^|]+\\|\\s*", "", labs)
    # Underscores → spaces so count_words can tokenise
    gsub("_", " ", labs, fixed = TRUE)
}

#' Pool significant GSEA pathways for one contrast across all collections
#'
#' @param by_coll Named list of GSEA result data frames (one per collection)
#'   as in `gsea_out$by_comparison[[comp]]`.
#' @param alpha Adjusted p-value threshold. Default `0.05`.
#' @return A data frame of significant pathways (possibly 0 rows) with
#'   columns `pathway`, `collection`, `NES`, `padj`, `size`, `leadingEdge`,
#'   `n_leadingEdge`, `label` (unique display key).
#' @export
pool_sig_gsea_pathways <- function(by_coll, alpha = 0.05) {
    if (is.null(by_coll) || !length(by_coll)) {
        return(data.frame(
            pathway = character(0), collection = character(0),
            NES = numeric(0), padj = numeric(0), size = integer(0),
            leadingEdge = character(0), n_leadingEdge = integer(0),
            label = character(0),
            stringsAsFactors = FALSE
        ))
    }
    rows <- list()
    for (coll in names(by_coll)) {
        df <- by_coll[[coll]]
        if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) next
        need <- c("pathway", "NES", "padj", "leadingEdge")
        missing <- setdiff(need, colnames(df))
        if (length(missing) > 0L) {
            stop("GSEA table for collection \"", coll, "\" missing column(s): ",
                 paste(missing, collapse = ", "), call. = FALSE)
        }
        keep <- !is.na(df$padj) & df$padj < alpha
        if (!any(keep)) next
        sub <- df[keep, , drop = FALSE]
        size_col <- if ("size" %in% colnames(sub)) as.integer(sub$size) else NA_integer_
        le <- as.character(sub$leadingEdge)
        n_le <- vapply(le, function(s) length(.parse_leading_edge(s)), integer(1))
        rows[[length(rows) + 1L]] <- data.frame(
            pathway = as.character(sub$pathway),
            collection = coll,
            NES = as.numeric(sub$NES),
            padj = as.numeric(sub$padj),
            size = size_col,
            leadingEdge = le,
            n_leadingEdge = n_le,
            stringsAsFactors = FALSE
        )
    }
    if (!length(rows)) {
        return(pool_sig_gsea_pathways(NULL, alpha = alpha))
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    raw_label <- paste(out$collection, out$pathway, sep = " | ")
    out$label <- make.unique(raw_label, sep = " #")
    out
}

#' Pairwise Jaccard index of leading-edge gene sets
#'
#' @param gene_sets Named list of character vectors (pathway → genes).
#' @return Square numeric matrix of Jaccard indices (dimnames = names).
#' @export
jaccard_matrix <- function(gene_sets) {
    if (is.null(gene_sets) || !length(gene_sets)) {
        return(matrix(numeric(0), nrow = 0L, ncol = 0L))
    }
    nms <- names(gene_sets)
    if (is.null(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
        stop("`gene_sets` must be a uniquely named list.", call. = FALSE)
    }
    sets <- lapply(gene_sets, function(g) unique(as.character(g)))
    n <- length(sets)
    mat <- matrix(0, nrow = n, ncol = n, dimnames = list(nms, nms))
    for (i in seq_len(n)) {
        mat[i, i] <- 1
        if (i == n) next
        for (j in seq.int(i + 1L, n)) {
            a <- sets[[i]]
            b <- sets[[j]]
            if (!length(a) && !length(b)) {
                jac <- 0
            } else {
                inter <- length(intersect(a, b))
                uni <- length(union(a, b))
                jac <- if (uni > 0L) inter / uni else 0
            }
            mat[i, j] <- jac
            mat[j, i] <- jac
        }
    }
    mat
}

#' Fallback pheatmap Jaccard plot (no pathway name labels)
#' @noRd
.plot_pathway_jaccard_pheatmap <- function(jac) {
    n <- nrow(jac)
    jac_dist <- stats::as.dist(1 - jac)
    can_cluster <- !anyNA(jac_dist) && n >= 2L
    display_mat <- jac
    diag(display_mat) <- NA
    colours <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
    lo <- min(display_mat, na.rm = TRUE)
    hi <- max(display_mat, na.rm = TRUE)
    if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
        lo <- 0
        hi <- 1
    }
    breaks <- seq(lo, hi, length.out = length(colours) + 1)
    cell_size <- pheatmap_cell_size(n, min_size = 6, max_size = 18, target_total = 360)
    # Match ComplexHeatmap footprint (~40% smaller).
    cell_size <- cell_size * 0.6
    p <- pheatmap::pheatmap(
        display_mat,
        cluster_rows = can_cluster,
        cluster_cols = can_cluster,
        clustering_distance_rows = if (can_cluster) jac_dist else NA,
        clustering_distance_cols = if (can_cluster) jac_dist else NA,
        treeheight_row = if (can_cluster) 25 else 0,
        treeheight_col = if (can_cluster) 25 else 0,
        color = colours,
        breaks = breaks,
        display_numbers = FALSE,
        show_rownames = FALSE,
        show_colnames = FALSE,
        na_col = "grey20",
        cellwidth = cell_size,
        cellheight = cell_size,
        silent = TRUE
    )
    attr(p, "jaccard_matrix") <- jac
    attr(p, "plot_engine") <- "pheatmap"
    attr(p, "clusters") <- NULL
    body_in <- (cell_size * n) / 72
    attr(p, "heatmap_body_in") <- body_in
    attr(p, "plot_size") <- list(
        width = body_in + 0.7,
        height = body_in + 0.5
    )
    p
}

#' simplifyEnrichment-style Jaccard heatmap + side word clouds
#'
#' Clusters pathways with `binary_cut` on the Jaccard similarity matrix,
#' draws a ComplexHeatmap without row/column names, and annotates clusters
#' with pathway-name word clouds (via `anno_word_cloud`).
#'
#' @param jac Square Jaccard matrix from [jaccard_matrix()].
#' @param term_labels Optional character vector (length = nrow(jac)) of
#'   pathway names used for word clouds. Default uses `rownames(jac)`.
#' @return A drawable object (`HeatmapList` when simplifyEnrichment is
#'   available, else a `pheatmap`), with attributes `jaccard_matrix`,
#'   `clusters`, `plot_engine`, and `plot_size`.
#' @export
plot_pathway_jaccard_heatmap <- function(jac, term_labels = NULL) {
    if (is.null(jac) || !is.matrix(jac) || nrow(jac) < 2L) {
        return(NULL)
    }
    n <- nrow(jac)
    if (is.null(term_labels)) {
        term_labels <- rownames(jac)
    }
    if (length(term_labels) != n) {
        stop("`term_labels` must have length nrow(jac).", call. = FALSE)
    }
    terms_clean <- .pathway_terms_for_cloud(term_labels)

    use_se <- requireNamespace("simplifyEnrichment", quietly = TRUE) &&
        requireNamespace("ComplexHeatmap", quietly = TRUE) &&
        requireNamespace("circlize", quietly = TRUE) &&
        requireNamespace("grid", quietly = TRUE)

    if (!use_se) {
        message("Packages simplifyEnrichment / ComplexHeatmap not available; ",
                "using pheatmap fallback (no side word clouds).")
        return(.plot_pathway_jaccard_pheatmap(jac))
    }

    cl_bin <- tryCatch(
        simplifyEnrichment::cluster_terms(jac, method = "binary_cut", verbose = FALSE),
        error = function(e) {
            warning("binary_cut clustering failed (", conditionMessage(e),
                    "); falling back to pheatmap.", call. = FALSE)
            NULL
        }
    )
    if (is.null(cl_bin) || length(cl_bin) != n) {
        return(.plot_pathway_jaccard_pheatmap(jac))
    }
    cl_bin <- as.vector(cl_bin)
    k_bin <- length(unique(cl_bin))

    # Order by binary_cut blocks (fallback if hclust fails)
    od <- unlist(lapply(sort(unique(cl_bin)), function(le) {
        idx <- which(cl_bin == le)
        if (length(idx) <= 1L) return(idx)
        mm <- jac[idx, idx, drop = FALSE]
        idx[stats::hclust(stats::dist(mm))$order]
    }), use.names = FALSE)

    pos <- jac[jac > 0 & is.finite(jac)]
    hi <- if (length(pos)) {
        as.numeric(stats::quantile(pos, 0.975, na.rm = TRUE))
    } else {
        1
    }
    if (!is.finite(hi) || hi <= 0) hi <- 1
    col_fun <- circlize::colorRamp2(c(0, hi), c("white", "red"))

    cell <- pheatmap_cell_size(n, min_size = 4, max_size = 12, target_total = 320)
    # Square heatmap body; dendrogram / legend / word-cloud links sit
    # *outside* that body. Device size must include those extras or the
    # PNG frame crops them (which looked like an HTML "frame").
    size_scale <- 0.6
    body_in <- size_scale * max(3, min(5.5, (cell * n) / 72))
    body_u <- grid::unit(body_in, "in")
    anno_in <- size_scale * 2.6
    legend_in <- 0.55
    dend_in <- 0.45
    link_in <- 0.7

    # Dendrogram + partition lines must share one clustering. Previously
    # binary_cut drew the splits while hclust drew the trees, so lines cut
    # through branches. Use hclust for the tree, take k from binary_cut,
    # then cutree() so splits / word clouds match the dendrogram.
    dmat <- 1 - jac
    dmat[dmat < 0] <- 0
    storage.mode(dmat) <- "double"
    jac_dist <- stats::as.dist(dmat)
    hc <- tryCatch(
        stats::hclust(jac_dist, method = "average"),
        error = function(e) NULL
    )
    has_dend <- !is.null(hc)
    if (has_dend) {
        k_use <- max(1L, min(as.integer(k_bin), n))
        cl <- as.integer(stats::cutree(hc, k = k_use))
        od <- hc$order
    } else {
        cl <- cl_bin
    }

    ht_args <- list(
        matrix = jac,
        name = "Jaccard",
        col = col_fun,
        show_row_names = FALSE,
        show_column_names = FALSE,
        border = "#404040",
        row_title = NULL,
        width = body_u,
        height = body_u,
        use_raster = TRUE,
        heatmap_legend_param = list(
            title = "Jaccard",
            legend_height = grid::unit(body_in * 0.55, "in")
        )
    )
    if (has_dend) {
        ht_args$cluster_rows <- hc
        ht_args$cluster_columns <- hc
        ht_args$show_row_dend <- TRUE
        ht_args$show_column_dend <- TRUE
        ht_args$row_dend_width <- grid::unit(dend_in, "in")
        ht_args$column_dend_height <- grid::unit(dend_in, "in")
    } else {
        ht_args$cluster_rows <- FALSE
        ht_args$cluster_columns <- FALSE
        ht_args$show_row_dend <- FALSE
        ht_args$show_column_dend <- FALSE
        ht_args$row_order <- od
        ht_args$column_order <- od
        dend_in <- 0
    }
    ht <- do.call(ComplexHeatmap::Heatmap, ht_args)

    # Side word clouds from the same clusters as the dendrogram splits
    cl_f <- factor(cl, levels = sort(unique(cl)))
    align_to <- split(seq_len(n), cl_f)
    term_list <- lapply(align_to, function(idx) terms_clean[idx])
    min_term <- max(2L, round(n * 0.01))
    keep_ann <- names(align_to)[vapply(align_to, length, integer(1)) >= min_term]
    has_anno <- length(keep_ann) > 0L
    if (has_anno) {
        ht <- ht + ComplexHeatmap::rowAnnotation(
            keywords = simplifyEnrichment::anno_word_cloud(
                align_to = align_to[keep_ann],
                term = term_list[keep_ann],
                max_words = 10,
                fontsize_range = c(3, 9),
                exclude_words = c(
                    "pathway", "of", "the", "and", "in", "to", "for",
                    "via", "by", "a", "an", "or", "go", "bp", "mf", "cc",
                    "hallmark", "reactome", "kegg", "pid", "gobp", "gocc", "gomf"
                )
            ),
            width = grid::unit(anno_in, "in"),
            show_annotation_name = FALSE
        )
    } else {
        anno_in <- 0
        link_in <- 0
    }

    plot_size <- list(
        # left: dend + legend; center: square body; right: link + clouds
        width = dend_in + legend_in + body_in + link_in + anno_in + 0.25,
        height = dend_in + body_in + 0.3
    )
    attr(ht, "jaccard_matrix") <- jac
    attr(ht, "plot_engine") <- "complexheatmap"
    attr(ht, "clusters") <- cl
    attr(ht, "row_order") <- od
    attr(ht, "heatmap_body_in") <- body_in
    attr(ht, "plot_size") <- plot_size
    ht
}

#' Draw a GSEA-overlap heatmap object to the current device
#'
#' @param heat Object from [plot_pathway_jaccard_heatmap()].
#' @return Invisibly `NULL`.
#' @export
draw_pathway_jaccard_heatmap <- function(heat) {
    if (is.null(heat)) return(invisible(NULL))
    engine <- attr(heat, "plot_engine")
    if (identical(engine, "complexheatmap") ||
            inherits(heat, "Heatmap") || inherits(heat, "HeatmapList")) {
        if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
            suppressMessages(ComplexHeatmap::ht_opt(message = FALSE))
        }
        ComplexHeatmap::draw(
            heat,
            newpage = TRUE,
            # Keep legend off the word-cloud column (was overlapping / crowding).
            heatmap_legend_side = "left",
            annotation_legend_side = "left",
            padding = grid::unit(c(2, 2, 2, 2), "mm")
        )
        # Partition lines
        cl <- attr(heat, "clusters")
        od <- attr(heat, "row_order")
        if (!is.null(cl) && !is.null(od)) {
            cl_ord <- factor(cl[od], levels = unique(cl[od]))
            tb <- table(cl_ord)
            ncl <- length(cl_ord)
            x <- cumsum(c(0, as.numeric(tb))) / ncl
            tryCatch(
                ComplexHeatmap::decorate_heatmap_body("Jaccard", {
                    grid::grid.rect(gp = grid::gpar(fill = NA, col = "#404040"))
                    grid::grid.segments(x, 0, x, 1, default.units = "npc",
                                       gp = grid::gpar(col = "#404040"))
                    grid::grid.segments(0, 1 - x, 1, 1 - x, default.units = "npc",
                                       gp = grid::gpar(col = "#404040"))
                }),
                error = function(e) invisible(NULL)
            )
        }
    } else if (!is.null(heat$gtable)) {
        grid::grid.newpage()
        grid::grid.draw(heat$gtable)
    } else {
        stop("Unrecognised heatmap object.", call. = FALSE)
    }
    invisible(NULL)
}

#' Save GSEA-overlap heatmap PDF
#' @noRd
.save_pathway_jaccard_heatmap <- function(heat, path) {
    if (is.null(heat)) return(invisible(FALSE))
    sz <- attr(heat, "plot_size")
    if (is.null(sz)) sz <- list(width = 8, height = 6)
    engine <- attr(heat, "plot_engine")
    if (identical(engine, "pheatmap")) {
        save_pheatmap(heat, path)
        return(invisible(TRUE))
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(path, width = sz$width, height = sz$height)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_pathway_jaccard_heatmap(heat)
    invisible(TRUE)
}

#' Build per-contrast GSEA leading-edge Jaccard overlap artifacts
#'
#' For each contrast in `gsea_out$by_comparison`, pools significant pathways
#' (`padj < alpha`) across all collections, computes pairwise Jaccard
#' similarity of leading-edge gene sets, clusters pathways with
#' simplifyEnrichment `binary_cut`, and writes a heatmap (no pathway name
#' labels; side word clouds of pathway-name tokens) plus CSVs under
#' `save_dir/<contrast>/gsea_overlap/`.
#'
#' The on-disk Jaccard CSV always covers **all** significant pathways. The
#' heatmap is limited to the top `max_heatmap_pathways` by `|NES|`.
#'
#' @param gsea_out Output of [run_gsea_pipeline()].
#' @param alpha Pathway significance threshold (`padj`). Default `0.05`.
#' @param max_heatmap_pathways Maximum pathways shown in the Jaccard
#'   heatmap (ranked by `|NES|`). Default `40`. Set `Inf` to plot all.
#' @param save_dir Directory for outputs. Default `getwd()`.
#' @param save_plot Logical. Default `TRUE`.
#' @return A list with `by_comparison` (named list of per-contrast results:
#'   `sig_pathways`, `jaccard`, `heatmap`, `clusters`, `n_sig`,
#'   `n_heatmap`) and `alpha`.
#' @export
run_gsea_overlap_pipeline <- function(
    gsea_out,
    alpha = 0.05,
    max_heatmap_pathways = 40L,
    save_dir = getwd(),
    save_plot = TRUE) {

    if (is.null(gsea_out) || is.null(gsea_out$by_comparison) ||
            !length(gsea_out$by_comparison)) {
        stop("`gsea_out` must contain a non-empty `by_comparison` list.",
             call. = FALSE)
    }
    alpha <- validate_alpha(alpha)
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)
    if (!is.numeric(max_heatmap_pathways) || length(max_heatmap_pathways) != 1L ||
            is.na(max_heatmap_pathways) || max_heatmap_pathways < 2) {
        stop("`max_heatmap_pathways` must be a number >= 2 (or Inf).", call. = FALSE)
    }

    by_comparison <- list()
    for (comp in names(gsea_out$by_comparison)) {
        message("GSEA leading-edge overlap: ", comp, " ...")
        sig <- pool_sig_gsea_pathways(gsea_out$by_comparison[[comp]], alpha = alpha)
        gene_sets <- stats::setNames(
            lapply(sig$leadingEdge, .parse_leading_edge),
            sig$label
        )
        jac_full <- if (length(gene_sets) >= 2L) jaccard_matrix(gene_sets) else NULL

        # Heatmap subset: top by |NES| to keep HTML embeds manageable.
        heat_idx <- seq_len(nrow(sig))
        if (nrow(sig) > 0L && is.finite(max_heatmap_pathways) &&
                nrow(sig) > max_heatmap_pathways) {
            ord <- order(abs(sig$NES), decreasing = TRUE)
            heat_idx <- ord[seq_len(as.integer(max_heatmap_pathways))]
        }
        heat_labels <- if (length(heat_idx)) sig$label[heat_idx] else character(0)
        heat_terms <- if (length(heat_idx)) sig$pathway[heat_idx] else character(0)
        jac_heat <- if (!is.null(jac_full) && length(heat_labels) >= 2L) {
            jac_full[heat_labels, heat_labels, drop = FALSE]
        } else {
            NULL
        }
        heat <- plot_pathway_jaccard_heatmap(jac_heat, term_labels = heat_terms)
        clusters <- if (!is.null(heat)) attr(heat, "clusters") else NULL

        # Alongside GSEA outputs: clusterprofiler/<comp>/gsea_overlap/
        comp_dir <- file.path(save_dir, comp, "gsea_overlap")
        dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

        write_output(
            sig[, c("pathway", "collection", "NES", "padj", "size",
                    "n_leadingEdge", "label", "leadingEdge"), drop = FALSE],
            file.path(comp_dir, paste0(comp, "_sig_pathways.csv")),
            format = "csv"
        )
        if (!is.null(jac_full)) {
            jac_df <- as.data.frame(jac_full)
            jac_df$pathway <- rownames(jac_df)
            jac_df <- jac_df[, c("pathway", setdiff(colnames(jac_df), "pathway"))]
            write_output(jac_df, file.path(comp_dir, paste0(comp, "_pathway_jaccard.csv")),
                         format = "csv")
        }
        if (!is.null(clusters) && length(heat_labels)) {
            cl_df <- data.frame(
                label = heat_labels,
                pathway = heat_terms,
                cluster = as.integer(clusters),
                stringsAsFactors = FALSE
            )
            write_output(cl_df, file.path(comp_dir, paste0(comp, "_pathway_clusters.csv")),
                         format = "csv")
        }

        heat_pdf <- file.path(comp_dir, paste0(comp, "_pathway_jaccard_heatmap.pdf"))
        if (save_plot && !is.null(heat)) {
            .save_pathway_jaccard_heatmap(heat, heat_pdf)
        }

        by_comparison[[comp]] <- list(
            sig_pathways = sig,
            jaccard = jac_full,
            # Keep heatmap for HTML; drop after page render via gc if needed.
            heatmap = heat,
            heatmap_pdf = if (file.exists(heat_pdf)) heat_pdf else NULL,
            clusters = clusters,
            n_sig = nrow(sig),
            n_heatmap = if (is.null(jac_heat)) 0L else nrow(jac_heat),
            output_dir = comp_dir
        )
    }

    list(
        by_comparison = by_comparison,
        alpha = alpha,
        max_heatmap_pathways = max_heatmap_pathways,
        save_dir = save_dir
    )
}
