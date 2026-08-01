#' Run the sample-similarity / PCA stage of the pipeline
#'
#' Variance-stabilizes the counts, optionally applies batch correction to the
#' VST matrix for *visualization only* (the batch term should instead be
#' included in the DESeq2 `design` for the actual differential expression
#' model), then separately renders:
#' \itemize{
#'   \item a **sample-similarity** heatmap (pairwise Spearman correlation
#'     between samples in VST/rlog space) -- a square (aspect ratio 1)
#'     heatmap coloured on a gradient scaled to the *actual* min/max of the
#'     correlation matrix, plus the underlying correlation matrix as CSV;
#'   \item an interactive (hover-for-sample-name) and a static **PCA**
#'     scatterplot (PC1 vs PC2, square panel), a scree plot of variance
#'     explained per PC, and a top-loadings plot/CSV showing the genes that
#'     contribute most to each of the leading PCs.
#' }
#' When `batch` is set and `batch_method != "none"`, the PCA/similarity
#' computations use the batch-corrected matrix (visualization only), and an
#' additional *uncorrected* PCA plot (interactive + static), scree plot, and
#' top-loadings plot/CSV are returned/saved for comparison.
#'
#' @param dds A `DESeqDataSet` object (typically post-[run_qc()]).
#' @param var Column name in `colData(dds)` used to colour points/rows by
#'   group.
#' @param batch Column name in `colData(dds)` to correct for, or `NULL`
#'   (default) to skip batch correction.
#' @param batch_method One of `"limma"` (`limma::removeBatchEffect()`),
#'   `"combat"` (`sva::ComBat()`), or `"none"`. Ignored if `batch` is `NULL`.
#'   Default `"limma"`.
#' @param transform One of `"vst"` (default, `DESeq2::vst()`) or `"rlog"`
#'   (`DESeq2::rlog()`) -- the variance-stabilizing transform used for PCA
#'   and the similarity heatmap.
#' @param pals Optional character vector of colours, one per level of `var`.
#'   Default `NULL` uses `ggplot2`'s default discrete palette.
#' @param shape_var Optional column name in `colData(dds)` to map to point
#'   shape in the PCA plot (e.g. a secondary grouping factor). Default `NULL`.
#' @param ntop Number of most-variable genes to use for PCA. Default `500`.
#' @param n_scree_pcs Number of leading PCs to show in the scree plot.
#'   Default `10` (capped to the number of PCs available).
#' @param n_loading_pcs Number of leading PCs to show top genes for in the
#'   loadings plot/CSV. Default `6` (capped to the number of PCs available).
#' @param n_top_genes Number of top **positive**-loading and top
#'   **negative**-loading genes to show per PC in the loadings plot/CSV
#'   (i.e. `n_top_genes` up + `n_top_genes` down = `2 * n_top_genes` genes
#'   per PC, ranked separately by sign). Default `10`.
#' @param save_dir Directory to write outputs to. Default `getwd()`.
#' @param save_plot Logical. If `TRUE`, saves plots as PDF. Default `TRUE`.
#' @return A list with elements `dds` (unmodified, returned for chaining),
#'   `pca` (data frame of PC scores + metadata), `pca_plot` (square PCA
#'   scatterplot on the corrected matrix if `batch` was corrected, else the
#'   plain PCA), `pca_plot_uncorrected` (the pre-correction PCA plot, or
#'   `NULL` if no batch correction was applied), `pca_plot_interactive` /
#'   `pca_plot_uncorrected_interactive` (`plotly` versions of the above,
#'   hoverable by sample name), `scree_plot` / `scree_plot_uncorrected`,
#'   `loadings` / `loadings_uncorrected` (data frames, genes x PCs),
#'   `top_genes` / `top_genes_uncorrected` (data frames of the top
#'   `n_top_genes` genes per PC), `top_loadings_plot` /
#'   `top_loadings_plot_uncorrected`, `similarity_plot` (square heatmap
#'   `grob`), `similarity_matrix` (data frame), and `percent_var`. The
#'   `*_uncorrected` elements are `NULL` unless batch correction was applied.
#' @examples
#' \dontrun{
#' out <- run_distance_pca(dds, var = "condition", batch = "batch")
#' }
#' @export
run_distance_pca <- function(
    dds,
    var,
    batch = NULL,
    batch_method = "limma",
    transform = "vst",
    pals = NULL,
    shape_var = NULL,
    ntop = 500,
    n_scree_pcs = 10,
    n_loading_pcs = 6,
    n_top_genes = 10,
    save_dir = getwd(),
    save_plot = TRUE) {

    dds <- validate_dds(dds)
    dds <- validate_var(dds, var)
    if (!is.null(batch)) dds <- validate_batch(dds, batch)
    batch_method <- validate_batch_method(batch_method)
    if (!is.null(pals)) pals <- validate_pals(dds, var, pals)
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)
    if (!transform %in% c("vst", "rlog")) {
        stop("`transform` must be \"vst\" or \"rlog\".", call. = FALSE)
    }

    vsd <- if (transform == "vst") {
        safe_vst(dds, blind = TRUE)
    } else {
        DESeq2::rlog(dds, blind = TRUE)
    }
    mat <- SummarizedExperiment::assay(vsd)

    batch_corrected <- !is.null(batch) && batch_method != "none"
    mat_primary <- mat
    if (batch_corrected) {
        mat_primary <- apply_batch_correction(mat, SummarizedExperiment::colData(dds), batch, batch_method)
    }

    coldata_df <- as.data.frame(SummarizedExperiment::colData(dds))
    coldata_df$sample <- rownames(coldata_df)

    # --- PCA ---------------------------------------------------------------
    pca_result <- compute_pca(mat_primary, ntop = ntop)
    pca_df <- pca_scores_df(pca_result, coldata_df)
    percent_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

    pca_plot <- plot_pca(pca_df, var = var, shape_var = shape_var, pals = pals, percent_var = percent_var)

    pca_plot_uncorrected <- NULL
    scree_plot_uncorrected <- NULL
    loadings_df_uncorrected <- NULL
    top_genes_df_uncorrected <- NULL
    top_loadings_plot_uncorrected <- NULL
    if (batch_corrected) {
        pca_result_raw <- compute_pca(mat, ntop = ntop)
        pca_df_raw <- pca_scores_df(pca_result_raw, coldata_df)
        percent_var_raw <- round(100 * pca_result_raw$sdev^2 / sum(pca_result_raw$sdev^2), 1)
        pca_plot_uncorrected <- plot_pca(pca_df_raw, var = var, shape_var = shape_var, pals = pals,
                                          percent_var = percent_var_raw) +
            ggplot2::ggtitle("Before batch correction")
        pca_plot <- pca_plot + ggplot2::ggtitle("After batch correction")

        scree_plot_uncorrected <- plot_pca_scree(percent_var_raw, n_pcs = n_scree_pcs) +
            ggplot2::ggtitle("Before batch correction")
        loadings_df_uncorrected <- compute_pca_loadings(pca_result_raw)
        top_genes_df_uncorrected <- top_genes_per_pc(loadings_df_uncorrected,
                                                       n_pcs = n_loading_pcs, n_top = n_top_genes)
        top_loadings_plot_uncorrected <- plot_pca_top_loadings(top_genes_df_uncorrected) +
            ggplot2::ggtitle("Before batch correction")
    }

    # Interactive (hover-for-sample-name) versions of the PCA scatterplots.
    # ~half of the original 550 px default, then +25% for readability
    # (550 * 0.5 * 1.25 ≈ 344).
    pca_plot_interactive <- make_interactive(pca_plot, base_size = 344)
    pca_plot_uncorrected_interactive <- if (!is.null(pca_plot_uncorrected)) {
        make_interactive(pca_plot_uncorrected, base_size = 344)
    } else {
        NULL
    }

    scree_plot <- plot_pca_scree(percent_var, n_pcs = n_scree_pcs)
    if (batch_corrected) scree_plot <- scree_plot + ggplot2::ggtitle("After batch correction")

    loadings_df <- compute_pca_loadings(pca_result)
    top_genes_df <- top_genes_per_pc(loadings_df, n_pcs = n_loading_pcs, n_top = n_top_genes)
    top_loadings_plot <- plot_pca_top_loadings(top_genes_df)
    if (batch_corrected) top_loadings_plot <- top_loadings_plot + ggplot2::ggtitle("After batch correction")

    # --- Sample similarity (Spearman correlation) ---------------------------
    similarity_plot <- plot_sample_similarity(mat_primary, dds, var = var)
    cor_mat <- attr(similarity_plot, "cor_matrix")
    cor_df <- as.data.frame(cor_mat)
    cor_df$sample <- rownames(cor_df)
    cor_df <- cor_df[, c("sample", setdiff(colnames(cor_df), "sample"))]

    if (save_plot) {
        save_ggplot(pca_plot, file.path(save_dir, "pca_plot.pdf"), width = 7.5, height = 7.5)
        if (!is.null(pca_plot_uncorrected)) {
            save_ggplot(pca_plot_uncorrected, file.path(save_dir, "pca_plot_uncorrected.pdf"),
                        width = 7.5, height = 7.5)
        }
        save_ggplot(scree_plot, file.path(save_dir, "pca_scree.pdf"), width = 6, height = 4.5)
        save_ggplot(top_loadings_plot, file.path(save_dir, "pca_top_loadings.pdf"),
                    width = 12, height = 6)
        if (!is.null(scree_plot_uncorrected)) {
            save_ggplot(scree_plot_uncorrected, file.path(save_dir, "pca_scree_uncorrected.pdf"),
                        width = 6, height = 4.5)
        }
        if (!is.null(top_loadings_plot_uncorrected)) {
            save_ggplot(top_loadings_plot_uncorrected,
                        file.path(save_dir, "pca_top_loadings_uncorrected.pdf"), width = 12, height = 6)
        }
        save_pheatmap(similarity_plot, file.path(save_dir, "sample_similarity_heatmap.pdf"))
    }
    write_output(pca_df, file.path(save_dir, "pca_coordinates.csv"), format = "csv")
    write_output(loadings_df, file.path(save_dir, "pca_loadings.csv"), format = "csv")
    write_output(top_genes_df, file.path(save_dir, "pca_top_genes_per_pc.csv"), format = "csv")
    if (!is.null(loadings_df_uncorrected)) {
        write_output(loadings_df_uncorrected, file.path(save_dir, "pca_loadings_uncorrected.csv"),
                     format = "csv")
        write_output(top_genes_df_uncorrected, file.path(save_dir, "pca_top_genes_per_pc_uncorrected.csv"),
                     format = "csv")
    }
    write_output(cor_df, file.path(save_dir, "sample_similarity_matrix.csv"), format = "csv")

    list(
        dds = dds,
        pca = pca_df,
        pca_plot = pca_plot,
        pca_plot_uncorrected = pca_plot_uncorrected,
        pca_plot_interactive = pca_plot_interactive,
        pca_plot_uncorrected_interactive = pca_plot_uncorrected_interactive,
        scree_plot = scree_plot,
        scree_plot_uncorrected = scree_plot_uncorrected,
        loadings = loadings_df,
        loadings_uncorrected = loadings_df_uncorrected,
        top_genes = top_genes_df,
        top_genes_uncorrected = top_genes_df_uncorrected,
        top_loadings_plot = top_loadings_plot,
        top_loadings_plot_uncorrected = top_loadings_plot_uncorrected,
        similarity_plot = similarity_plot,
        similarity_matrix = cor_df,
        percent_var = percent_var
    )
}

#' Build the PC-scores data frame (PC scores merged with sample metadata)
#' @param pca_result Result of `compute_pca()` (a `stats::prcomp` object).
#' @param coldata_df Data frame of sample metadata with a `sample` column.
#' @return A data frame with one row per sample: PC score columns plus
#'   the merged metadata columns.
#' @noRd
pca_scores_df <- function(pca_result, coldata_df) {
    pca_df <- as.data.frame(pca_result$x)
    pca_df$sample <- rownames(pca_df)
    merge(pca_df, coldata_df, by = "sample")
}

#' Apply batch correction to a VST/rlog expression matrix for visualization
#' @param mat Numeric matrix, genes x samples.
#' @param coldata `DataFrame`/data.frame of sample metadata.
#' @param batch Batch column name.
#' @param method `"limma"` or `"combat"`.
#' @return The batch-corrected matrix.
#' @noRd
apply_batch_correction <- function(mat, coldata, batch, method) {
    batch_factor <- factor(as.data.frame(coldata)[[batch]])
    if (method == "limma") {
        limma::removeBatchEffect(mat, batch = batch_factor)
    } else if (method == "combat") {
        sva::ComBat(dat = mat, batch = batch_factor, par.prior = TRUE, mean.only = FALSE)
    } else {
        mat
    }
}

#' Compute PCA on the top-variance genes of an expression matrix
#' @param mat Numeric matrix, genes x samples.
#' @param ntop Number of most-variable genes to use.
#' @return The result of `stats::prcomp()` on the transposed, top-variance
#'   subset of `mat`.
#' @noRd
compute_pca <- function(mat, ntop = 500) {
    row_vars <- matrixStats_rowVars(mat)
    ntop <- min(ntop, length(row_vars))
    select <- order(row_vars, decreasing = TRUE)[seq_len(ntop)]
    stats::prcomp(t(mat[select, , drop = FALSE]), center = TRUE, scale. = FALSE)
}

#' Row-wise variance without requiring the matrixStats package
#' @param mat Numeric matrix.
#' @return Numeric vector of per-row variances.
#' @noRd
matrixStats_rowVars <- function(mat) {
    apply(mat, 1, stats::var)
}

#' Plot a PCA scatterplot coloured by group (and optionally shaped by a
#' second variable)
#'
#' Rendered as a square panel (`aspect.ratio = 1`) regardless of the data
#' range on each axis, so PC1/PC2 plots are visually comparable across runs.
#' Axis expansion is tightened to 3% of the data range on each side (vs.
#' ggplot2's 5% default) so points sit closer to the panel edges.
#'
#' @param pca_df Data frame with columns `PC1`, `PC2`, `var`, and optionally
#'   `shape_var`.
#' @param var Grouping column name (colour).
#' @param shape_var Optional second grouping column name (shape).
#' @param pals Optional colour palette vector.
#' @param percent_var Numeric vector of % variance explained per PC.
#' @return A `ggplot` object.
#' @export
plot_pca <- function(pca_df, var, shape_var = NULL, pals = NULL, percent_var = NULL) {
    xlab_txt <- if (!is.null(percent_var)) paste0("PC1: ", percent_var[1], "% variance") else "PC1"
    ylab_txt <- if (!is.null(percent_var)) paste0("PC2: ", percent_var[2], "% variance") else "PC2"

    mapping <- if (!is.null(shape_var)) {
        ggplot2::aes(
            x = .data$PC1, y = .data$PC2, colour = .data[[var]], shape = .data[[shape_var]],
            text = paste0("Sample: ", .data$sample, "<br>", var, ": ", .data[[var]])
        )
    } else {
        ggplot2::aes(
            x = .data$PC1, y = .data$PC2, colour = .data[[var]],
            text = paste0("Sample: ", .data$sample, "<br>", var, ": ", .data[[var]])
        )
    }

    p <- ggplot2::ggplot(pca_df, mapping) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
        ggplot2::geom_point(size = 3, alpha = 0.85) +
        # 10% headroom on both axes so points near the extremes aren't
        # crowded against the panel border.
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
        unclipped_coord() +
        ggplot2::theme_bw(base_size = 12) +
        square_panel_theme() +
        ggplot2::xlab(xlab_txt) +
        ggplot2::ylab(ylab_txt)

    if (!is.null(pals)) {
        p <- p + ggplot2::scale_colour_manual(values = pals)
    }
    p
}

#' Estimate the on-screen pixel width a `ggplotly()` legend will occupy
#'
#' `plotly::plotly_build()`'s reported `layout$margin` (used by
#' [.square_plotly_canvas()] to size axis-label whitespace) never
#' includes the legend -- the legend is *not* part of the margin box at
#' all; plotly lays it out as a separate flex column beside the plotting
#' area, and that layout is resolved client-side, in the browser's own
#' JavaScript rendering of Plotly.js, which this R session cannot invoke.
#' Any square-canvas solver that only reads `layout$margin` is therefore
#' blind to the legend entirely, and a plot with a colour/shape legend
#' silently renders with a non-square *interior* data panel once the
#' legend's real width is carved out of the canvas by the browser -- the
#' panel narrows by exactly the legend's width, which nothing server-side
#' ever measured or reserved for.
#'
#' This estimates that width without a browser: legend entries are static
#' text (trace name / discrete level) at a known font size, so their
#' rendered width can be approximated with the same font-metrics engine
#' `graphics::strwidth()` uses for any other text sizing in R -- accurate
#' enough to reserve *at least* the true width (padded generously below)
#' without needing to actually lay the legend out.
#'
#' @param built The result of `plotly::plotly_build(plotly::ggplotly(...))`.
#' @return Numeric, estimated legend width in px (`0` if the plot has no
#'   legend).
#' @noRd
.estimate_legend_width_px <- function(built) {
    layout <- built$x$layout
    legend <- layout$legend
    if (is.null(legend) || isFALSE(layout$showlegend)) {
        return(0)
    }
    title_text <- legend$title$text
    if (is.null(title_text)) title_text <- ""
    title_text <- gsub("<br\\s*/?>", " ", title_text)

    font_size <- legend$font$size
    if (is.null(font_size) || !is.numeric(font_size)) font_size <- 12
    title_font_size <- legend$title$font$size
    if (is.null(title_font_size) || !is.numeric(title_font_size)) title_font_size <- font_size * 1.1

    labels <- vapply(built$x$data, function(tr) {
        sl <- tr$showlegend
        if (is.null(sl) || isTRUE(sl)) {
            nm <- tr$name
            if (is.null(nm)) "" else nm
        } else {
            ""
        }
    }, character(1))
    labels <- labels[nzchar(labels)]
    if (!length(labels) && !nzchar(title_text)) {
        return(0)
    }

    px_per_in <- 96
    swatch_and_padding_px <- 45 # colour/shape swatch + inter-element gaps
    outer_margin_px <- 10 # gap between panel and legend, and legend right edge

    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    label_w_in <- if (length(labels)) {
        max(graphics::strwidth(labels, units = "inches", cex = font_size / 12))
    } else {
        0
    }
    title_w_in <- if (nzchar(title_text)) {
        graphics::strwidth(title_text, units = "inches", cex = title_font_size / 12)
    } else {
        0
    }
    text_w_px <- max(label_w_in, title_w_in) * px_per_in
    # Round up generously (ceiling + a flat safety pad) -- reserving
    # slightly *too much* legend width still yields an exactly square
    # data panel (just a little extra whitespace next to the legend);
    # reserving too little is the failure mode this exists to avoid.
    ceiling(text_w_px + swatch_and_padding_px + outer_margin_px)
}

#' Solve for a `ggplotly()` canvas size and margins such that the *data
#' panel* is an exact square and the legend sits outside it on the right
#'
#' `ggplot2::theme(aspect.ratio = 1)` -- set by [plot_pca()], [plot_ma()],
#' and [plot_volcano()] -- is silently dropped by [plotly::ggplotly()]
#' (`plotly` warns "Aspect ratios aren't yet implemented..." and does not
#' carry the setting into the built widget). `scaleanchor`/`scaleratio`
#' looks like the natural substitute but is not: it locks the y-axis to
#' the x-axis at a 1:1 *data-unit* ratio (like `ggplot2::coord_equal()`),
#' which only yields a square panel by coincidence when the x- and y-data
#' ranges happen to match -- for the MA/volcano/PCA plots here they
#' generally don't, so `scaleanchor` alone does not reproduce the square
#' panel the static ggplot theme requests.
#'
#' A first, legend-*less* pass renders `gg_plot` at an arbitrary square
#' canvas to read back `plotly_build()`'s resolved axis margins (tick
#' labels, axis titles) -- these are what `layout$margin` actually
#' measures (confirmed identical whether or not the plot has a legend --
#' plotly's margin solver does not budget for the legend at all). A
#' second, legend-*bearing* pass at the same canvas feeds
#' [.estimate_legend_width_px()].
#'
#' The square target is the margin-subtracted *paper* area (axis-label
#' whitespace excluded; legend excluded):
#' `interior = height - top_margin - bottom_margin`. The right margin is
#' then expanded by the estimated legend width so the legend can sit at
#' `legend.x > 1` in paper coordinates -- i.e. in the right margin,
#' clearly outside the axes -- while the axes themselves keep
#' `domain = c(0, 1)` and therefore fill an `interior x interior` square.
#' (An earlier approach that shrunk `xaxis.domain` and pinned
#' `legend.x = domain[2]` looked correct in `plotly_build()`'s layout
#' object but still rendered the legend over the data panel in the
#' browser; putting the legend in the margin avoids that.)
#'
#' @param gg_plot A `ggplot` object.
#' @param tooltip Passed to `plotly::ggplotly(tooltip = )`.
#' @param base_size Canvas height (px); also the starting canvas width
#'   for the margin/legend discovery passes. Default `550`.
#' @return A list: `width`, `height` (px, canvas size -- fixed, not
#'   responsive), `margin` (named list `l`/`r`/`t`/`b` to apply via
#'   `plotly::layout()`), and `legend_x` (paper-coordinate left edge of
#'   the legend, `> 1` so it sits in the right margin; `NULL` when there
#'   is no legend).
#' @noRd
.square_plotly_canvas <- function(gg_plot, tooltip = "text", base_size = 550) {
    pass_nolegend <- plotly::ggplotly(
        gg_plot + ggplot2::theme(legend.position = "none"),
        tooltip = tooltip, width = base_size, height = base_size
    )
    m <- plotly::plotly_build(pass_nolegend)$x$layout$margin
    l <- if (is.null(m$l)) 0 else m$l
    r_axis <- if (is.null(m$r)) 0 else m$r
    t <- if (is.null(m$t)) 0 else m$t
    b <- if (is.null(m$b)) 0 else m$b

    pass_full <- plotly::ggplotly(gg_plot, tooltip = tooltip, width = base_size, height = base_size)
    built_full <- plotly::plotly_build(pass_full)
    legend_w <- .estimate_legend_width_px(built_full)

    # Height-driven square target for the paper area (axes only -- axis
    # title/tick margins already stripped via t/b; legend lives in the
    # expanded right margin below, not inside this square).
    interior <- base_size - t - b
    # Keep plotly's small right-axis pad and add the legend strip beyond it.
    r_total <- r_axis + legend_w
    total_width <- l + interior + r_total

    legend_x <- if (legend_w > 0) 1.02 else NULL

    list(
        width = total_width,
        height = base_size,
        margin = list(l = l, r = r_total, t = t, b = b),
        legend_x = legend_x
    )
}

#' Convert a `ggplot` object into an interactive `plotly` widget
#'
#' Thin wrapper around [plotly::ggplotly()] that shows only the mapped
#' `text` aesthetic in the hover tooltip (rather than plotly's default of
#' every mapped aesthetic, which tends to be noisy/redundant for plots that
#' already carry a purpose-built tooltip string), and drops the modebar
#' logo/some default buttons for a cleaner embed in the HTML report.
#'
#' When `square = TRUE` (the default), the *data panel* -- the axes'
#' plotting area, excluding axis-label whitespace AND the legend -- is
#' squared: an explicit canvas width/height and right margin are solved
#' (see [.square_plotly_canvas()]) so the margin-subtracted paper area is
#' an exact square and the legend (if any) is pinned just outside it at
#' `legend.x > 1` (in the right margin), rather than left to
#' `ggplotly()`'s own (legend-width-oblivious) automatic layout -- see
#' `.square_plotly_canvas()`'s documentation for why plotly's
#' `aspect.ratio`/`scaleanchor` cannot be used directly, and why the
#' legend must live in the margin rather than inside a shrunk axis
#' domain. The returned widget always carries an explicit fixed pixel
#' `width`/`height` (both from the solved canvas size and via
#' `config(responsive = FALSE)`) -- it does not resize when the browser
#' window is resized.
#'
#' When `print()`-ing the returned widget from inside a knitr
#' `results='asis'` loop (as the bundled report template's per-comparison
#' loop does), use [emit_widget()] rather than `print()`/`tagList()`
#' directly -- see that function's documentation for why.
#'
#' @param gg_plot A `ggplot` object whose mapping includes a `text`
#'   aesthetic (e.g. as added by [plot_pca()], [plot_ma()], [plot_volcano()]).
#' @param tooltip Passed to `plotly::ggplotly(tooltip = )`. Default
#'   `"text"` (show only the custom tooltip string).
#' @param square If `TRUE` (default), constrain the rendered widget's data
#'   panel (excluding the legend) to a square interior plotting area (see
#'   Details).
#' @param base_size Canvas height (px) for the square-solving pass; also
#'   the final canvas size when `square = FALSE`. Default `550`.
#' @return A `plotly` htmlwidget object, always with a fixed (non-responsive)
#'   pixel size.
#' @export
make_interactive <- function(gg_plot, tooltip = "text", square = TRUE, base_size = 550) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Package \"plotly\" is required for interactive plots. Install it with ",
             "`install.packages(\"plotly\")`.", call. = FALSE)
    }
    if (isTRUE(square)) {
        dims <- .square_plotly_canvas(gg_plot, tooltip = tooltip, base_size = base_size)
        p <- plotly::ggplotly(gg_plot, tooltip = tooltip, width = dims$width, height = dims$height)
        layout_args <- list(
            p = p,
            margin = dims$margin,
            autosize = FALSE,
            # Paper area is already the square target; keep axes filling it
            # and stop automargin from re-expanding into the legend strip.
            xaxis = list(domain = c(0, 1), automargin = FALSE),
            yaxis = list(domain = c(0, 1), automargin = FALSE)
        )
        if (!is.null(dims$legend_x)) {
            layout_args$legend <- list(
                x = dims$legend_x,
                xanchor = "left",
                xref = "paper",
                y = 0.5,
                yanchor = "middle",
                yref = "paper"
            )
        }
        p <- do.call(plotly::layout, layout_args)
    } else {
        p <- plotly::ggplotly(gg_plot, tooltip = tooltip, width = base_size, height = base_size)
    }
    # ggplotly drops ggplot2::annotate(); re-apply via plotly::layout so
    # they survive plotly_build() (stored in layoutAttrs until build).
    anns <- attr(gg_plot, "plotly_annotations")
    if (!is.null(anns) && length(anns) > 0) {
        p <- plotly::layout(p, annotations = anns)
    }
    plotly::config(p, displaylogo = FALSE, responsive = FALSE)
}

#' Emit an htmlwidget from inside a knitr `results='asis'` chunk/loop
#'
#' `print()`-ing an htmlwidget (a `plotly` or `DT::datatable` object, even
#' wrapped in `htmltools::tagList()`) from inside a `results='asis'` chunk
#' that runs a loop does not reliably register that widget's JavaScript
#' library (Plotly.js, DataTables.js, ...) as a knitr dependency -- the
#' widget's container `<div>` and JSON data get written into the output
#' HTML, but the library that actually turns that div into a live chart or
#' interactive table never gets embedded, so every affected widget renders
#' as a blank/dead placeholder in the browser. This affects both `plotly`
#' and `DT::datatable()` widgets, and is specific to the print-inside-a-loop
#' pattern (a single widget in its own non-looped, non-`asis` chunk embeds
#' its dependencies correctly via the normal knitr/rmarkdown print method).
#'
#' The fix is to bypass the automatic print-time dependency registration
#' entirely: extract the widget's HTML dependencies explicitly via
#' [htmltools::findDependencies()], register them directly with
#' [knitr::knit_meta_add()], and then emit only the tag HTML (as a plain
#' character string via `cat()`) into the `asis` output stream. Use this
#' in place of `print(htmltools::tagList(widget))` anywhere an htmlwidget
#' (from [make_interactive()], [make_deg_datatable()], or
#' [make_gsea_datatable()]) is emitted inside a `results='asis'` loop, as
#' the bundled report template does for the per-comparison MA/volcano
#' plots and DEG/GSEA tables.
#'
#' The tag HTML is rendered with `indent = FALSE` (via
#' [htmltools::doRenderTags()]) rather than `as.character()`'s default
#' pretty-printing, which indents each nested tag by two spaces per
#' nesting level. Pandoc's Markdown parser treats any line indented by
#' four or more spaces as a fenced/indented code block, so content nested
#' three or more `<div>`s deep (e.g. [wrap_with_corner_button()]'s
#' centering wrapper around its corner-button holder around a plot/table)
#' would otherwise render as a literal, HTML-escaped code block instead of
#' live markup -- exactly the failure mode this avoids.
#'
#' @param widget An htmlwidget object (e.g. a `plotly` object or a
#'   `DT::datatable`).
#' @return Invisibly, `NULL`. Called for its side effect of `cat()`-ing the
#'   widget's tag HTML and registering its dependencies with
#'   [knitr::knit_meta_add()].
#' @export
emit_widget <- function(widget) {
    if (!requireNamespace("htmltools", quietly = TRUE)) {
        stop("Package \"htmltools\" is required to emit interactive widgets.", call. = FALSE)
    }
    tags <- htmltools::as.tags(widget)
    knitr::knit_meta_add(htmltools::findDependencies(tags))
    cat(as.character(htmltools::doRenderTags(htmltools::tagList(tags), indent = FALSE)))
    invisible(NULL)
}

#' Compute a fixed pheatmap cell size (in points) so an n x n data matrix
#' renders with aspect ratio 1
#'
#' `pheatmap()`'s default auto-sizing derives the width and height of the
#' *data matrix* panel independently from the surrounding row/column label
#' text, dendrograms, and legend -- so a square overall device canvas does
#' not guarantee a square matrix panel whenever row and column label text
#' differs in extent (e.g. long sample names). Passing an explicit, equal
#' `cellwidth`/`cellheight` to `pheatmap()` instead fixes the matrix
#' panel's absolute size to `n * cell_size` on both axes, independent of
#' label/dendrogram overhead, which is exactly what a symmetric
#' sample-by-sample matrix (distance/correlation) needs for aspect ratio 1.
#'
#' @param n Number of rows/columns (samples) in the (square) matrix.
#' @param min_size,max_size Minimum/maximum per-cell size, in points.
#'   Defaults `12`/`40`.
#' @param target_total Target total matrix side length, in points, before
#'   clamping to `min_size`/`max_size` per cell. Default `320`.
#' @return A single numeric cell size, in points, usable as both
#'   `cellwidth` and `cellheight`.
#' @noRd
pheatmap_cell_size <- function(n, min_size = 12, max_size = 40, target_total = 320) {
    if (n <= 0) return(max_size)
    max(min_size, min(max_size, target_total / n))
}

#' Plot a sample-to-sample Spearman correlation (similarity) heatmap
#'
#' A heatmap whose **data-matrix panel** (excluding row/column name labels
#' and dendrograms) has aspect ratio 1, showing pairwise Spearman
#' correlation between samples in VST/rlog expression space. The colour
#' gradient is scaled to the *actual* min/max of the correlation matrix
#' (rather than a fixed range), so the visually-relevant contrast in
#' same-experiment correlations -- which are typically all high and
#' clustered in a narrow range -- is not compressed by unused scale.
#'
#' @param mat Numeric matrix, genes x samples (VST/rlog, possibly batch-corrected).
#' @param dds A `DESeqDataSet`, used to annotate rows/columns by `var`.
#' @param var Grouping column name for row/column annotation.
#' @return A `pheatmap` grob (invisibly printable/savable object), with the
#'   underlying correlation matrix attached as `attr(., "cor_matrix")`.
#' @export
plot_sample_similarity <- function(mat, dds, var) {
    cor_mat <- stats::cor(mat, method = "spearman")
    cor_dist <- stats::as.dist(1 - cor_mat)
    can_cluster <- !anyNA(cor_dist)

    annotation <- data.frame(group = SummarizedExperiment::colData(dds)[[var]])
    rownames(annotation) <- colnames(dds)
    colnames(cor_mat) <- colnames(dds)
    rownames(cor_mat) <- colnames(dds)

    # Self-correlation (diagonal == 1) carries no information and would
    # otherwise dominate the colour scale; blank it out to NA for display
    # only. The full matrix (diagonal included) is still attached below via
    # `attr(., "cor_matrix")` and is what gets written to
    # sample_similarity_matrix.csv.
    display_mat <- cor_mat
    diag(display_mat) <- NA

    # Similarity: higher correlation (more similar) = red, lower = blue.
    colours <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
    # Scale breaks to the actual off-diagonal data range rather than a fixed
    # [-1, 1] or [0, 1], so the colour scale uses its full range across the
    # values actually present in this dataset's correlation matrix.
    lo <- min(display_mat, na.rm = TRUE)
    hi <- max(display_mat, na.rm = TRUE)
    if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
        # Degenerate (e.g. single sample, or all-identical correlations):
        # fall back to a small padded range so pheatmap's cut() doesn't error.
        lo <- (if (is.finite(lo)) lo else 0) - 0.01
        hi <- (if (is.finite(hi)) hi else 1) + 0.01
    }
    breaks <- seq(lo, hi, length.out = length(colours) + 1)
    cell_size <- pheatmap_cell_size(ncol(cor_mat))

    number_labels <- matrix(
        sprintf("%.2f", display_mat),
        nrow = nrow(display_mat), dimnames = dimnames(display_mat)
    )
    number_labels[is.na(display_mat)] <- ""

    p <- pheatmap::pheatmap(
        display_mat,
        cluster_rows = can_cluster,
        cluster_cols = can_cluster,
        clustering_distance_rows = if (can_cluster) cor_dist else NA,
        clustering_distance_cols = if (can_cluster) cor_dist else NA,
        annotation_row = annotation,
        color = colours,
        breaks = breaks,
        display_numbers = number_labels,
        na_col = "grey20",
        cellwidth = cell_size,
        cellheight = cell_size,
        silent = TRUE
    )
    attr(p, "cor_matrix") <- cor_mat
    p
}

#' Plot a PCA scree plot (percent variance explained per PC)
#' @param percent_var Numeric vector of % variance explained per PC (as
#'   returned in `run_distance_pca()`'s output).
#' @param n_pcs Number of leading PCs to display. Default `10`.
#' @return A `ggplot` object: bar plot of % variance per PC with a
#'   cumulative-variance line overlaid.
#' @export
plot_pca_scree <- function(percent_var, n_pcs = 10) {
    n_pcs <- min(n_pcs, length(percent_var))
    df <- data.frame(
        pc = factor(paste0("PC", seq_len(n_pcs)), levels = paste0("PC", seq_len(n_pcs))),
        percent_var = percent_var[seq_len(n_pcs)]
    )
    df$cumulative <- cumsum(df$percent_var)
    # scale cumulative line onto the same axis as the bars
    scale_factor <- max(df$percent_var) / max(df$cumulative)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$pc)) +
        ggplot2::geom_col(ggplot2::aes(y = .data$percent_var), fill = "#3182bd") +
        ggplot2::geom_line(ggplot2::aes(y = .data$cumulative * scale_factor, group = 1),
                            colour = "#e6550d", linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes(y = .data$cumulative * scale_factor),
                             colour = "#e6550d", size = 2) +
        ggplot2::scale_y_continuous(
            name = "% variance explained",
            sec.axis = ggplot2::sec_axis(~ . / scale_factor, name = "Cumulative % variance"),
            expand = zero_expansion()
        ) +
        # Zero data-padding (above) can put the cumulative-variance line's
        # point markers exactly on the panel border; without disabling
        # clipping, their drawn radius would be sliced off there.
        unclipped_coord() +
        ggplot2::theme_bw(base_size = 12) +
        square_panel_theme() +
        ggplot2::xlab(NULL)
}

#' Extract per-gene PCA loadings (rotation matrix) as a tidy data frame
#' @param pca_result Result of `compute_pca()` (a `stats::prcomp` object).
#' @return A data frame with column `gene` plus one column per PC
#'   (`PC1`, `PC2`, ...) giving that gene's loading on each PC.
#' @export
compute_pca_loadings <- function(pca_result) {
    loadings <- as.data.frame(pca_result$rotation)
    loadings$gene <- rownames(loadings)
    loadings <- loadings[, c("gene", setdiff(colnames(loadings), "gene")), drop = FALSE]
    rownames(loadings) <- NULL
    loadings
}

#' Extract the top up- and down-loading genes for each of the leading PCs
#'
#' For each PC, separately ranks genes by their **positive** loadings and by
#' their **negative** loadings and keeps the top `n_top` of each, so that
#' both directions of the PC axis are represented regardless of any
#' imbalance in how many genes load one way vs. the other.
#'
#' @param loadings_df Result of [compute_pca_loadings()].
#' @param n_pcs Number of leading PCs to extract top genes for. Default `6`.
#' @param n_top Number of top **positive**-loading genes and top
#'   **negative**-loading genes to keep per PC (`2 * n_top` genes per PC
#'   total). Default `10`.
#' @return A long data frame with columns `pc`, `gene`, `loading`,
#'   `direction` (`"up"`/`"down"`), and `rank` (rank within `direction`,
#'   `1` = largest magnitude).
#' @export
top_genes_per_pc <- function(loadings_df, n_pcs = 6, n_top = 10) {
    pc_cols <- setdiff(colnames(loadings_df), "gene")
    n_pcs <- min(n_pcs, length(pc_cols))
    pc_cols <- pc_cols[seq_len(n_pcs)]

    result <- lapply(pc_cols, function(pc) {
        loading <- loadings_df[[pc]]
        gene <- loadings_df$gene

        up_ord <- order(loading, decreasing = TRUE)
        up_ord <- up_ord[loading[up_ord] > 0]
        up_ord <- up_ord[seq_len(min(n_top, length(up_ord)))]

        down_ord <- order(loading, decreasing = FALSE)
        down_ord <- down_ord[loading[down_ord] < 0]
        down_ord <- down_ord[seq_len(min(n_top, length(down_ord)))]

        rbind(
            data.frame(pc = pc, gene = gene[up_ord], loading = loading[up_ord],
                       direction = "up", rank = seq_along(up_ord)),
            data.frame(pc = pc, gene = gene[down_ord], loading = loading[down_ord],
                       direction = "down", rank = seq_along(down_ord))
        )
    })
    do.call(rbind, result)
}

#' Plot top up/down-loading genes per PC as a faceted bar chart
#' @param top_genes_df Result of [top_genes_per_pc()].
#' @return A `ggplot` object: one facet per PC, genes on the x-axis
#'   (ordered by loading), loading value on the y-axis, coloured by
#'   direction (up/down).
#' @export
plot_pca_top_loadings <- function(top_genes_df) {
    top_genes_df$sign <- ifelse(top_genes_df$loading >= 0, "positive", "negative")
    # order genes within each facet by loading (ggplot2 facet-local reordering)
    top_genes_df$gene_order <- paste(top_genes_df$pc, top_genes_df$gene, sep = "___")
    top_genes_df$gene_order <- stats::reorder(top_genes_df$gene_order, top_genes_df$loading)

    ggplot2::ggplot(top_genes_df, ggplot2::aes(x = .data$gene_order, y = .data$loading, fill = .data$sign)) +
        ggplot2::geom_col() +
        ggplot2::scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
        ggplot2::scale_y_continuous(expand = zero_expansion()) +
        ggplot2::scale_fill_manual(values = c(positive = "#de2d26", negative = "#3182bd")) +
        ggplot2::facet_wrap(~pc, scales = "free") +
        ggplot2::theme_bw(base_size = 11) +
        square_panel_theme() +
        ggplot2::theme(legend.position = "none",
                        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab("Loading")
}

#' Save a pheatmap grob to PDF
#'
#' By default, sizes the PDF device to the gtable's own physical size (as
#' computed by `pheatmap`/`grid`), so panels built with a fixed
#' `cellwidth`/`cellheight` (e.g. via `pheatmap_cell_size()`) render at
#' their intended size without clipping or added whitespace. Pass explicit
#' `width`/`height` to override.
#'
#' @param heatmap_grob A `pheatmap` object (as returned by [plot_sample_similarity()]).
#' @param path Output file path.
#' @param width,height Plot dimensions in inches, or `NULL` (default) to
#'   auto-size from the gtable.
#' @return Invisibly, `path`.
#' @noRd
save_pheatmap <- function(heatmap_grob, path, width = NULL, height = NULL) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    gt <- heatmap_grob$gtable
    if (is.null(width) || is.null(height)) {
        grDevices::pdf(NULL)
        auto_w <- grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE)
        auto_h <- grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
        grDevices::dev.off()
        if (is.null(width)) width <- auto_w
        if (is.null(height)) height <- auto_h
    }
    grDevices::pdf(path, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(gt)
    grDevices::dev.off()
    invisible(path)
}
