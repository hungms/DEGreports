# Package-wide geometry conventions, checked directly against the shared
# helpers in R/plot_helpers.R (zero_expansion(), square_panel_theme()):
#
#   1. A *narrow, explicitly-scoped* set of ggplot figures render their
#      data panel as an exact square via `ggplot2::theme(aspect.ratio = 1)`
#      (`square_panel_theme()`): plot_pca, plot_pca_scree,
#      plot_pca_top_loadings, plot_ma, plot_volcano, plot_library_size,
#      plot_count_distribution, and the sample-similarity
#      pheatmaps (test-report-helpers.R). Squaring these -- PCA/QC/DE
#      overview plots users compare side-by-side, plus the two
#      symmetric sample x sample matrices -- makes cross-plot comparison
#      easier.
#   2. GSEA (plot_gsea_barplot) and every WGCNA ggplot function
#      (plot_soft_threshold, plot_module_sizes, plot_module_trait_heatmap,
#      plot_module_ora_dotplot) are explicitly EXEMPT from squaring: their
#      natural aspect ratio (one bar per pathway/module, a lopsided
#      module x trait grid, ...) is more informative than a forced
#      square, and forcing it produced distorted-looking panels for
#      typical gene-set/module counts. plot_module_eigengene_heatmap
#      (a pheatmap) is similarly exempt: its cellwidth/cellheight are
#      sized independently per axis (module count vs. sample count),
#      so the overall panel is rectangular whenever those counts differ.
#   3. Continuous axes on the squared plots carry zero data-padding
#      (`expand = zero_expansion()`, i.e. c(0, 0, 0, 0)) so plotted data
#      touches the panel edges -- EXCEPT both axes of plot_pca /
#      plot_volcano and the y-axis of plot_ma, which instead use a
#      10%-both-sides multiplicative expansion (`expansion(mult = c(0.1,
#      0.1))`) so points near the extremes aren't crowded against the
#      border (see test-distance-pca.R / test-diffexp-plots.R for those).
#      Discrete/categorical axes deliberately keep ggplot2's small built-in
#      discrete expansion, to avoid clipping the outermost bar/box/tile.
#
#   4. Zero data-padding (point 3) puts the outermost data value's
#      *coordinate* exactly on the panel border. Any geom with rendered
#      extent beyond that coordinate -- a point marker's radius, an
#      outlined bar's stroke, a boxplot's whisker/outlier -- would
#      otherwise have that extent hard-clipped by ggplot2's default
#      `coord_cartesian(clip = "on")`, visibly cutting off data at the
#      panel edges. Every such plot instead uses `unclipped_coord()` /
#      `unclipped_coord_flip()` (`clip = "off"`) so the marker/stroke
#      renders in full while the zero-padded axis range itself is
#      unaffected. Plots whose geoms have no rendered extent beyond their
#      data coordinate (a bare, unstroked `geom_col()`; a `colour = NA`
#      `geom_density()` fill) correctly keep the default `clip = "on"` --
#      checked below alongside the `clip = "off"` cases so a regression
#      in either direction is caught.
#
# Functions covered elsewhere (not repeated here): plot_pca, plot_ma,
# plot_volcano (test-distance-pca.R / test-diffexp-plots.R), and the
# square-count sample-similarity pheatmaps
# (test-report-helpers.R).

test_that("zero_expansion() returns exactly zero multiplicative and additive expansion", {
    ze <- zero_expansion()
    expect_equal(as.numeric(ze), c(0, 0, 0, 0))
})

test_that("unclipped_coord()/unclipped_coord_flip() disable clipping and pass extra args through", {
    cc <- unclipped_coord()
    expect_s3_class(cc, "CoordCartesian")
    expect_equal(cc$clip, "off")

    cc_ylim <- unclipped_coord(ylim = c(0, 1.05))
    expect_equal(cc_ylim$clip, "off")
    expect_equal(cc_ylim$limits$y, c(0, 1.05))

    cf <- unclipped_coord_flip()
    expect_s3_class(cf, "CoordFlip")
    expect_equal(cf$clip, "off")
})

test_that("boundary-touching point markers are not clipped under zero expansion + unclipped_coord() (regression test)", {
    # Reproduces the reported defect directly: with expand = zero_expansion(),
    # the outermost point's *coordinate* sits exactly on the panel border, so
    # under ggplot2's default coord_cartesian(clip = "on") its drawn radius
    # is sliced off there. unclipped_coord() (clip = "off") is the fix.
    df <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
    p_clipped <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
        ggplot2::geom_point(size = 8) +
        ggplot2::scale_x_continuous(expand = zero_expansion()) +
        ggplot2::scale_y_continuous(expand = zero_expansion())
    expect_equal(p_clipped$coordinates$clip, "on")

    p_fixed <- p_clipped + unclipped_coord()
    expect_equal(p_fixed$coordinates$clip, "off")

    # Confirm at the grob level: rendering the clip="off" version produces
    # a viewport/grob tree without the panel clip rectangle being applied
    # to the point layer (grid::grid.draw is not invoked here -- checking
    # the built ggplot's own coord object, which is what every renderer
    # (svglite/ragg/plotly) reads to decide whether to install a clip path,
    # is the stable, renderer-independent assertion).
    built <- ggplot2::ggplot_build(p_fixed)
    expect_equal(built$plot$coordinates$clip, "off")
})

test_that("square_panel_theme() sets aspect.ratio = 1 and nothing else", {
    th <- square_panel_theme()
    expect_s3_class(th, "theme")
    expect_equal(th$aspect.ratio, 1)
})

test_that("pheatmap_cell_size sizes a single dimension independently (no cross-axis squaring helper exists)", {
    # pheatmap_square_cell_size() (which solved cellwidth/cellheight so an
    # unequal-row/column matrix rendered as an overall square) was removed
    # together with the WGCNA/GSEA squaring convention; only the
    # single-dimension pheatmap_cell_size() remains, applied independently
    # per axis (see plot_module_eigengene_heatmap below).
    expect_false(exists("pheatmap_square_cell_size", where = asNamespace("deseq2reports"), inherits = FALSE))

    cell_w <- pheatmap_cell_size(30)
    cell_h <- pheatmap_cell_size(3)
    # applying the function to different n's independently need not (and,
    # for a lopsided n, does not) yield a square overall panel
    expect_false(isTRUE(all.equal(cell_w * 30, cell_h * 3)))
})

test_that("plot_pca_scree renders a square panel with zero y-axis expansion and clip = \"off\" (the cumulative-variance line's point markers sit on the zero-padded panel border)", {
    p <- plot_pca_scree(percent_var = c(40, 20, 15, 10, 8, 4, 2, 1), n_pcs = 5)
    expect_s3_class(p, "ggplot")
    expect_equal(p$theme$aspect.ratio, 1)
    y_scale <- p$scales$get_scales("y")
    expect_equal(y_scale$expand, c(0, 0, 0, 0))
    expect_equal(p$coordinates$clip, "off")
})

test_that("plot_pca_top_loadings renders square facet panels with zero y-axis expansion and default clip (bare geom_col() has no stroke beyond its data extent)", {
    top_genes_df <- data.frame(
        gene_order = factor(paste0("PC1___g", 1:6), levels = paste0("PC1___g", 1:6)),
        gene = paste0("g", 1:6),
        pc = "PC1",
        loading = c(0.5, 0.4, 0.3, -0.3, -0.4, -0.5),
        sign = rep(c("positive", "negative"), each = 3)
    )
    p <- plot_pca_top_loadings(top_genes_df)
    expect_s3_class(p, "ggplot")
    expect_equal(p$theme$aspect.ratio, 1)
    y_scale <- p$scales$get_scales("y")
    expect_equal(y_scale$expand, c(0, 0, 0, 0))
    expect_equal(p$coordinates$clip, "on")
})

test_that("plot_count_distribution renders a square panel with zero expansion on both continuous axes and default clip (colour = NA fill has no stroke beyond its data extent)", {
    gene_stat <- c(1, 2, 3, 4, 5, 6, 7, 8)
    keep <- c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
    p <- plot_count_distribution(gene_stat, cutoff_log2 = 2.5, keep = keep, xlab = "log2 mean CPM")
    expect_equal(p$theme$aspect.ratio, 1)
    expect_equal(p$scales$get_scales("x")$expand, c(0, 0, 0, 0))
    expect_equal(p$scales$get_scales("y")$expand, c(0, 0, 0, 0))
    expect_equal(p$coordinates$clip, "on")
})

test_that("plot_library_size renders a square panel: zero expansion on the continuous y-axis (discrete x-axis unaffected), clip = \"off\" so boxplot whiskers/outliers at the y-extremes aren't clipped", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    dds <- DESeq2::estimateSizeFactors(dds)
    p <- plot_library_size(dds)
    expect_s3_class(p, "ggplot")
    expect_equal(p$theme$aspect.ratio, 1)
    expect_equal(p$scales$get_scales("y")$expand, c(0, 0, 0, 0))
    expect_equal(p$coordinates$clip, "off")
})

test_that("plot_soft_threshold is exempt from square-panel theming (no aspect.ratio) but keeps zero expansion and clip = \"off\" on both sub-plots (point markers and power labels sit at the zero-padded axis extremes)", {
    fit_indices <- data.frame(
        Power = 1:10,
        SFT.R.sq = c(0.1, 0.3, 0.5, 0.7, 0.82, 0.88, 0.9, 0.91, 0.92, 0.93),
        mean.k. = c(500, 300, 180, 100, 60, 35, 20, 12, 7, 4)
    )
    combined <- plot_soft_threshold(fit_indices, power = 6, rsq_cutoff = 0.85)
    expect_s3_class(combined, "patchwork")

    # patchwork::wrap_plots(p_fit, p_conn) stores all but the last
    # sub-plot in `$patches$plots`; the combined object itself *is* the
    # last sub-plot (p_conn) with the patchwork layout attached.
    p_fit <- combined$patches$plots[[1]]
    p_conn <- combined
    expect_s3_class(p_fit, "ggplot")
    # WGCNA plots are explicitly exempt from square_panel_theme()
    expect_null(p_fit$theme$aspect.ratio)
    expect_null(p_conn$theme$aspect.ratio)
    expect_equal(p_fit$scales$get_scales("x")$expand, c(0, 0, 0, 0))
    expect_equal(p_conn$scales$get_scales("x")$expand, c(0, 0, 0, 0))
    expect_equal(p_fit$coordinates$clip, "off")
    expect_equal(p_conn$coordinates$clip, "off")
})

test_that("plot_module_sizes is exempt from square-panel theming (no aspect.ratio) but keeps zero expansion on the continuous x-axis and clip = \"off\" (outlined bars' stroke sits on the zero-padded panel border)", {
    module_df <- data.frame(
        gene = paste0("g", 1:12),
        module = rep(c("Module 1", "Module 2", "Unassigned"), c(6, 4, 2)),
        module_color = rep(c("turquoise", "blue", "grey"), c(6, 4, 2))
    )
    p <- plot_module_sizes(module_df)
    expect_null(p$theme$aspect.ratio)
    expect_equal(p$scales$get_scales("x")$expand, c(0, 0, 0, 0))
    expect_equal(p$coordinates$clip, "off")
})

test_that("plot_module_ora_dotplot is exempt from square-panel theming (no aspect.ratio) but keeps zero x-axis expansion and clip = \"off\" (point markers at the largest gene-set size sit on the zero-padded panel border)", {
    ora_res <- data.frame(
        pathway = paste0("PATHWAY_", 1:6),
        padj = c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04),
        overlap = 1:6,
        size = c(20, 40, 60, 80, 100, 120)
    )
    p <- plot_module_ora_dotplot(ora_res, top_n = 6, title = "Module 1")
    expect_s3_class(p, "ggplot")
    expect_null(p$theme$aspect.ratio)
    expect_equal(p$scales$get_scales("x")$expand, c(0, 0, 0, 0))
    expect_equal(p$coordinates$clip, "off")
})

test_that("plot_module_trait_heatmap is exempt from square-panel theming (no aspect.ratio), independent of module/trait-level count", {
    # deliberately unequal module count (3) vs. trait-level count (2)
    module_trait_cor <- data.frame(
        module = rep(c("Module 1", "Module 2", "Module 3"), each = 2),
        level = rep(c("Treated", "Control"), 3),
        cor = c(0.8, -0.8, 0.3, -0.3, 0.1, -0.1),
        pvalue = c(0.001, 0.001, 0.2, 0.2, 0.5, 0.5)
    )
    p <- plot_module_trait_heatmap(module_trait_cor)
    expect_null(p$theme$aspect.ratio)
})

test_that("plot_module_eigengene_heatmap sizes cellwidth/cellheight independently per axis, so a lopsided (non-square-count) panel is NOT forced square", {
    # 3 modules x 12 samples: a deliberately lopsided matrix.
    n_samples <- 12
    MEs <- data.frame(
        `Module 1` = stats::rnorm(n_samples),
        `Module 2` = stats::rnorm(n_samples),
        Unassigned = stats::rnorm(n_samples),
        check.names = FALSE,
        row.names = paste0("sample", seq_len(n_samples))
    )
    p <- plot_module_eigengene_heatmap(MEs)
    expect_s3_class(p, "pheatmap")

    gt <- p$gtable
    idx <- which(gt$layout$name == "matrix")
    w <- grid::convertWidth(gt$widths[gt$layout$l[idx]], "in", valueOnly = TRUE)
    h <- grid::convertHeight(gt$heights[gt$layout$t[idx]], "in", valueOnly = TRUE)
    # cellwidth/cellheight are solved independently off ncol/nrow
    # (pheatmap_cell_size() applied per-axis) -- for this deliberately
    # lopsided 3 x 12 matrix the resulting panel must NOT be square.
    expect_false(isTRUE(all.equal(w, h, tolerance = 1e-6)))
})

test_that("plot_gsea_barplot is exempt from square-panel theming (no aspect.ratio), uses 25% NES-axis expansion, and clip = \"off\" via unclipped_coord_flip()", {
    gsea_res <- data.frame(
        pathway = paste0("PATH", 1:6),
        NES = c(2.5, 2, 1.5, -1.5, -2, -2.5),
        padj = c(0.001, 0.9, 0.95, 0.95, 0.9, 0.001),
        pval = c(0.001, 0.9, 0.95, 0.95, 0.9, 0.001)
    )
    p <- plot_gsea_barplot(gsea_res, top_n = 3)
    expect_null(p$theme$aspect.ratio)
    # the NES axis is `y` pre-coord_flip() (visual x after flip)
    y_scale <- p$scales$get_scales("y")
    expect_equal(y_scale$expand, c(0.25, 0, 0.25, 0))
    expect_s3_class(p$coordinates, "CoordFlip")
    expect_equal(p$coordinates$clip, "off")
})

test_that("plot_module_ora_dotplot is exempt from square-panel theming (no aspect.ratio) but keeps zero expansion on the continuous x-axis", {
    ora_res <- data.frame(
        pathway = paste0("HALLMARK_SET_", 1:5),
        padj = c(0.001, 0.01, 0.02, 0.03, 0.5),
        size = c(20, 15, 10, 8, 2)
    )
    p <- plot_module_ora_dotplot(ora_res, top_n = 3, title = "Module: Module 1")
    expect_null(p$theme$aspect.ratio)
    expect_equal(p$scales$get_scales("x")$expand, c(0, 0, 0, 0))
})
