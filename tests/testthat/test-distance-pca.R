# Tests for run_distance_pca() (sample distance / PCA stage).

test_that("run_distance_pca returns the expected structure and writes outputs", {
    dds <- simulate_dds(n_genes = 100, n_per_group = 4)
    save_dir <- withr::local_tempdir()

    out <- run_distance_pca(dds, var = "condition", save_dir = save_dir)

    expect_s4_class(out$dds, "DESeqDataSet")
    expect_true(all(c("PC1", "PC2", "sample", "condition") %in% colnames(out$pca)))
    expect_s3_class(out$pca_plot, "ggplot")
    expect_null(out$pca_plot_uncorrected)
    expect_null(out$pca_plot_uncorrected_interactive)
    expect_null(out$scree_plot_uncorrected)
    expect_null(out$loadings_uncorrected)
    expect_null(out$top_genes_uncorrected)
    expect_null(out$top_loadings_plot_uncorrected)
    expect_null(out$pca_assoc_uncorrected)
    expect_null(out$pca_assoc_plot_uncorrected)
    expect_s3_class(out$pca_plot_interactive, "plotly")
    expect_s3_class(out$scree_plot, "ggplot")
    expect_true(all(c("gene", "PC1", "PC2") %in% colnames(out$loadings)))
    expect_true(all(c("pc", "gene", "loading", "direction", "rank") %in% colnames(out$top_genes)))
    expect_s3_class(out$top_loadings_plot, "ggplot")
    expect_true(is.data.frame(out$pca_assoc))
    expect_true(all(c("pc", "covariate", "p_value", "neglog10_p") %in% colnames(out$pca_assoc)))
    expect_false(is.null(out$pca_assoc_plot))
    expect_true(is.data.frame(out$similarity_matrix))
    expect_null(out$similarity_plot_corrected)
    expect_null(out$similarity_matrix_corrected)
    expect_equal(out$batch_correction, "none")
    expect_true(is.numeric(out$percent_var))

    expect_true(file.exists(file.path(save_dir, "pca_plot.pdf")))
    expect_true(file.exists(file.path(save_dir, "pca_scree.pdf")))
    expect_true(file.exists(file.path(save_dir, "pca_top_loadings.pdf")))
    expect_true(file.exists(file.path(save_dir, "pca_covariate_association.pdf")))
    expect_true(file.exists(file.path(save_dir, "pca_covariate_association.csv")))
    expect_true(file.exists(file.path(save_dir, "sample_similarity_heatmap.pdf")))
    expect_true(file.exists(file.path(save_dir, "pca_coordinates.csv")))
    expect_true(file.exists(file.path(save_dir, "pca_loadings.csv")))
    expect_true(file.exists(file.path(save_dir, "pca_top_genes_per_pc.csv")))
    expect_true(file.exists(file.path(save_dir, "sample_similarity_matrix.csv")))
    expect_false(file.exists(file.path(save_dir, "pca_plot_uncorrected.pdf")))
    expect_false(file.exists(file.path(save_dir, "pca_scree_uncorrected.pdf")))
    expect_false(file.exists(file.path(save_dir, "pca_top_loadings_uncorrected.pdf")))
})

test_that("run_distance_pca with batch correction also returns/saves the uncorrected PCA/scree/loadings", {
    dds <- simulate_dds(n_genes = 60, n_per_group = 4, batch = TRUE)
    save_dir <- withr::local_tempdir()

    out <- run_distance_pca(dds, var = "condition", batch = "batch",
                             save_dir = save_dir)

    expect_s3_class(out$pca_plot_uncorrected, "ggplot")
    expect_s3_class(out$pca_plot_uncorrected_interactive, "plotly")
    expect_s3_class(out$scree_plot_uncorrected, "ggplot")
    expect_s3_class(out$top_loadings_plot_uncorrected, "ggplot")
    expect_true(is.data.frame(out$loadings_uncorrected))
    expect_true(is.data.frame(out$top_genes_uncorrected))
    expect_true(file.exists(file.path(save_dir, "pca_plot_uncorrected.pdf")))
    expect_true(file.exists(file.path(save_dir, "pca_scree_uncorrected.pdf")))
    expect_true(file.exists(file.path(save_dir, "pca_top_loadings_uncorrected.pdf")))
    expect_true(file.exists(file.path(save_dir, "pca_loadings_uncorrected.csv")))
    expect_true(file.exists(file.path(save_dir, "pca_top_genes_per_pc_uncorrected.csv")))
    expect_true(is.data.frame(out$pca_assoc_uncorrected))
    expect_false(is.null(out$pca_assoc_plot_uncorrected))
    expect_true(file.exists(file.path(save_dir, "pca_covariate_association_uncorrected.pdf")))
    expect_true(file.exists(file.path(save_dir, "pca_covariate_association_uncorrected.csv")))
    # design is untouched by visualization-only batch correction
    expect_identical(as.character(DESeq2::design(out$dds))[2], "condition")
    expect_false(is.null(out$similarity_plot_corrected))
    expect_true(is.data.frame(out$similarity_matrix_corrected))
    expect_equal(out$batch_correction, "removeBatchEffect")
    expect_true(file.exists(file.path(save_dir, "sample_similarity_heatmap_batch_corrected.pdf")))
    expect_true(file.exists(file.path(save_dir, "sample_similarity_matrix_batch_corrected.csv")))
    # VST (before) vs corrected (after) correlation matrices should differ
    expect_false(isTRUE(all.equal(
        as.matrix(out$similarity_matrix[, -1]),
        as.matrix(out$similarity_matrix_corrected[, -1])
    )))
})

test_that("run_distance_pca removes several covariates at once (batch + pair)", {
    dds <- simulate_dds(n_genes = 60, n_per_group = 4, batch = TRUE)
    SummarizedExperiment::colData(dds)$donor <- factor(rep(paste0("d", 1:4), times = 2))
    save_dir <- withr::local_tempdir()

    out <- run_distance_pca(dds, var = "condition", batch = c("batch", "donor"),
                            save_plot = FALSE, save_dir = save_dir)

    expect_equal(out$batch_correction, "removeBatchEffect")
    expect_true(is.data.frame(out$similarity_matrix_corrected))
})

test_that("plot_pca renders a square panel (aspect.ratio = 1)", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    expect_equal(out$pca_plot$theme$aspect.ratio, 1)
})

test_that("plot_pca draws dotted zero-reference lines at PC1 = 0 and PC2 = 0", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    geom_classes <- vapply(out$pca_plot$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomHline" %in% geom_classes)
    expect_true("GeomVline" %in% geom_classes)
    hline_layer <- out$pca_plot$layers[[which(geom_classes == "GeomHline")]]
    vline_layer <- out$pca_plot$layers[[which(geom_classes == "GeomVline")]]
    expect_equal(hline_layer$aes_params$linetype, "dotted")
    expect_equal(vline_layer$aes_params$linetype, "dotted")
})

test_that("plot_pca uses 10%-both-sides multiplicative expansion on both axes", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    x_scale <- out$pca_plot$scales$get_scales("x")
    y_scale <- out$pca_plot$scales$get_scales("y")
    expect_equal(x_scale$expand, c(0.1, 0, 0.1, 0))
    expect_equal(y_scale$expand, c(0.1, 0, 0.1, 0))
})

test_that("plot_pca sets clip = \"off\" so sample point markers at the zero-padded axis extremes are not clipped", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    expect_equal(out$pca_plot$coordinates$clip, "off")
})

test_that("make_interactive converts a ggplot into a hoverable plotly widget", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    pl <- make_interactive(out$pca_plot)
    expect_s3_class(pl, "plotly")
    # tooltip text should mention the sample names present in the PCA data
    tooltip_text <- unlist(lapply(pl$x$data, function(d) d$text))
    expect_true(any(grepl(out$pca$sample[1], tooltip_text, fixed = TRUE)))
})

test_that("make_interactive(square = TRUE) is legend-aware: the *data panel* (margin-subtracted paper area) is an exact square, and the legend sits outside it in the right margin", {
    # Regression #1: ggplot2::theme(aspect.ratio = 1), set by plot_pca()/
    # plot_ma()/plot_volcano(), is silently dropped by plotly::ggplotly()
    # (no error, just not carried into the built widget). scaleanchor/
    # scaleratio (the old mechanism here) only locks a 1:1 *data-unit*
    # ratio -- it does not reproduce a square *panel* unless the x/y data
    # ranges happen to match.
    #
    # Regression #2: a naive "canvas minus margins is square" fix (an
    # earlier iteration of this helper) is *also* wrong whenever the plot
    # carries a legend -- plotly's own margin solver does not budget any
    # space for the legend at all, so "canvas minus margins" silently
    # includes the legend's width inside what should be the square data
    # panel, squashing it into a non-square rectangle. plot_pca() always
    # maps a discrete grouping variable to colour, so this dataset
    # exercises exactly that legend-bearing case.
    #
    # Regression #3: shrinking xaxis.domain and pinning legend.x =
    # domain[2] looked correct in plotly_build()'s layout object but
    # still rendered the legend over the data panel in the browser.
    # The legend now lives in the expanded right margin at legend.x > 1
    # (paper coordinates), while the axes keep domain = c(0, 1) over a
    # square paper area.
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    pl <- make_interactive(out$pca_plot, square = TRUE)
    built <- plotly::plotly_build(pl)
    m <- built$x$layout$margin
    xaxis <- built$x$layout$xaxis
    legend <- built$x$layout$legend

    expect_true(!is.null(legend))
    # axes fill the (square) paper area; legend is outside it
    expect_equal(xaxis$domain, c(0, 1))
    expect_true(legend$x > 1)
    expect_equal(legend$xanchor, "left")

    panel_w_px <- built$x$layout$width - m$l - m$r
    panel_h_px <- built$x$layout$height - m$t - m$b
    expect_equal(panel_w_px, panel_h_px, tolerance = 0.01)

    # the old data-unit-locking mechanism must not be present
    expect_null(xaxis$scaleanchor)

    # fixed pixel size, not window-responsive
    expect_false(pl$x$config$responsive)
})

test_that("make_interactive(square = TRUE) on a legend-less plot fills the full inner width with the squared panel (no legend strip reserved)", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    # strip the colour legend so this exercises the no-legend branch
    p_no_legend <- out$pca_plot + ggplot2::theme(legend.position = "none")
    pl <- make_interactive(p_no_legend, square = TRUE)
    built <- plotly::plotly_build(pl)
    m <- built$x$layout$margin
    xaxis <- built$x$layout$xaxis
    legend <- built$x$layout$legend

    expect_equal(xaxis$domain, c(0, 1))
    # no outside-margin legend placement when there is nothing to show
    expect_true(is.null(legend$x) || isTRUE(legend$x <= 1))
    inner_width <- built$x$layout$width - m$l - m$r
    panel_h_px <- built$x$layout$height - m$t - m$b
    expect_equal(inner_width, panel_h_px, tolerance = 0.01)
    expect_null(xaxis$scaleanchor)
})

test_that("make_interactive(square = FALSE) leaves the canvas size and axis layout unconstrained, but still fixed (non-responsive)", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    pl <- make_interactive(out$pca_plot, square = FALSE)
    built <- plotly::plotly_build(pl)
    expect_null(built$x$layout$xaxis$scaleanchor)
    expect_false(pl$x$config$responsive)
})

test_that("plot_sample_similarity colours high correlation red / low correlation blue, scaled to data range", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    cor_mat <- attr(out$similarity_plot, "cor_matrix")
    expect_true(is.matrix(cor_mat))
    expect_equal(dim(cor_mat), c(ncol(dds), ncol(dds)))
    expect_true(all(diag(cor_mat) > 0.999))  # self-correlation ~= 1

    # recompute the palette the same way plot_sample_similarity does and
    # check its endpoints directly (blue at low end, red at high end)
    ramp <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
    start_rgb <- grDevices::col2rgb(ramp[1])
    end_rgb <- grDevices::col2rgb(ramp[100])
    expect_true(start_rgb["blue", 1] > start_rgb["red", 1])
    expect_true(end_rgb["red", 1] > end_rgb["blue", 1])
})

test_that("compute_pca_loadings and top_genes_per_pc identify the correct top up/down genes", {
    dds <- simulate_dds(n_genes = 50, n_per_group = 4)
    out <- run_distance_pca(dds, var = "condition", n_loading_pcs = 2, n_top_genes = 5,
                             save_plot = FALSE, save_dir = withr::local_tempdir())

    expect_setequal(unique(out$top_genes$pc), c("PC1", "PC2"))
    # 5 "up" + 5 "down" per PC = 10 rows for PC1
    expect_equal(sum(out$top_genes$pc == "PC1"), 10)
    expect_setequal(unique(out$top_genes$direction), c("up", "down"))

    pc1_loadings <- out$loadings$PC1
    names(pc1_loadings) <- out$loadings$gene
    expected_up5 <- names(sort(pc1_loadings, decreasing = TRUE))[1:5]
    expected_down5 <- names(sort(pc1_loadings, decreasing = FALSE))[1:5]
    actual_up5 <- out$top_genes$gene[out$top_genes$pc == "PC1" & out$top_genes$direction == "up"]
    actual_down5 <- out$top_genes$gene[out$top_genes$pc == "PC1" & out$top_genes$direction == "down"]
    expect_setequal(actual_up5, expected_up5)
    expect_setequal(actual_down5, expected_down5)
})

test_that("plot_pca_scree caps n_pcs to the number of PCs available", {
    dds <- simulate_dds(n_genes = 20, n_per_group = 2)
    out <- run_distance_pca(dds, var = "condition", n_scree_pcs = 100,
                             save_plot = FALSE, save_dir = withr::local_tempdir())
    expect_s3_class(out$scree_plot, "ggplot")
    n_pcs_available <- length(out$percent_var)
    expect_lte(nrow(out$scree_plot$data), n_pcs_available)
})

