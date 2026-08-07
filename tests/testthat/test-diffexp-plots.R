# Tests for the direction-coloured MA/volcano plots (classify_de_direction,
# plot_ma, plot_volcano) and the clusterProfiler GSEA add-on
# (run_clusterprofiler_gsea, import_msigdbr_term2gene).

test_that("classify_de_direction labels up/down/ns from an existing `sig` column", {
    res_df <- data.frame(
        log2FoldChange = c(2, -2, 0.1, 3, -3),
        padj = c(0.001, 0.001, 0.9, NA, 0.02),
        sig = c(TRUE, TRUE, FALSE, NA, TRUE)
    )
    direction <- classify_de_direction(res_df)
    expect_equal(as.character(direction), c("up", "down", "ns", "ns", "down"))
})

test_that("classify_de_direction recomputes significance from alpha/lfc_threshold when `sig` is absent", {
    res_df <- data.frame(
        log2FoldChange = c(2, -2, 0.1, 0.5),
        padj = c(0.001, 0.001, 0.9, 0.001)
    )
    direction <- classify_de_direction(res_df, alpha = 0.05, lfc_threshold = 1)
    # gene 4: padj significant but |LFC| = 0.5 <= lfc_threshold = 1 -> not sig
    expect_equal(as.character(direction), c("up", "down", "ns", "ns"))
})

test_that("plot_ma colours points by direction and annotates up/down counts", {
    res_df <- mock_res_df()
    res_df$sig <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    res_df$comparison <- "A_vs_B"

    p <- plot_ma(res_df, alpha = 0.05)
    expect_s3_class(p, "ggplot")
    expect_setequal(levels(p$data$direction), c("up", "down", "ns"))
    expect_equal(sum(p$data$direction == "up"), 1)
    expect_equal(sum(p$data$direction == "down"), 1)

    # count annotations present as geom_text/annotate layers
    geom_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomText" %in% geom_classes)

    # significant points are drawn larger than non-significant ones
    expect_true(.de_direction_sizes[["up"]] > .de_direction_sizes[["ns"]])
    expect_true(.de_direction_sizes[["down"]] > .de_direction_sizes[["ns"]])
})

test_that("plot_ma and plot_volcano render square panels (aspect.ratio = 1) with the legend outside on the right", {
    res_df <- mock_res_df()
    res_df$sig <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    res_df$comparison <- "A_vs_B"
    res_df$rank_score <- res_df$log2FoldChange

    p_ma <- plot_ma(res_df, alpha = 0.05)
    expect_equal(p_ma$theme$aspect.ratio, 1)
    expect_equal(p_ma$theme$legend.position, "right")

    p_volcano <- plot_volcano(res_df, alpha = 0.05, top_label_n = 2)
    expect_equal(p_volcano$theme$aspect.ratio, 1)
    expect_equal(p_volcano$theme$legend.position, "right")
})

test_that("plot_ma keeps zero x-axis expansion but uses 10%-both-sides multiplicative expansion on the y-axis (log2FoldChange)", {
    res_df <- mock_res_df()
    res_df$sig <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    res_df$comparison <- "A_vs_B"
    res_df$rank_score <- res_df$log2FoldChange

    p_ma <- plot_ma(res_df, alpha = 0.05)
    x_scale_ma <- p_ma$scales$get_scales("x")
    y_scale_ma <- p_ma$scales$get_scales("y")
    # x (log2 mean CPM) keeps exact zero data-padding.
    expect_equal(x_scale_ma$expand, c(0, 0, 0, 0))
    # y (log2FoldChange, the axis the up/down count and direction-arrow
    # annotations are anchored near) instead uses expansion(mult = c(0.1,
    # 0.1)): 10% headroom on each side so points/annotations near the
    # extremes aren't crowded against the panel border.
    expect_equal(y_scale_ma$expand, c(0.1, 0, 0.1, 0))
})

test_that("plot_volcano uses 10%-both-sides multiplicative expansion on both axes", {
    res_df <- mock_res_df()
    res_df$sig <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    res_df$comparison <- "A_vs_B"
    res_df$rank_score <- res_df$log2FoldChange

    p_volcano <- plot_volcano(res_df, alpha = 0.05, top_label_n = 2)
    x_scale_volcano <- p_volcano$scales$get_scales("x")
    y_scale_volcano <- p_volcano$scales$get_scales("y")
    expect_equal(x_scale_volcano$expand, c(0.1, 0, 0.1, 0))
    expect_equal(y_scale_volcano$expand, c(0.1, 0, 0.1, 0))
})

test_that("plot_ma annotates up/down direction arrows on the left side, naming the test/reference comparison levels", {
    res_df <- mock_res_df()
    res_df$sig <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    res_df$comparison <- "TreatmentA_vs_Control"

    p_ma <- plot_ma(res_df, alpha = 0.05)
    text_layers <- Filter(function(l) inherits(l$geom, "GeomText"), p_ma$layers)
    labels <- vapply(text_layers, function(l) l$aes_params$label, character(1))
    xpos <- vapply(text_layers, function(l) l$data$x, numeric(1))
    # up-direction arrow names the numerator ("test") level, down-direction
    # arrow names the denominator ("reference") level
    expect_true(any(grepl("\u2191", labels) & grepl("TreatmentA", labels)))
    expect_true(any(grepl("\u2193", labels) & grepl("Control", labels)))
    # anchored on the left side of the panel (x = -Inf), distinct from the
    # right-anchored up/down count labels (x = Inf)
    left_layers <- text_layers[is.infinite(xpos) & xpos < 0]
    expect_true(length(left_layers) >= 2)
})

test_that("plot_volcano annotates up/down direction arrows just under the title (top of panel), naming the test/reference comparison levels", {
    res_df <- data.frame(
        gene = paste0("g", 1:10),
        log2FoldChange = c(3, 2.5, 2, -3, -2.5, -2, 0.1, -0.1, 0.05, 0.2),
        pvalue = c(0.0001, 0.0002, 0.0003, 0.0001, 0.0002, 0.0003, 0.9, 0.8, 0.95, 0.5),
        padj = c(0.001, 0.002, 0.003, 0.001, 0.002, 0.003, 0.9, 0.8, 0.95, 0.5),
        stat = c(6, 5, 4, -6, -5, -4, 0.1, -0.1, 0.05, 0.2),
        rank_score = c(6, 5, 4, -6, -5, -4, 0.1, -0.1, 0.05, 0.2),
        comparison = "TreatmentA_vs_Control"
    )
    p <- plot_volcano(res_df, alpha = 0.05, lfc_threshold = 0, top_label_n = 2)
    text_layers <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)
    labels <- vapply(text_layers, function(l) l$aes_params$label, character(1))
    ypos <- vapply(text_layers, function(l) l$data$y, numeric(1))
    expect_true(any(grepl("\u2192", labels) & grepl("TreatmentA", labels)))
    expect_true(any(grepl("\u2190", labels) & grepl("Control", labels)))
    # all top-of-panel annotations (count labels and direction arrows) are
    # anchored at y = Inf, i.e. at the top under the title, not mid-panel
    top_layers <- text_layers[is.infinite(ypos) & ypos > 0]
    expect_true(length(top_layers) >= 4) # 2 count labels + 2 direction arrows
})

test_that("plot_ma and plot_volcano set clip = \"off\" so point markers at the zero-padded axis extremes are not clipped", {
    res_df <- mock_res_df()
    res_df$sig <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    res_df$comparison <- "A_vs_B"
    res_df$rank_score <- res_df$log2FoldChange

    p_ma <- plot_ma(res_df, alpha = 0.05)
    expect_equal(p_ma$coordinates$clip, "off")

    p_volcano <- plot_volcano(res_df, alpha = 0.05, top_label_n = 2)
    expect_equal(p_volcano$coordinates$clip, "off")
})

test_that("plot_volcano colours points by direction, labels top up/down genes, and annotates counts", {
    res_df <- data.frame(
        gene = paste0("g", 1:10),
        log2FoldChange = c(3, 2.5, 2, -3, -2.5, -2, 0.1, -0.1, 0.05, 0.2),
        pvalue = c(0.0001, 0.0002, 0.0003, 0.0001, 0.0002, 0.0003, 0.9, 0.8, 0.95, 0.5),
        padj = c(0.001, 0.002, 0.003, 0.001, 0.002, 0.003, 0.9, 0.8, 0.95, 0.5),
        stat = c(6, 5, 4, -6, -5, -4, 0.1, -0.1, 0.05, 0.2),
        rank_score = c(6, 5, 4, -6, -5, -4, 0.1, -0.1, 0.05, 0.2),
        comparison = "A_vs_B"
    )

    p <- plot_volcano(res_df, alpha = 0.05, lfc_threshold = 0, top_label_n = 2)
    expect_s3_class(p, "ggplot")
    expect_equal(sum(p$data$direction == "up"), 3)
    expect_equal(sum(p$data$direction == "down"), 3)

    geom_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomText" %in% geom_classes)
    # ggrepel labels: exactly top_label_n up + top_label_n down = 4 labelled genes
    label_layer <- p$layers[[which(vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1)))]]
    expect_equal(nrow(label_layer$data), 4)
    expect_setequal(label_layer$data$gene, c("g1", "g2", "g4", "g5"))
})

test_that("plot_ma and plot_volcano carry a gene/log2FC/padj tooltip usable by make_interactive", {
    res_df <- mock_res_df()
    res_df$sig <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    res_df$comparison <- "A_vs_B"

    p_ma <- plot_ma(res_df, alpha = 0.05)
    expect_true(".tooltip" %in% colnames(p_ma$data))
    expect_true(all(grepl("^Gene: .*log2FC: .*padj: ", p_ma$data$.tooltip)))

    volcano_df <- data.frame(
        gene = paste0("g", 1:10),
        log2FoldChange = c(3, 2.5, 2, -3, -2.5, -2, 0.1, -0.1, 0.05, 0.2),
        pvalue = c(0.0001, 0.0002, 0.0003, 0.0001, 0.0002, 0.0003, 0.9, 0.8, 0.95, 0.5),
        padj = c(0.001, 0.002, 0.003, 0.001, 0.002, 0.003, 0.9, 0.8, 0.95, 0.5),
        stat = c(6, 5, 4, -6, -5, -4, 0.1, -0.1, 0.05, 0.2),
        rank_score = c(6, 5, 4, -6, -5, -4, 0.1, -0.1, 0.05, 0.2),
        comparison = "A_vs_B"
    )
    p_volcano <- plot_volcano(volcano_df, alpha = 0.05, top_label_n = 2)
    expect_true(".tooltip" %in% colnames(p_volcano$data))
    expect_true(all(grepl("^Gene: .*log2FC: .*padj: ", p_volcano$data$.tooltip)))
})

test_that("plot_volcano sets a symmetric x-axis limit from the larger absolute log2FoldChange bound", {
    res_df <- data.frame(
        gene = paste0("g", 1:4),
        log2FoldChange = c(5, -1, 2, -0.5),
        pvalue = c(0.001, 0.5, 0.2, 0.9),
        padj = c(0.01, 0.5, 0.3, 0.9),
        stat = c(5, -1, 2, -0.5),
        rank_score = c(5, -1, 2, -0.5),
        comparison = "A_vs_B"
    )
    p <- plot_volcano(res_df, alpha = 0.05, top_label_n = 1)
    x_limits <- ggplot2::layer_scales(p)$x$get_limits()
    expect_equal(x_limits, c(-5, 5))
})

test_that("make_interactive wraps a ggplot into a plotly widget with the modebar logo suppressed", {
    skip_if_not_installed("plotly")
    df <- data.frame(x = 1:5, y = 1:5, lab = paste0("s", 1:5))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, text = lab)) + ggplot2::geom_point()
    pl <- make_interactive(p)
    expect_s3_class(pl, "plotly")
    expect_false(pl$x$config$displaylogo)
})

test_that("make_deg_datatable returns a DT::datatable over the result table", {
    skip_if_not_installed("DT")
    res_df <- mock_res_df()
    res_df$sig <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    res_df$comparison <- "A_vs_B"

    dt <- make_deg_datatable(res_df)
    expect_s3_class(dt, "datatables")
    # unfiltered: every row of res_df is present, not just top/bottom N
    expect_equal(nrow(dt$x$data), nrow(res_df))
})

test_that("import_msigdbr_term2gene returns a term/gene data frame with the same content as import_msigdbr", {
    skip_if_not_installed("msigdbr")
    list_form <- import_msigdbr("human", "HALLMARK")
    term2gene <- import_msigdbr_term2gene("human", "HALLMARK")

    expect_s3_class(term2gene, "data.frame")
    expect_setequal(colnames(term2gene), c("term", "gene"))
    expect_setequal(unique(term2gene$term), names(list_form))
    one_term <- names(list_form)[1]
    expect_setequal(term2gene$gene[term2gene$term == one_term], list_form[[one_term]])
})

test_that("run_clusterprofiler_gsea returns an unfiltered gseaResult object", {
    skip_if_not_installed("clusterProfiler")
    set.seed(42)
    genes <- paste0("G", 1:300)
    ranks <- sort(stats::rnorm(300), decreasing = TRUE)
    names(ranks) <- genes
    term2gene <- data.frame(
        term = rep(paste0("SET", 1:10), each = 25),
        gene = sample(genes, 250, replace = TRUE)
    )

    res <- run_clusterprofiler_gsea(ranks, term2gene, min_size = 5, max_size = 500)
    expect_s4_class(res, "gseaResult")
    res_df <- as.data.frame(res)
    # pvalueCutoff = 1 internally -> every tested gene set is retained,
    # not just significant ones
    expect_true(any(res_df$p.adjust > 0.05))
})

test_that("run_clusterprofiler_gsea warns and returns NULL when no gene set clears min_size", {
    ranks <- stats::setNames(sort(stats::rnorm(20), decreasing = TRUE), paste0("G", 1:20))
    term2gene <- data.frame(term = "BIGSET", gene = paste0("OTHER", 1:100))

    expect_warning(
        res <- run_clusterprofiler_gsea(ranks, term2gene, min_size = 50, max_size = 500),
        "no gene sets"
    )
    expect_null(res)
})
