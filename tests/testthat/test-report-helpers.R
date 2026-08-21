# Tests for report-support helpers: the before/after gene-distribution
# plot, the fixed-aspect-ratio heatmap sizing, the GSEA top/bottom-by-NES
# table + significance stars, and the DE top/bottom-by-rank_score table.

test_that("plot_count_distribution overlays before/after densities on the same axes", {
    gene_stat <- c(1, 2, 3, 4, 5, 6, 7, 8)
    keep <- c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
    p <- plot_count_distribution(gene_stat, cutoff_log2 = 2.5, keep = keep, xlab = "log2 mean CPM")

    expect_s3_class(p, "ggplot")
    # one row per gene for "Before filtering", plus one row per *retained*
    # gene for "After filtering" -- i.e. more rows than genes when some are
    # removed, all on the same `value`/`stage` columns.
    expect_setequal(unique(p$data$stage), c("Before filtering", "After filtering"))
    expect_equal(sum(p$data$stage == "Before filtering"), length(gene_stat))
    expect_equal(sum(p$data$stage == "After filtering"), sum(keep))
    # cutoff line present
    expect_true(length(p$layers) >= 2)
})

test_that("plot_count_distribution omits the cutoff line when cutoff_log2 is NA", {
    gene_stat <- c(1, 2, 3, 4)
    p <- plot_count_distribution(gene_stat, cutoff_log2 = NA_real_, keep = rep(TRUE, 4), xlab = "x")
    geom_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_false("GeomVline" %in% geom_classes)
})

test_that("pheatmap_cell_size returns a single value clamped between floor and ceiling", {
    expect_equal(pheatmap_cell_size(1000), 12)   # clamped to min_size for a huge n
    expect_equal(pheatmap_cell_size(1), 40)      # clamped to max_size for n = 1
    expect_true(pheatmap_cell_size(8) >= 12 && pheatmap_cell_size(8) <= 40)
})

test_that("plot_sample_similarity renders an aspect-ratio-1 data matrix", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())

    gt <- out$similarity_plot$gtable
    idx <- which(gt$layout$name == "matrix")
    w <- grid::convertWidth(gt$widths[gt$layout$l[idx]], "in", valueOnly = TRUE)
    h <- grid::convertHeight(gt$heights[gt$layout$t[idx]], "in", valueOnly = TRUE)
    expect_equal(w, h, tolerance = 1e-6)
})

test_that("plot_sample_similarity blanks the self-correlation diagonal (rendered grey, not on the correlation scale) but keeps the full matrix attached", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())

    # the attached matrix (used for the CSV export / other tests) keeps the
    # real diagonal (self-correlation ~= 1)
    cor_mat <- attr(out$similarity_plot, "cor_matrix")
    n <- nrow(cor_mat)
    expect_true(all(diag(cor_mat) > 0.999))

    # the matrix panel actually drawn by pheatmap() renders the diagonal
    # cells in the dedicated na_col (grey20), not a colour off the
    # correlation-scale gradient -- i.e. self-correlation is genuinely
    # blanked out, not just numerically close to a real gradient colour.
    idx <- which(out$similarity_plot$gtable$layout$name == "matrix")
    matrix_grob <- out$similarity_plot$gtable$grobs[[idx]]
    rect_grob <- matrix_grob$children[[1]]
    fills <- matrix(rect_grob$gp$fill, nrow = n, ncol = n)
    expect_true(all(diag(fills) == "grey20"))
    expect_true(all(fills[upper.tri(fills)] != "grey20"))
    expect_true(all(fills[lower.tri(fills)] != "grey20"))
})

test_that("save_pheatmap auto-sizes the PDF device from the gtable's physical size", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 3)
    out <- run_distance_pca(dds, var = "condition", save_plot = FALSE, save_dir = withr::local_tempdir())
    path <- tempfile(fileext = ".pdf")
    save_pheatmap(out$similarity_plot, path)
    expect_true(file.exists(path))
})

test_that("padj_to_stars maps significance thresholds to conventional stars", {
    expect_equal(padj_to_stars(c(0.0001, 0.005, 0.03, 0.5, NA)),
                 c("***", "**", "*", "", ""))
})

test_that("gsea_top_bottom_table returns top_n per direction ranked by NES, not filtered by significance", {
    gsea_res <- data.frame(
        pathway = paste0("PATH", 1:10),
        NES = c(3, 2.5, 2, 1.5, 1, -1, -1.5, -2, -2.5, -3),
        padj = c(0.001, 0.2, 0.5, 0.9, 0.99, 0.99, 0.9, 0.5, 0.2, 0.001),
        leadingEdge = paste0("geneA", 1:10, ",geneB", 1:10)
    )
    out <- gsea_top_bottom_table(gsea_res, top_n = 3)

    expect_equal(nrow(out), 6)
    expect_true("leadingEdge" %in% colnames(out))
    # top 3 by NES (most positive) and bottom 3 (most negative) present,
    # regardless of padj
    expect_setequal(out$pathway, c("PATH1", "PATH2", "PATH3", "PATH8", "PATH9", "PATH10"))
})

test_that("plot_gsea_barplot shows top/bottom pathways by NES regardless of significance by default", {
    gsea_res <- data.frame(
        pathway = paste0("PATH", 1:6),
        NES = c(2.5, 2, 1.5, -1.5, -2, -2.5),
        padj = c(0.8, 0.9, 0.95, 0.95, 0.9, 0.8),  # none significant
        pval = c(0.8, 0.9, 0.95, 0.95, 0.9, 0.8)
    )
    p <- plot_gsea_barplot(gsea_res, top_n = 3, signif_only = FALSE)
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 6)
    expect_true("stars" %in% colnames(p$data))
    expect_true(all(p$data$stars == ""))  # none significant, so no stars

    expect_error(
        plot_gsea_barplot(gsea_res, top_n = 3, signif_only = TRUE),
        "No pathways to plot"
    )
})

test_that("plot_gsea_barplot draws black-outlined bars, an optional bold title, and pathway names as in-plot text", {
    gsea_res <- data.frame(
        pathway = paste0("PATH", 1:6),
        NES = c(2.5, 2, 1.5, -1.5, -2, -2.5),
        padj = c(0.001, 0.9, 0.95, 0.95, 0.9, 0.001),
        pval = c(0.001, 0.9, 0.95, 0.95, 0.9, 0.001)
    )
    p <- plot_gsea_barplot(gsea_res, top_n = 3, title = "My Comparison")

    bar_layer <- p$layers[[which(vapply(p$layers, function(l) inherits(l$geom, "GeomCol"), logical(1)))]]
    expect_equal(bar_layer$aes_params$colour, "black")

    # pathway names are drawn via geom_text (using the derived
    # `pathway_label` column) rather than as axis labels
    text_layers <- p$layers[vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1))]
    pathway_label_layer <- text_layers[[which(vapply(text_layers, function(l) {
        "label" %in% names(l$mapping) &&
            identical(rlang::quo_get_expr(l$mapping$label), quote(.data$pathway_label))
    }, logical(1)))]]
    expect_true(all(p$data$pathway %in% gsea_res$pathway))
    expect_true(all(p$data$pathway_label == p$data$pathway))  # no HALLMARK_ prefix to strip here

    expect_equal(p$labels$title, "My Comparison")
    expect_true(inherits(p$theme$axis.text.y, "element_blank"))

    # a y = 0 reference line is drawn (geom_hline), and no panel border /
    # pathway (y-)axis line remains -- only the bottom NES-value axis line
    expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomHline"), logical(1))))
    expect_true(inherits(p$theme$panel.border, "element_blank"))
    expect_true(inherits(p$theme$axis.line.y, "element_blank"))
})

test_that("plot_gsea_barplot uses symmetric x-axis limits (|xmin| == |xmax|, from 0)", {
    gsea_res <- data.frame(
        pathway = paste0("PATH", 1:6),
        NES = c(3.5, 2, 1.5, -1.5, -2, -2.2),  # asymmetric NES range: max 3.5, min -2.2
        padj = c(0.001, 0.9, 0.95, 0.95, 0.9, 0.001),
        pval = c(0.001, 0.9, 0.95, 0.95, 0.9, 0.001)
    )
    p <- plot_gsea_barplot(gsea_res, top_n = 3)
    y_limits <- ggplot2::layer_scales(p)$y$get_limits()
    # scale_y_continuous(limits=) is set pre-coord_flip on the NES axis;
    # after coord_flip the panel's x-axis is this same scale
    expect_equal(abs(y_limits[1]), abs(y_limits[2]), tolerance = 1e-6)
})

test_that("plot_gsea_barplot hjusts significance-star labels away from the bar tip based on direction", {
    gsea_res <- data.frame(
        pathway = paste0("PATH", 1:4),
        NES = c(2.5, 1.5, -1.5, -2.5),
        padj = c(0.001, 0.9, 0.9, 0.001),
        pval = c(0.001, 0.9, 0.9, 0.001)
    )
    p <- plot_gsea_barplot(gsea_res, top_n = 2)
    expect_true("star_hjust" %in% colnames(p$data))
    up <- p$data[p$data$NES > 0, ]
    down <- p$data[p$data$NES < 0, ]
    expect_true(all(up$star_hjust == 0))
    expect_true(all(down$star_hjust == 1))
})

test_that("plot_gsea_barplot strips the HALLMARK_ prefix from displayed pathway names without altering the data", {
    gsea_res <- data.frame(
        pathway = c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_APOPTOSIS",
                    "HALLMARK_ANGIOGENESIS", "HALLMARK_COAGULATION"),
        NES = c(3, 2, -2, -3),
        padj = c(0.001, 0.01, 0.01, 0.001)
    )
    p <- plot_gsea_barplot(gsea_res, top_n = 2)

    expect_setequal(p$data$pathway, gsea_res$pathway)  # underlying data untouched
    expect_false(any(grepl("^HALLMARK_", p$data$pathway_label)))
    expect_setequal(p$data$pathway_label,
                     c("INTERFERON_GAMMA_RESPONSE", "APOPTOSIS", "ANGIOGENESIS", "COAGULATION"))
})

test_that("make_gsea_datatable returns a DT::datatable, styled like make_deg_datatable, over the full result table", {
    skip_if_not_installed("DT")
    gsea_res <- data.frame(
        pathway = paste0("HALLMARK_PATH", 1:6),
        NES = c(2.5, 2, 1.5, -1.5, -2, -2.5),
        padj = c(0.001, 0.9, 0.95, 0.95, 0.9, 0.001),
        leadingEdge = paste0("g", 1:6)
    )
    dt <- make_gsea_datatable(gsea_res)
    expect_s3_class(dt, "datatables")
    # unfiltered: every row is present, not just top/bottom N
    expect_equal(nrow(dt$x$data), nrow(gsea_res))
    expect_true("csv" %in% unlist(dt$x$options$buttons))
})

test_that("emit_widget cat()s the widget's tag HTML and registers its dependencies with knitr", {
    skip_if_not_installed("DT")
    skip_if_not_installed("knitr")
    df <- data.frame(x = 1:5, y = rnorm(5))
    w <- DT::datatable(df)

    knitr::knit_meta(class = NULL, clean = TRUE)  # reset any accumulated meta from other tests
    out <- capture.output(emit_widget(w))
    html <- paste(out, collapse = "\n")

    # the widget's container div/JSON payload was written to the asis stream
    expect_true(grepl("htmlwidget", html, fixed = TRUE))
    expect_true(grepl("datatables", html, fixed = TRUE))

    # and its HTML dependencies (the DataTables JS library itself) were
    # registered via knit_meta_add(), not silently dropped -- this is the
    # mechanism that keeps DT/plotly widgets emitted from a results='asis'
    # loop (as the bundled report template's per-comparison sections do)
    # from rendering as dead, unstyled placeholders in the browser.
    meta <- knitr::knit_meta(class = NULL, clean = TRUE)
    dep_names <- vapply(meta, function(d) if (inherits(d, "html_dependency")) d$name else NA_character_,
                         character(1))
    expect_true(any(grepl("datatables", dep_names, ignore.case = TRUE), na.rm = TRUE))
})

test_that("emit_widget works for a plotly widget too (generic htmlwidget dependency fix, not DT-specific)", {
    skip_if_not_installed("plotly")
    skip_if_not_installed("knitr")
    df <- data.frame(x = 1:5, y = rnorm(5))
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point()
    w <- plotly::ggplotly(p)

    knitr::knit_meta(class = NULL, clean = TRUE)
    out <- capture.output(emit_widget(w))
    html <- paste(out, collapse = "\n")
    expect_true(grepl("htmlwidget", html, fixed = TRUE))

    meta <- knitr::knit_meta(class = NULL, clean = TRUE)
    dep_names <- vapply(meta, function(d) if (inherits(d, "html_dependency")) d$name else NA_character_,
                         character(1))
    expect_true(any(grepl("plotly", dep_names, ignore.case = TRUE), na.rm = TRUE))
})

test_that("emit_contrast_nav heads each contrast class and links its own page", {
    specs <- data.frame(
        label = c("C_vs_noC", "C_vs_noC/CB__noK", "C_vs_noC_x_K_vs_noK_in_CB"),
        class = c("main", "within_B", "interaction"),
        stringsAsFactors = FALSE
    )
    md <- paste(capture.output(
        emit_contrast_nav(specs, c("celltype", "C_mutation", "K_mutation"))
    ), collapse = "\n")

    expect_true(grepl("#### Main effects (across all other factors)", md, fixed = TRUE))
    expect_true(grepl("#### C_vs_noC — by stratum", md, fixed = TRUE))
    expect_true(grepl("#### Interactions (difference of differences)", md, fixed = TRUE))
    expect_false(grepl("#### Other", md, fixed = TRUE))
    expect_true(grepl("(C_vs_noC/CB__noK/CB__noK.html)", md, fixed = TRUE))
    # Underscores in the link text are escaped so pandoc does not read a
    # stem like CB__noK as emphasis.
    expect_true(grepl("[**CB\\_\\_noK**]", md, fixed = TRUE))
})

test_that("emit_contrast_nav leaves a single-factor grid as a flat list", {
    specs <- data.frame(label = c("B_vs_A", "C_vs_A"), class = "single",
                        stringsAsFactors = FALSE)
    md <- paste(capture.output(emit_contrast_nav(specs, "condition")), collapse = "\n")
    expect_false(grepl("####", md, fixed = TRUE))
    expect_true(grepl("(B_vs_A/B_vs_A.html)", md, fixed = TRUE))
})

test_that("interaction_question_table states each interaction as a question", {
    specs <- data.frame(
        label = c("C_vs_noC", "C_vs_noC_x_K_vs_noK_in_CB",
                  "CA_vs_CB+CC_x_K_vs_noK_in_C",
                  "CA_vs_CB+CC_x_C_vs_noC_x_K_vs_noK"),
        class = c("main", "interaction", "interaction", "interaction"),
        stratum_var = c(NA, "celltype", "C_mutation", NA),
        stratum_level = c(NA, "CB", "C", NA),
        ix_terms = c(
            NA,
            "C_mutation:C_vs_noC|K_mutation:K_vs_noK",
            "celltype:CA_vs_CB+CC|K_mutation:K_vs_noK",
            "celltype:CA_vs_CB+CC|C_mutation:C_vs_noC|K_mutation:K_vs_noK"
        ),
        stringsAsFactors = FALSE
    )
    tbl <- interaction_question_table(specs, link = FALSE)

    # Main effects are not interactions and drop out.
    expect_equal(tbl$contrast, specs$label[-1])
    expect_equal(
        tbl$question[[1]],
        paste0("Does the K_mutation effect (`K_vs_noK`) differ between ",
               "C_mutation C and noC, within celltype = CB?")
    )
    # A pooled side reads as an average, not as a merged group.
    expect_true(grepl("between celltype CA and the average of CB and CC, within C_mutation = C?",
                      tbl$question[[2]], fixed = TRUE))
    # The three-way row asks about the interaction itself, not about an effect.
    expect_true(grepl("C_mutation \u00d7 K_mutation interaction", tbl$question[[3]],
                      fixed = TRUE))
    expect_false(grepl("within", tbl$question[[3]], fixed = TRUE))
    # Both readings of a positive log2FC, since the pair is symmetric.
    expect_true(grepl("equivalently", tbl$positive[[1]], fixed = TRUE))

    linked <- interaction_question_table(specs)
    expect_equal(
        linked$contrast[[2]],
        "[`CA_vs_CB+CC_x_K_vs_noK_in_C`](CA_vs_CB+CC_x_K_vs_noK_in_C/CA_vs_CB+CC_x_K_vs_noK_in_C.html)"
    )
})

test_that("interaction_question_table is empty without interaction rows", {
    flat <- data.frame(label = c("B_vs_A", "C_vs_A"), class = "single",
                       stringsAsFactors = FALSE)
    expect_equal(nrow(interaction_question_table(flat)), 0L)
    expect_named(interaction_question_table(flat),
                 c("contrast", "question", "positive"))
    expect_equal(nrow(interaction_question_table(NULL)), 0L)
})

test_that("diffexp_top_bottom_table returns top_n genes per direction ranked by rank_score", {
    res_df <- data.frame(
        gene = paste0("g", 1:10),
        rank_score = c(5, 4, 3, 2, 1, -1, -2, -3, -4, -5),
        comparison = "A_vs_B"
    )
    out <- diffexp_top_bottom_table(res_df, top_n = 3)

    expect_equal(nrow(out), 6)
    expect_setequal(out$gene, c("g1", "g2", "g3", "g8", "g9", "g10"))
    # ordered by descending rank_score
    expect_true(all(diff(out$rank_score) <= 0))
})
