# Unit tests for PPI connectivity helpers (no live STRINGdb download).

test_that(".build_ppi_subgraph ranks hubs by degree and caps plot size", {
    net <- data.frame(
        from = c("A", "A", "A", "B", "C", "X"),
        to   = c("B", "C", "D", "C", "D", "Y"),
        stringsAsFactors = FALSE
    )
    sub <- deseq2reports:::.build_ppi_subgraph(
        genes = c("A", "B", "C", "D", "Z"),
        net = net,
        top_n_hubs = 2L,
        max_plot_nodes = 80L
    )
    expect_false(is.null(sub))
    expect_equal(sub$n_seed, 5L)
    expect_equal(sub$n_with_edge, 4L) # Z has no edges
    expect_equal(sub$hubs[1], "A")   # degree 3
    expect_true(all(c("A", "B", "C", "D") %in% sub$nodes$gene))
    expect_false("Z" %in% sub$nodes$gene)
    expect_equal(sum(sub$nodes$hub), 2L)
})

test_that(".build_ppi_subgraph returns NULL when no edges overlap the gene set", {
    net <- data.frame(from = "A", to = "B", stringsAsFactors = FALSE)
    expect_null(deseq2reports:::.build_ppi_subgraph(c("X", "Y"), net))
})

test_that(".ppi_seed_genes_de prefers significant genes", {
    res <- data.frame(
        gene = paste0("g", 1:6),
        sig = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
        rank_score = c(1, 2, 10, 9, 8, 7),
        stringsAsFactors = FALSE
    )
    seeds <- deseq2reports:::.ppi_seed_genes_de(res, max_genes = 10L)
    expect_setequal(seeds, c("g1", "g2"))
})

test_that("plot_ppi_connectivity and make_interactive_ppi work for a small graph", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggrepel")
    skip_if_not_installed("plotly")
    nodes <- data.frame(
        gene = c("A", "B", "C", "D"),
        degree = c(3L, 2L, 2L, 1L),
        hub = c(TRUE, TRUE, FALSE, FALSE),
        stringsAsFactors = FALSE
    )
    edges <- data.frame(
        from = c("A", "A", "A", "B"),
        to = c("B", "C", "D", "C"),
        combined_score = c(900L, 800L, 700L, 600L),
        stringsAsFactors = FALSE
    )
    p <- plot_ppi_connectivity(nodes, edges, title = "test", top_n_hubs = 2L)
    expect_s3_class(p, "ggplot")
    w <- make_interactive_ppi(p, top_n_hubs = 2L)
    expect_s3_class(w, "plotly")
})

test_that(".ppi_seed_genes_de keeps all significant genes when max_genes is Inf", {
    res <- data.frame(
        gene = paste0("g", 1:5),
        sig = TRUE,
        rank_score = 5:1,
        stringsAsFactors = FALSE
    )
    expect_equal(deseq2reports:::.ppi_seed_genes_de(res, max_genes = Inf), paste0("g", 1:5))
    expect_equal(deseq2reports:::.ppi_seed_genes_de(res, max_genes = 2L), c("g1", "g2"))
})
