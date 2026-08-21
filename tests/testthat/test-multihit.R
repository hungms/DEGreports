test_that("rank_multihit_genes_contrast scores DEG + GSEA + TF layers", {
    de <- data.frame(
        gene = c("GeneA", "GeneB", "GeneC", "GeneD"),
        comparison = "A_vs_B",
        sig = c(TRUE, TRUE, FALSE, FALSE),
        padj = c(0.01, 0.02, 0.5, 0.8),
        log2FoldChange = c(2, -1.5, 0.1, 0),
        stat = c(5, -4, 0.2, 0.1),
        stringsAsFactors = FALSE
    )
    diffexp_out <- list(by_comparison = list(A_vs_B = de), results = de)

    gsea_out <- list(
        by_comparison = list(
            A_vs_B = list(
                HALLMARK = data.frame(
                    pathway = c("P1", "P2"),
                    NES = c(2, 1),
                    padj = c(0.01, 0.02),
                    leadingEdge = c("GeneA,GeneX", "GeneA,GeneB"),
                    stringsAsFactors = FALSE
                )
            )
        ),
        org = "human"
    )

    viper_out <- list(
        contrast_acts = list(
            A_vs_B = data.frame(
                tf = c("GeneA", "GeneZ"),
                padj = c(0.01, 0.9),
                stringsAsFactors = FALSE
            )
        ),
        org = "human"
    )

    decoupler_tf_out <- list(
        contrast_acts = list(
            A_vs_B = data.frame(
                source = c("GeneB", "GeneY"),
                p_value = c(0.01, 0.5),
                stringsAsFactors = FALSE
            )
        ),
        org = "human"
    )

    ranked <- rank_multihit_genes_contrast(
        "A_vs_B", diffexp_out, gsea_out, viper_out, decoupler_tf_out,
        alpha = 0.05, top_n = 10, org = "human"
    )
    expect_true(nrow(ranked) >= 1L)
    expect_true("direction" %in% colnames(ranked))
    # GeneA: DEG + GSEA (once) + VIPER = 3; up by LFC
    row_a <- ranked[ranked$gene == "GeneA", , drop = FALSE]
    expect_equal(nrow(row_a), 1L)
    expect_equal(row_a$deg[[1]], 1L)
    expect_equal(row_a$gsea[[1]], 1L)
    expect_equal(row_a$viper[[1]], 1L)
    expect_equal(row_a$cspa[[1]], 0L)
    expect_equal(row_a$score[[1]], 3L)
    expect_equal(as.character(row_a$direction[[1]]), "up")
    # GeneB: down by LFC
    row_b <- ranked[ranked$gene == "GeneB", , drop = FALSE]
    expect_equal(as.character(row_b$direction[[1]]), "down")
    expect_equal(ranked$gene[[1]], "GeneA")
})

test_that("rank_multihit_genes_contrast returns top_n per DEG direction", {
    de <- data.frame(
        gene = paste0("G", 1:8),
        comparison = "A_vs_B",
        sig = rep(TRUE, 8),
        padj = rep(0.01, 8),
        log2FoldChange = c(3, 2, 1, 0.5, -0.5, -1, -2, -3),
        # scores will all be DEG-only (=1); |stat| breaks ties
        stat = c(8, 7, 6, 5, -5, -6, -7, -8),
        stringsAsFactors = FALSE
    )
    diffexp_out <- list(by_comparison = list(A_vs_B = de), results = de)
    ranked <- rank_multihit_genes_contrast(
        "A_vs_B", diffexp_out, alpha = 0.05, top_n = 2, org = "human"
    )
    expect_equal(nrow(ranked), 4L)
    expect_equal(as.character(ranked$direction), c("up", "up", "down", "down"))
    expect_equal(ranked$gene, c("G1", "G2", "G8", "G7"))
})

test_that("rank_multihit_genes_contrast counts CSPA markers", {
    de <- data.frame(
        gene = c("Surf1", "Surf2", "Other"),
        comparison = "A_vs_B",
        sig = c(TRUE, FALSE, TRUE),
        padj = c(0.01, 0.2, 0.02),
        log2FoldChange = c(2, 0.1, -1),
        stat = c(5, 0.2, -3),
        stringsAsFactors = FALSE
    )
    diffexp_out <- list(by_comparison = list(A_vs_B = de), results = de)
    cspa_out <- list(
        markers_by_comparison = list(
            A_vs_B = data.frame(gene = c("Surf1", "Surf2"), stringsAsFactors = FALSE)
        )
    )
    ranked <- rank_multihit_genes_contrast(
        "A_vs_B", diffexp_out, cspa_out = cspa_out,
        alpha = 0.05, top_n = 10, org = "human"
    )
    row1 <- ranked[ranked$gene == "Surf1", , drop = FALSE]
    expect_equal(row1$cspa[[1]], 1L)
    expect_equal(row1$deg[[1]], 1L)
    expect_equal(row1$score[[1]], 2L)
    expect_equal(as.character(row1$direction[[1]]), "up")
    row2 <- ranked[ranked$gene == "Surf2", , drop = FALSE]
    expect_equal(row2$cspa[[1]], 1L)
    expect_equal(row2$deg[[1]], 0L)
    expect_equal(as.character(row2$direction[[1]]), "up")
})

test_that("run_multihit_pipeline writes outputs and plots", {
    de <- data.frame(
        gene = paste0("G", 1:5),
        comparison = "A_vs_B",
        sig = c(TRUE, TRUE, TRUE, FALSE, FALSE),
        padj = c(0.01, 0.02, 0.03, 0.5, 0.6),
        log2FoldChange = c(2, 1, -1, 0.1, 0),
        stat = c(4, 3, -2, 0.1, 0),
        stringsAsFactors = FALSE
    )
    diffexp_out <- list(by_comparison = list(A_vs_B = de), results = de)
    save_dir <- withr::local_tempdir()
    out <- run_multihit_pipeline(
        diffexp_out, alpha = 0.05, top_n = 3,
        save_dir = save_dir, save_plot = TRUE
    )
    expect_equal(out$layers_used, "DEG")
    # top_n=3 per direction: 2 up DEGs + 1 down DEG
    expect_equal(nrow(out$by_comparison$A_vs_B), 3L)
    expect_true(all(c("up", "down") %in% as.character(out$by_comparison$A_vs_B$direction)))
    expect_s3_class(out$plots$A_vs_B, "ggplot")
    expect_true(file.exists(file.path(save_dir, "A_vs_B", "A_vs_B_top_genes.csv")))
})
