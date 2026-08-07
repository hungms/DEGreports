test_that("jaccard_matrix returns known values and is symmetric", {
    sets <- list(
        A = c("g1", "g2", "g3"),
        B = c("g2", "g3", "g4"),
        C = c("g9")
    )
    jac <- jaccard_matrix(sets)
    expect_equal(dim(jac), c(3L, 3L))
    expect_equal(jac["A", "A"], 1)
    expect_equal(jac["A", "B"], 2 / 4)  # {g2,g3} / {g1..g4}
    expect_equal(jac["A", "C"], 0)
    expect_equal(jac["A", "B"], jac["B", "A"])
})

test_that("pool_sig_gsea_pathways filters by padj and pools collections", {
    by_coll <- list(
        HALLMARK = data.frame(
            pathway = c("P1", "P2"),
            NES = c(2, 1),
            padj = c(0.01, 0.2),
            size = c(10L, 20L),
            leadingEdge = c("a,b,c", "d,e"),
            stringsAsFactors = FALSE
        ),
        GOBP = data.frame(
            pathway = "P3",
            NES = -1.5,
            padj = 0.04,
            size = 12L,
            leadingEdge = "a,x",
            stringsAsFactors = FALSE
        )
    )
    sig <- pool_sig_gsea_pathways(by_coll, alpha = 0.05)
    expect_equal(nrow(sig), 2L)
    expect_true(all(sig$padj < 0.05))
    expect_setequal(sig$pathway, c("P1", "P3"))
    expect_true(all(c("label", "n_leadingEdge", "collection") %in% colnames(sig)))
})

test_that("run_gsea_overlap_pipeline soft-handles <2 sig pathways", {
    gsea_out <- list(
        by_comparison = list(
            A_vs_B = list(
                HALLMARK = data.frame(
                    pathway = "ONLY",
                    NES = 2,
                    padj = 0.01,
                    size = 5L,
                    leadingEdge = "g1,g2",
                    stringsAsFactors = FALSE
                )
            )
        )
    )
    out <- run_gsea_overlap_pipeline(gsea_out, save_dir = withr::local_tempdir(), save_plot = FALSE)
    expect_equal(out$by_comparison$A_vs_B$n_sig, 1L)
    expect_null(out$by_comparison$A_vs_B$jaccard)
    expect_null(out$by_comparison$A_vs_B$heatmap)
})

test_that("run_gsea_overlap_pipeline builds heatmap for >=2 pathways", {
    gsea_out <- list(
        by_comparison = list(
            A_vs_B = list(
                HALLMARK = data.frame(
                    pathway = c("P1", "P2", "P3"),
                    NES = c(2, 1.5, -1),
                    padj = c(0.01, 0.02, 0.03),
                    size = c(5L, 6L, 7L),
                    leadingEdge = c("g1,g2", "g2,g3", "g9"),
                    stringsAsFactors = FALSE
                )
            )
        )
    )
    save_dir <- withr::local_tempdir()
    out <- run_gsea_overlap_pipeline(gsea_out, save_dir = save_dir, save_plot = TRUE)
    expect_equal(nrow(out$by_comparison$A_vs_B$jaccard), 3L)
    expect_equal(out$by_comparison$A_vs_B$n_heatmap, 3L)
    expect_false(is.null(out$by_comparison$A_vs_B$heatmap))
    expect_true(file.exists(file.path(
        save_dir, "A_vs_B", "gsea_overlap", "A_vs_B_pathway_jaccard.csv"
    )))
    expect_true(file.exists(file.path(
        save_dir, "A_vs_B", "gsea_overlap", "A_vs_B_sig_pathways.csv"
    )))
    expect_true(file.exists(file.path(
        save_dir, "A_vs_B", "gsea_overlap", "A_vs_B_pathway_jaccard_heatmap.pdf"
    )))
})

test_that("run_gsea_overlap_pipeline caps heatmap pathways by |NES|", {
    pathways <- paste0("P", 1:60)
    gsea_out <- list(
        by_comparison = list(
            A_vs_B = list(
                HALLMARK = data.frame(
                    pathway = pathways,
                    NES = seq(3, 0.1, length.out = 60),
                    padj = rep(0.01, 60),
                    size = rep(10L, 60),
                    leadingEdge = paste0("g", seq_len(60), ",gX"),
                    stringsAsFactors = FALSE
                )
            )
        )
    )
    out <- run_gsea_overlap_pipeline(
        gsea_out, max_heatmap_pathways = 40,
        save_dir = withr::local_tempdir(), save_plot = FALSE
    )
    expect_equal(out$by_comparison$A_vs_B$n_sig, 60L)
    expect_equal(out$by_comparison$A_vs_B$n_heatmap, 40L)
    expect_equal(nrow(out$by_comparison$A_vs_B$jaccard), 60L)
    expect_equal(nrow(attr(out$by_comparison$A_vs_B$heatmap, "jaccard_matrix")), 40L)
})

test_that("plot_pathway_jaccard_heatmap returns clustered heatmap without name labels", {
    sets <- list(
        HALLMARK_A = c("g1", "g2"),
        HALLMARK_B = c("g2", "g3"),
        REACTOME_C = c("g1", "g3")
    )
    jac <- jaccard_matrix(sets)
    heat <- plot_pathway_jaccard_heatmap(jac, term_labels = names(sets))
    skip_if(is.null(heat))
    engine <- attr(heat, "plot_engine")
    expect_true(engine %in% c("pheatmap", "complexheatmap"))
    expect_equal(nrow(attr(heat, "jaccard_matrix")), 3L)
    if (identical(engine, "complexheatmap")) {
        expect_false(is.null(attr(heat, "clusters")))
        expect_length(attr(heat, "clusters"), 3L)
    }
})
