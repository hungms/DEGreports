skip_if_not_installed("readxl")
skip_if_not_installed("pheatmap")

test_that("load_cspa_genes returns human and mouse symbols from bundled workbook", {
    h <- load_cspa_genes("human")
    m <- load_cspa_genes("mouse")
    expect_true(length(h) > 1000)
    expect_true(length(m) > 1000)
    expect_true("MME" %in% toupper(h) || "MME" %in% h)
    expect_true(any(toupper(m) == "MME"))
})

test_that("select_cspa_markers prefers significant CSPA genes by |rank_score|", {
    cspa <- c("GeneA", "GeneB", "GeneC", "GeneD")
    res <- data.frame(
        gene = c("GeneA", "GeneB", "GeneX", "GeneC", "GeneD"),
        rank_score = c(5, -4, 10, 0.5, -8),
        log2FoldChange = c(2, -1.5, 3, 0.2, -2),
        sig = c(TRUE, TRUE, TRUE, FALSE, TRUE),
        stringsAsFactors = FALSE
    )
    sel <- select_cspa_markers(res, cspa, n_markers = 2, significant_only = TRUE)
    expect_equal(nrow(sel), 2L)
    expect_equal(sel$gene, c("GeneD", "GeneA"))
})

test_that("run_cspa_marker_pipeline writes per-contrast and combined outputs", {
    skip_if_not_installed("ashr")
    dds <- simulate_dds_real_genes(n_genes = 80, n_per_group = 3, seed = 21)
    # Ensure some CSPA genes exist among rownames when possible
    cspa <- tryCatch(load_cspa_genes(detect_organism(rownames(dds))), error = function(e) character(0))
    if (length(cspa)) {
        rn <- rownames(dds)
        # Force-overlap a few genes for a deterministic smoke test
        take <- utils::head(cspa[toupper(cspa) %in% toupper(rn)], 5)
        if (!length(take)) {
            # Rename a few rows so they intersect CSPA (case-matched)
            n_force <- min(5L, length(cspa), nrow(dds))
            rownames(dds)[seq_len(n_force)] <- cspa[seq_len(n_force)]
        }
    }
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    tmp <- withr::local_tempdir()
    de <- run_diffexp_pipeline(dds, var = "condition", save_dir = file.path(tmp, "deseq2"),
                               save_plot = FALSE)
    out <- run_cspa_marker_pipeline(
        dds, de, var = "condition", n_markers = 5,
        save_dir = file.path(tmp, "deseq2"), save_plot = TRUE
    )
    expect_true(out$cspa_n > 0)
    comps <- names(de$by_comparison)
    expect_true(file.exists(file.path(tmp, "deseq2", "cspa_markers_combined.csv")) ||
                    length(out$combined_markers) == 0 || nrow(out$combined_markers) == 0)
    if (length(comps) && nrow(out$markers_by_comparison[[comps[[1]]]])) {
        expect_true(file.exists(file.path(
            tmp, "deseq2", comps[[1]], "cspa", "cspa_markers.csv"
        )))
    }
    expect_true("combined" %in% names(out$heatmaps))
    expect_equal(out$expr_correction, "vst")
})
