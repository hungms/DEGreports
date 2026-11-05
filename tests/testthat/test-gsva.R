test_that("GSVA pathway labels wrap at 40 characters", {
    short <- "HALLMARK INTERFERON GAMMA RESPONSE"
    expect_equal(.gsva_wrap_labels(short), short)

    long <- paste(
        "REGULATION OF TRANSCRIPTION FROM RNA POLYMERASE II PROMOTER",
        "IN RESPONSE TO STRESS"
    )
    wrapped <- .gsva_wrap_labels(long, width = 40L)
    expect_true(grepl("\n", wrapped, fixed = TRUE))
    lines <- strsplit(wrapped, "\n", fixed = TRUE)[[1]]
    expect_true(all(nchar(lines) <= 40L))
    expect_equal(gsub("\\s+", " ", gsub("\n", " ", wrapped)),
                 gsub("\\s+", " ", long))
})

test_that("GSVA heatmap keeps all HALLMARK sets and top-n for large collections", {
    stats_df <- data.frame(
        pathway = paste0("P", 1:8),
        variance = c(0.1, 0.9, 0.2, 0.3, 0.4, 0.5, 0.05, 0.8),
        p_group = c(0.4, 0.001, 0.2, 0.05, 0.8, 0.01, 0.9, 0.002),
        stringsAsFactors = FALSE
    )
    all_hm <- .gsva_heatmap_pathways(stats_df, n = 3, plot_all = TRUE)
    expect_equal(all_hm, stats_df$pathway)

    top3 <- .gsva_heatmap_pathways(stats_df, n = 3, plot_all = FALSE)
    expect_equal(top3, c("P2", "P8", "P6"))
})

test_that("GSVA pathway stats rank by Kruskal-Wallis vs grouping factor", {
    set.seed(1)
    scores <- rbind(
        A = c(1, 1.1, 1.2, 5, 5.1, 5.2),
        B = c(3, 3, 3, 3.1, 2.9, 3),
        C = rnorm(6)
    )
    colnames(scores) <- paste0("s", 1:6)
    group <- factor(rep(c("x", "y"), each = 3))
    st <- .gsva_pathway_stats(scores, group)
    expect_equal(st$pathway[[1]], "A")
    expect_true(st$p_group[st$pathway == "A"] < st$p_group[st$pathway == "B"])
})

test_that("run_gsva_pipeline scores supplied gene sets and heatmaps HALLMARK fully", {
    skip_if_not_installed("GSVA")
    skip_if_not_installed("DESeq2")
    dds <- simulate_dds(n_genes = 40, n_per_group = 3, n_de = 10, seed = 4)
    dds <- DESeq2::estimateSizeFactors(dds)
    genes <- rownames(dds)
    gene_sets <- list(
        HALLMARK = list(
            HALLMARK_SET1 = genes[1:12],
            HALLMARK_SET2 = genes[5:20],
            HALLMARK_SET3 = genes[15:30]
        ),
        GOBP = setNames(
            lapply(1:12, function(i) genes[seq(i, min(i + 8, length(genes)))]),
            paste0("GOBP_SET", 1:12)
        )
    )
    save_dir <- withr::local_tempdir()
    out <- run_gsva_pipeline(
        dds, var = "condition",
        gene_sets = gene_sets,
        n_heatmap = 4,
        min_size = 5,
        save_dir = save_dir,
        save_plot = FALSE
    )
    expect_true("HALLMARK" %in% out$collections)
    expect_true("GOBP" %in% out$collections)
    expect_equal(nrow(out$scores$HALLMARK), 3L)
    expect_true(nrow(out$scores$GOBP) >= 4L)
    expect_true(file.exists(file.path(save_dir, "HALLMARK", "gsva_scores.csv")))
    expect_true(file.exists(file.path(save_dir, "GOBP", "gsva_pathway_stats.csv")))
    expect_equal(out$expr_correction, "vst")
})

test_that("sample_level_expr_matrix applies removeBatchEffect on VST when batch is set", {
    dds <- simulate_dds(n_genes = 40, n_per_group = 4, batch = TRUE)
    raw <- sample_level_expr_matrix(dds, var = "condition", batch = NULL)
    cor <- sample_level_expr_matrix(dds, var = "condition", batch = "batch")
    expect_equal(dim(raw), dim(cor))
    expect_equal(colnames(raw), colnames(cor))
    expect_equal(attr(raw, "expr_correction"), "vst")
    expect_equal(attr(cor, "expr_correction"), "rbe_vst")
    expect_false(isTRUE(all.equal(unname(raw), unname(cor), check.attributes = FALSE)))
})

test_that("expression-correction labels describe each correction", {
    expect_equal(batch_correction_label("rbe_vst"), "removeBatchEffect")
    expect_equal(batch_correction_label("rbe_log2expr"), "removeBatchEffect")
    expect_equal(batch_correction_label("vst"), "none")
    expect_equal(batch_correction_label(NA_character_), "none")
    expect_equal(expr_correction_label("rbe_vst"), "removeBatchEffect on VST")
    expect_equal(expr_correction_label("rbe_log2expr"), "removeBatchEffect on log2 expression")
    expect_equal(expr_correction_label("vst"), "blind VST")
})
