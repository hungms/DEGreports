test_that("compute_ranking_score('pxfc') matches the documented formula", {
    df <- mock_res_df()
    score <- compute_ranking_score(df, "pxfc")

    p <- df$pvalue; p[is.na(p)] <- 1
    fc <- df$log2FoldChange; fc[is.na(fc)] <- 0
    expected <- -log10(p) * sign(fc) * abs(fc)

    expect_equal(score, expected, tolerance = 1e-8)
})

test_that("compute_ranking_score('padj') matches the documented formula", {
    df <- mock_res_df()
    score <- compute_ranking_score(df, "padj")

    p <- df$padj; p[is.na(p)] <- 1
    fc <- df$log2FoldChange; fc[is.na(fc)] <- 0
    expected <- -log10(p) * sign(fc)

    expect_equal(score, expected, tolerance = 1e-8)
})

test_that("compute_ranking_score('log2FoldChange') returns the fold-change column with NAs zeroed", {
    df <- mock_res_df()
    score <- compute_ranking_score(df, "log2FoldChange")
    expected <- df$log2FoldChange
    expected[is.na(expected)] <- 0
    expect_equal(score, expected)
})

test_that("compute_ranking_score('stat') returns the Wald statistic column with NAs zeroed", {
    df <- mock_res_df()
    score <- compute_ranking_score(df, "stat")
    expected <- df$stat
    expected[is.na(expected)] <- 0
    expect_equal(score, expected)
})

test_that("compute_ranking_score rejects an unknown ranking method", {
    df <- mock_res_df()
    expect_error(compute_ranking_score(df, "not_a_real_method"))
})

test_that("compute_ranking_score errors when required columns are missing", {
    df <- mock_res_df()[, c("gene", "log2FoldChange")]
    expect_error(compute_ranking_score(df, "pxfc"), "missing required column")
})

test_that("compute_ranking_score handles pvalue == 0 without producing Inf", {
    df <- mock_res_df()
    df$pvalue[1] <- 0
    df$padj[1] <- 0
    score <- compute_ranking_score(df, "pxfc")
    expect_true(all(is.finite(score)))
})

test_that("compute_ranking_score never returns NA even with all-NA input rows", {
    df <- mock_res_df()
    score <- compute_ranking_score(df, "pxfc")
    expect_false(anyNA(score))
})

test_that("build_gsea_ranks returns a sorted, named, deduplicated numeric vector", {
    df <- mock_res_df()
    ranks <- build_gsea_ranks(df, gene_col = "gene", ranking = "pxfc")

    expect_true(is.numeric(ranks))
    expect_false(is.null(names(ranks)))
    expect_true(all(diff(ranks) <= 0))  # sorted decreasing
    expect_false(anyNA(ranks))
})

test_that("build_gsea_ranks drops rows with NA gene identifiers", {
    df <- mock_res_df()
    df$gene[2] <- NA
    ranks <- build_gsea_ranks(df, gene_col = "gene", ranking = "pxfc")
    expect_false("g2" %in% names(ranks))
    expect_equal(length(ranks), nrow(df) - 1)
})

test_that("build_gsea_ranks keeps the highest-magnitude score on duplicate gene names", {
    df <- data.frame(
        gene = c("dup", "dup", "unique"),
        log2FoldChange = c(1, 5, 2),
        pvalue = c(0.01, 0.001, 0.02),
        padj = c(0.02, 0.005, 0.03),
        stat = c(2, 6, 3),
        baseMean = c(10, 10, 10)
    )
    ranks <- build_gsea_ranks(df, gene_col = "gene", ranking = "stat")
    expect_equal(sum(names(ranks) == "dup"), 1)
    expect_equal(unname(ranks["dup"]), 6)
})

test_that("build_gsea_ranks errors if gene_col is not present", {
    df <- mock_res_df()
    expect_error(build_gsea_ranks(df, gene_col = "not_a_column"), "not a column")
})
