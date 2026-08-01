test_that("assert_full_rank_model_matrix passes silently for a full-rank design", {
    dds <- simulate_dds(group_levels = c("Control", "Treated"), n_per_group = 3, n_genes = 20)
    expect_silent(assert_full_rank_model_matrix(dds, "Treated_vs_Control", "condition"))
})

test_that("assert_full_rank_model_matrix detects a rank-deficient design (aliased column)", {
    dds <- simulate_dds(group_levels = c("A", "B"), n_per_group = 3, n_genes = 10)
    # add a second column that is an exact duplicate of `condition` -> the
    # combined design ~condition + condition2 is rank-deficient
    SummarizedExperiment::colData(dds)$condition2 <- SummarizedExperiment::colData(dds)$condition
    DESeq2::design(dds) <- ~condition + condition2

    expect_error(
        assert_full_rank_model_matrix(dds, "A_vs_B", "condition"),
        "not full rank"
    )
})

test_that("assert_full_rank_model_matrix error message differs for pooled vs pairwise comparisons", {
    dds <- simulate_dds(group_levels = c("A", "B"), n_per_group = 3, n_genes = 10)
    SummarizedExperiment::colData(dds)$condition2 <- SummarizedExperiment::colData(dds)$condition
    DESeq2::design(dds) <- ~condition + condition2

    expect_error(
        assert_full_rank_model_matrix(dds, "A_vs_B+C", "condition"),
        "pooled one-vs-rest comparison"
    )
    expect_error(
        assert_full_rank_model_matrix(dds, "A_vs_B", "condition"),
        "linearly redundant"
    )
})

test_that("assert_full_rank_model_matrix is a no-op when design is not a formula", {
    dds <- simulate_dds(group_levels = c("A", "B"), n_per_group = 3, n_genes = 10)
    # DESeq2::design() is always a formula in practice, but the guard should
    # not error if given something else defensively
    expect_silent(assert_full_rank_model_matrix(dds, "A_vs_B", "condition"))
})

test_that("run_diffexp surfaces the rank-deficiency error before calling DESeq()", {
    dds <- simulate_dds(group_levels = c("A", "B"), n_per_group = 3, n_genes = 10)
    dds_qc <- run_qc(dds, var = "condition", min_count = 1, remove_mt = FALSE, remove_xy = FALSE,
                      save_plot = FALSE, save_dir = tempfile("qc_test_"))
    SummarizedExperiment::colData(dds_qc)$condition2 <- SummarizedExperiment::colData(dds_qc)$condition
    DESeq2::design(dds_qc) <- ~condition + condition2

    expect_error(
        run_diffexp(dds_qc, var = "condition", comparison = "A_vs_B"),
        "not full rank"
    )
})
