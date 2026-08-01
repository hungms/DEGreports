test_that("generate_comparisons returns a single sorted comparison for 2 levels", {
    out <- generate_comparisons(c("Treated", "Control"))
    expect_equal(out, "Treated_vs_Control")

    # order of input levels should not matter -- alphabetical winner is
    # always the "test" (first) side
    out2 <- generate_comparisons(c("Control", "Treated"))
    expect_equal(out2, "Treated_vs_Control")
})

test_that("generate_comparisons never generates both A_vs_B and B_vs_A", {
    out <- generate_comparisons(c("A", "B", "C"))
    pairwise <- out[!grepl("+", out, fixed = TRUE)]
    expect_length(pairwise, choose(3, 2))
    # no reversed duplicate present
    reversed <- vapply(strsplit(pairwise, "_vs_", fixed = TRUE), function(p) paste0(p[2], "_vs_", p[1]), character(1))
    expect_true(length(intersect(pairwise, reversed)) == 0)
})

test_that("generate_comparisons produces correct one-vs-rest comparisons for >2 levels", {
    out <- generate_comparisons(c("A", "B", "C"))
    onetoall <- out[grepl("+", out, fixed = TRUE)]
    expect_length(onetoall, 3)
    expect_true("A_vs_B+C" %in% onetoall)
    expect_true("B_vs_A+C" %in% onetoall)
    expect_true("C_vs_A+B" %in% onetoall)
})

test_that("generate_comparisons produces no one-vs-rest comparisons for exactly 2 levels", {
    out <- generate_comparisons(c("A", "B"))
    expect_false(any(grepl("+", out, fixed = TRUE)))
    expect_length(out, 1)
})

test_that("generate_comparisons deduplicates input levels", {
    out <- generate_comparisons(c("A", "B", "A"))
    expect_equal(out, generate_comparisons(c("A", "B")))
})

test_that("generate_comparisons errors on fewer than 2 unique levels", {
    expect_error(generate_comparisons("A"), "at least 2 unique levels")
    expect_error(generate_comparisons(character(0)), "at least 2 unique levels")
})

test_that("generate_comparisons scales correctly for 4 levels", {
    out <- generate_comparisons(c("W", "X", "Y", "Z"))
    pairwise <- out[!grepl("+", out, fixed = TRUE)]
    onetoall <- out[grepl("+", out, fixed = TRUE)]
    expect_length(pairwise, choose(4, 2))
    expect_length(onetoall, 4)
    # every one-vs-rest pools exactly the 3 other levels
    for (comp in onetoall) {
        parts <- strsplit(comp, "_vs_", fixed = TRUE)[[1]]
        pooled <- strsplit(parts[2], "+", fixed = TRUE)[[1]]
        expect_length(pooled, 3)
    }
})

test_that("validate_comparisons(NULL) auto-generates via generate_comparisons", {
    dds <- simulate_dds(group_levels = c("A", "B", "C"), n_per_group = 3, n_genes = 20)
    out <- validate_comparisons(dds, "condition", NULL)
    expect_equal(sort(out), sort(generate_comparisons(c("A", "B", "C"))))
})

test_that("validate_comparisons accepts and passes through valid user comparisons", {
    dds <- simulate_dds(group_levels = c("A", "B", "C"), n_per_group = 3, n_genes = 20)
    out <- validate_comparisons(dds, "condition", c("A_vs_B"))
    expect_equal(out, "A_vs_B")
})

test_that("validate_comparisons rejects a comparison referencing an unknown level", {
    dds <- simulate_dds(group_levels = c("A", "B"), n_per_group = 3, n_genes = 20)
    expect_error(validate_comparisons(dds, "condition", "A_vs_Z"), "not found in the data")
})

test_that("validate_comparison rejects malformed strings", {
    expect_error(validate_comparison("A-vs-B", c("A", "B")), "not in .A_vs_B. format")
    expect_error(validate_comparison("A_vs_B_vs_C", c("A", "B", "C")), "exactly one")
})

test_that("merge_comparison_levels merges levels by exact match, not substring", {
    # deliberately include a level name that is a substring of another
    # ("B" is a substring of "B2") to catch any regex/substring-matching
    # regression in the merge logic
    dds <- simulate_dds(group_levels = c("A", "B", "B2"), n_per_group = 3, n_genes = 10)
    merged <- merge_comparison_levels(dds, "condition", "A_vs_B+B2")

    grp_before <- as.character(SummarizedExperiment::colData(dds)$condition)
    grp_after <- as.character(SummarizedExperiment::colData(merged$dds)[[merged$merged_var]])

    # all samples (both "B" and "B2") should be present and merged into the
    # pooled reference level "B+B2" -- confirming both are matched exactly,
    # neither dropped nor double counted by a substring slip
    expect_true(all(grp_after[grp_before == "B"] == "B+B2"))
    expect_true(all(grp_after[grp_before == "B2"] == "B+B2"))
    expect_true(all(grp_after[grp_before == "A"] == "A"))
    expect_equal(nrow(SummarizedExperiment::colData(merged$dds)), ncol(dds))
})

test_that("classify_comparison_type identifies pairwise vs one-vs-all by the '+' marker", {
    expect_equal(classify_comparison_type("B_vs_A"), "pairwise")
    expect_equal(classify_comparison_type("A_vs_B+C"), "one_vs_all")
    expect_equal(
        classify_comparison_type(c("B_vs_A", "A_vs_B+C", "C_vs_A+B", "A_vs_B")),
        c("pairwise", "one_vs_all", "one_vs_all", "pairwise")
    )
})

test_that("classify_comparison_type agrees with generate_comparisons's own pairwise/one-vs-rest split", {
    comps <- generate_comparisons(c("A", "B", "C"))
    pairwise_expected <- comps[!grepl("+", comps, fixed = TRUE)]
    onetoall_expected <- comps[grepl("+", comps, fixed = TRUE)]
    types <- classify_comparison_type(comps)
    expect_setequal(comps[types == "pairwise"], pairwise_expected)
    expect_setequal(comps[types == "one_vs_all"], onetoall_expected)
})

test_that("merge_comparison_levels drops samples from levels outside the comparison", {
    dds <- simulate_dds(group_levels = c("A", "B", "C"), n_per_group = 3, n_genes = 10)
    merged <- merge_comparison_levels(dds, "condition", "A_vs_B")
    grp_before <- as.character(SummarizedExperiment::colData(dds)$condition)
    expect_equal(ncol(merged$dds), sum(grp_before %in% c("A", "B")))
})
