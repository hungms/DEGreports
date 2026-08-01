test_that("resolve_design auto-builds from var and batch", {
    expect_equal(deparse(resolve_design("Group2")), deparse(~Group2))
    expect_equal(deparse(resolve_design("Group2", batch = "Batch")), deparse(~Batch + Group2))
    expect_equal(deparse(resolve_design(c("Tissue", "Isotype"))), deparse(~Tissue * Isotype))
    expect_equal(deparse(resolve_design(c("Tissue", "Isotype"), batch = "Batch")),
                 deparse(~Batch + Tissue * Isotype))
    expect_equal(deparse(resolve_design(c("Tissue", "Isotype"), batch = c("Batch", "Sex"))),
                 deparse(~Batch + Sex + Tissue * Isotype))
    expect_equal(deparse(resolve_design("Group2", design = ~Batch + Group2)),
                 deparse(~Batch + Group2))
})

test_that("expand_contrast_specs single-factor matches generate_comparisons", {
    dds <- simulate_dds(n_genes = 20, group_levels = c("A", "B", "C"), n_per_group = 2)
    specs <- expand_contrast_specs(dds, "condition")
    expect_equal(sort(specs$label), sort(generate_comparisons(levels(dds$condition))))
    expect_true(all(specs$class == "single"))
})

test_that("expand_contrast_specs 2x3 yields 17 contrasts", {
    set.seed(1)
    n_genes <- 30
    counts <- matrix(rnbinom(n_genes * 12, mu = 50, size = 10), nrow = n_genes)
    rownames(counts) <- paste0("g", seq_len(n_genes))
    colnames(counts) <- paste0("s", seq_len(12))
    coldata <- data.frame(
        Tissue = factor(rep(c("BM", "Spleen"), each = 6), levels = c("BM", "Spleen")),
        Isotype = factor(rep(c("IgM", "IgG", "IgA"), times = 4), levels = c("IgM", "IgG", "IgA")),
        row.names = colnames(counts)
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~Tissue * Isotype)

    specs <- expand_contrast_specs(dds, c("Tissue", "Isotype"))
    expect_equal(nrow(specs), 17L)
    expect_equal(sum(specs$class == "within_A" & specs$type == "pairwise"), 6L)
    expect_equal(sum(specs$class == "within_A" & specs$type == "one_vs_all"), 6L)
    expect_equal(sum(specs$class == "within_B" & specs$type == "pairwise"), 3L)
    expect_equal(sum(specs$class == "within_B" & specs$type == "one_vs_all"), 0L)
    expect_equal(sum(specs$class == "interaction"), 2L)
    expect_true(all(c("BM__IgG_vs_IgM", "Spleen__IgA_vs_IgM+IgG",
                      "IgM__Spleen_vs_BM",
                      "interaction__TissueSpleen_x_IsotypeIgG") %in% specs$label))
})

test_that("contrast_class_title uses factor names", {
    expect_equal(contrast_class_title("within_A", c("Tissue", "Isotype")),
                 "Within Tissue \u2014 Isotype effects")
    expect_equal(contrast_class_title("within_B", c("Tissue", "Isotype")),
                 "Within Isotype \u2014 Tissue effects")
    expect_equal(contrast_class_title("interaction", c("Tissue", "Isotype")), "Interaction")
})
