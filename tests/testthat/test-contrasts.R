test_that("resolve_design auto-builds from var and batch", {
    expect_equal(deparse(resolve_design("Group2")), deparse(~Group2))
    expect_equal(deparse(resolve_design("Group2", batch = "Batch")), deparse(~Batch + Group2))
    # Default two-factor: additive (no interactions)
    expect_equal(deparse(resolve_design(c("Tissue", "Isotype"))), deparse(~Tissue + Isotype))
    expect_equal(deparse(resolve_design(c("Tissue", "Isotype"), batch = "Batch")),
                 deparse(~Batch + Tissue + Isotype))
    expect_equal(
        deparse(resolve_design(c("Tissue", "Isotype"), include_interactions = TRUE)),
        deparse(~Tissue * Isotype)
    )
    expect_equal(
        deparse(resolve_design(c("Tissue", "Isotype"), batch = "Batch",
                               include_interactions = TRUE)),
        deparse(~Batch + Tissue * Isotype)
    )
    expect_equal(
        deparse(resolve_design(c("Tissue", "Isotype"), batch = c("Batch", "Sex"),
                               include_interactions = TRUE)),
        deparse(~Batch + Sex + Tissue * Isotype)
    )
    expect_equal(deparse(resolve_design("Group2", design = ~Batch + Group2)),
                 deparse(~Batch + Group2))
})

test_that("expand_contrast_specs single-factor matches generate_comparisons", {
    dds <- simulate_dds(n_genes = 20, group_levels = c("A", "B", "C"), n_per_group = 2)
    specs <- expand_contrast_specs(dds, "condition")
    expect_equal(sort(specs$label), sort(generate_comparisons(levels(dds$condition))))
    expect_true(all(specs$class == "single"))
    expect_true("coef_subtract" %in% names(specs))
})

test_that("expand_contrast_specs 2x3 default omits interactions (15 contrasts)", {
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
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~Tissue + Isotype)

    specs <- expand_contrast_specs(dds, c("Tissue", "Isotype"))
    expect_equal(nrow(specs), 15L)
    expect_equal(sum(specs$class == "within_A" & specs$type == "pairwise"), 6L)
    expect_equal(sum(specs$class == "within_A" & specs$type == "one_vs_all"), 6L)
    expect_equal(sum(specs$class == "within_B" & specs$type == "pairwise"), 3L)
    expect_equal(sum(specs$class == "within_B" & specs$type == "one_vs_all"), 0L)
    expect_equal(sum(specs$class == "interaction"), 0L)
})

test_that("expand_contrast_specs 2x3 with interactions adds ref + pairwise", {
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

    specs <- expand_contrast_specs(dds, c("Tissue", "Isotype"),
                                  include_interactions = TRUE)
    # 15 within + 2 ref interactions + 1 pairwise among Isotype non-ref at Spleen
    expect_equal(nrow(specs), 18L)
    expect_equal(sum(specs$class == "interaction"), 3L)
    expect_true(all(c("BM__IgG_vs_IgM", "Spleen__IgA_vs_IgM+IgG",
                      "IgM__Spleen_vs_BM",
                      "interaction__TissueSpleen_x_IsotypeIgG",
                      "interaction__TissueSpleen_x_IsotypeIgA_vs_IsotypeIgG") %in% specs$label))

    ix <- specs[specs$class == "interaction", , drop = FALSE]
    ref_ix <- ix[is.na(ix$coef_subtract) | !nzchar(ix$coef_subtract), , drop = FALSE]
    pw_ix <- ix[!is.na(ix$coef_subtract) & nzchar(ix$coef_subtract), , drop = FALSE]
    expect_equal(nrow(ref_ix), 2L)
    expect_equal(nrow(pw_ix), 1L)
    expect_equal(pw_ix$comparison[[1]], "TissueSpleen.IsotypeIgA")
    expect_equal(pw_ix$coef_subtract[[1]], "TissueSpleen.IsotypeIgG")
})

test_that("expand_contrast_specs genotype x celltype includes mutant–mutant interactions", {
    set.seed(2)
    n_genes <- 20
    # 4 genotypes x 2 celltypes x 2 reps
    n <- 16
    counts <- matrix(rnbinom(n_genes * n, mu = 40, size = 8), nrow = n_genes)
    rownames(counts) <- paste0("g", seq_len(n_genes))
    colnames(counts) <- paste0("s", seq_len(n))
    coldata <- data.frame(
        genotype = factor(rep(c("WT", "C", "K", "CK"), each = 4),
                          levels = c("WT", "C", "K", "CK")),
        celltype = factor(rep(rep(c("CB", "CC"), each = 2), times = 4),
                          levels = c("CB", "CC")),
        row.names = colnames(counts)
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~genotype * celltype)

    specs <- expand_contrast_specs(dds, c("genotype", "celltype"),
                                  include_interactions = TRUE)
    ix <- specs[specs$class == "interaction", , drop = FALSE]
    # 3 ref (C/K/CK x CC) + combn(3,2)=3 pairwise among genotypes at CC
    expect_equal(nrow(ix), 6L)
    expect_true("interaction__genotypeCK_vs_genotypeC_x_celltypeCC" %in% ix$label)
    pw <- ix[ix$label == "interaction__genotypeCK_vs_genotypeC_x_celltypeCC", , drop = FALSE]
    expect_equal(pw$comparison[[1]], "genotypeCK.celltypeCC")
    expect_equal(pw$coef_subtract[[1]], "genotypeC.celltypeCC")
})

test_that("contrast_class_title uses factor names", {
    expect_equal(contrast_class_title("within_A", c("Tissue", "Isotype")),
                 "Within Tissue \u2014 Isotype effects")
    expect_equal(contrast_class_title("within_B", c("Tissue", "Isotype")),
                 "Within Isotype \u2014 Tissue effects")
    expect_equal(contrast_class_title("interaction", c("Tissue", "Isotype")), "Interaction")
})
