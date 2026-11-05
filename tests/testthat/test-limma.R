# limma-trend engine: carrier object, log2 transform, limma results.

simulate_tpm_se <- function(levels = c("Control", "Treated"), n_per_group = 4,
                            seed = 1, ...) {
    dds <- simulate_dds(group_levels = levels, n_per_group = n_per_group,
                        seed = seed, ...)
    x <- DESeq2::counts(dds) * stats::runif(nrow(dds), 0.5, 1.5)
    tpm <- sweep(x, 2, colSums(x), "/") * 1e6
    SummarizedExperiment::SummarizedExperiment(
        assays = list(tpm = tpm), colData = SummarizedExperiment::colData(dds))
}

test_that("expr_to_dds keeps exact values beside an integer stand-in", {
    se <- simulate_tpm_se()
    dds <- expr_to_dds(se, assay = "tpm")
    expect_s4_class(dds, "DESeqDataSet")
    expect_true(.is_expr_dds(dds))
    expect_equal(.raw_count_matrix(dds), SummarizedExperiment::assay(se))
    sub <- dds[1:50, 1:6]
    DESeq2::design(sub) <- ~condition
    expect_equal(.raw_count_matrix(sub), SummarizedExperiment::assay(se)[1:50, 1:6])

    y <- list(counts = SummarizedExperiment::assay(se),
              samples = data.frame(as.data.frame(SummarizedExperiment::colData(se)),
                                   lib.size = 1, norm.factors = 1))
    class(y) <- "DGEList"
    d2 <- expr_to_dds(y)
    expect_false("lib.size" %in% colnames(SummarizedExperiment::colData(d2)))
})

test_that("expr_to_dds rejects bad input", {
    se <- simulate_tpm_se()
    expect_error(expr_to_dds(SummarizedExperiment::assay(se)), "SummarizedExperiment")
    SummarizedExperiment::assay(se)[1, 1] <- -1
    expect_error(expr_to_dds(se), "non-negative")
})

test_that("safe_vst / safe_rlog return log2(x + prior_count) for a carrier", {
    se <- simulate_tpm_se()
    dds <- expr_to_dds(se, prior_count = 0.5)
    expected <- log2(SummarizedExperiment::assay(se) + 0.5)
    expect_equal(SummarizedExperiment::assay(safe_vst(dds)), expected)
    expect_equal(SummarizedExperiment::assay(safe_rlog(dds)), expected)
    m <- sample_level_expr_matrix(dds, var = "condition")
    expect_identical(attr(m, "expr_correction"), "log2expr")
})

test_that("run_limma_diffexp matches limma-trend and keeps the DESeq2 schema", {
    se <- simulate_tpm_se()
    dds <- expr_to_dds(se)
    DESeq2::design(dds) <- ~condition
    res <- run_limma_diffexp(dds, var = "condition",
                             comparison = "Treated_vs_Control", shrinkage = "none")
    expect_named(res, c("gene", "comparison", "baseMean", "mean_cpm", "log2FoldChange",
                        "log2FoldChange_unshrunken", "lfcSE", "stat",
                        "pvalue", "padj", "rank_score", "sig"))

    mm <- stats::model.matrix(~condition, data = as.data.frame(SummarizedExperiment::colData(se)))
    fit <- limma::eBayes(limma::lmFit(log2(SummarizedExperiment::assay(se) + 1), mm),
                         trend = TRUE, robust = TRUE)
    ref <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
    expect_equal(res$log2FoldChange, ref$logFC)
    expect_equal(res$stat, ref$t)
    expect_equal(res$pvalue, ref$P.Value)
    expect_equal(res$padj, ref$adj.P.Val)
    expect_equal(res$log2FoldChange / res$lfcSE, ref$t)
    # The simulated up-shift lands in the first 20 genes.
    expect_gt(mean(res$sig[1:20]), 0.5)

    rev <- run_limma_diffexp(dds, var = "condition",
                             comparison = "Control_vs_Treated", shrinkage = "none")
    expect_equal(rev$log2FoldChange, -res$log2FoldChange)

    shrunk <- run_limma_diffexp(dds, var = "condition", comparison = "Treated_vs_Control")
    expect_equal(shrunk$log2FoldChange_unshrunken, ref$logFC)
    expect_lte(mean(abs(shrunk$log2FoldChange)), mean(abs(ref$logFC)))
})

test_that("limma pipeline handles pooled and interaction contrasts", {
    dds <- expr_to_dds(simulate_tpm_se(levels = c("A", "B", "C"), n_per_group = 3))
    DESeq2::design(dds) <- ~condition
    out <- run_diffexp_pipeline(dds, var = "condition", shrinkage = "none",
                                save_dir = withr::local_tempdir(), save_plot = FALSE)
    expect_true(all(c("A_vs_B", "A_vs_C", "B_vs_C", "A_vs_B+C", "C_vs_A+B") %in%
                        names(out$by_comparison)))
    expect_true(all(vapply(out$by_comparison, nrow, 1L) == nrow(dds)))

    se2 <- simulate_tpm_se(levels = c("A", "B"), n_per_group = 8, seed = 2)
    SummarizedExperiment::colData(se2)$geno <- factor(rep(c("wt", "ko"), 8),
                                                      levels = c("wt", "ko"))
    d2 <- expr_to_dds(se2)
    DESeq2::design(d2) <- ~condition * geno
    ix <- run_diffexp_pipeline(d2, var = c("condition", "geno"),
                               include_interactions = TRUE, shrinkage = "none",
                               save_dir = withr::local_tempdir(), save_plot = FALSE)
    ix_labels <- ix$contrast_specs$label[ix$contrast_specs$fit_method == "coefs"]
    expect_gt(length(ix_labels), 0L)
    expect_true(all(ix_labels %in% names(ix$by_comparison)))
})

test_that(".deseq2_coef_names follows DESeq2 resultsNames()", {
    dds <- simulate_dds(group_levels = c("A", "B", "C"), n_per_group = 4)
    SummarizedExperiment::colData(dds)$geno <- factor(rep(c("wt", "ko"), 6),
                                                      levels = c("wt", "ko"))
    DESeq2::design(dds) <- ~condition * geno
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    mm <- stats::model.matrix(~condition * geno, data = cd)
    fit <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
    expect_identical(.deseq2_coef_names(mm, ~condition * geno, cd),
                     DESeq2::resultsNames(fit))
})
