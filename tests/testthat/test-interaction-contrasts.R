var_3f <- c("celltype", "C_mutation", "K_mutation")

# One product-design fit reused by the tests below, since DESeq() on the
# 48-sample fixture is the slow part.
fit_interaction_fixture <- function() {
    dds <- simulate_interaction_dds()
    DESeq2::design(dds) <- resolve_design(var_3f, include_interactions = TRUE)
    list(dds = dds, fit = DESeq2::DESeq(dds, quiet = TRUE))
}

interaction_result <- function(fixture, label, shrinkage = "none") {
    specs <- expand_contrast_specs(fixture$dds, var_3f, include_interactions = TRUE)
    row <- specs[specs$label == label, , drop = FALSE]
    expect_equal(nrow(row), 1L)
    cw <- .parse_coef_weights(row$coef_weights[[1]])
    run_diffexp_coef(
        fixture$dds,
        coef = if (is.null(cw)) row$comparison[[1]] else NULL,
        coef_subtract = if (is.null(cw)) row$coef_subtract[[1]] else NULL,
        coef_weights = cw,
        shrinkage = shrinkage, comparison = label, prefit = fixture$fit
    )
}

test_that("a summed list contrast equals the sum of the model's coefficients", {
    fixture <- fit_interaction_fixture()
    res <- interaction_result(fixture, "C_vs_noC_x_K_vs_noK_in_CC")

    betas <- stats::coef(fixture$fit)
    expected <- betas[, "C_mutationC.K_mutationK"] +
        betas[, "celltypeCC.C_mutationC.K_mutationK"]
    expect_equal(res$log2FoldChange, unname(expected[res$gene]), tolerance = 1e-6)

    # A single coefficient with nothing subtracted still goes through
    # results(name = ), which must agree with the same coefficient's beta.
    single <- interaction_result(fixture, "C_vs_noC_x_K_vs_noK_in_CB")
    expect_equal(single$log2FoldChange,
                 unname(betas[single$gene, "C_mutationC.K_mutationK"]),
                 tolerance = 1e-6)
})

test_that("run_diffexp_coef rejects a coefficient the fit does not have", {
    fixture <- fit_interaction_fixture()
    expect_error(
        run_diffexp_coef(fixture$dds, coef = "celltypeCX.K_mutationK",
                         prefit = fixture$fit, comparison = "bogus"),
        "not in resultsNames"
    )
})

test_that("an interaction contrast recovers a K effect confined to one stratum", {
    fixture <- fit_interaction_fixture()
    ix_genes <- S4Vectors::metadata(fixture$dds)$ix_genes
    main_genes <- S4Vectors::metadata(fixture$dds)$main_genes

    in_cb <- interaction_result(fixture, "C_vs_noC_x_K_vs_noK_in_CB", shrinkage = "ashr")
    in_cc <- interaction_result(fixture, "C_vs_noC_x_K_vs_noK_in_CC", shrinkage = "ashr")
    lfc <- function(res, genes) res$log2FoldChange[match(genes, res$gene)]
    sig <- function(res, genes) {
        padj <- res$padj[match(genes, res$gene)]
        !is.na(padj) & padj < 0.05
    }

    # Implanted: K matters only in C, and only within CB (log2(8) = 3).
    expect_gt(stats::median(lfc(in_cb, ix_genes)), 2)
    expect_true(mean(sig(in_cb, ix_genes)) > 0.8)
    # Same question inside CC, where nothing was implanted.
    expect_lt(max(abs(lfc(in_cc, ix_genes))), 1)
    expect_equal(sum(sig(in_cc, ix_genes)), 0L)
    # A plain additive K effect is not an interaction anywhere.
    expect_lt(max(abs(lfc(in_cb, main_genes))), 1)
    expect_equal(sum(sig(in_cb, main_genes)), 0L)
    expect_equal(sum(sig(in_cc, main_genes)), 0L)
})

test_that("a weighted contrast averages the pooled levels' coefficients", {
    dds <- simulate_pooled_interaction_dds()
    DESeq2::design(dds) <- resolve_design(var_3f, include_interactions = TRUE)
    fit <- DESeq2::DESeq(dds, quiet = TRUE)
    fixture <- list(dds = dds, fit = fit)

    # CA is the reference celltype, so its own corners are structurally zero
    # and the whole contrast is minus the mean of the other two.
    res <- interaction_result(fixture, "CA_vs_CB+CC_x_K_vs_noK_in_C")
    betas <- stats::coef(fit)
    expected <- -0.5 * (betas[, "celltypeCB.K_mutationK"] +
                            betas[, "celltypeCB.C_mutationC.K_mutationK"] +
                            betas[, "celltypeCC.K_mutationK"] +
                            betas[, "celltypeCC.C_mutationC.K_mutationK"])
    expect_equal(res$log2FoldChange, unname(expected[res$gene]), tolerance = 1e-6)
    # Each pooled level carries -1/2: an equal weight per level, not per
    # sample, so the contrast is a mean of coefficients and not a fit on
    # pooled counts.
    specs <- expand_contrast_specs(dds, var_3f, include_interactions = TRUE)
    w <- .parse_coef_weights(
        specs$coef_weights[specs$label == "CA_vs_CB+CC_x_K_vs_noK_in_C"][[1]]
    )
    expect_length(w, 4L)
    expect_true(all(w == -0.5))
})

test_that("a pooled interaction recovers an effect confined to one level", {
    dds <- simulate_pooled_interaction_dds()
    DESeq2::design(dds) <- resolve_design(var_3f, include_interactions = TRUE)
    fixture <- list(dds = dds, fit = DESeq2::DESeq(dds, quiet = TRUE))
    ix_genes <- S4Vectors::metadata(dds)$ix_genes
    main_genes <- S4Vectors::metadata(dds)$main_genes

    lfc <- function(res, genes) res$log2FoldChange[match(genes, res$gene)]
    # The pipeline's own flag: padj < alpha *and* past the LFC threshold.
    # Averaging four coefficients over a low-dispersion fixture makes tiny
    # differences formally significant, which is what the threshold is for.
    sig <- function(res, genes) res$sig[match(genes, res$gene)]
    in_ca <- interaction_result(fixture, "CA_vs_CB+CC_x_K_vs_noK_in_C",
                               shrinkage = "ashr")
    in_cb <- interaction_result(fixture, "CB_vs_CA+CC_x_K_vs_noK_in_C",
                               shrinkage = "ashr")

    # Implanted at log2(8) = 3 inside CA only. CA against the mean of the
    # other two recovers it whole; CB against the mean of CA and CC sees
    # 0 - (3 + 0) / 2, i.e. half of it with the sign flipped.
    expect_gt(stats::median(lfc(in_ca, ix_genes)), 2.5)
    expect_true(mean(sig(in_ca, ix_genes)) > 0.8)
    expect_lt(stats::median(lfc(in_cb, ix_genes)), -1.2)
    expect_gt(stats::median(lfc(in_cb, ix_genes)), -1.8)
    # A plain additive K effect is not an interaction in any celltype.
    expect_lt(max(abs(lfc(in_ca, main_genes))), 1)
    expect_equal(sum(sig(in_ca, main_genes)), 0L)
    expect_equal(sum(sig(in_cb, main_genes)), 0L)
})

test_that("run_diffexp_pipeline writes a pooled interaction to its own folder", {
    dds <- simulate_pooled_interaction_dds()
    save_dir <- withr::local_tempdir()
    label <- "CA_vs_CB+CC_x_K_vs_noK_in_C"
    out <- run_diffexp_pipeline(
        dds, var = var_3f, comparisons = label, include_interactions = TRUE,
        save_dir = save_dir, save_plot = FALSE
    )
    expect_equal(out$contrast_specs$label, label)
    # The pooled label keeps its "+" on disk, as one-vs-rest main effects do.
    expect_true(file.exists(file.path(save_dir, label,
                                      paste0(label, "_diffexp.csv"))))
    res <- out$by_comparison[[label]]
    expect_true(all(res$comparison == label))
    expect_gt(stats::median(res$log2FoldChange[
        res$gene %in% S4Vectors::metadata(dds)$ix_genes
    ]), 2)
})

test_that("run_diffexp_coef refuses a half-specified weighted contrast", {
    fixture <- fit_interaction_fixture()
    expect_error(
        run_diffexp_coef(fixture$dds, coef = "C_mutationC.K_mutationK",
                         coef_weights = c(C_mutationC.K_mutationK = 1),
                         prefit = fixture$fit, comparison = "both"),
        "not both", fixed = TRUE
    )
    expect_error(
        run_diffexp_coef(fixture$dds,
                         coef_weights = c(celltypeCX.K_mutationK = -0.5),
                         prefit = fixture$fit, comparison = "bogus"),
        "not in resultsNames"
    )
    expect_error(
        run_diffexp_coef(fixture$dds, coef_weights = c(0.5, -0.5),
                         prefit = fixture$fit, comparison = "unnamed"),
        "fully named", fixed = TRUE
    )
})

test_that("interaction labels are not read as group names on the plots", {
    parse <- .parse_contrast_levels
    expect_equal(parse("C_vs_noC_x_K_vs_noK_in_CB"),
                 list(test = "positive LFC", ref = "negative LFC"))
    expect_equal(parse("CC_vs_CB_x_C_vs_noC_x_K_vs_noK"),
                 list(test = "positive LFC", ref = "negative LFC"))
    # Plain and stratum-prefixed labels keep their group names.
    expect_equal(parse("C_vs_noC"), list(test = "C", ref = "noC"))
    expect_equal(parse("CB__IgG_vs_IgM"), list(test = "IgG", ref = "IgM"))
})
