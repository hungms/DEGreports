# Paired designs: `pair` validation, design formulas, and removing the
# pairing effect from the sample-level matrix.

# 4 donors, each measured under Control and Treated (8 samples), with a
# strong donor effect on every gene so pairing has something to remove.
simulate_paired_dds <- function(n_genes = 200, n_donors = 4, donor_shift = 3, seed = 1) {
    dds <- simulate_dds(n_genes = n_genes, n_per_group = n_donors, seed = seed)
    donor <- factor(rep(paste0("d", seq_len(n_donors)), times = 2))
    cnt <- DESeq2::counts(dds)
    set.seed(seed)
    mult <- 2^(stats::rnorm(n_donors, sd = donor_shift / 2))
    cnt <- round(sweep(cnt, 2, mult[as.integer(donor)], "*"))
    storage.mode(cnt) <- "integer"
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    cd$donor <- donor
    DESeq2::DESeqDataSetFromMatrix(cnt, cd, design = ~condition)
}

test_that("validate_pair passes NULL through and factorises the column", {
    dds <- simulate_paired_dds()
    expect_identical(validate_pair(dds, NULL, var = "condition"), dds)
    SummarizedExperiment::colData(dds)$donor <- as.character(SummarizedExperiment::colData(dds)$donor)
    out <- validate_pair(dds, "donor", var = "condition")
    expect_true(is.factor(SummarizedExperiment::colData(out)$donor))
})

test_that("validate_pair rejects missing, overlapping and single-level columns", {
    dds <- simulate_paired_dds()
    expect_error(validate_pair(dds, "nope", var = "condition"), "not a column")
    expect_error(validate_pair(dds, "condition", var = "condition"), "must not also be given")
    expect_error(validate_pair(dds, c("donor", "condition"), var = "condition"), "single column")
    SummarizedExperiment::colData(dds)$one <- factor("x")
    expect_error(validate_pair(dds, "one", var = "condition"), "at least 2 levels")
})

test_that("validate_pair errors when the tested factor is constant within every unit", {
    dds <- simulate_paired_dds()
    # each 'subject' sees only one condition: a between-subject factor
    SummarizedExperiment::colData(dds)$subject <- factor(paste0("s", seq_len(ncol(dds))))
    expect_error(validate_pair(dds, "subject", var = "condition"), "nested within `condition`")
})

test_that("validate_pair warns about units seen in only one condition", {
    dds <- simulate_paired_dds()
    d <- as.character(SummarizedExperiment::colData(dds)$donor)
    d[8] <- "d9"   # last Treated sample from a new, unpaired donor
    SummarizedExperiment::colData(dds)$donor <- factor(d)
    expect_warning(validate_pair(dds, "donor", var = "condition"), "only one level")
})

test_that("validate_pair rejects a batch that is constant within every unit", {
    dds <- simulate_paired_dds()
    d <- SummarizedExperiment::colData(dds)$donor
    SummarizedExperiment::colData(dds)$site <- factor(ifelse(d %in% c("d1", "d2"), "A", "B"))
    expect_error(validate_pair(dds, "donor", var = "condition", batch = "site"),
                 "already accounts for it")
})

test_that("pair enters the design as a blocking term next to batch", {
    expect_equal(deparse(resolve_design("condition", batch = c("batch", "donor"))),
                 "~batch + donor + condition")
    expect_equal(deparse(resolve_design(c("A", "B"), batch = "donor",
                                        include_interactions = TRUE)),
                 "~donor + A * B")
})

test_that("removing the pair effect makes samples cluster by condition", {
    dds <- simulate_paired_dds(donor_shift = 4)
    raw <- sample_level_expr_matrix(dds, var = "condition", batch = NULL)
    cor <- sample_level_expr_matrix(dds, var = "condition", batch = "donor")
    expect_equal(attr(cor, "expr_correction"), "rbe_vst")
    donor <- SummarizedExperiment::colData(dds)$donor
    # Donor explains much less of the per-sample mean after correction.
    r2 <- function(m) summary(stats::lm(colMeans(m) ~ donor))$r.squared
    expect_lt(r2(cor), r2(raw))
    expect_lt(r2(cor), 0.2)
})

test_that("run_diffexp_pipeline fits the paired design for every contrast", {
    dds <- simulate_paired_dds()
    out <- run_diffexp_pipeline(dds, var = "condition", batch = "donor",
                                save_plot = FALSE, save_dir = withr::local_tempdir())
    expect_match(out$manifest$design, "donor", fixed = TRUE)
    res <- out$by_comparison[[1]]
    expect_true(sum(res$padj < 0.05, na.rm = TRUE) > 0)
})

test_that("run_deseq2 accepts pair and validates it before fitting", {
    dds <- simulate_paired_dds()
    dir <- withr::local_tempdir()
    SummarizedExperiment::colData(dds)$subject <- factor(paste0("s", seq_len(ncol(dds))))
    expect_error(run_deseq2(dds, var = "condition", pair = "subject",
                            render_html = FALSE, save_dir = dir),
                 "nested within `condition`")
})

test_that("platform is no longer an argument", {
    expect_false("platform" %in% names(formals(run_deseq2)))
    expect_false("platform" %in% names(formals(run_limma)))
    expect_false("platform" %in% names(formals(run_distance_pca)))
    expect_true("pair" %in% names(formals(run_deseq2)))
    expect_true("pair" %in% names(formals(run_limma)))
})
