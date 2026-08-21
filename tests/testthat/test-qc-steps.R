# Tests for the per-step gene-removal tracking added to run_qc(): the
# steps_summary data.frame, the qc_steps_summary.csv it writes, and the
# metadata attached to `dds` for render_report() to consume.

test_that("run_qc() tracks genes removed at each step and adds up correctly", {
    dds <- simulate_dds(n_genes = 100, group_levels = c("A", "B"), n_per_group = 4)
    # rename a handful of genes to look mitochondrial, and mark a disjoint
    # set of genes as X/Y in gene_chr, so the two removal steps don't
    # overlap and each removes exactly the count it claims to.
    genes <- rownames(dds)
    genes[1:3] <- paste0("MT-", genes[1:3])
    rownames(dds) <- genes
    chr <- rep("1", length(genes))
    chr[4:7] <- "X"    # 4 X genes, disjoint from the 3 MT genes above
    chr[8:9] <- "Y"    # 2 Y genes
    gene_chr <- data.frame(gene = genes, chr = chr)

    save_dir <- withr::local_tempdir()
    dds_qc <- run_qc(dds, var = "condition", filter_method = "count",
                      min_count = 1, min_samples_per_group = 1,
                      remove_mt = TRUE, mt_pattern = "^MT-",
                      remove_xy = TRUE, gene_chr = gene_chr,
                      save_plot = FALSE, save_dir = save_dir)

    qc_meta <- S4Vectors::metadata(dds_qc)$deseq2reports_qc
    expect_false(is.null(qc_meta))
    steps <- qc_meta$steps_summary

    expect_true(all(c("step", "n_before", "n_removed", "n_after", "pct_removed") %in% colnames(steps)))
    expect_identical(steps$step,
                      c("Input", "Mitochondrial gene removal", "X/Y chromosome gene removal",
                        "Low-expression gene removal (raw count, per-group)"))
    # n_removed = n_before - n_after at every step
    expect_equal(steps$n_removed, steps$n_before - steps$n_after)
    # each step's n_before matches the previous step's n_after (chained correctly)
    expect_equal(steps$n_before[-1], steps$n_after[-nrow(steps)])
    # first step starts at the pre-filter gene count, last step ends at nrow(dds_qc)
    expect_equal(steps$n_before[1], 100)
    expect_equal(steps$n_after[nrow(steps)], nrow(dds_qc))
    # exactly 3 MT genes removed, exactly the 6 X/Y genes removed
    expect_equal(steps$n_removed[steps$step == "Mitochondrial gene removal"], 3)
    expect_equal(steps$n_removed[steps$step == "X/Y chromosome gene removal"], 6)

    expect_true(!is.null(qc_meta$sample_qc_table))
    expect_equal(nrow(qc_meta$sample_qc_table), ncol(dds))
    expect_true(all(c("sample", "n_detected_genes", "n_total_transcripts") %in%
                        colnames(qc_meta$sample_qc_table)))
    expect_true(file.exists(file.path(save_dir, "qc_sample_summary.csv")))

    # CSV written and matches the in-memory summary
    csv_path <- file.path(save_dir, "qc_steps_summary.csv")
    expect_true(file.exists(csv_path))
    steps_csv <- utils::read.csv(csv_path)
    expect_equal(steps_csv$n_before, steps$n_before)
    expect_equal(steps_csv$n_after, steps$n_after)
})

test_that("run_qc() steps_summary has one row when only the expression filter runs", {
    dds <- simulate_dds(n_genes = 30, group_levels = c("A", "B"), n_per_group = 3)
    save_dir <- withr::local_tempdir()
    dds_qc <- run_qc(dds, var = "condition", filter_method = "count",
                      min_count = 1, min_samples_per_group = 1,
                      remove_mt = FALSE, remove_xy = FALSE,
                      save_plot = FALSE, save_dir = save_dir)

    steps <- S4Vectors::metadata(dds_qc)$deseq2reports_qc$steps_summary
    expect_identical(steps$step, c("Input", "Low-expression gene removal (raw count, per-group)"))
    expect_equal(nrow(steps), 2)
})

test_that("run_qc(filter_method = NULL) applies no expression filtering", {
    dds <- simulate_dds(n_genes = 25, group_levels = c("A", "B"), n_per_group = 3)
    save_dir <- withr::local_tempdir()
    dds_qc <- run_qc(dds, var = "condition", filter_method = NULL,
                      remove_mt = FALSE, remove_xy = FALSE,
                      save_plot = FALSE, save_dir = save_dir)

    expect_equal(nrow(dds_qc), 25)
    qc_meta <- S4Vectors::metadata(dds_qc)$deseq2reports_qc
    expect_null(qc_meta$filter_method)
    # no expression-filtering step recorded when filter_method = NULL
    expect_identical(qc_meta$steps_summary$step, "Input")
})

test_that("run_qc() metadata records filter_method/cpm_percentile for the CPM percentile path", {
    dds <- simulate_dds(n_genes = 200, group_levels = c("A", "B"), n_per_group = 5)
    save_dir <- withr::local_tempdir()
    dds_qc <- run_qc(dds, var = "condition", filter_method = "cpm",
                      cpm_percentile = 0.1,
                      remove_mt = FALSE, remove_xy = FALSE,
                      save_plot = FALSE, save_dir = save_dir)

    qc_meta <- S4Vectors::metadata(dds_qc)$deseq2reports_qc
    expect_identical(qc_meta$filter_method, "cpm")
    expect_equal(qc_meta$cpm_percentile, 0.1)
    expect_s3_class(qc_meta$count_dist_plot, "ggplot")
    expect_s3_class(qc_meta$library_size_plot, "ggplot")
    expect_true(any(grepl("Low-expression gene removal \\(CPM percentile\\)", qc_meta$steps_summary$step)))
    # no unexpressed genes in this simulation, so the 10% percentile is the
    # whole filter: ~20 of 200 genes
    expect_equal(nrow(dds_qc), 200 - 20, tolerance = 0.05)
})

# A percentile taken over the raw mean-CPM distribution silently removes
# *nothing* on real bulk RNA-seq: a full annotation leaves half the genes zero
# in every sample, the quantile lands inside that tie at zero, and every gene
# satisfies log2(mean CPM + 1) >= 0. GSE224465 lost 0 of 54,860 genes this way.

test_that("run_qc() CPM filter removes unexpressed genes even when they outnumber cpm_percentile", {
    dds <- simulate_dds(n_genes = 100, group_levels = c("A", "B"), n_per_group = 3)
    # 60% of genes zero in every sample -- far more than the 5% percentile.
    counts <- DESeq2::counts(dds)
    counts[1:60, ] <- 0L
    SummarizedExperiment::assay(dds, "counts") <- counts
    save_dir <- withr::local_tempdir()

    dds_qc <- run_qc(dds, var = "condition", filter_method = "cpm",
                      cpm_percentile = 0.05, remove_mt = FALSE, remove_xy = FALSE,
                      save_plot = FALSE, save_dir = save_dir)

    expect_false(any(rownames(dds_qc) %in% rownames(dds)[1:60]))
    # all 60 unexpressed, plus 5% of the 40 that are left
    expect_equal(nrow(dds_qc), 38, tolerance = 0.1)
    steps <- S4Vectors::metadata(dds_qc)$deseq2reports_qc$steps_summary
    expect_gte(steps$n_removed[steps$step == "Low-expression gene removal (CPM percentile)"], 60)
})

test_that("run_qc() CPM filter with cpm_percentile = 0 removes only the unexpressed genes", {
    dds <- simulate_dds(n_genes = 50, group_levels = c("A", "B"), n_per_group = 3)
    counts <- DESeq2::counts(dds)
    counts[1:10, ] <- 0L
    SummarizedExperiment::assay(dds, "counts") <- counts
    save_dir <- withr::local_tempdir()

    dds_qc <- run_qc(dds, var = "condition", filter_method = "cpm",
                      cpm_percentile = 0, remove_mt = FALSE, remove_xy = FALSE,
                      save_plot = FALSE, save_dir = save_dir)

    expect_equal(nrow(dds_qc), 40)
})

test_that("run_qc() metadata is NULL-safe: render_report() tolerates a dds without it", {
    # Simulate a dds that was never passed through run_qc() (metadata absent).
    dds <- simulate_dds(n_genes = 20, group_levels = c("A", "B"), n_per_group = 3)
    expect_null(S4Vectors::metadata(dds)$deseq2reports_qc)
})
