test_that("run_counts_export writes raw/vst/rlog/cpm tables with matching gene sets", {
    dds <- simulate_dds(n_genes = 30, n_per_group = 3)
    out_dir <- tempfile("counts_")
    out <- run_counts_export(dds, save_dir = out_dir)

    expect_setequal(names(out), c("raw", "vst", "rlog", "cpm"))
    expect_setequal(
        list.files(out_dir),
        c("counts_raw.csv", "counts_vst.csv", "counts_rlog.csv", "counts_cpm.csv")
    )
    # all four tables cover the same genes, in the same order
    expect_equal(out$raw$gene, out$vst$gene)
    expect_equal(out$raw$gene, out$rlog$gene)
    expect_equal(out$raw$gene, out$cpm$gene)
    expect_equal(nrow(out$raw), nrow(dds))

    # CPM is computed correctly relative to raw counts and library size
    raw_mat <- as.matrix(out$raw[, -1])
    cpm_mat <- as.matrix(out$cpm[, -1])
    expected_cpm <- sweep(raw_mat, 2, colSums(raw_mat), "/") * 1e6
    expect_equal(unname(cpm_mat), unname(expected_cpm), tolerance = 1e-6)
})

test_that("run_counts_export requires `var` when `batch` is supplied", {
    dds <- simulate_dds(n_genes = 20, n_per_group = 3, batch = TRUE)
    expect_error(
        run_counts_export(dds, batch = "batch", save_dir = tempfile("counts_")),
        "`var` must be supplied"
    )
})

test_that("run_counts_export writes a ComBat_seq-corrected table when `batch` is supplied", {
    dds <- simulate_dds(n_genes = 20, n_per_group = 3, batch = TRUE)
    out_dir <- tempfile("counts_batch_")
    out <- run_counts_export(dds, batch = "batch", var = "condition", save_dir = out_dir)

    expect_true("batch_corrected" %in% names(out))
    expect_true(file.exists(file.path(out_dir, "counts_combat_seq.csv")))
    expect_equal(dim(out$batch_corrected), dim(out$raw))
    expect_equal(out$batch_corrected$gene, out$raw$gene)
})

test_that("run_counts_export writes a ComBat-corrected table for microarray platform", {
    dds <- simulate_dds(n_genes = 20, n_per_group = 3, batch = TRUE)
    out_dir <- tempfile("counts_batch_array_")
    out <- run_counts_export(dds, batch = "batch", var = "condition",
                             platform = "microarray", save_dir = out_dir)

    expect_true("batch_corrected" %in% names(out))
    expect_true(file.exists(file.path(out_dir, "counts_combat_vst.csv")))
    expect_false(file.exists(file.path(out_dir, "counts_combat_seq.csv")))
    expect_equal(out$batch_corrected$gene, out$raw$gene)
})

test_that("run_qc filter_method = \"cpm\" retains fewer genes at a stricter (higher) percentile", {
    dds <- simulate_dds(n_genes = 200, n_per_group = 4)
    out_dir1 <- tempfile("qc_cpm_")
    dds_loose <- run_qc(dds, var = "condition", filter_method = "cpm", cpm_percentile = 0.05,
                         save_plot = FALSE, save_dir = out_dir1)
    out_dir2 <- tempfile("qc_cpm_")
    dds_strict <- run_qc(dds, var = "condition", filter_method = "cpm", cpm_percentile = 0.3,
                          save_plot = FALSE, save_dir = out_dir2)

    expect_gte(nrow(dds_loose), nrow(dds_strict))
    summary_csv <- read.csv(file.path(out_dir1, "qc_filter_summary.csv"))
    expect_true(any(summary_csv$metric == "filter_method" & summary_csv$value == "cpm"))
})

test_that("run_qc filter_method = \"cpm\" rejects an invalid cpm_percentile", {
    dds <- simulate_dds(n_genes = 20, n_per_group = 3)
    expect_error(
        run_qc(dds, var = "condition", filter_method = "cpm", cpm_percentile = -1,
               save_plot = FALSE, save_dir = tempfile("qc_")),
        "single number in"
    )
})

test_that("run_qc filter_method = \"count\" keeps a gene expressed in only one group", {
    # 4 genes, 2 groups of 3 samples each. Gene 1 expressed in all of group A
    # only (0 in group B) -- should be retained because *some* group clears
    # the per-group sample-count threshold.
    counts <- matrix(0L, nrow = 4, ncol = 6)
    rownames(counts) <- paste0("g", 1:4)
    colnames(counts) <- paste0("s", 1:6)
    counts[1, 1:3] <- 20L         # gene1: high in group A (cols 1-3), zero in B
    counts[2, ] <- 20L            # gene2: high everywhere
    counts[3, ] <- 0L             # gene3: never expressed
    counts[4, 1:2] <- 20L         # gene4: only 2/3 samples in group A -- below min_samples_per_group=3

    coldata <- data.frame(condition = factor(rep(c("A", "B"), each = 3)),
                           row.names = colnames(counts))
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~condition)

    dds_qc <- run_qc(dds, var = "condition", filter_method = "count",
                      min_count = 10, min_samples_per_group = 3,
                      save_plot = FALSE, save_dir = tempfile("qc_count_"))

    kept <- rownames(dds_qc)
    expect_true("g1" %in% kept)
    expect_true("g2" %in% kept)
    expect_false("g3" %in% kept)
    expect_false("g4" %in% kept)
})
