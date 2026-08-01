# Shared helper: build a small simulated DESeqDataSet for tests.
#
# `n_genes` genes x samples across the levels of `group_levels`
# (`n_per_group` samples each), with the first `n_de` genes given a
# deliberate mean-count shift in the last group so that differential
# expression tests have true positives to find. Counts are simulated from
# a negative binomial so DESeq2's dispersion estimation has something
# realistic to work with.
simulate_dds <- function(
    n_genes = 200,
    group_levels = c("Control", "Treated"),
    n_per_group = 4,
    n_de = 20,
    batch = FALSE,
    seed = 1) {

    set.seed(seed)
    n_samples <- length(group_levels) * n_per_group
    condition <- factor(rep(group_levels, each = n_per_group), levels = group_levels)

    base_mean <- 2^stats::rnorm(n_genes, mean = 6, sd = 2)
    counts <- matrix(NA_integer_, nrow = n_genes, ncol = n_samples)
    for (j in seq_len(n_samples)) {
        lambda <- base_mean
        if (j > (n_samples - n_per_group) && n_de > 0) {
            # last group gets an upregulation shift on the first n_de genes
            lambda[seq_len(n_de)] <- lambda[seq_len(n_de)] * 4
        }
        counts[, j] <- stats::rnbinom(n_genes, mu = pmax(lambda, 1), size = 5)
    }
    rownames(counts) <- paste0("gene", seq_len(n_genes))
    colnames(counts) <- paste0("sample", seq_len(n_samples))

    coldata <- data.frame(
        row.names = colnames(counts),
        condition = condition
    )
    if (batch) {
        coldata$batch <- factor(rep(c("b1", "b2"), length.out = n_samples))
    }

    DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        design = ~condition
    )
}

# Like simulate_dds(), but with real HGNC Hallmark-pathway gene symbols as
# rownames instead of generic "gene1", "gene2", ... names. Required by any
# test that exercises organism auto-detection (detect_organism()) or
# MSigDB-backed analysis (run_gsea_pipeline(), run_wgcna_pipeline()'s ORA
# step) -- those need real gene symbols to have anything to detect/match
# against. Requires the (Suggests-only) msigdbr package.
simulate_dds_real_genes <- function(
    n_genes = 300,
    group_levels = c("Control", "Treated"),
    n_per_group = 5,
    n_de = 30,
    batch = FALSE,
    seed = 1) {

    dds <- simulate_dds(
        n_genes = n_genes, group_levels = group_levels, n_per_group = n_per_group,
        n_de = n_de, batch = batch, seed = seed
    )
    set.seed(seed)
    hallmark_genes <- unique(msigdbr::msigdbr(species = "Homo sapiens", collection = "H")$gene_symbol)
    rownames(dds) <- sample(hallmark_genes, n_genes)
    dds
}

# Minimal results-shaped data frame for ranking-metric tests.
mock_res_df <- function() {
    data.frame(
        gene = paste0("g", 1:6),
        log2FoldChange = c(3, -2, 0.05, NA, 1.5, -0.2),
        pvalue = c(0.0001, 0.001, 0.9, 0.5, NA, 0.04),
        padj = c(0.001, 0.01, 0.95, 0.6, NA, 0.08),
        stat = c(6, -4, 0.1, NA, 2.5, -1.8),
        baseMean = c(100, 200, 50, 10, 30, 80),
        mean_cpm = c(120, 210, 55, 12, 33, 88)
    )
}
