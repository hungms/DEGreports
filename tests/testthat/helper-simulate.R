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
# against. Requires the (Imports) msigdbr package.
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

# Balanced 2 x 2 x 2 design (celltype x C_mutation x K_mutation), the shape
# the three-factor contrast grid is built for. `n_per_cell` samples per cell,
# with the first `n_de` genes shifted in the C_mutation = "C" samples.
simulate_three_factor_dds <- function(n_genes = 60, n_per_cell = 3, n_de = 10,
                                      seed = 1) {
    set.seed(seed)
    grid <- expand.grid(
        rep = seq_len(n_per_cell),
        K_mutation = c("noK", "K"),
        C_mutation = c("noC", "C"),
        celltype = c("CB", "CC"),
        stringsAsFactors = FALSE
    )
    n_samples <- nrow(grid)
    base_mean <- 2^stats::rnorm(n_genes, mean = 6, sd = 1.5)
    counts <- vapply(seq_len(n_samples), function(j) {
        lambda <- base_mean
        if (n_de > 0 && identical(grid$C_mutation[j], "C")) {
            lambda[seq_len(n_de)] <- lambda[seq_len(n_de)] * 4
        }
        stats::rnbinom(n_genes, mu = pmax(lambda, 1), size = 5)
    }, numeric(n_genes))
    rownames(counts) <- paste0("gene", seq_len(n_genes))
    colnames(counts) <- paste0("sample", seq_len(n_samples))

    coldata <- data.frame(
        row.names = colnames(counts),
        celltype = factor(grid$celltype, levels = c("CB", "CC")),
        C_mutation = factor(grid$C_mutation, levels = c("noC", "C")),
        K_mutation = factor(grid$K_mutation, levels = c("noK", "K"))
    )
    DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        design = ~ celltype * C_mutation * K_mutation
    )
}

# Balanced 2 x 2 x 2 design carrying a known three-way structure, for the
# difference-of-differences contrasts. Genes fall into three blocks:
#
#   1..n_ix                     K effect present only in C_mutation = "C"
#                               samples of celltype = "CB" -- i.e. a C x K
#                               interaction inside CB and none inside CC
#   n_ix+1 .. n_ix+n_main       plain additive K effect in every cell, so
#                               every interaction contrast should be null
#   the rest                    no effect at all
#
# The block boundaries are returned in the object's metadata (`ix_genes`,
# `main_genes`) so tests need not re-derive them.
simulate_interaction_dds <- function(n_genes = 60, n_per_cell = 6, n_ix = 10,
                                    n_main = 10, effect = 8, seed = 42) {
    set.seed(seed)
    grid <- expand.grid(
        rep = seq_len(n_per_cell),
        K_mutation = c("noK", "K"),
        C_mutation = c("noC", "C"),
        celltype = c("CB", "CC"),
        stringsAsFactors = FALSE
    )
    ix <- seq_len(n_ix)
    main <- if (n_main > 0) n_ix + seq_len(n_main) else integer(0)
    base_mean <- 2^stats::rnorm(n_genes, mean = 7, sd = 1)
    counts <- vapply(seq_len(nrow(grid)), function(j) {
        lambda <- base_mean
        in_cb_c_k <- identical(grid$celltype[j], "CB") &&
            identical(grid$C_mutation[j], "C") &&
            identical(grid$K_mutation[j], "K")
        if (length(ix) && in_cb_c_k) lambda[ix] <- lambda[ix] * effect
        if (length(main) && identical(grid$K_mutation[j], "K")) {
            lambda[main] <- lambda[main] * effect
        }
        # size = 50 keeps dispersion low enough for a 48-sample fixture to
        # resolve the implanted interaction.
        stats::rnbinom(n_genes, mu = pmax(lambda, 1), size = 50)
    }, numeric(n_genes))
    rownames(counts) <- paste0("gene", seq_len(n_genes))
    colnames(counts) <- paste0("sample", seq_len(nrow(grid)))

    coldata <- data.frame(
        row.names = colnames(counts),
        celltype = factor(grid$celltype, levels = c("CB", "CC")),
        C_mutation = factor(grid$C_mutation, levels = c("noC", "C")),
        K_mutation = factor(grid$K_mutation, levels = c("noK", "K"))
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        design = ~ celltype * C_mutation * K_mutation
    )
    S4Vectors::metadata(dds)$ix_genes <- rownames(counts)[ix]
    S4Vectors::metadata(dds)$main_genes <- rownames(counts)[main]
    dds
}

# Same idea over a 3 x 2 x 2 design, for the one-vs-rest interactions: the
# K-inside-C effect is confined to celltype "CA", so the pooled comparison
# `CA_vs_CB+CC` should recover it at full size while `CB_vs_CA+CC` sees the
# mirror image, halved (0 - (effect + 0) / 2).
simulate_pooled_interaction_dds <- function(n_genes = 60, n_per_cell = 4,
                                           n_ix = 10, n_main = 10, effect = 8,
                                           seed = 11) {
    set.seed(seed)
    grid <- expand.grid(
        rep = seq_len(n_per_cell),
        K_mutation = c("noK", "K"),
        C_mutation = c("noC", "C"),
        celltype = c("CA", "CB", "CC"),
        stringsAsFactors = FALSE
    )
    ix <- seq_len(n_ix)
    main <- if (n_main > 0) n_ix + seq_len(n_main) else integer(0)
    # A pooled contrast averages four coefficients, so it needs deeper counts
    # and tighter dispersion than the 2 x 2 x 2 fixture for the null genes to
    # stay flat: a low-count gene's noise alone would clear the LFC threshold.
    base_mean <- 2^stats::rnorm(n_genes, mean = 9, sd = 0.5)
    counts <- vapply(seq_len(nrow(grid)), function(j) {
        lambda <- base_mean
        in_ca_c_k <- identical(grid$celltype[j], "CA") &&
            identical(grid$C_mutation[j], "C") &&
            identical(grid$K_mutation[j], "K")
        if (length(ix) && in_ca_c_k) lambda[ix] <- lambda[ix] * effect
        if (length(main) && identical(grid$K_mutation[j], "K")) {
            lambda[main] <- lambda[main] * effect
        }
        stats::rnbinom(n_genes, mu = pmax(lambda, 1), size = 500)
    }, numeric(n_genes))
    rownames(counts) <- paste0("gene", seq_len(n_genes))
    colnames(counts) <- paste0("sample", seq_len(nrow(grid)))

    coldata <- data.frame(
        row.names = colnames(counts),
        celltype = factor(grid$celltype, levels = c("CA", "CB", "CC")),
        C_mutation = factor(grid$C_mutation, levels = c("noC", "C")),
        K_mutation = factor(grid$K_mutation, levels = c("noK", "K"))
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        design = ~ celltype * C_mutation * K_mutation
    )
    S4Vectors::metadata(dds)$ix_genes <- rownames(counts)[ix]
    S4Vectors::metadata(dds)$main_genes <- rownames(counts)[main]
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
