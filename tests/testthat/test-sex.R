# Sex imputation in QC, and imputed sex as a model covariate.

# Mouse-style symbols plus real sex markers. `sex` sets which samples are
# female (high Xist, no Y-gene reads); `condition` is independent of it.
simulate_sexed_dds <- function(sex = rep(c("female", "male"), 4),
                               condition = rep(c("Control", "Treated"), each = 4),
                               n_genes = 200, seed = 1) {
    set.seed(seed)
    markers_f <- c("Xist", "Kdm6a", "Kdm5c", "Eif2s3x", "Ddx3x")
    markers_m <- c("Ddx3y", "Kdm5d", "Uty", "Eif2s3y")
    genes <- c(paste0("Gene", seq_len(n_genes)), markers_f, markers_m)
    n <- length(sex)
    mu <- matrix(2^stats::rnorm(length(genes), 6, 1.5), length(genes), n)
    fem <- sex == "female"
    mu[genes == "Xist", fem] <- 4000
    mu[genes == "Xist", !fem] <- 2
    mu[genes %in% markers_m, fem] <- 0.5
    mu[genes %in% markers_m, !fem] <- 300
    cnt <- matrix(stats::rnbinom(length(mu), mu = mu, size = 10), nrow = length(genes),
                  dimnames = list(genes, paste0("s", seq_len(n))))
    cd <- data.frame(condition = factor(condition), row.names = colnames(cnt))
    DESeq2::DESeqDataSetFromMatrix(cnt, cd, design = ~condition)
}

test_that("impute_sample_sex recovers a two-sex cohort", {
    dds <- simulate_sexed_dds()
    out <- impute_sample_sex(dds)
    expect_equal(out$status, "imputed")
    expect_equal(out$org, "mouse")
    expect_equal(out$table$imputed_sex, rep(c("female", "male"), 4))
    expect_true(all(out$table$female_score[c(TRUE, FALSE)] > out$table$female_score[c(FALSE, TRUE)]))
})

test_that("impute_sample_sex leaves a single-sex cohort undetermined", {
    dds <- simulate_sexed_dds(sex = rep("female", 8))
    out <- impute_sample_sex(dds)
    expect_equal(out$status, "undetermined")
    expect_true(all(is.na(out$table$imputed_sex)))
    expect_match(out$message, "single-sex")
})

test_that("impute_sample_sex fails softly on non-symbol gene IDs", {
    dds <- simulate_dds()   # rownames gene1, gene2, ...
    out <- expect_silent(impute_sample_sex(dds))
    expect_equal(out$status, "failed")
    expect_true(all(is.na(out$table$imputed_sex)))
})

test_that("run_qc always imputes sex, before X/Y genes are removed", {
    dds <- simulate_sexed_dds()
    gene_chr <- data.frame(gene = c("Xist", "Kdm6a", "Kdm5c", "Eif2s3x", "Ddx3x",
                                    "Ddx3y", "Kdm5d", "Uty", "Eif2s3y"),
                           chr = c(rep("X", 5), rep("Y", 4)))
    qc <- run_qc(dds, var = "condition", remove_xy = TRUE, gene_chr = gene_chr,
                 save_plot = FALSE, save_dir = withr::local_tempdir())
    expect_false(any(gene_chr$gene %in% rownames(qc)))
    expect_equal(as.character(SummarizedExperiment::colData(qc)$imputed_sex),
                 rep(c("female", "male"), 4))
    meta <- S4Vectors::metadata(qc)$degreports_qc
    expect_equal(meta$sex$status, "imputed")
    steps <- meta$steps_summary$step
    expect_lt(which(grepl("Low-expression", steps)), which(grepl("X/Y", steps)))
})

test_that("run_qc does not fail when sex cannot be imputed", {
    qc <- run_qc(simulate_dds(), var = "condition", save_plot = FALSE,
                 save_dir = withr::local_tempdir())
    expect_true(all(is.na(SummarizedExperiment::colData(qc)$imputed_sex)))
})

test_that("sex joins the model only when it is not confounded", {
    dds <- simulate_sexed_dds()
    SummarizedExperiment::colData(dds)$imputed_sex <- factor(rep(c("female", "male"), 4))
    used <- .model_covariates(dds, var = "condition")
    expect_true(used$sex_used)
    expect_equal(used$covars, "imputed_sex")

    # sex identical to condition: confounded
    SummarizedExperiment::colData(dds)$imputed_sex <-
        factor(ifelse(SummarizedExperiment::colData(dds)$condition == "Control", "female", "male"))
    conf <- .model_covariates(dds, var = "condition")
    expect_false(conf$sex_used)
    expect_null(conf$covars)
    expect_match(conf$sex_reason, "confounded")

    # each donor has one sex: absorbed by pair
    SummarizedExperiment::colData(dds)$donor <- factor(rep(c("d1", "d2", "d3", "d4"), times = 2))
    SummarizedExperiment::colData(dds)$imputed_sex <- factor(c("female", "male")[
        as.integer(SummarizedExperiment::colData(dds)$donor) %% 2 + 1])
    nested <- .model_covariates(dds, var = "condition", covars = "donor")
    expect_false(nested$sex_used)
    expect_equal(nested$covars, "donor")

    # not imputed
    SummarizedExperiment::colData(dds)$imputed_sex <- factor(NA, levels = c("female", "male"))
    expect_false(.model_covariates(dds, var = "condition")$sex_used)
})

test_that(".estimable_covariates drops covariates constant within a subset", {
    dds <- simulate_sexed_dds()
    SummarizedExperiment::colData(dds)$imputed_sex <- factor(rep(c("female", "male"), 4))
    sub <- dds[, SummarizedExperiment::colData(dds)$imputed_sex == "female"]
    expect_message(out <- .estimable_covariates(sub, "imputed_sex", "condition"), "dropped")
    expect_null(out)
})

test_that("run_distance_pca draws one PCA per colour variable, before and after", {
    dds <- simulate_sexed_dds()
    SummarizedExperiment::colData(dds)$imputed_sex <- factor(rep(c("female", "male"), 4))
    SummarizedExperiment::colData(dds)$donor <- factor(rep(paste0("d", 1:4), times = 2))
    out <- run_distance_pca(dds, var = "condition", batch = c("donor", "imputed_sex"),
                            colour_vars = c(Group = "condition", Pair = "donor",
                                            Sex = "imputed_sex"),
                            save_plot = FALSE, save_dir = withr::local_tempdir())
    expect_named(out$pca_by_uncorrected, c("Group", "Pair", "Sex"))
    expect_named(out$pca_by_corrected, c("Group", "Pair", "Sex"))
    expect_s3_class(out$pca_by_corrected$Sex$static, "ggplot")

    plain <- run_distance_pca(dds, var = "condition", save_plot = FALSE,
                              save_dir = withr::local_tempdir())
    expect_named(plain$pca_by_uncorrected, "Group")
    expect_null(plain$pca_by_corrected)
})

test_that(".pca_colour_vars lists group, pair, sex and batch when present", {
    dds <- simulate_sexed_dds()
    SummarizedExperiment::colData(dds)$imputed_sex <- factor(rep(c("female", "male"), 4))
    expect_equal(.pca_colour_vars(dds, "condition", pair = "donor", batch = "run"),
                 c(Group = "condition", Pair = "donor", Sex = "imputed_sex", Batch = "run"))
    SummarizedExperiment::colData(dds)$imputed_sex <- NULL
    expect_equal(.pca_colour_vars(dds, "condition"), c(Group = "condition"))
})

test_that("remove_xy without gene_chr uses biowrenchR's X/Y list, after sex imputation", {
    dds <- simulate_sexed_dds()
    qc <- run_qc(dds, var = "condition", remove_xy = TRUE, save_plot = FALSE,
                 save_dir = withr::local_tempdir())
    expect_false(any(c("Xist", "Ddx3y", "Kdm5d") %in% rownames(qc)))
    expect_equal(as.character(SummarizedExperiment::colData(qc)$imputed_sex),
                 rep(c("female", "male"), 4))
})
