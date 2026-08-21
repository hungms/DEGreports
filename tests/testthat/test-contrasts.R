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
    expect_equal(
        deparse(resolve_design(c("celltype", "C_mutation", "K_mutation"),
                               include_interactions = TRUE)),
        deparse(~ celltype * C_mutation * K_mutation)
    )
    expect_equal(
        deparse(resolve_design(c("celltype", "C_mutation", "K_mutation"))),
        deparse(~ celltype + C_mutation + K_mutation)
    )
    expect_equal(deparse(resolve_design("Group2", design = ~Batch + Group2)),
                 deparse(~Batch + Group2))
})

test_that("expand_contrast_specs three-factor main effects ignore the other factors", {
    dds <- simulate_three_factor_dds()
    var <- c("celltype", "C_mutation", "K_mutation")

    # include_interactions = FALSE: main effects only, no strata.
    specs <- expand_contrast_specs(dds, var)
    expect_true(all(specs$class == "main"))
    expect_equal(specs$label, c("CC_vs_CB", "C_vs_noC", "K_vs_noK"))
    expect_equal(specs$test_var, var)
    expect_true(all(is.na(specs$stratum_var)))
    expect_true(all(is.na(specs$stratum_level)))
    expect_false(any(grepl("/", specs$label, fixed = TRUE)))

    # include_interactions = TRUE keeps them, and adds the interactions.
    both <- expand_contrast_specs(dds, var, include_interactions = TRUE)
    expect_equal(sum(both$class == "main"), 3L)
    expect_equal(both$label[both$class == "main"], specs$label)
    expect_equal(sum(both$class == "interaction"), 7L)
    expect_false(any(both$class %in% c("within_A", "within_B", "within_C")))
    expect_false(anyDuplicated(both$label) > 0)

    # stratified = TRUE is what brings the within-stratum grid back.
    strat <- expand_contrast_specs(dds, var, stratified = TRUE)
    expect_equal(sum(strat$class == "main"), 3L)
    expect_equal(sum(strat$class != "main"), 12L)
    expect_equal(sum(strat$class == "interaction"), 0L)
})

test_that("interaction specs are difference-of-differences over one full model", {
    dds <- simulate_three_factor_dds()
    var <- c("celltype", "C_mutation", "K_mutation")
    ix <- expand_contrast_specs(dds, var, include_interactions = TRUE)
    ix <- ix[ix$class == "interaction", , drop = FALSE]

    # One row per factor pair per level of the third factor (3 x 2), plus the
    # three-way row: 7, replacing the 12 within-stratum contrasts.
    expect_equal(nrow(ix), 7L)
    expect_setequal(ix$label, c(
        "C_vs_noC_x_K_vs_noK_in_CB", "C_vs_noC_x_K_vs_noK_in_CC",
        "CC_vs_CB_x_K_vs_noK_in_noC", "CC_vs_CB_x_K_vs_noK_in_C",
        "CC_vs_CB_x_C_vs_noC_in_noK", "CC_vs_CB_x_C_vs_noC_in_K",
        "CC_vs_CB_x_C_vs_noC_x_K_vs_noK"
    ))
    expect_true(all(ix$type == "interaction"))

    coefs <- function(label) {
        row <- ix[ix$label == label, , drop = FALSE]
        list(num = strsplit(row$comparison[[1]], "+", fixed = TRUE)[[1]],
             den = if (is.na(row$coef_subtract[[1]])) character(0) else
                 strsplit(row$coef_subtract[[1]], "+", fixed = TRUE)[[1]])
    }
    # In the reference stratum the interaction is a single coefficient; in
    # the other stratum the three-way term is added to it.
    expect_equal(coefs("C_vs_noC_x_K_vs_noK_in_CB"),
                 list(num = "C_mutationC.K_mutationK", den = character(0)))
    expect_equal(coefs("C_vs_noC_x_K_vs_noK_in_CC"),
                 list(num = c("C_mutationC.K_mutationK",
                              "celltypeCC.C_mutationC.K_mutationK"),
                      den = character(0)))
    expect_equal(coefs("CC_vs_CB_x_K_vs_noK_in_C"),
                 list(num = c("celltypeCC.K_mutationK",
                              "celltypeCC.C_mutationC.K_mutationK"),
                      den = character(0)))
    expect_equal(coefs("CC_vs_CB_x_C_vs_noC_x_K_vs_noK"),
                 list(num = "celltypeCC.C_mutationC.K_mutationK",
                      den = character(0)))

    # The held-fixed factor is recorded, except for the three-way row.
    held <- contrast_stratum_labels(ix)
    expect_equal(sum(is.na(held)), 1L)
    expect_true("celltype = CB" %in% held)
    expect_equal(ix$test_var[ix$label == "C_vs_noC_x_K_vs_noK_in_CB"],
                 "C_mutation+K_mutation")

    # Asking a pair in the other order must not add a second, mirrored row:
    # "does K differ between C and noC" is the same contrast as "does C
    # differ between K and noK".
    flipped <- expand_contrast_specs(
        dds, c("celltype", "K_mutation", "C_mutation"), include_interactions = TRUE
    )
    expect_equal(sum(flipped$class == "interaction"), 7L)
})

test_that("interaction specs pool the rest when a factor has three levels", {
    dds <- simulate_three_factor_dds()
    SummarizedExperiment::colData(dds)$celltype <- factor(
        rep(c("CB", "CC", "CX"), length.out = ncol(dds)), levels = c("CB", "CC", "CX")
    )
    ix <- expand_contrast_specs(dds, c("celltype", "C_mutation", "K_mutation"),
                                include_interactions = TRUE)
    ix <- ix[ix$class == "interaction", , drop = FALSE]
    # celltype one-vs-rest (3) x other factor (1) x strata (2), twice over,
    # plus C x K in each of 3 celltypes, plus 3 three-way rows.
    expect_equal(nrow(ix), 18L)
    expect_false(anyDuplicated(ix$label) > 0)
    # One-vs-rest, not pairwise: no CX_vs_CC anywhere.
    expect_true("CX_vs_CB+CC_x_K_vs_noK_in_C" %in% ix$label)
    expect_false(any(grepl("CX_vs_CC|CC_vs_CB_x", ix$label)))
    expect_setequal(
        grep("_x_K_vs_noK_in_noC$", ix$label, value = TRUE),
        c("CB_vs_CC+CX_x_K_vs_noK_in_noC", "CC_vs_CB+CX_x_K_vs_noK_in_noC",
          "CX_vs_CB+CC_x_K_vs_noK_in_noC")
    )

    # A pooled side averages, so the row carries weights instead of the
    # `+`-joined names a summed list contrast would use.
    weights_of <- function(label) {
        row <- ix[ix$label == label, , drop = FALSE]
        expect_true(is.na(row$comparison[[1]]))
        expect_true(is.na(row$coef_subtract[[1]]))
        .parse_coef_weights(row$coef_weights[[1]])
    }
    # CB is the reference level, so its corner is structurally zero and only
    # the other pooled level survives.
    expect_equal(weights_of("CX_vs_CB+CC_x_K_vs_noK_in_noC"),
                 c(celltypeCX.K_mutationK = 1, celltypeCC.K_mutationK = -0.5))
    expect_equal(weights_of("CB_vs_CC+CX_x_K_vs_noK_in_noC"),
                 c(celltypeCC.K_mutationK = -0.5, celltypeCX.K_mutationK = -0.5))
    # Another stratum adds each corner's three-way term at the same weight.
    expect_equal(
        weights_of("CX_vs_CB+CC_x_K_vs_noK_in_C"),
        c(celltypeCX.K_mutationK = 1, celltypeCX.C_mutationC.K_mutationK = 1,
          celltypeCC.K_mutationK = -0.5, celltypeCC.C_mutationC.K_mutationK = -0.5)
    )
    # The three-way row: the two-way interaction common to every celltype
    # cancels between the two sides and must not be left in.
    expect_equal(
        weights_of("CX_vs_CB+CC_x_C_vs_noC_x_K_vs_noK"),
        c(celltypeCX.C_mutationC.K_mutationK = 1,
          celltypeCC.C_mutationC.K_mutationK = -0.5)
    )

    # The two-level pair is untouched: still a plain sum of coefficients.
    summed <- ix[ix$label == "C_vs_noC_x_K_vs_noK_in_CX", , drop = FALSE]
    expect_true(is.na(summed$coef_weights[[1]]))
    expect_equal(summed$comparison[[1]],
                 "C_mutationC.K_mutationK+celltypeCX.C_mutationC.K_mutationK")
})

test_that("expand_contrast_specs main effects with >2 levels include one-vs-rest", {
    dds <- simulate_dds(n_genes = 20, group_levels = c("A", "B", "C"), n_per_group = 4)
    SummarizedExperiment::colData(dds)$treatment <-
        factor(rep(c("ctrl", "drug"), length.out = ncol(dds)))
    specs <- expand_contrast_specs(dds, c("condition", "treatment"))

    one_vs_rest <- specs$label[specs$test_var == "condition" &
                                   specs$type == "one_vs_all"]
    expect_setequal(one_vs_rest, c("A_vs_B+C", "B_vs_A+C", "C_vs_A+B"))
    expect_equal(specs$label[specs$test_var == "treatment"], "drug_vs_ctrl")
})

test_that("expand_contrast_specs disambiguates main effects sharing level names", {
    dds <- simulate_dds(n_genes = 20, group_levels = c("yes", "no"), n_per_group = 4)
    SummarizedExperiment::colData(dds)$other <-
        factor(rep(c("yes", "no"), length.out = ncol(dds)))
    specs <- expand_contrast_specs(dds, c("condition", "other"))
    expect_equal(specs$label, c("condition_yes_vs_no", "other_yes_vs_no"))
    expect_equal(specs$comparison, c("yes_vs_no", "yes_vs_no"))
})

test_that("expand_contrast_specs three-factor yields unique within-stratum labels", {
    set.seed(1)
    n_genes <- 20
    n_samp <- 24
    counts <- matrix(rnbinom(n_genes * n_samp, mu = 50, size = 10), nrow = n_genes)
    rownames(counts) <- paste0("g", seq_len(n_genes))
    colnames(counts) <- paste0("s", seq_len(n_samp))
    coldata <- data.frame(
        celltype = factor(rep(c("CB", "CC"), each = 12), levels = c("CB", "CC")),
        C_mutation = factor(rep(rep(c("noC", "C"), each = 6), times = 2),
                            levels = c("noC", "C")),
        K_mutation = factor(rep(c("noK", "K"), times = 12), levels = c("noK", "K")),
        row.names = colnames(counts)
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(
        counts, coldata, design = ~ celltype * C_mutation * K_mutation
    )
    specs <- expand_contrast_specs(
        dds, c("celltype", "C_mutation", "K_mutation"),
        stratified = TRUE
    )
    specs <- specs[specs$class != "main", , drop = FALSE]
    expect_equal(nrow(specs), 12L)
    expect_equal(sum(specs$class == "within_A"), 4L)
    expect_equal(sum(specs$class == "within_B"), 4L)
    expect_equal(sum(specs$class == "within_C"), 4L)
    expect_equal(sum(specs$class == "interaction"), 0L)
    expect_false(anyDuplicated(specs$label) > 0)
    expect_true(all(grepl("^CC_vs_CB/", specs$label[specs$class == "within_A"])))
    expect_true(all(grepl("^C_vs_noC/", specs$label[specs$class == "within_B"])))
    expect_true(all(grepl("^K_vs_noK/", specs$label[specs$class == "within_C"])))
    expect_true("C_vs_noC/CB__noK" %in% specs$label)
    expect_true("CC_vs_CB/noC__noK" %in% specs$label)
    expect_true("K_vs_noK/CB__noC" %in% specs$label)
    expect_equal(specs$comparison[specs$label == "C_vs_noC/CB__noK"], "C_vs_noC")
    expect_equal(specs$stratum_tag[specs$label == "C_vs_noC/CB__noK"], "CB__noK")
    # Mutant vs WT direction (not alphabetical noC_vs_C); no long flat prefixes.
    expect_false(any(grepl("noC_vs_C|noK_vs_K|CB_vs_CC", specs$label)))
    expect_false(any(grepl("C_mutation_|K_mutation_|celltype_", specs$label)))
    expect_equal(contrast_file_stem("C_vs_noC/CB__noK"), "CB__noK")
    expect_equal(contrast_file_stem("IgG_vs_IgM"), "IgG_vs_IgM")
})

test_that("expand_contrast_specs single-factor matches generate_comparisons", {
    dds <- simulate_dds(n_genes = 20, group_levels = c("A", "B", "C"), n_per_group = 2)
    specs <- expand_contrast_specs(dds, "condition")
    expect_equal(sort(specs$label), sort(generate_comparisons(levels(dds$condition))))
    expect_true(all(specs$class == "single"))
    expect_true("coef_subtract" %in% names(specs))
})

test_that("expand_contrast_specs 2x3 default gives main effects only", {
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
    # Tissue: 1 pairwise. Isotype: 3 pairwise + 3 one-vs-rest.
    expect_equal(nrow(specs), 7L)
    expect_true(all(specs$class == "main"))
    expect_equal(sum(specs$type == "pairwise"), 4L)
    expect_equal(sum(specs$type == "one_vs_all"), 3L)
    expect_true(all(c("Spleen_vs_BM", "IgG_vs_IgM", "IgA_vs_IgM+IgG") %in% specs$label))
})

test_that("expand_contrast_specs 2x3 interactions cover every level pair", {
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
    # 7 main effects + one interaction per Isotype level, each against the
    # mean of the other two, with no third factor to hold fixed.
    expect_equal(nrow(specs), 10L)
    expect_equal(sum(specs$class == "main"), 7L)

    ix <- specs[specs$class == "interaction", , drop = FALSE]
    expect_setequal(ix$label, c("Spleen_vs_BM_x_IgM_vs_IgG+IgA",
                                "Spleen_vs_BM_x_IgG_vs_IgM+IgA",
                                "Spleen_vs_BM_x_IgA_vs_IgM+IgG"))
    expect_true(all(is.na(ix$stratum_var)))
    # Every row averages, so none of them is a summable list contrast.
    expect_true(all(is.na(ix$comparison)))
    expect_true(all(is.na(ix$coef_subtract)))
    weights_of <- function(label) {
        .parse_coef_weights(ix$coef_weights[ix$label == label][[1]])
    }
    expect_equal(weights_of("Spleen_vs_BM_x_IgG_vs_IgM+IgA"),
                 c(TissueSpleen.IsotypeIgG = 1, TissueSpleen.IsotypeIgA = -0.5))
    # IgM is the reference level: its own corner drops out and only the two
    # pooled levels are left, both negative.
    expect_equal(weights_of("Spleen_vs_BM_x_IgM_vs_IgG+IgA"),
                 c(TissueSpleen.IsotypeIgG = -0.5, TissueSpleen.IsotypeIgA = -0.5))

    # The within-stratum grid is now opt-in and unaffected by interactions.
    strat <- expand_contrast_specs(dds, c("Tissue", "Isotype"), stratified = TRUE)
    expect_equal(sum(strat$class %in% c("within_A", "within_B")), 15L)
    expect_true(all(c("BM__IgG_vs_IgM", "Spleen__IgA_vs_IgM+IgG",
                      "IgM__Spleen_vs_BM") %in% strat$label))
    expect_equal(sum(strat$class == "interaction"), 0L)
})

test_that("expand_contrast_specs genotype x celltype pools the other genotypes", {
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
    # One row per genotype against the mean of the other three (4), not one
    # per combn(4, 2) pair, each against the one celltype pair.
    expect_equal(nrow(ix), 4L)
    expect_true("CK_vs_WT+C+K_x_CC_vs_CB" %in% ix$label)
    pooled <- ix[ix$label == "CK_vs_WT+C+K_x_CC_vs_CB", , drop = FALSE]
    expect_true(is.na(pooled$comparison[[1]]))
    expect_equal(
        .parse_coef_weights(pooled$coef_weights[[1]]),
        c(genotypeCK.celltypeCC = 1, genotypeC.celltypeCC = -1 / 3,
          genotypeK.celltypeCC = -1 / 3)
    )
})

test_that("interaction labels stay unique when factors share level names", {
    dds <- simulate_dds(n_genes = 20, group_levels = c("yes", "no"), n_per_group = 4)
    SummarizedExperiment::colData(dds)$other <-
        factor(rep(c("yes", "no"), length.out = ncol(dds)), levels = c("yes", "no"))
    specs <- expand_contrast_specs(dds, c("condition", "other"),
                                   include_interactions = TRUE)
    ix <- specs[specs$class == "interaction", , drop = FALSE]
    expect_equal(nrow(ix), 1L)
    expect_equal(ix$label, "conditionno_vs_conditionyes_x_otherno_vs_otheryes")
    expect_false(anyDuplicated(specs$label) > 0)
})

test_that("run_diffexp_pipeline fits main effects on all samples, additively", {
    dds <- simulate_three_factor_dds()
    var <- c("celltype", "C_mutation", "K_mutation")
    seen <- list()
    fake_res <- function(comparison) {
        data.frame(
            gene = rownames(dds), comparison = comparison,
            baseMean = 1, mean_cpm = 1, log2FoldChange = 0,
            log2FoldChange_unshrunken = 0, lfcSE = 1, stat = 0,
            pvalue = 1, padj = 1, rank_score = 0, sig = FALSE,
            stringsAsFactors = FALSE
        )
    }
    local_mocked_bindings(
        run_diffexp = function(dds, var, comparison, ...) {
            seen[[length(seen) + 1L]] <<- list(
                design = paste(deparse(DESeq2::design(dds)), collapse = ""),
                n_samples = ncol(dds), var = var, comparison = comparison
            )
            fake_res(comparison)
        }
    )

    out <- run_diffexp_pipeline(
        dds, var = var, stratified = TRUE,
        save_dir = withr::local_tempdir(), save_plot = FALSE
    )
    specs <- out$contrast_specs
    seen <- do.call(rbind, lapply(seen, as.data.frame))

    main <- seen[seen$comparison %in% specs$comparison[specs$class == "main"] &
                     seen$n_samples == ncol(dds), , drop = FALSE]
    expect_equal(nrow(main), 3L)
    # Every factor in the design, tested one additively, on every sample.
    expect_setequal(main$design, c(
        "~C_mutation + K_mutation + celltype",
        "~celltype + K_mutation + C_mutation",
        "~celltype + C_mutation + K_mutation"
    ))
    # Strata still get their own subset and a single-term design.
    strata <- seen[seen$n_samples < ncol(dds), , drop = FALSE]
    expect_equal(nrow(strata), 12L)
    expect_setequal(strata$design, c("~celltype", "~C_mutation", "~K_mutation"))
    expect_true(all(strata$n_samples == 6L))
})

test_that("a main effect and its strata write to distinct paths", {
    dds <- simulate_three_factor_dds()
    save_dir <- withr::local_tempdir()
    local_mocked_bindings(
        run_diffexp = function(dds, var, comparison, ...) {
            data.frame(gene = rownames(dds), comparison = comparison,
                       log2FoldChange = 0, padj = 1, stringsAsFactors = FALSE)
        }
    )
    run_diffexp_pipeline(
        dds, var = c("celltype", "C_mutation", "K_mutation"),
        stratified = TRUE, save_dir = save_dir, save_plot = FALSE
    )
    expect_true(file.exists(file.path(save_dir, "C_vs_noC", "C_vs_noC_diffexp.csv")))
    expect_true(file.exists(file.path(save_dir, "C_vs_noC", "CB__noK",
                                      "CB__noK_diffexp.csv")))
})

test_that("run_diffexp_pipeline fits one product model for every interaction row", {
    dds <- simulate_three_factor_dds()
    var <- c("celltype", "C_mutation", "K_mutation")
    seen <- list()
    local_mocked_bindings(
        run_diffexp_coef = function(dds, coef, coef_subtract = NULL, ...,
                                    comparison, prefit = NULL) {
            seen[[length(seen) + 1L]] <<- list(
                comparison = comparison, coef = list(coef),
                coef_subtract = list(coef_subtract),
                design = paste(deparse(DESeq2::design(prefit)), collapse = ""),
                resultsNames = list(DESeq2::resultsNames(prefit))
            )
            data.frame(gene = rownames(dds), comparison = comparison,
                       log2FoldChange = 0, padj = 1, stringsAsFactors = FALSE)
        },
        run_diffexp = function(dds, var, comparison, ...) {
            data.frame(gene = rownames(dds), comparison = comparison,
                       log2FoldChange = 0, padj = 1, stringsAsFactors = FALSE)
        }
    )

    run_diffexp_pipeline(
        dds, var = var, include_interactions = TRUE,
        save_dir = withr::local_tempdir(), save_plot = FALSE
    )
    expect_length(seen, 7L)
    # One prefit of the full product design, shared by every row, and the
    # coefficient names reach the fitter already split.
    expect_true(all(vapply(seen, function(x) x$design, character(1)) ==
                        "~celltype * C_mutation * K_mutation"))
    in_cc <- Filter(function(x) x$comparison == "C_vs_noC_x_K_vs_noK_in_CC", seen)[[1]]
    expect_equal(in_cc$coef[[1]], c("C_mutationC.K_mutationK",
                                    "celltypeCC.C_mutationC.K_mutationK"))
    expect_length(in_cc$coef_subtract[[1]], 0L)
    expect_true(all(unlist(lapply(seen, function(x) x$coef[[1]])) %in%
                        in_cc$resultsNames[[1]]))
})

test_that("interaction contrasts refuse an incomplete factorial by name", {
    dds <- simulate_three_factor_dds()
    keep <- !(dds$celltype == "CC" & dds$C_mutation == "C" & dds$K_mutation == "K")
    dds <- dds[, keep]
    expect_error(
        run_diffexp_pipeline(
            dds, var = c("celltype", "C_mutation", "K_mutation"),
            include_interactions = TRUE,
            save_dir = withr::local_tempdir(), save_plot = FALSE
        ),
        "celltype=CC, C_mutation=C, K_mutation=K",
        fixed = TRUE
    )
})

test_that("expand_lrt_specs prefixes batch onto every full and reduced model", {
    specs <- expand_lrt_specs(var = NULL, time = "isotype", batch = "batch")
    expect_true(nrow(specs) >= 1L)
    expect_true(all(grepl("~ batch", specs$full, fixed = TRUE)))
    expect_true(all(grepl("~ batch", specs$reduced, fixed = TRUE)))
    expect_equal(specs$full[[1]], "~ batch + isotype")
    expect_equal(specs$reduced[[1]], "~ batch")
})

test_that("contrast_class_title uses factor names", {
    expect_equal(contrast_class_title("within_A", c("Tissue", "Isotype")),
                 "Within Tissue \u2014 Isotype effects")
    expect_equal(contrast_class_title("within_B", c("Tissue", "Isotype")),
                 "Within Isotype \u2014 Tissue effects")
    expect_equal(contrast_class_title("interaction", c("Tissue", "Isotype")),
                 "Interactions (difference of differences)")
    expect_equal(contrast_class_title("main", c("Tissue", "Isotype")),
                 "Main effects (across all other factors)")
})
