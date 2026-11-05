test_that("plot_lr_contrast_barplot draws a GSEA-style ggplot", {
    acts <- data.frame(
        source = c("TP53&MDM2", "IL6&IL6R", "TNF&TNFRSF1A",
                   "WNT3A&FZD1", "BMP2&BMPR1A", "VEGFA&FLT1"),
        score = c(4.2, 2.1, 0.8, -0.4, -1.8, -3.5),
        p_value = c(1e-4, 0.008, 0.2, 0.4, 0.03, 2e-4),
        stringsAsFactors = FALSE
    )
    p <- plot_lr_contrast_barplot(acts, n_lr = 3, title = "A_vs_B")
    expect_s3_class(p, "ggplot")
    sz <- attr(p, "plot_size")
    expect_true(is.list(sz) && sz$width > 0 && sz$height > 0)
    expect_equal(formals(plot_lr_contrast_barplot)$n_lr, 10)
    expect_equal(formals(run_decoupler_ccc_pipeline)$n_lr, 10)
})

test_that(".ccc_sig_genes expands significant pairs to ligand and receptor", {
    ccc_out <- list(
        net = data.frame(
            source = c("TP53&MDM2", "TP53&MDM2", "IL6&IL6R", "IL6&IL6R"),
            target = c("TP53", "MDM2", "IL6", "IL6R"),
            mor = 1,
            stringsAsFactors = FALSE
        ),
        contrast_acts = list(
            A_vs_B = data.frame(
                source = c("TP53&MDM2", "IL6&IL6R"),
                p_value = c(0.01, 0.8),
                padj = c(0.02, 0.9),
                stringsAsFactors = FALSE
            )
        )
    )
    genes <- DEGreports:::.ccc_sig_genes(ccc_out, "A_vs_B", alpha = 0.05)
    expect_setequal(genes, c("TP53", "MDM2"))
})

test_that("rank_multihit_genes_contrast scores CCC separately from TF decoupleR", {
    de <- data.frame(
        gene = c("TP53", "MDM2", "IL6", "MYC"),
        comparison = "A_vs_B",
        sig = c(TRUE, TRUE, FALSE, FALSE),
        padj = c(0.01, 0.02, 0.4, 0.5),
        log2FoldChange = c(2, -1, 0.1, 0.2),
        stat = c(5, -4, 0.2, 0.1),
        stringsAsFactors = FALSE
    )
    diffexp_out <- list(by_comparison = list(A_vs_B = de), results = de)
    decoupler_tf_out <- list(
        contrast_acts = list(
            A_vs_B = data.frame(
                source = c("MYC"),
                p_value = c(0.01),
                score = c(2),
                stringsAsFactors = FALSE
            )
        )
    )
    ccc_out <- list(
        net = data.frame(
            source = c("TP53&MDM2", "TP53&MDM2"),
            target = c("TP53", "MDM2"),
            mor = 1,
            stringsAsFactors = FALSE
        ),
        contrast_acts = list(
            A_vs_B = data.frame(
                source = "TP53&MDM2",
                score = 2,
                p_value = 0.001,
                padj = 0.001,
                stringsAsFactors = FALSE
            )
        )
    )
    ranked <- rank_multihit_genes_contrast(
        "A_vs_B", diffexp_out,
        decoupler_tf_out = decoupler_tf_out,
        decoupler_ccc_out = ccc_out,
        alpha = 0.05, top_n = 10, org = "human"
    )
    row_p53 <- ranked[ranked$gene == "TP53", , drop = FALSE]
    expect_equal(row_p53[["decoupler-ccc"]][[1]], 1L)
    expect_equal(row_p53[["decoupler-tf"]][[1]], 0L)
    row_myc <- ranked[ranked$gene == "MYC", , drop = FALSE]
    expect_equal(row_myc[["decoupler-tf"]][[1]], 1L)
    expect_equal(row_myc[["decoupler-ccc"]][[1]], 0L)
})

test_that(".ccc_cluster_location splits a cluster label into test and folder", {
    loc <- DEGreports:::.ccc_cluster_location("lrt__full_a_vs_reduced_b cluster 3")
    expect_equal(loc$test, "lrt__full_a_vs_reduced_b")
    expect_equal(loc$dir, "cluster_3")
    plain <- DEGreports:::.ccc_cluster_location("A_vs_B")
    expect_equal(plain$test, "A_vs_B")
})

test_that("run_decoupler_ccc_pipeline scores each degPatterns cluster separately", {
    skip_if_not_installed("decoupleR")
    net <- data.frame(
        source = c("TP53&MDM2", "TP53&MDM2", "BRCA1&IL6", "BRCA1&IL6"),
        target = c("TP53", "MDM2", "BRCA1", "IL6"),
        mor = 1,
        stringsAsFactors = FALSE
    )
    filler <- paste0("G", 1:8)
    test_label <- "lrt__full_time_vs_reduced_intercept"
    de <- data.frame(
        gene = c("TP53", "MDM2", "BRCA1", "IL6", filler),
        stat = c(30, 28, 25, 24, seq_along(filler)),
        comparison = test_label,
        stringsAsFactors = FALSE
    )
    diffexp_out <- list(
        by_comparison = stats::setNames(list(de), test_label),
        results = de
    )
    labels <- paste(test_label, "cluster", 1:3)
    # Cluster 1 holds one complete pair, cluster 2 the other, cluster 3 only
    # a single partner (IL6) and so cannot score any pair.
    cluster_genes <- list(
        c("TP53", "MDM2", filler[1:4]),
        c("BRCA1", "IL6", filler[5:8]),
        c("IL6", filler[1:4])
    )
    by_cluster <- stats::setNames(
        Map(function(genes, label) {
            sub <- de[de$gene %in% genes, , drop = FALSE]
            sub$comparison <- label
            sub
        }, cluster_genes, labels),
        labels
    )
    save_dir <- withr::local_tempdir()
    out <- run_decoupler_ccc_pipeline(
        diffexp_out, cluster_diffexp = list(by_comparison = by_cluster),
        lr_net = net, n_lr = 5, save_dir = save_dir, save_plot = FALSE
    )

    expect_setequal(names(out$cluster_acts), labels[1:2])
    expect_setequal(out$cluster_acts[[labels[1]]]$source, "TP53&MDM2")
    expect_setequal(out$cluster_acts[[labels[2]]]$source, "BRCA1&IL6")
    # Every cluster attempted is listed, including the one that scored nothing.
    expect_equal(out$clusters_by_test[[test_label]], labels)
    expect_true(file.exists(file.path(save_dir, test_label, "cluster_1",
                                      "ulm_lr.csv")))
    expect_true(file.exists(file.path(save_dir, test_label, "cluster_2",
                                      "ulm_lr.csv")))
    expect_false(dir.exists(file.path(save_dir, test_label, "cluster_3")))
    # Whole-test scoring is untouched and still sees both pairs.
    expect_setequal(out$contrast_acts[[test_label]]$source,
                    c("TP53&MDM2", "BRCA1&IL6"))
})

test_that("run_decoupler_ccc_pipeline warns on cluster input with no signatures", {
    skip_if_not_installed("decoupleR")
    de <- data.frame(
        gene = c("TP53", "MDM2", "BRCA1", "IL6"),
        stat = c(4, 3.5, -2, 0.1),
        comparison = "A_vs_B",
        stringsAsFactors = FALSE
    )
    net <- data.frame(
        source = rep("TP53&MDM2", 2),
        target = c("TP53", "MDM2"),
        mor = 1,
        stringsAsFactors = FALSE
    )
    expect_warning(
        run_decoupler_ccc_pipeline(
            list(by_comparison = list(A_vs_B = de), results = de),
            cluster_diffexp = list(by_comparison = list()),
            lr_net = net, save_dir = withr::local_tempdir(), save_plot = FALSE
        ),
        "no `by_comparison` signatures"
    )
})

test_that("run_decoupler_ccc_pipeline writes per-contrast ULM tables", {
    skip_if_not_installed("decoupleR")
    genes <- c("TP53", "MDM2", "BRCA1", "IL6", "MYC")
    de <- data.frame(
        gene = genes,
        comparison = "A_vs_B",
        stat = c(4, 3.5, -2, 0.1, -0.2),
        stringsAsFactors = FALSE
    )
    diffexp_out <- list(by_comparison = list(A_vs_B = de), results = de)
    net <- data.frame(
        source = c("TP53&MDM2", "TP53&MDM2", "BRCA1&IL6", "BRCA1&IL6"),
        target = c("TP53", "MDM2", "BRCA1", "IL6"),
        mor = 1,
        stringsAsFactors = FALSE
    )
    save_dir <- withr::local_tempdir()
    out <- run_decoupler_ccc_pipeline(
        diffexp_out, lr_net = net, n_lr = 5,
        save_dir = save_dir, save_plot = FALSE
    )
    expect_true("A_vs_B" %in% names(out$contrast_acts))
    expect_true(file.exists(file.path(save_dir, "A_vs_B", "ulm_lr.csv")))
    expect_true(nrow(out$contrast_acts$A_vs_B) >= 1L)
})
