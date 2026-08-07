# Tests for the WGCNA module: run_wgcna_pipeline() end-to-end, and its
# individual plotting helpers (plot_soft_threshold, plot_module_sizes,
# plot_module_trait_heatmap, plot_module_eigengene_heatmap,
# plot_module_ora_dotplot).

skip_if_not_installed("WGCNA")
skip_if_not_installed("withr")
skip_if_not_installed("msigdbr")

test_that("run_wgcna_pipeline returns a fully-populated, consistently-keyed result list", {
    dds <- simulate_dds_real_genes(n_genes = 300, n_per_group = 5, seed = 7)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)

    tmp <- withr::local_tempdir()
    out <- run_wgcna_pipeline(dds, var = "condition", save_dir = tmp, save_plot = TRUE,
                               min_module_size = 15)

    expect_type(out, "list")
    expect_true(is.numeric(out$power))
    expect_true(out$power_source %in% c("estimated", "fallback_best_fit"))

    expect_s3_class(out$module_df, "data.frame")
    expect_true(all(c("gene", "module", "module_color") %in% colnames(out$module_df)))

    # Modules are named "Module 1", "Module 2", ... ranked by descending
    # size, with the leftover/unassigned co-expression bucket (WGCNA's
    # "grey") relabelled "Unassigned" and excluded from the numbering.
    module_names <- setdiff(unique(out$module_df$module), "Unassigned")
    expect_true(all(grepl("^Module [0-9]+$", module_names)))
    sizes_by_rank <- as.numeric(table(out$module_df$module)[module_names][
        order(as.numeric(sub("Module ", "", module_names)))
    ])
    expect_true(all(diff(sizes_by_rank) <= 0))

    # Every module-keyed output structure must use the same set of
    # human-readable module names.
    expect_true(all(module_names %in% colnames(out$MEs)) || nrow(out$MEs) == 0)
    expect_true("Unassigned" %in% colnames(out$MEs))
    # hub_genes is keyed by every detected module label, including
    # "Unassigned" (with an empty data frame for it);
    # ora_by_module/ora_plots_by_module deliberately exclude "Unassigned"
    # since running over-representation analysis on co-expression
    # leftovers isn't meaningful, and are keyed first by MSigDB collection
    # then by module (ora_by_module[[collection]][[module]]) since ORA now
    # runs against every requested collection independently.
    all_module_labels <- unique(out$module_df$module)
    expect_setequal(names(out$hub_genes), all_module_labels)
    expect_true(all(sapply(out$ora_by_module, function(x) setequal(names(x), module_names))))
    expect_true(all(sapply(out$ora_plots_by_module, function(x) setequal(names(x), module_names))))
    expect_true(all(out$module_trait_cor$module %in% c(module_names, "Unassigned")))

    # Plot objects are always built (not just inside the save_plot branch),
    # so callers who skip disk-saving still get usable objects back.
    expect_s3_class(out$soft_threshold_plot, "ggplot")
    expect_s3_class(out$module_size_plot, "ggplot")
    expect_s3_class(out$module_trait_heatmap, "ggplot")
    expect_s3_class(out$module_eigengene_heatmap, "pheatmap")

    # save_plot = TRUE writes the expected run-level files to disk.
    expect_true(file.exists(file.path(tmp, "wgcna_soft_threshold.pdf")))
    expect_true(file.exists(file.path(tmp, "wgcna_module_sizes.pdf")))
    expect_true(file.exists(file.path(tmp, "wgcna_module_trait_heatmap.pdf")))
    expect_true(file.exists(file.path(tmp, "wgcna_module_eigengene_heatmap.pdf")))
    expect_true(file.exists(file.path(tmp, "wgcna_dendrogram.pdf")))

    # Per-module artifacts land in sanitized subfolders (Module 1 -> Module_1/).
    real_mods <- setdiff(unique(as.character(out$module_df$module)), "Unassigned")
    for (mod in real_mods) {
        mod_dir <- file.path(tmp, gsub("[^A-Za-z0-9]+", "_", mod))
        expect_true(dir.exists(mod_dir))
        expect_true(file.exists(file.path(mod_dir, "hub_genes.csv")))
    }
})

test_that("run_wgcna_pipeline still returns plot objects when save_plot = FALSE", {
    dds <- simulate_dds_real_genes(n_genes = 200, n_per_group = 4, seed = 11)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)

    tmp <- withr::local_tempdir()
    out <- run_wgcna_pipeline(dds, var = "condition", save_dir = tmp, save_plot = FALSE,
                               min_module_size = 15)

    expect_s3_class(out$soft_threshold_plot, "ggplot")
    expect_s3_class(out$module_size_plot, "ggplot")
    expect_s3_class(out$module_eigengene_heatmap, "pheatmap")
    expect_false(file.exists(file.path(tmp, "wgcna_soft_threshold.pdf")))
})

test_that("plot_module_sizes renders a bar chart with one bar per module, filled by original WGCNA colour", {
    module_df <- data.frame(
        gene = paste0("g", 1:12),
        module = rep(c("Module 1", "Module 2", "Unassigned"), c(6, 4, 2)),
        module_color = rep(c("turquoise", "blue", "grey"), c(6, 4, 2))
    )
    p <- plot_module_sizes(module_df)
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), length(unique(module_df$module)))
})

test_that("plot_module_trait_heatmap draws a tile heatmap annotated with r and p", {
    module_trait_cor <- data.frame(
        module = rep(c("Module 1", "Module 2"), each = 2),
        level = rep(c("Treated", "Control"), 2),
        cor = c(0.8, -0.8, 0.3, -0.3),
        pvalue = c(0.001, 0.001, 0.2, 0.2)
    )
    p <- plot_module_trait_heatmap(module_trait_cor)
    expect_s3_class(p, "ggplot")
    geom_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomTile" %in% geom_classes)
    expect_true("GeomText" %in% geom_classes)
})

test_that("plot_module_eigengene_heatmap draws a clustered pheatmap of modules x samples", {
    MEs <- data.frame(
        `Module 1` = c(0.5, 0.3, -0.4, -0.2),
        `Module 2` = c(-0.1, 0.2, 0.6, -0.3),
        Unassigned = c(0.05, -0.05, 0.1, -0.1),
        check.names = FALSE,
        row.names = paste0("sample", 1:4)
    )
    p <- plot_module_eigengene_heatmap(MEs)
    expect_s3_class(p, "pheatmap")
    expect_true(!is.null(p$gtable))
})

test_that("plot_module_ora_dotplot renders top_n enriched terms for a module, sized by gene-set size", {
    ora_res <- data.frame(
        pathway = paste0("HALLMARK_SET_", 1:5),
        padj = c(0.001, 0.01, 0.02, 0.03, 0.5),
        size = c(20, 15, 10, 8, 2)
    )
    p <- plot_module_ora_dotplot(ora_res, top_n = 3, title = "Module: Module 1")
    expect_s3_class(p, "ggplot")
    expect_lte(nrow(p$data), 3)
    geom_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomPoint" %in% geom_classes)
    expect_equal(rlang::as_label(p$mapping$x), "size")
})

test_that("run_wgcna_pipeline runs per-module ORA against multiple MSigDB collections by default", {
    dds <- simulate_dds_real_genes(n_genes = 300, n_per_group = 5, seed = 7)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)

    tmp <- withr::local_tempdir()
    out <- run_wgcna_pipeline(dds, var = "condition", save_dir = tmp, save_plot = FALSE,
                               min_module_size = 15, ora_collection = c("HALLMARK", "GOBP"))

    expect_setequal(names(out$ora_by_module), c("HALLMARK", "GOBP"))
    expect_setequal(names(out$ora_plots_by_module), c("HALLMARK", "GOBP"))
    module_names <- setdiff(unique(out$module_df$module), "Unassigned")
    expect_setequal(names(out$ora_by_module[["HALLMARK"]]), module_names)
    expect_setequal(names(out$ora_by_module[["GOBP"]]), module_names)
})
