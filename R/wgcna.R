#' Pick a soft-thresholding power for WGCNA network construction
#'
#' Thin wrapper around `WGCNA::pickSoftThreshold()` that also resolves a
#' single recommended power: `WGCNA::pickSoftThreshold()`'s own
#' `powerEstimate` (the smallest power whose scale-free-topology fit
#' `SFT.R.sq` clears `rsq_cutoff`) when available, falling back to the
#' power with the single highest `SFT.R.sq` across `power_vector` when no
#' power clears the cutoff (with a warning, since the resulting network
#' may not be very scale-free).
#'
#' @param datExpr Samples-by-genes numeric matrix (e.g. `t(vst_matrix)`).
#' @param power_vector Candidate powers to evaluate. Default `1:10` plus
#'   every even power from `12` to `20`.
#' @param rsq_cutoff Minimum scale-free-topology model fit (`SFT.R.sq`) to
#'   accept. Default `0.85`.
#' @param network_type One of `"signed"`, `"unsigned"`, `"signed hybrid"`,
#'   passed to `WGCNA::pickSoftThreshold()`. Default `"signed"`.
#' @return A list with elements `fit_indices` (the full per-power table
#'   from `WGCNA::pickSoftThreshold()`), `power` (the resolved
#'   recommended power), and `power_source` (`"estimated"` if
#'   `powerEstimate` cleared the cutoff, `"fallback_best_fit"` otherwise).
#' @noRd
pick_wgcna_power <- function(datExpr, power_vector = c(1:10, seq(12, 20, 2)),
                              rsq_cutoff = 0.85, network_type = "signed") {
    sft <- quiet(WGCNA::pickSoftThreshold(
        datExpr, powerVector = power_vector, networkType = network_type,
        RsquaredCut = rsq_cutoff, verbose = 0
    ))
    fit <- sft$fitIndices

    if (!is.na(sft$powerEstimate)) {
        return(list(fit_indices = fit, power = sft$powerEstimate, power_source = "estimated"))
    }

    best <- fit$Power[which.max(fit$SFT.R.sq)]
    warning(
        "No candidate power reached the scale-free-topology fit cutoff (R.sq >= ", rsq_cutoff,
        "); falling back to power = ", best, " (the best-fitting candidate, R.sq = ",
        sprintf("%.2f", max(fit$SFT.R.sq)), "). Consider widening `power_vector` or relaxing `rsq_cutoff`.",
        call. = FALSE
    )
    list(fit_indices = fit, power = best, power_source = "fallback_best_fit")
}

#' Plot scale-free-topology fit and mean connectivity across candidate powers
#'
#' Renders the two diagnostic curves conventionally used to justify a
#' WGCNA soft-thresholding power choice: scale-free-topology model fit
#' (`SFT.R.sq`) and mean connectivity, both against candidate power, with
#' the chosen power marked by a vertical dashed line.
#'
#' @param fit_indices The `fit_indices` element of `pick_wgcna_power()`'s
#'   return value (a data frame with columns `Power`, `SFT.R.sq`, `mean.k.`).
#' @param power The chosen power (as returned by `pick_wgcna_power()`), used
#'   to draw the reference line on both panels.
#' @param rsq_cutoff The fit cutoff used, drawn as a horizontal reference
#'   line on the scale-free-topology panel. Default `0.85`.
#' @return A `patchwork`-combined `ggplot` object (two panels side by side).
#' @export
plot_soft_threshold <- function(fit_indices, power, rsq_cutoff = 0.85) {
    # Both panels: x is continuous (candidate power) with zero-padding via
    # scale_x_continuous(); the y-range on the fit panel is set explicitly
    # via coord_cartesian() (to accommodate the fixed 1.05 ceiling regardless
    # of data), so its expansion is zeroed through the underlying continuous
    # scale rather than coord_cartesian()'s own (non-expansion) argument.
    p_fit <- ggplot2::ggplot(fit_indices, ggplot2::aes(x = .data$Power, y = .data$SFT.R.sq)) +
        ggplot2::geom_hline(yintercept = rsq_cutoff, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_vline(xintercept = power, linetype = "dashed", colour = "#B2182B") +
        ggplot2::geom_line(colour = "grey40") +
        ggplot2::geom_point(size = 2, colour = "#2166AC") +
        ggplot2::geom_text(ggplot2::aes(label = .data$Power), vjust = -0.8, size = 3) +
        ggplot2::labs(x = "Soft-thresholding power", y = expression(paste("Scale-free topology fit, ", R^2)),
                      title = "Scale independence") +
        ggplot2::scale_x_continuous(expand = zero_expansion()) +
        ggplot2::scale_y_continuous(expand = zero_expansion()) +
        # unclipped_coord() (rather than a bare coord_cartesian()) keeps
        # this panel's explicit ylim while still disabling clipping, so
        # the zero-padded axes don't slice off the point markers / power
        # labels sitting at the data extrema.
        unclipped_coord(ylim = c(min(0, min(fit_indices$SFT.R.sq)), 1.05)) +
        # Intentionally NOT squared (no square_panel_theme()) -- WGCNA
        # diagnostic plots keep their natural wide/rectangular shape.
        ggplot2::theme_bw()

    p_conn <- ggplot2::ggplot(fit_indices, ggplot2::aes(x = .data$Power, y = .data$mean.k.)) +
        ggplot2::geom_vline(xintercept = power, linetype = "dashed", colour = "#B2182B") +
        ggplot2::geom_line(colour = "grey40") +
        ggplot2::geom_point(size = 2, colour = "#2166AC") +
        ggplot2::geom_text(ggplot2::aes(label = .data$Power), vjust = -0.8, size = 3) +
        ggplot2::labs(x = "Soft-thresholding power", y = "Mean connectivity", title = "Mean connectivity") +
        ggplot2::scale_x_continuous(expand = zero_expansion()) +
        ggplot2::scale_y_continuous(expand = zero_expansion()) +
        unclipped_coord() +
        # Intentionally NOT squared (no square_panel_theme()), matching
        # p_fit above.
        ggplot2::theme_bw()

    if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p_fit, p_conn, nrow = 1)
    } else {
        gridExtra_ok <- requireNamespace("gridExtra", quietly = TRUE)
        if (gridExtra_ok) gridExtra::grid.arrange(p_fit, p_conn, nrow = 1) else p_fit
    }
}

#' Run the WGCNA co-expression module detection pipeline
#'
#' Builds a signed weighted gene co-expression network from
#' variance-stabilized counts, picks a soft-thresholding power (see
#' `pick_wgcna_power()`), detects modules via `WGCNA::blockwiseModules()`,
#' correlates module eigengenes with the grouping variable, and runs a
#' Hallmark over-representation analysis (`fgsea::fora()`) for every
#' non-grey module. This is a heavier, exploratory-analysis companion to
#' [run_diffexp_pipeline()]/[run_gsea_pipeline()]: rather than testing one
#' gene at a time against a design, it groups genes into co-expression
#' modules using every sample's expression profile regardless of group
#' labels, then asks which modules' summary profile (eigengene) tracks the
#' grouping variable and which pathways each such module is enriched for.
#'
#' @param dds A `DESeqDataSet` object (typically post-[run_qc()]).
#' @param var Column name in `colData(dds)` identifying the biological
#'   group of interest; module eigengenes are correlated against this.
#' @param power Soft-thresholding power. Default `NULL` auto-selects via
#'   `pick_wgcna_power()`.
#' @param power_vector Candidate powers evaluated when `power = NULL`.
#'   Default `c(1:10, seq(12, 20, 2))`.
#' @param rsq_cutoff Minimum scale-free-topology fit accepted when
#'   auto-selecting `power`. Default `0.85`.
#' @param network_type One of `"signed"` (default), `"unsigned"`, `"signed
#'   hybrid"`.
#' @param min_module_size Minimum module size, passed to
#'   `WGCNA::blockwiseModules()`. Default `30`.
#' @param merge_cut_height Dendrogram height below which similar modules
#'   are merged, passed to `WGCNA::blockwiseModules()`. Default `0.25`.
#' @param max_modules Soft upper bound on the number of non-grey modules
#'   kept after detection. Default `15`. When WGCNA finds more modules,
#'   closely related modules are iteratively merged (`mergeCloseModules`)
#'   by raising the cut height until the count is at most `max_modules`
#'   (or the cut height hits `0.95`). When the natural module count is
#'   already ≤ `max_modules`, that count is kept unchanged. Set to `NULL`
#'   to disable post-hoc merging.
#' @param cor_type One of `"pearson"` (default) or `"bicor"` (robust
#'   biweight midcorrelation), passed to `WGCNA::blockwiseModules()`.
#' @param top_n_hub_genes Number of highest-|kME| hub genes reported per
#'   module. Default `10`.
#' @param ora_collection MSigDB collection alias(es) for the per-module
#'   over-representation analysis (see `validate_msigdb_collections()`):
#'   one or more of `"HALLMARK"`, `"GOBP"`, `"GOCC"`, `"GOMF"`,
#'   `"KEGG_MEDICUS"`, `"KEGG_LEGACY"`, `"REACTOME"`, `"PID"`. Default
#'   matches the GSEA default
#'   (`c("GOBP", "GOCC", "GOMF", "KEGG_MEDICUS", "KEGG_LEGACY", "REACTOME", "PID")`).
#'   [run_deseq2_pipeline()] passes its `collections` argument through so
#'   WGCNA ORA and GSEA use the same set. The leftover/unassigned "grey"
#'   module (renamed "Unassigned") is always excluded (see `ora_by_module`
#'   below and the package vignette).
#' @param top_n_ora Number of top (by adjusted p-value) pathways per
#'   module shown in the ORA dot plot. Default `10`.
#' @param gene_col Name of the gene-identifier column in output tables.
#'   Default `"gene"`.
#' @param save_dir Directory to write outputs to. Default `getwd()`.
#' @param save_plot Logical. If `TRUE` (default), saves the soft-threshold
#'   diagnostic plot, module dendrogram, module-trait heatmap, and one ORA
#'   barplot per module to `save_dir`.
#' @return A list with elements: `power` (the power used), `power_source`
#'   (`"user"`, `"estimated"`, or `"fallback_best_fit"`), `fit_indices`
#'   (the soft-threshold candidate table), `module_df` (data frame: gene,
#'   module color, module label), `MEs` (module eigengene data frame, one
#'   row per sample), `module_trait_cor` (data frame: module, correlation
#'   and p-value against `var`), `hub_genes` (named list of per-module top
#'   hub-gene data frames), `ora_by_module` (named list, keyed first by
#'   MSigDB collection alias and then by module, of `fgsea::fora()` result
#'   data frames -- i.e. `ora_by_module[[collection]][[module]]`; the
#'   "grey"/unassigned module is excluded from every collection, since it
#'   is a heterogeneous leftover bucket of genes with no coherent
#'   co-expression signal, not a real module, so pathway enrichment on it
#'   would not be meaningful), `net` (the raw
#'   `WGCNA::blockwiseModules()` result, for advanced use), and the
#'   pre-built plot objects `soft_threshold_plot`, `module_size_plot`,
#'   `module_trait_heatmap` (condition-level, `ggplot`),
#'   `module_eigengene_heatmap` (sample-level, `pheatmap`, from
#'   [plot_module_eigengene_heatmap()]), and `ora_plots_by_module` (named
#'   list with the same `[[collection]][[module]]` keying as
#'   `ora_by_module`; `NULL` entries for modules with no
#'   significant/overlapping pathways) -- these are the same objects
#'   written to `save_dir` when `save_plot = TRUE`, exposed directly so
#'   callers (e.g. the bundled report) don't need to re-read the PDFs
#'   from disk.
#' @examples
#' \dontrun{
#' wgcna_out <- run_wgcna_pipeline(dds_qc, var = "condition", save_dir = "results")
#' }
#' @export
run_wgcna_pipeline <- function(
    dds,
    var,
    power = NULL,
    power_vector = c(1:10, seq(12, 20, 2)),
    rsq_cutoff = 0.85,
    network_type = "signed",
    min_module_size = 30,
    merge_cut_height = 0.25,
    max_modules = 15L,
    cor_type = "pearson",
    top_n_hub_genes = 10,
    ora_collection = c("GOBP", "GOCC", "GOMF", "KEGG_MEDICUS", "KEGG_LEGACY", "REACTOME", "PID"),
    top_n_ora = 10,
    gene_col = "gene",
    save_dir = getwd(),
    save_plot = TRUE) {

    if (!requireNamespace("WGCNA", quietly = TRUE)) {
        stop("Package \"WGCNA\" is required for run_wgcna_pipeline() but is not installed.", call. = FALSE)
    }
    if (!requireNamespace("withr", quietly = TRUE)) {
        stop("Package \"withr\" is required for run_wgcna_pipeline() but is not installed.", call. = FALSE)
    }
    dds <- validate_dds(dds)
    dds <- validate_var(dds, var)
    ora_collection <- validate_msigdb_collections(ora_collection)
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)

    vst_mat <- SummarizedExperiment::assay(quiet(safe_vst(dds, blind = TRUE)))
    datExpr <- t(vst_mat)

    message("Checking for genes/samples with excessive missing data...")
    gsg <- quiet(WGCNA::goodSamplesGenes(datExpr, verbose = 0))
    if (!(all(gsg$goodSamples) && all(gsg$goodGenes))) {
        n_bad_genes <- sum(!gsg$goodGenes)
        n_bad_samples <- sum(!gsg$goodSamples)
        message("Removing ", n_bad_genes, " gene(s) and ", n_bad_samples,
                " sample(s) flagged by WGCNA::goodSamplesGenes().")
        datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    }

    if (is.null(power)) {
        message("Selecting soft-thresholding power...")
        power_pick <- pick_wgcna_power(datExpr, power_vector = power_vector,
                                        rsq_cutoff = rsq_cutoff, network_type = network_type)
        power <- power_pick$power
        power_source <- power_pick$power_source
        fit_indices <- power_pick$fit_indices
    } else {
        power_source <- "user"
        fit_indices <- quiet(WGCNA::pickSoftThreshold(
            datExpr, powerVector = sort(unique(c(power, power_vector))),
            networkType = network_type, RsquaredCut = rsq_cutoff, verbose = 0
        ))$fitIndices
    }
    message("Using soft-thresholding power = ", power, " (", power_source, ").")

    message("Detecting co-expression modules (this can take a while)...")
    # WGCNA::blockwiseModules() internally calls a bare `cor()` (not
    # `stats::cor()`) expecting to resolve WGCNA's own `cor()` override
    # (which accepts the extra `weights.x`/`weights.y`/`cosine` arguments
    # blockwiseModules() passes). That override is only found via the
    # search path, i.e. only if the WGCNA namespace is *attached*
    # (library(WGCNA)), not merely loaded via `WGCNA::` qualification --
    # so this call must run with WGCNA temporarily attached.
    net <- withr::with_package("WGCNA", quiet(WGCNA::blockwiseModules(
        datExpr, power = power, networkType = network_type,
        TOMType = if (network_type == "unsigned") "unsigned" else "signed",
        minModuleSize = min_module_size, mergeCutHeight = merge_cut_height,
        corType = cor_type, numericLabels = TRUE, verbose = 0
    )))

    # Prefer WGCNA's natural module count when it is already within
    # max_modules; otherwise merge closest modules until the cap is met.
    n_mod <- length(unique(net$colors[net$colors != 0]))
    if (!is.null(max_modules)) {
        max_modules <- as.integer(max_modules)[1]
        if (is.na(max_modules) || max_modules < 1) {
            stop("`max_modules` must be a positive integer or NULL.", call. = FALSE)
        }
        if (n_mod > max_modules) {
            message("Detected ", n_mod, " modules; merging closest until <= ",
                    max_modules, " ...")
            colors <- net$colors
            cut_h <- merge_cut_height
            merged <- NULL
            while (length(unique(colors[colors != 0])) > max_modules && cut_h < 0.95) {
                cut_h <- min(cut_h + 0.05, 0.95)
                merged <- withr::with_package("WGCNA", quiet(WGCNA::mergeCloseModules(
                    datExpr, colors, cutHeight = cut_h, verbose = 0
                )))
                colors <- merged$colors
            }
            n_after <- length(unique(colors[colors != 0]))
            net$colors <- colors
            if (!is.null(merged)) {
                net$MEs <- merged$newMEs
            }
            if (n_after > max_modules) {
                warning("Could not reduce modules to <= ", max_modules,
                        " (still ", n_after, " after cutHeight=", cut_h,
                        "); keeping the merged set. The HTML report will still ",
                        "cap displayed modules at ", max_modules, ".",
                        call. = FALSE)
            } else {
                message("Merged to ", n_after, " modules (cutHeight=", cut_h, ").")
            }
        } else {
            message("Keeping natural module count (", n_mod,
                    " <= max_modules=", max_modules, ").")
        }
    }

    module_colors <- WGCNA::labels2colors(net$colors)

    # WGCNA's own module labels are arbitrary colour names tied to its
    # internal palette order, not to anything about the module itself.
    # Rename them to "Module 1", "Module 2", ... ranked by descending gene
    # count, so identity is stable and readable across runs/organisms.
    # "grey" -- WGCNA's catch-all for genes it could not assign to any
    # real co-expression module -- is excluded from the numbering entirely
    # and relabelled "Unassigned", since it isn't a real module.
    module_sizes_desc <- sort(table(module_colors[module_colors != "grey"]), decreasing = TRUE)
    module_name_map <- stats::setNames(paste0("Module ", seq_along(module_sizes_desc)), names(module_sizes_desc))
    if ("grey" %in% module_colors) module_name_map <- c(module_name_map, c(grey = "Unassigned"))
    module_names <- unname(module_name_map[module_colors])

    module_df <- data.frame(
        gene = colnames(datExpr),
        module_label = net$colors,
        module_color = module_colors,
        module = module_names,
        stringsAsFactors = FALSE
    )
    colnames(module_df)[1] <- gene_col

    # `net$MEs` columns are named by numeric module *label* ("ME0", "ME1",
    # ...), not by the module *name* used everywhere else in this
    # function's output (module_df$module, the ORA results, plot titles).
    # Resolve each numeric label to its colour and then to its "Module N"
    # / "Unassigned" name via the correspondences already computed above,
    # and rename MEs columns to those names so every output object in this
    # function is keyed the same way, matching module_df/ora_by_module.
    numeric_labels <- sub("^ME", "", colnames(net$MEs))
    label_to_color <- stats::setNames(module_colors, net$colors)
    me_colors <- unname(label_to_color[numeric_labels])
    me_names <- unname(module_name_map[me_colors])

    MEs <- net$MEs
    colnames(MEs) <- me_names
    MEs_out <- cbind(sample = rownames(datExpr), MEs)
    write_output(MEs_out, file.path(save_dir, "wgcna_module_eigengenes.csv"))

    message("Correlating module eigengenes against `", var, "`...")
    trait <- SummarizedExperiment::colData(dds)[rownames(datExpr), var, drop = TRUE]
    trait_mat <- stats::model.matrix(~ 0 + factor(trait))
    colnames(trait_mat) <- levels(factor(trait))

    module_trait_cor <- do.call(rbind, lapply(colnames(MEs), function(mod) {
        do.call(rbind, lapply(colnames(trait_mat), function(lvl) {
            r <- stats::cor(MEs[[mod]], trait_mat[, lvl])
            p <- WGCNA::corPvalueStudent(r, nrow(MEs))
            data.frame(module = mod, level = lvl, cor = r, pvalue = p, stringsAsFactors = FALSE)
        }))
    }))
    write_output(module_trait_cor, file.path(save_dir, "wgcna_module_trait_correlation.csv"))

    message("Identifying hub genes per module...")
    # WGCNA::signedKME() strips a hardcoded "ME" prefix off `datME`'s
    # column names to build its own "kME<module>" output names; since
    # `MEs`'s columns were already de-prefixed above, restore the "ME"
    # prefix just for this call so signedKME's own stripping recovers the
    # correct module names (rather than dropping the first two characters
    # of a name that no longer starts with "ME").
    kME <- WGCNA::signedKME(datExpr, stats::setNames(MEs, paste0("ME", colnames(MEs))))
    colnames(kME) <- sub("^kME", "", colnames(kME))
    # Attach each gene's kME to its own module (membership strength).
    module_df$kME <- vapply(seq_len(nrow(module_df)), function(i) {
        mod <- module_df$module[i]
        if (!mod %in% colnames(kME)) return(NA_real_)
        kME[[mod]][i]
    }, numeric(1))
    write_output(module_df, file.path(save_dir, "wgcna_modules.csv"))

    safe_file_tag <- function(mod) gsub("[^A-Za-z0-9]+", "_", mod)
    hub_genes <- stats::setNames(lapply(colnames(MEs), function(mod) {
        if (mod == "Unassigned") {
            empty <- data.frame(gene = character(0), kME = numeric(0), stringsAsFactors = FALSE)
            colnames(empty)[1] <- gene_col
            return(empty)
        }
        ord <- order(-abs(kME[[mod]]))
        top <- utils::head(ord, top_n_hub_genes)
        df <- data.frame(gene = colnames(datExpr)[top], kME = kME[[mod]][top], stringsAsFactors = FALSE)
        colnames(df)[1] <- gene_col
        write_output(df, file.path(save_dir, paste0("wgcna_hub_genes_", safe_file_tag(mod), ".csv")))
        df
    }), colnames(MEs))

    # "Unassigned" (WGCNA's "grey" label) groups genes WGCNA could not
    # assign to any real co-expression module -- a heterogeneous leftover
    # bucket, not a module with a coherent expression profile -- so
    # pathway over-representation on it would not be biologically
    # meaningful and it is excluded from every collection below.
    real_modules <- setdiff(unique(module_df$module), "Unassigned")
    universe <- colnames(datExpr)
    org_for_ora <- detect_organism(universe)

    message("Running per-module over-representation analysis against ",
            length(ora_collection), " collection(s): ", paste(ora_collection, collapse = ", "), " ...")
    ora_by_module <- stats::setNames(lapply(ora_collection, function(coll) {
        pathways <- import_msigdbr(org_for_ora, coll)
        stats::setNames(lapply(real_modules, function(mod) {
            mod_genes <- module_df[[gene_col]][module_df$module == mod]
            res <- as.data.frame(fgsea::fora(pathways, genes = mod_genes, universe = universe))
            res <- res[order(res$padj), ]
            res$overlapGenes <- vapply(res$overlapGenes, paste, character(1), collapse = ",")
            write_output(res, file.path(save_dir, paste0("wgcna_ora_", coll, "_", safe_file_tag(mod), ".csv")))
            res
        }), real_modules)
    }), ora_collection)

    p_power <- plot_soft_threshold(fit_indices, power, rsq_cutoff = rsq_cutoff)
    p_size <- plot_module_sizes(module_df)
    p_heat <- plot_module_trait_heatmap(module_trait_cor)
    p_me_heat <- plot_module_eigengene_heatmap(MEs)
    p_ora_by_module <- stats::setNames(lapply(ora_collection, function(coll) {
        stats::setNames(lapply(real_modules, function(mod) {
            res <- ora_by_module[[coll]][[mod]]
            if (nrow(res) == 0) return(NULL)
            plot_module_ora_dotplot(res, top_n = top_n_ora, title = paste0("Module: ", mod))
        }), real_modules)
    }), ora_collection)

    if (save_plot) {
        save_ggplot(p_power, file.path(save_dir, "wgcna_soft_threshold.pdf"), width = 10, height = 5)

        save_ggplot(p_size, file.path(save_dir, "wgcna_module_sizes.pdf"),
                    width = 7, height = max(4, 1 + length(unique(module_df$module)) * 0.35))

        grDevices::pdf(file.path(save_dir, "wgcna_dendrogram.pdf"), width = 10, height = 6)
        WGCNA::plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                                    groupLabels = "Module", dendroLabels = FALSE, addGuide = TRUE)
        grDevices::dev.off()

        save_ggplot(p_heat, file.path(save_dir, "wgcna_module_trait_heatmap.pdf"),
                    width = max(6, 1 + ncol(trait_mat) * 1.2), height = max(4, 1 + length(unique(module_trait_cor$module)) * 0.4))

        me_heat_w <- max(6, 1 + nrow(MEs) * 0.5)
        me_heat_h <- max(4, 1 + ncol(MEs) * 0.4)
        grDevices::pdf(file.path(save_dir, "wgcna_module_eigengene_heatmap.pdf"), width = me_heat_w, height = me_heat_h)
        grid::grid.newpage()
        grid::grid.draw(p_me_heat$gtable)
        grDevices::dev.off()

        for (coll in names(p_ora_by_module)) {
            for (mod in names(p_ora_by_module[[coll]])) {
                p <- p_ora_by_module[[coll]][[mod]]
                if (is.null(p)) next
                save_ggplot(p, file.path(save_dir, paste0("wgcna_ora_dotplot_", coll, "_", safe_file_tag(mod), ".pdf")),
                            width = 9, height = 4.8)
            }
        }
    }

    list(
        power = power, power_source = power_source, fit_indices = fit_indices,
        module_df = module_df, MEs = MEs_out, module_trait_cor = module_trait_cor,
        hub_genes = hub_genes, ora_by_module = ora_by_module, net = net,
        soft_threshold_plot = p_power, module_size_plot = p_size,
        module_trait_heatmap = p_heat, module_eigengene_heatmap = p_me_heat,
        ora_plots_by_module = p_ora_by_module
    )
}

#' Plot a module-eigengene vs. trait-level correlation heatmap
#'
#' @param module_trait_cor A data frame as returned in
#'   [run_wgcna_pipeline()]'s output (columns `module`, `level`, `cor`,
#'   `pvalue`).
#' @return A `ggplot` object: tile heatmap of correlation, coloured on a
#'   red-blue diverging scale, with `r` (top) / `p` (bottom) annotated in
#'   each cell.
#' @export
plot_module_trait_heatmap <- function(module_trait_cor) {
    df <- module_trait_cor
    df$label <- sprintf("%.2f\n(%.1g)", df$cor, df$pvalue)
    # Intentionally NOT squared (no square_panel_theme()) -- like the
    # other WGCNA plots, this heatmap keeps its natural panel shape:
    # with discrete x/y and geom_tile(), tile width/height is uniform and
    # follows however many trait levels (columns) vs. modules (rows)
    # there are, rather than being forced into a square overall panel.
    ggplot2::ggplot(df, ggplot2::aes(x = .data$level, y = .data$module, fill = .data$cor)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.8) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
        ggplot2::labs(x = NULL, y = "Module", title = "Module-trait relationships") +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Plot a sample-level module-eigengene heatmap
#'
#' A `pheatmap` with one row per module (module eigengene, "ME") and one
#' column per sample, cell value = that module's eigengene value for that
#' sample, with both rows and columns hierarchically clustered. This is
#' the sample-resolution companion to [plot_module_trait_heatmap()]'s
#' condition-level view, styled the same way as the existing
#' sample-similarity heatmap (diverging blue-white-red colour scale,
#' scaled to the actual data range, clustered on both axes).
#'
#' @param MEs Module eigengene data frame (samples x modules), as returned
#'   in [run_wgcna_pipeline()]'s `MEs` output (column names are module
#'   names, e.g. `"Module 1"`, `"Unassigned"`; row names/order are
#'   samples).
#' @return A `pheatmap` object (a `list` with a `$gtable` element,
#'   printable directly and usable with [download_grid_plot_button()]).
#' @export
plot_module_eigengene_heatmap <- function(MEs) {
    me_mat <- t(as.matrix(MEs))

    can_cluster_rows <- nrow(me_mat) > 1 && !anyNA(me_mat)
    can_cluster_cols <- ncol(me_mat) > 1 && !anyNA(me_mat)

    colours <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
    lo <- min(me_mat, na.rm = TRUE)
    hi <- max(me_mat, na.rm = TRUE)
    if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
        lo <- (if (is.finite(lo)) lo else -1) - 0.01
        hi <- (if (is.finite(hi)) hi else 1) + 0.01
    }
    # Centre the scale on zero (eigengenes, like correlations, are signed
    # and zero is a meaningful reference point) while still spanning the
    # actual data range.
    bound <- max(abs(lo), abs(hi))
    breaks <- seq(-bound, bound, length.out = length(colours) + 1)

    # Intentionally NOT squared -- unlike the sample-similarity/correlation
    # heatmaps, this panel is left at its natural size: cellwidth/cellheight
    # are sized independently off each dimension's own count (modules vs.
    # samples, typically very lopsided), so the overall data-matrix panel
    # is rectangular whenever nrow != ncol, matching the "no forced square
    # aspect ratio for any WGCNA plot" convention.
    cell_w <- pheatmap_cell_size(ncol(me_mat))
    cell_h <- pheatmap_cell_size(nrow(me_mat))

    pheatmap::pheatmap(
        me_mat,
        cluster_rows = can_cluster_rows,
        cluster_cols = can_cluster_cols,
        color = colours,
        breaks = breaks,
        display_numbers = matrix(sprintf("%.2f", me_mat), nrow = nrow(me_mat), dimnames = dimnames(me_mat)),
        number_color = "black",
        fontsize_number = 6,
        cellwidth = cell_w,
        cellheight = cell_h,
        main = "Module eigengenes by sample",
        silent = TRUE
    )
}

#' Plot module sizes (gene counts per module)
#'
#' @param module_df A data frame as returned in [run_wgcna_pipeline()]'s
#'   `module_df` output (must contain `module` and `module_color`
#'   columns -- `module` is the "Module N" / "Unassigned" display name,
#'   `module_color` is WGCNA's original colour label used to fill each
#'   bar so the plot stays visually tied to the underlying co-expression
#'   colour even though modules are labelled/ordered by name).
#' @return A `ggplot` object: horizontal barplot of gene counts per
#'   module, ordered by descending size ("Unassigned" always last), bars
#'   filled with their original WGCNA module colour.
#' @export
plot_module_sizes <- function(module_df) {
    sizes <- stats::aggregate(list(n_genes = module_df$module), by = list(module = module_df$module, module_color = module_df$module_color), FUN = length)
    sizes <- sizes[order(sizes$module == "Unassigned", -sizes$n_genes), ]
    sizes$module <- factor(sizes$module, levels = rev(sizes$module))

    # The gene-count label is drawn just past the tip of each bar (`hjust =
    # -0.2`); with zero scale expansion that headroom has to come from the
    # explicit x-limit itself instead, so the label on the longest bar isn't
    # clipped by the (now flush-to-data) panel edge.
    x_max <- max(sizes$n_genes, na.rm = TRUE)

    ggplot2::ggplot(sizes, ggplot2::aes(x = .data$n_genes, y = .data$module, fill = .data$module_color)) +
        ggplot2::geom_col(colour = "black", linewidth = 0.2) +
        ggplot2::geom_text(ggplot2::aes(label = .data$n_genes), hjust = -0.2, size = 3) +
        ggplot2::scale_fill_identity() +
        # y (module) is discrete: keep ggplot2's minimal built-in discrete
        # expansion so the outermost bar isn't sliced by the panel border.
        ggplot2::scale_x_continuous(limits = c(0, x_max * 1.15), expand = zero_expansion()) +
        # The outlined bars (colour = "black" above) have their x = 0 edge
        # sitting exactly on the panel border under zero data-padding;
        # without disabling clipping, half of that outline stroke would be
        # sliced off there.
        unclipped_coord() +
        ggplot2::labs(x = "Number of genes", y = "Module", title = "Module sizes") +
        # Intentionally NOT squared (no square_panel_theme()) -- matches
        # the other WGCNA plots' natural (unforced) panel shape.
        ggplot2::theme_bw()
}

#' Plot a per-module over-representation analysis dot plot
#'
#' @param ora_res A data frame as returned by `fgsea::fora()` (columns
#'   `pathway`, `padj`, `overlap`, `size`), as stored in
#'   [run_wgcna_pipeline()]'s `ora_by_module[[collection]][[module]]`
#'   output.
#' @param top_n Number of top pathways (by adjusted p-value) to plot.
#'   Default `10`.
#' @param title Optional plot title (e.g. the module name). Default `NULL`.
#' @return A `ggplot` object: dot plot with `-log10(padj)` on the x-axis,
#'   pathway name on the y-axis (most to least significant, top to bottom),
#'   and point size mapped to gene-set ratio (`overlap / size`).
#' @export
plot_module_ora_dotplot <- function(ora_res, top_n = 10, title = NULL) {
    df <- ora_res[order(ora_res$padj), ]
    df <- utils::head(df, top_n)
    df$pathway_label <- as.character(df$pathway)
    for (spec in .msigdb_collection_map) {
        df$pathway_label <- sub(paste0("^", spec$prefix), "", df$pathway_label)
    }
    # Most significant (smallest padj / highest -log10) at the top.
    df$pathway_label <- factor(df$pathway_label, levels = rev(df$pathway_label))
    df$neg_log10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))
    df$gene_ratio <- df$overlap / pmax(df$size, 1)

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$neg_log10_padj, y = .data$pathway_label, size = .data$gene_ratio)
    ) +
        ggplot2::geom_point(colour = "#B2182B", alpha = 0.85) +
        ggplot2::scale_size_continuous(name = "Gene ratio", range = c(2.5, 9)) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
        unclipped_coord() +
        ggplot2::labs(x = expression(-log[10](p[adj])), y = NULL, title = title) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    p
}
