#' Regex matching BCR and/or TCR VDJ gene symbols (human HGNC + mouse MGI)
#' @param kinds Character vector containing `"BCR"` and/or `"TCR"`.
#' @return A single alternation regex suitable for `grepl()` on gene symbols.
#' @noRd
.vdj_gene_pattern <- function(kinds) {
    pats <- character(0)
    if ("BCR" %in% kinds) {
        pats <- c(
            pats,
            # mouse (title-case) immunoglobulin loci + Jchain
            "^Igh[vdj]", "^Igh[amge]$", "^Ighg[0-9]", "^Igha",
            "^Igk[vdjc]", "^Igl[vdjc]", "^Jchain$",
            # human (UPPER) immunoglobulin loci + JCHAIN
            "^IGH[VDJ]", "^IGH[AMGE]$", "^IGHG[0-9]", "^IGHA",
            "^IGK[VDJC]", "^IGL[VDJC]", "^JCHAIN$"
        )
    }
    if ("TCR" %in% kinds) {
        pats <- c(
            pats,
            "^Tr[abgd][vdjc]",  # mouse Trav/Trbv/...
            "^TR[ABGD][VDJC]"   # human TRAV/TRBV/...
        )
    }
    if (!length(pats)) {
        stop("`.vdj_gene_pattern()` requires at least one of \"BCR\", \"TCR\".", call. = FALSE)
    }
    paste(pats, collapse = "|")
}

#' Run the quality-control stage of the pipeline
#'
#' Filters lowly-expressed genes, imputes each sample's sex from the
#' sex-chromosome marker genes ([impute_sample_sex()], stored as
#' `colData(dds)$imputed_sex`), then optionally removes mitochondrial,
#' sex-chromosome and VDJ genes -- in that order, so sex is imputed before
#' its marker genes can be removed -- and writes a library-size
#' distribution plot, a gene-filtering diagnostic plot, and a filtering
#' summary CSV. Returns the filtered `dds` object so it can be fed into
#' [run_distance_pca()] and [run_diffexp_pipeline()].
#'
#' Low-expression filtering is on by default (`filter_method = "cpm"`); pass
#' `filter_method = NULL` to remove no genes for expression. Two filtering
#' strategies are available via `filter_method`:
#' \describe{
#'   \item{`"cpm"`}{drop genes with no reads in any sample, then rank what
#'     remains by mean counts-per-million (CPM) across all samples and remove
#'     the bottom `cpm_percentile` fraction of *those* (default `0.05`). This
#'     adapts automatically to sequencing depth (via CPM) without requiring an
#'     absolute count threshold to be chosen up front. The unexpressed genes
#'     are taken out first because the percentile is otherwise meaningless: a
#'     full annotation typically leaves half the genes zero everywhere, so the
#'     quantile lands inside that tie and the filter removes nothing at all.}
#'   \item{`"count"`}{a group-size-aware absolute threshold: keep a gene if,
#'     in at least one level of `var`, at least `min_samples_per_group`
#'     samples of that group have a raw count >= `min_count` (default `10`).
#'     `min_samples_per_group` defaults to the size of the smallest group in
#'     `var`. This is similar in spirit to `edgeR::filterByExpr()` but
#'     operates on raw counts rather than CPM.}
#' }
#' In both cases, [plot_count_distribution()] renders the per-gene
#' expression distribution actually used for the decision (log2 mean CPM,
#' or log2 total count) with the effective cutoff overlaid, so the chosen
#' threshold can be checked against the real data rather than trusted
#' blindly, and is displayed in the HTML report *before* the low-expression
#' filtering step is applied.
#'
#' @param dds A `DESeqDataSet` object with raw counts.
#' @param var Column name in `colData(dds)` used for grouping (used to
#'   validate sample groups, to determine the default
#'   `min_samples_per_group` for `filter_method = "count"`, and to label the
#'   library-size plot).
#' @param filter_method Low-expression filtering strategy: `"cpm"` (default,
#'   percentile-based), `"count"`, or `NULL` (no expression filtering).
#'   See Details.
#' @param cpm_percentile Fraction (0-1) of the *expressed* genes (by mean CPM
#'   across samples) to remove, on top of the genes with no reads at all. Used
#'   only when `filter_method = "cpm"`. Default `0.05` (bottom 5%); `0` removes
#'   only the unexpressed genes.
#' @param min_count Minimum raw count for a sample to count as "expressing"
#'   a gene. Used only when `filter_method = "count"`. Default `10`.
#' @param min_samples_per_group Minimum number of samples *within a single
#'   group* that must clear `min_count` for a gene to be retained. Used only
#'   when `filter_method = "count"`. Default `NULL` uses the size of the
#'   smallest group in `var`.
#' @param remove_mt Logical. If `TRUE`, removes genes matching `mt_pattern`
#'   against `rownames(dds)`. Default `FALSE`.
#' @param mt_pattern Regular expression used to identify mitochondrial genes
#'   by name when `remove_mt = TRUE`. Default matches common human/mouse
#'   mitochondrial gene symbol prefixes (`"^MT-"` / `"^mt-"`).
#' @param remove_xy Logical. If `TRUE`, removes genes on X/Y chromosomes,
#'   using `gene_chr` to determine chromosome membership. Default `FALSE`.
#' @param gene_chr Optional data frame with columns `gene` and `chr` mapping
#'   `rownames(dds)` to chromosome names (e.g. from `biomaRt` or a GTF).
#'   Required if `remove_xy = TRUE`; ignored otherwise.
#' @param remove_vdj VDJ gene removal: `NULL`/`FALSE` (default, off),
#'   `TRUE`/`"both"` (BCR + TCR), or one or more of `"BCR"` / `"TCR"`.
#'   Matches immunoglobulin / T-cell receptor gene-symbol conventions for
#'   human (HGNC) and mouse (MGI) against `rownames(dds)`.
#' @param save_dir Directory to write outputs to. Default `getwd()`.
#' @param save_plot Logical. If `TRUE`, saves the library-size boxplot and
#'   the count-distribution diagnostic plot as PDFs in `save_dir`. Default
#'   `TRUE`.
#' @return The filtered `DESeqDataSet` object. The number/percent of genes
#'   removed at each step, the count-distribution and library-size plots,
#'   and the filtering parameters used are attached as
#'   `S4Vectors::metadata(dds)$degreports_qc` so that [render_report()]
#'   can display them without recomputation; also written to
#'   `qc_steps_summary.csv` in `save_dir`.
#' @examples
#' \dontrun{
#' dds <- run_qc(dds, var = "condition")  # CPM percentile filter (default)
#' dds <- run_qc(dds, var = "condition", filter_method = "count",
#'               min_count = 10, remove_mt = TRUE)
#' dds <- run_qc(dds, var = "condition", filter_method = NULL)  # no filter
#' }
#' @export
run_qc <- function(
    dds,
    var,
    filter_method = "cpm",
    cpm_percentile = 0.05,
    min_count = 10,
    min_samples_per_group = NULL,
    remove_mt = FALSE,
    mt_pattern = "^MT-|^mt-",
    remove_xy = FALSE,
    gene_chr = NULL,
    remove_vdj = NULL,
    save_dir = getwd(),
    save_plot = TRUE) {

    dds <- validate_dds(dds)
    dds <- validate_var(dds, var)
    if (!is.null(filter_method)) {
        filter_method <- validate_choice(filter_method, c("cpm", "count"), "filter_method")
    }
    min_count <- validate_min_count(min_count)
    validate_logical(remove_mt, "remove_mt")
    validate_logical(remove_xy, "remove_xy")
    vdj_kinds <- validate_remove_vdj(remove_vdj)
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)

    if (!is.null(filter_method) && filter_method == "cpm") {
        if (!is.numeric(cpm_percentile) || length(cpm_percentile) != 1 || is.na(cpm_percentile) ||
                cpm_percentile < 0 || cpm_percentile >= 1) {
            stop("`cpm_percentile` must be a single number in [0, 1).", call. = FALSE)
        }
    }
    if (!is.null(filter_method) && filter_method == "count") {
        if (is.null(min_samples_per_group)) {
            min_samples_per_group <- min(table(SummarizedExperiment::colData(dds)[[var]]))
        } else if (!is.numeric(min_samples_per_group) || length(min_samples_per_group) != 1 ||
                min_samples_per_group < 1) {
            stop("`min_samples_per_group` must be `NULL` or a single positive number.", call. = FALSE)
        }
    }

    n_genes_start <- nrow(dds)
    step_names <- "Input"
    step_n_before <- n_genes_start
    step_n_after <- n_genes_start

    # Per-sample detection / depth on the input count matrix (before gene
    # filtering). Reported under Sample QC ahead of Library QC.
    counts_input <- .raw_count_matrix(dds)
    sample_qc_table <- data.frame(
        sample = colnames(counts_input),
        n_detected_genes = as.integer(colSums(counts_input > 0)),
        n_total_transcripts = as.numeric(colSums(counts_input)),
        stringsAsFactors = FALSE
    )

    raw_counts <- .raw_count_matrix(dds)
    n_before_expr_filter <- nrow(dds)

    # Always compute the mean-CPM distribution for the diagnostic plot, even
    # when `filter_method = NULL`, so the report can show "what would the
    # distribution look like" prior to any filtering decision.
    lib_sizes <- colSums(raw_counts)
    cpm_mat <- sweep(raw_counts, 2, lib_sizes, "/") * 1e6
    gene_stat <- rowMeans(cpm_mat)
    gene_stat_log2 <- log2(gene_stat + 1)
    xlab <- "log2(mean CPM across samples + 1)"
    cutoff_log2 <- NA_real_
    keep_expr <- rep(TRUE, nrow(dds))
    filter_desc <- "none"
    filter_step_name <- NULL

    if (is.null(filter_method)) {
        message("No low-expression filtering applied (filter_method = NULL).")
    } else if (filter_method == "cpm") {
        # Genes with no reads at all go first, then the percentile applies to
        # what is left. Taking the percentile over the whole distribution makes
        # it meaningless on typical bulk RNA-seq, where a full annotation leaves
        # half the genes zero in every sample: the quantile then lands inside
        # that tie at zero, and since log2(0 + 1) == 0, a `>= cutoff` rule keeps
        # every gene in the dataset -- the filter silently removed nothing.
        expressed <- gene_stat > 0
        n_zero <- sum(!expressed)
        if (any(expressed) && cpm_percentile > 0) {
            cutoff_log2 <- stats::quantile(gene_stat_log2[expressed],
                                           probs = cpm_percentile,
                                           names = FALSE, type = 7)
            keep_expr <- expressed & gene_stat_log2 > cutoff_log2
        } else {
            cutoff_log2 <- 0
            keep_expr <- expressed
        }
        message("Filtering low-expression genes: ", n_zero,
                " gene(s) with no reads in any sample, plus the bottom ",
                round(cpm_percentile * 100, 2), "% of the ", sum(expressed),
                " expressed gene(s) (below log2(CPM+1) = ", round(cutoff_log2, 4),
                ") -- ", sum(!keep_expr), " gene(s) removed in total...")
        filter_desc <- paste0("unexpressed genes plus the bottom ",
                              round(cpm_percentile * 100, 2),
                              "% of expressed genes by mean CPM removed")
        filter_step_name <- "Low-expression gene removal (CPM percentile)"
    } else {
        groups <- SummarizedExperiment::colData(dds)[[var]]
        group_levels <- levels(groups)
        count_ok <- raw_counts >= min_count
        per_group_n <- vapply(group_levels, function(lv) {
            cols <- which(groups == lv)
            if (length(cols) == 1) as.integer(count_ok[, cols]) else rowSums(count_ok[, cols, drop = FALSE])
        }, numeric(nrow(raw_counts)))
        keep_expr <- apply(per_group_n >= min_samples_per_group, 1, any)
        # The actual per-gene decision is per-group (>= min_samples_per_group
        # samples within *some* group must clear min_count), which cannot be
        # reduced to a single point on the total-count axis plotted here; the
        # line marks min_count itself as a visual reference for where
        # "expressed" begins on a per-sample basis, not an exact boundary.
        total_count_log2 <- log2(rowSums(raw_counts) + 1)
        gene_stat_log2 <- total_count_log2
        xlab <- "log2(total count across all samples + 1)"
        cutoff_log2 <- log2(min_count + 1)
        message("Filtering genes with count < ", min_count, " in fewer than ", min_samples_per_group,
                " sample(s) within every group (", sum(!keep_expr), " gene(s) removed)...")
        filter_desc <- paste0("count >= ", min_count, " in >= ", min_samples_per_group,
                               " sample(s) of at least one group")
        filter_step_name <- "Low-expression gene removal (raw count, per-group)"
    }

    count_dist_plot <- plot_count_distribution(gene_stat_log2, cutoff_log2, keep_expr, xlab)

    dds <- dds[keep_expr, ]
    if (!is.null(filter_step_name)) {
        step_names <- c(step_names, filter_step_name)
        step_n_before <- c(step_n_before, n_before_expr_filter)
        step_n_after <- c(step_n_after, nrow(dds))
    }

    # Sex is imputed on every run: after the expression filter (so the marker
    # genes are on the same footing as everything else) and before X/Y genes
    # can be removed (the markers live there).
    sex_out <- impute_sample_sex(dds)
    message(sex_out$message)
    SummarizedExperiment::colData(dds)[[.sex_col]] <- factor(
        sex_out$table$imputed_sex, levels = c("female", "male"))
    write_output(sex_out$table, file.path(save_dir, "qc_imputed_sex.csv"), format = "csv")
    sex_plot <- plot_sex_markers(sex_out$table)
    if (save_plot && !is.null(sex_plot)) {
        save_ggplot(sex_plot, file.path(save_dir, "qc_imputed_sex.pdf"), width = 4, height = 3.5)
    }

    if (remove_mt) {
        n_before <- nrow(dds)
        keep_mt <- !grepl(mt_pattern, rownames(dds))
        message("Removing ", sum(!keep_mt), " mitochondrial gene(s) (pattern: \"", mt_pattern, "\")...")
        dds <- dds[keep_mt, ]
        step_names <- c(step_names, "Mitochondrial gene removal")
        step_n_before <- c(step_n_before, n_before)
        step_n_after <- c(step_n_after, nrow(dds))
    }

    if (remove_xy) {
        if (is.null(gene_chr)) {
            warning("`remove_xy = TRUE` but `gene_chr` was not supplied; skipping X/Y gene removal. ",
                    "Provide a data.frame(gene, chr) via `gene_chr` to enable this step.", call. = FALSE)
        } else {
            if (!all(c("gene", "chr") %in% colnames(gene_chr))) {
                stop("`gene_chr` must be a data.frame with columns `gene` and `chr`.", call. = FALSE)
            }
            n_before <- nrow(dds)
            xy_genes <- gene_chr$gene[gene_chr$chr %in% c("X", "Y", "chrX", "chrY")]
            keep_xy <- !rownames(dds) %in% xy_genes
            message("Removing ", sum(!keep_xy), " X/Y chromosome gene(s)...")
            dds <- dds[keep_xy, ]
            step_names <- c(step_names, "X/Y chromosome gene removal")
            step_n_before <- c(step_n_before, n_before)
            step_n_after <- c(step_n_after, nrow(dds))
        }
    }

    if (length(vdj_kinds) > 0) {
        n_before <- nrow(dds)
        vdj_pat <- .vdj_gene_pattern(vdj_kinds)
        keep_vdj <- !grepl(vdj_pat, rownames(dds))
        message("Removing ", sum(!keep_vdj), " VDJ gene(s) (",
                paste(vdj_kinds, collapse = "+"), ")...")
        dds <- dds[keep_vdj, ]
        step_names <- c(step_names, paste0("VDJ gene removal (", paste(vdj_kinds, collapse = "+"), ")"))
        step_n_before <- c(step_n_before, n_before)
        step_n_after <- c(step_n_after, nrow(dds))
    }

    message(nrow(dds), " / ", n_genes_start, " genes retained after QC filtering.")

    step_n_removed <- step_n_before - step_n_after
    step_pct_removed <- ifelse(step_n_before > 0, 100 * step_n_removed / step_n_before, NA_real_)
    steps_summary <- data.frame(
        step = step_names,
        n_before = step_n_before,
        n_removed = step_n_removed,
        pct_removed = round(step_pct_removed, 2),
        n_after = step_n_after
    )
    write_output(steps_summary, file.path(save_dir, "qc_steps_summary.csv"), format = "csv")
    write_output(sample_qc_table, file.path(save_dir, "qc_sample_summary.csv"), format = "csv")

    filter_summary <- data.frame(
        metric = c("genes_input", "genes_after_qc", "filter_method", "filter_criterion", "samples"),
        value = c(n_genes_start, nrow(dds), if (is.null(filter_method)) "none" else filter_method,
                  filter_desc, ncol(dds))
    )
    write_output(filter_summary, file.path(save_dir, "qc_filter_summary.csv"), format = "csv")

    p <- plot_library_size(dds)
    # Half the previous default height; width still scales with sample count.
    lib_w <- max(2.5, 2 + 0.15 * ncol(dds))
    lib_h <- 2.5
    if (save_plot) {
        save_ggplot(p, file.path(save_dir, "qc_library_size.pdf"),
                    width = lib_w, height = lib_h)
        save_ggplot(count_dist_plot, file.path(save_dir, "qc_count_distribution.pdf"),
                    width = 3, height = 2.5)
    }

    # Attach everything render_report() needs to build a QC section, so the
    # report doesn't need to recompute (or be separately passed) QC detail.
    qc_meta <- list(
        steps_summary = steps_summary,
        sample_qc_table = sample_qc_table,
        filter_method = filter_method,
        cpm_percentile = if (!is.null(filter_method) && filter_method == "cpm") cpm_percentile else NA_real_,
        min_samples_per_group = if (!is.null(filter_method) && filter_method == "count") min_samples_per_group else NA_integer_,
        min_count = min_count,
        remove_vdj = if (length(vdj_kinds)) vdj_kinds else NULL,
        n_samples = ncol(dds),
        count_dist_plot = count_dist_plot,
        library_size_plot = p,
        sex = sex_out,
        sex_plot = sex_plot
    )
    S4Vectors::metadata(dds)$degreports_qc <- qc_meta

    dds
}

#' Plot the per-gene expression distribution before/after QC filtering
#'
#' A single overlaid density plot of a per-gene summary statistic (log2 mean
#' CPM, or log2 total count): one density for *all* genes ("before
#' filtering") and one for the subset retained after low-expression removal
#' ("after filtering"), distinguished by colour on the same axes so the
#' shift in the distribution is directly visible. The cutoff is overlaid as
#' a vertical dashed line. For `filter_method = "cpm"` this line is the exact
#' boundary: it is the `cpm_percentile` cutoff taken over the expressed genes,
#' and everything at or below it -- including the unexpressed genes piled up at
#' zero -- is removed. For
#' `filter_method = "count"` the underlying decision is per-group (whether
#' *any* group has `min_samples_per_group` samples clearing `min_count`),
#' which cannot be reduced to a single point on the total-count axis, so the
#' line marks `min_count` itself as a reference for where "expressed"
#' begins on a per-sample basis. When `filter_method = NULL`, `keep` is all
#' `TRUE` (the two densities coincide) and no cutoff line is drawn.
#'
#' @param gene_stat_log2 Numeric vector, one value per gene (log2 scale).
#' @param cutoff_log2 The cutoff value on the same log2 scale, or `NA` to
#'   omit the vertical line.
#' @param keep Logical vector, same length as `gene_stat_log2`: whether
#'   each gene was retained.
#' @param xlab X-axis label describing the statistic plotted.
#' @return A `ggplot` object.
#' @export
plot_count_distribution <- function(gene_stat_log2, cutoff_log2, keep, xlab) {
    df <- rbind(
        data.frame(value = gene_stat_log2, stage = "Before filtering"),
        data.frame(value = gene_stat_log2[keep], stage = "After filtering")
    )
    df$stage <- factor(df$stage, levels = c("Before filtering", "After filtering"))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, fill = .data$stage)) +
        ggplot2::geom_density(alpha = 0.5, colour = NA) +
        ggplot2::scale_x_continuous(expand = zero_expansion()) +
        ggplot2::scale_y_continuous(expand = zero_expansion()) +
        ggplot2::scale_fill_manual(values = c("Before filtering" = "#bdbdbd", "After filtering" = "#3182bd")) +
        ggplot2::theme_bw(base_size = 12) +
        square_panel_theme() +
        ggplot2::xlab(xlab) +
        ggplot2::ylab("Density") +
        ggplot2::labs(fill = NULL)

    if (!is.na(cutoff_log2)) {
        p <- p + ggplot2::geom_vline(xintercept = cutoff_log2, linetype = "dashed", colour = "black")
    }

    p
}

#' Plot library-size (VST expression) distribution across samples
#'
#' @param dds A `DESeqDataSet` object.
#' @return A `ggplot` object: per-sample boxplots of VST-transformed expression.
#' @export
plot_library_size <- function(dds) {
    dds <- validate_dds(dds)
    vsd <- safe_vst(dds, blind = FALSE)
    assay_df <- as.data.frame(SummarizedExperiment::assay(vsd))
    assay_long <- tidyr::pivot_longer(assay_df, cols = dplyr::everything(),
                                       names_to = "sample", values_to = "expression")
    assay_long$sample <- factor(assay_long$sample, levels = colnames(dds))

    ggplot2::ggplot(assay_long, ggplot2::aes(x = .data$sample, y = .data$expression)) +
        ggplot2::geom_boxplot(fill = "grey85", width = 0.7) +
        # x is discrete (per-sample): keep ggplot2's minimal built-in discrete
        # expansion so the outermost boxplot isn't sliced in half by the
        # panel border. y is continuous: zero data-padding.
        ggplot2::scale_y_continuous(expand = zero_expansion()) +
        # Zero y-padding (above) can put a boxplot's whisker end or an
        # outlier point exactly on the panel border; without disabling
        # clipping, its drawn extent would be sliced off there.
        unclipped_coord() +
        ggplot2::theme_bw(base_size = 12) +
        square_panel_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab(if (.is_expr_dds(dds)) "log2 expression" else "VST expression")
}
