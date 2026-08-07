#' Map human gene symbols to the analysis organism (human or mouse)
#'
#' Human input is returned unchanged. Mouse mapping uses
#' `babelgene::orthologs(..., species = "mouse", human = TRUE)`.
#'
#' @param symbols Character vector of human symbols.
#' @param org `"human"` or `"mouse"`.
#' @return Character vector of organism-native symbols (possibly empty).
#' @noRd
.map_human_symbols_to_org <- function(symbols, org) {
    symbols <- unique(as.character(symbols))
    symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
    if (!length(symbols)) return(character(0))
    if (is.null(org) || is.na(org) || identical(org, "human")) {
        return(symbols)
    }
    if (!identical(org, "mouse")) {
        return(symbols)
    }
    if (!requireNamespace("babelgene", quietly = TRUE)) {
        warning("Package \"babelgene\" not installed; cannot map VIPER human TFs ",
                "to mouse for multi-hit scoring.", call. = FALSE)
        return(character(0))
    }
    ortho <- tryCatch(
        babelgene::orthologs(symbols, species = "mouse", human = TRUE),
        error = function(e) {
            warning("babelgene human→mouse mapping failed: ", conditionMessage(e),
                    call. = FALSE)
            NULL
        }
    )
    if (is.null(ortho) || !nrow(ortho)) return(character(0))
    # babelgene typically returns human_symbol + symbol (mouse)
    human_col <- intersect(c("human_symbol", "human_ensembl"), colnames(ortho))
    mouse_col <- intersect(c("symbol", "ensembl"), colnames(ortho))
    if (!length(mouse_col)) return(character(0))
    unique(as.character(ortho[[mouse_col[[1]]]]))
}

#' Detect organism from available pipeline outputs
#' @noRd
.multihit_org <- function(gsea_out = NULL, viper_out = NULL, decoupler_out = NULL,
                          diffexp_out = NULL) {
    org <- tryCatch(gsea_out$org, error = function(e) NULL)
    if (is.null(org) || is.na(org)) {
        org <- tryCatch(gsea_out$results$org[[1]], error = function(e) NULL)
    }
    if (is.null(org) || is.na(org)) {
        org <- tryCatch(viper_out$org, error = function(e) NULL)
    }
    if (is.null(org) || is.na(org)) {
        org <- tryCatch(decoupler_out$org, error = function(e) NULL)
    }
    if (is.null(org) || is.na(org)) {
        genes <- tryCatch({
            df <- diffexp_out$results
            head(df[[intersect(c("gene"), colnames(df))[[1]]]], 200)
        }, error = function(e) character(0))
        if (length(genes)) {
            org <- tryCatch(detect_organism(genes), error = function(e) NA_character_)
        }
    }
    if (is.null(org) || is.na(org)) org <- "human"
    org
}

#' Rank genes by cross-analysis hit frequency for one contrast
#'
#' Scores each gene as:
#' \itemize{
#'   \item DEG: +1 if `sig`
#'   \item GSEA: +1 if the gene appears in **any** significant pathway
#'     leading edge (`padj < alpha`; all collections pooled)
#'   \item VIPER: +1 if the TF has `padj < alpha` (human TFs mapped to
#'     organism symbols when needed)
#'   \item decoupleR: +1 if `p_value < alpha`
#'   \item CSPA: +1 if the gene was selected as a CSPA surface marker for
#'     this contrast
#' }
#' Genes are ranked by total score, then `|stat|`, then `|log2FoldChange|`.
#' Direction (`up` / `down`) comes from the sign of DEG `log2FoldChange`
#' (positive = up, negative = down). The top `top_n` genes are returned
#' **per direction**.
#'
#' @param contrast Contrast label.
#' @param diffexp_out Output of [run_diffexp_pipeline()].
#' @param gsea_out Output of [run_gsea_pipeline()], or `NULL`.
#' @param viper_out Output of [run_viper_pipeline()], or `NULL`.
#' @param decoupler_out Output of [run_decoupler_pipeline()], or `NULL`.
#' @param cspa_out Output of [run_cspa_marker_pipeline()], or `NULL`.
#' @param alpha Significance threshold. Default `0.05`.
#' @param top_n Number of genes to return **per direction** (up and down).
#'   Default `10` (up to `2 * top_n` genes total).
#' @param org Organism (`"human"` / `"mouse"`). Default `NULL` (auto-detect).
#' @return A data frame with columns `gene`, `comparison`, `direction`,
#'   `deg`, `gsea`, `viper`, `decoupler`, `cspa`, `score`, plus DEG
#'   tie-break columns when present. Rows are ordered up-then-down, each
#'   block ranked by score.
#' @export
rank_multihit_genes_contrast <- function(
    contrast,
    diffexp_out,
    gsea_out = NULL,
    viper_out = NULL,
    decoupler_out = NULL,
    cspa_out = NULL,
    alpha = 0.05,
    top_n = 10,
    org = NULL) {

    alpha <- validate_alpha(alpha)
    top_n <- max(1L, as.integer(top_n)[1])
    if (is.null(org) || is.na(org)) {
        org <- .multihit_org(gsea_out, viper_out, decoupler_out, diffexp_out)
    }

    empty <- data.frame(
        gene = character(0), comparison = character(0),
        direction = character(0),
        deg = integer(0), gsea = integer(0), viper = integer(0),
        decoupler = integer(0), cspa = integer(0), score = integer(0),
        stringsAsFactors = FALSE
    )
    de <- diffexp_out$by_comparison[[contrast]]
    if (is.null(de) || !nrow(de)) {
        return(empty)
    }
    gene_col <- if ("gene" %in% colnames(de)) "gene" else colnames(de)[[1]]
    genes <- as.character(de[[gene_col]])
    score <- data.frame(
        gene = genes,
        comparison = contrast,
        deg = as.integer(!is.na(de$sig) & de$sig),
        gsea = 0L,
        viper = 0L,
        decoupler = 0L,
        cspa = 0L,
        stringsAsFactors = FALSE
    )
    if ("stat" %in% colnames(de)) {
        score$stat <- as.numeric(de$stat)
    } else {
        score$stat <- NA_real_
    }
    if ("log2FoldChange" %in% colnames(de)) {
        score$log2FoldChange <- as.numeric(de$log2FoldChange)
    } else {
        score$log2FoldChange <- NA_real_
    }

    # GSEA leading-edge membership (at most +1 total per gene)
    if (!is.null(gsea_out) && !is.null(gsea_out$by_comparison[[contrast]])) {
        sig <- pool_sig_gsea_pathways(gsea_out$by_comparison[[contrast]], alpha = alpha)
        if (nrow(sig) > 0L) {
            le_genes <- unique(unlist(
                lapply(sig$leadingEdge, .parse_leading_edge),
                use.names = FALSE
            ))
            score$gsea <- as.integer(score$gene %in% le_genes)
        }
    }

    # VIPER TFs
    if (!is.null(viper_out) && !is.null(viper_out$contrast_acts[[contrast]])) {
        va <- viper_out$contrast_acts[[contrast]]
        if (!is.null(va) && nrow(va) > 0L && all(c("tf", "padj") %in% colnames(va))) {
            tf_sig <- as.character(va$tf[!is.na(va$padj) & va$padj < alpha])
            tf_org <- .map_human_symbols_to_org(tf_sig, org)
            score$viper <- as.integer(score$gene %in% tf_org)
        }
    }

    # decoupleR sources
    if (!is.null(decoupler_out) && !is.null(decoupler_out$contrast_acts[[contrast]])) {
        da <- decoupler_out$contrast_acts[[contrast]]
        if (!is.null(da) && nrow(da) > 0L &&
                all(c("source", "p_value") %in% colnames(da))) {
            src_sig <- as.character(
                da$source[!is.na(da$p_value) & da$p_value < alpha]
            )
            score$decoupler <- as.integer(score$gene %in% unique(src_sig))
        }
    }

    # CSPA surface markers selected for this contrast
    if (!is.null(cspa_out) && !is.null(cspa_out$markers_by_comparison[[contrast]])) {
        mk <- cspa_out$markers_by_comparison[[contrast]]
        if (!is.null(mk) && nrow(mk) > 0L) {
            mk_col <- if ("gene" %in% colnames(mk)) "gene" else colnames(mk)[[1]]
            score$cspa <- as.integer(score$gene %in% as.character(mk[[mk_col]]))
        }
    }

    score$score <- score$deg + score$gsea + score$viper + score$decoupler + score$cspa
    # DEG direction: sign of log2FoldChange (fallback: sign of Wald stat)
    lfc <- score$log2FoldChange
    direction <- ifelse(!is.na(lfc) & lfc > 0, "up",
                 ifelse(!is.na(lfc) & lfc < 0, "down", NA_character_))
    miss_dir <- is.na(direction)
    if (any(miss_dir)) {
        st <- score$stat
        direction[miss_dir] <- ifelse(
            !is.na(st[miss_dir]) & st[miss_dir] > 0, "up",
            ifelse(!is.na(st[miss_dir]) & st[miss_dir] < 0, "down", NA_character_)
        )
    }
    score$direction <- direction

    # Prefer genes with any hit; then score, |stat|, |LFC|
    keep <- score$score > 0L
    if (!any(keep)) {
        # Fall back to top DEG by |stat| so the plot is not empty
        keep <- !is.na(score$stat)
        if (!any(keep)) keep <- rep(TRUE, nrow(score))
    }
    ranked_all <- score[keep & !is.na(score$direction), , drop = FALSE]
    if (!nrow(ranked_all)) {
        return(empty)
    }

    pick_direction <- function(d) {
        sub <- ranked_all[ranked_all$direction == d, , drop = FALSE]
        if (!nrow(sub)) {
            return(sub[0, , drop = FALSE])
        }
        sub <- sub[order(
            -sub$score,
            -abs(sub$stat),
            -abs(sub$log2FoldChange),
            sub$gene
        ), , drop = FALSE]
        sub[seq_len(min(top_n, nrow(sub))), , drop = FALSE]
    }
    ranked <- rbind(pick_direction("up"), pick_direction("down"))
    ranked$direction <- factor(as.character(ranked$direction), levels = c("up", "down"))
    col_order <- c(
        "gene", "comparison", "direction",
        "deg", "gsea", "viper", "decoupler", "cspa", "score",
        "stat", "log2FoldChange"
    )
    ranked <- ranked[, intersect(col_order, colnames(ranked)), drop = FALSE]
    rownames(ranked) <- NULL
    ranked
}

#' Plot a bidirectional multi-hit gene barplot (GSEA-style)
#'
#' Up genes extend right of zero and down genes left of zero (direction from
#' DEG `log2FoldChange`). Bars are stacked by analysis layer. Gene names sit
#' in the empty half opposite each bar, matching [plot_gsea_barplot()].
#'
#' @param ranked Data frame from [rank_multihit_genes_contrast()].
#' @param title Optional title (also used to label up/down contrast sides).
#' @return A `ggplot` object with `attr(,"plot_size")` width/height inches,
#'   or `NULL` if empty.
#' @export
plot_multihit_genes <- function(ranked, title = NULL) {
    if (is.null(ranked) || !nrow(ranked)) return(NULL)
    layers <- c("deg", "gsea", "viper", "decoupler", "cspa")
    layers <- layers[layers %in% colnames(ranked)]
    has_dir <- "direction" %in% colnames(ranked) &&
        any(!is.na(ranked$direction) & nzchar(as.character(ranked$direction)))

    dir_chr <- if (has_dir) {
        as.character(ranked$direction)
    } else {
        rep("up", nrow(ranked))
    }
    dir_chr[is.na(dir_chr) | !nzchar(dir_chr)] <- "up"
    signed_score <- ifelse(dir_chr == "down", -as.numeric(ranked$score),
                           as.numeric(ranked$score))

    long <- do.call(rbind, lapply(layers, function(ly) {
        hits <- as.numeric(ranked[[ly]])
        data.frame(
            gene = as.character(ranked$gene),
            layer = ly,
            hits = hits,
            hits_signed = ifelse(dir_chr == "down", -hits, hits),
            stringsAsFactors = FALSE
        )
    }))
    long <- long[long$hits > 0, , drop = FALSE]
    if (!nrow(long)) return(NULL)

    # Same ordering convention as GSEA NES bars: most negative at bottom
    ord_idx <- order(signed_score, as.character(ranked$gene))
    gene_levels <- as.character(ranked$gene)[ord_idx]
    long$gene <- factor(long$gene, levels = gene_levels)
    long$layer <- factor(long$layer, levels = layers)

    score_limit <- max(abs(signed_score), na.rm = TRUE)
    if (!is.finite(score_limit) || score_limit <= 0) score_limit <- 1
    label_nudge <- score_limit * 0.06

    lab_df <- data.frame(
        gene = factor(as.character(ranked$gene), levels = gene_levels),
        gene_label = as.character(ranked$gene),
        signed_score = signed_score,
        stringsAsFactors = FALSE
    )
    lab_df <- lab_df[!duplicated(lab_df$gene), , drop = FALSE]
    lab_df$name_hjust <- ifelse(lab_df$signed_score >= 0, 1, 0)
    lab_df$name_x <- ifelse(lab_df$signed_score >= 0, -label_nudge, label_nudge)

    lvls <- .parse_contrast_levels(
        if (is.null(title) || !nzchar(title)) "up_vs_down" else title
    )
    up_label <- paste0(lvls$test, " \u2192")
    down_label <- paste0("\u2190 ", lvls$ref)

    pal <- c(
        deg = "#4C72B0", gsea = "#55A868",
        viper = "#C44E52", decoupler = "#8172B2", cspa = "#CCB974"
    )
    lab <- c(
        deg = "DEG", gsea = "GSEA LE", viper = "VIPER",
        decoupler = "decoupleR", cspa = "CSPA"
    )

    p <- ggplot2::ggplot(
        long,
        ggplot2::aes(x = .data$gene, y = .data$hits_signed, fill = .data$layer)
    ) +
        ggplot2::geom_col(width = 0.86, colour = "black", linewidth = 0.3) +
        ggplot2::geom_text(
            data = lab_df,
            ggplot2::aes(
                x = .data$gene, y = .data$name_x,
                label = .data$gene_label, hjust = .data$name_hjust
            ),
            inherit.aes = FALSE, size = 3.2
        ) +
        ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
        unclipped_coord_flip() +
        ggplot2::scale_y_continuous(
            limits = c(-score_limit, score_limit),
            expand = ggplot2::expansion(mult = c(0.25, 0.25))
        ) +
        ggplot2::scale_fill_manual(
            values = pal[layers],
            breaks = layers,
            labels = unname(lab[layers]),
            name = NULL
        ) +
        ggplot2::annotate(
            "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
            label = up_label, colour = "#B2182B",
            size = 5, fontface = "bold"
        ) +
        ggplot2::annotate(
            "text", x = Inf, y = -Inf, hjust = -0.05, vjust = 1.3,
            label = down_label, colour = "#2166AC",
            size = 5, fontface = "bold"
        ) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            legend.position = "bottom",
            plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
        ) +
        ggplot2::labs(
            x = NULL, y = "Hits",
            title = if (!is.null(title) && nzchar(title)) title else "Cross-analysis gene hits"
        )

    attr(p, "plot_size") <- .multihit_barplot_size(nrow(ranked))
    p
}

#' Suggested PDF / HTML size for a multihit barplot (~30% smaller footprint)
#' @noRd
.multihit_barplot_size <- function(n_genes) {
    n <- max(1L, as.integer(n_genes))
    list(
        width = 0.7 * max(5.5, min(8, 5 + 0.05 * n)),
        height = 0.7 * max(2.8, min(5.5, 0.28 * n + 1.6))
    )
}

#' Rank multi-hit genes for every contrast and write artifacts
#'
#' @inheritParams rank_multihit_genes_contrast
#' @param contrasts Optional character vector of contrast labels. Default
#'   `NULL` uses `names(diffexp_out$by_comparison)`.
#' @param save_dir Directory for CSVs/plots. Default `getwd()`.
#' @param save_plot Logical. Default `TRUE`.
#' @return A list with `by_comparison` (named list of data frames),
#'   `plots` (named list of ggplots), `layers_used`, `org`, `alpha`.
#' @export
run_multihit_pipeline <- function(
    diffexp_out,
    gsea_out = NULL,
    viper_out = NULL,
    decoupler_out = NULL,
    cspa_out = NULL,
    contrasts = NULL,
    alpha = 0.05,
    top_n = 10,
    org = NULL,
    save_dir = getwd(),
    save_plot = TRUE) {

    if (is.null(diffexp_out) || is.null(diffexp_out$by_comparison) ||
            !length(diffexp_out$by_comparison)) {
        stop("`diffexp_out` must contain a non-empty `by_comparison` list.",
             call. = FALSE)
    }
    validate_logical(save_plot, "save_plot")
    save_dir <- validate_save_dir(save_dir)
    if (is.null(org) || is.na(org)) {
        org <- .multihit_org(gsea_out, viper_out, decoupler_out, diffexp_out)
    }
    if (is.null(contrasts)) {
        contrasts <- names(diffexp_out$by_comparison)
    }

    layers_used <- c("DEG")
    if (!is.null(gsea_out)) layers_used <- c(layers_used, "GSEA")
    if (!is.null(viper_out)) layers_used <- c(layers_used, "VIPER")
    if (!is.null(decoupler_out)) layers_used <- c(layers_used, "decoupleR")
    if (!is.null(cspa_out)) layers_used <- c(layers_used, "CSPA")

    by_comparison <- list()
    plots <- list()
    for (comp in contrasts) {
        message("Multi-hit gene ranking: ", comp, " ...")
        ranked <- rank_multihit_genes_contrast(
            contrast = comp,
            diffexp_out = diffexp_out,
            gsea_out = gsea_out,
            viper_out = viper_out,
            decoupler_out = decoupler_out,
            cspa_out = cspa_out,
            alpha = alpha,
            top_n = top_n,
            org = org
        )
        p <- plot_multihit_genes(ranked, title = comp)
        by_comparison[[comp]] <- ranked
        plots[[comp]] <- p

        comp_dir <- file.path(save_dir, comp)
        dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
        write_output(ranked, file.path(comp_dir, paste0(comp, "_top_genes.csv")),
                     format = "csv")
        if (save_plot && !is.null(p)) {
            sz <- attr(p, "plot_size")
            w <- if (!is.null(sz)) sz$width else 4.9
            h <- if (!is.null(sz)) sz$height else max(2.8, 0.2 * nrow(ranked) + 1.5)
            save_ggplot(p, file.path(comp_dir, paste0(comp, "_multihit_barplot.pdf")),
                        width = w, height = h)
        }
    }

    list(
        by_comparison = by_comparison,
        plots = plots,
        layers_used = layers_used,
        org = org,
        alpha = alpha,
        top_n = top_n,
        save_dir = save_dir
    )
}

#' Collect significant TF symbols from a contrast/cluster acts table
#' @noRd
.sig_tf_symbols <- function(acts, alpha, org, kind = c("viper", "decoupler")) {
    kind <- match.arg(kind)
    if (is.null(acts) || !nrow(acts)) return(character(0))
    if (identical(kind, "viper")) {
        if (!all(c("tf", "padj") %in% colnames(acts))) return(character(0))
        tf_sig <- as.character(acts$tf[!is.na(acts$padj) & acts$padj < alpha])
        return(.map_human_symbols_to_org(tf_sig, org))
    }
    if (!all(c("source", "p_value") %in% colnames(acts))) return(character(0))
    unique(as.character(acts$source[!is.na(acts$p_value) & acts$p_value < alpha]))
}

#' Rank genes within one LRT temporal cluster by multi-analysis hits
#'
#' Scores each cluster gene as:
#' \itemize{
#'   \item LRT DEG: +1 (always, for genes in the cluster of significant LRT genes)
#'   \item WGCNA: +1 if the gene is assigned to any WGCNA module
#'   \item VIPER: +1 if the gene is a significant TF for this cluster
#'   \item decoupleR: +1 if the gene is a significant TF for this cluster
#' }
#' Returns the top `top_n` genes by score (then by LRT `stat` when available).
#'
#' @param module Module / cluster label (`"<lrt> cluster <n>"`).
#' @param module_genes Character vector of genes in the cluster.
#' @param lrt_stats Optional named numeric vector of LRT statistics.
#' @param wgcna_genes Character vector of genes in any WGCNA module.
#' @param viper_tfs Character vector of significant VIPER TFs (organism symbols).
#' @param decoupler_tfs Character vector of significant decoupleR TFs.
#' @param top_n Maximum genes to return. Default `10`.
#' @return A data frame with hit columns and `score`.
#' @export
rank_multihit_genes_module <- function(
    module,
    module_genes,
    lrt_stats = NULL,
    wgcna_genes = character(0),
    viper_tfs = character(0),
    decoupler_tfs = character(0),
    top_n = 10) {

    top_n <- max(1L, as.integer(top_n)[1])
    genes <- unique(as.character(module_genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    empty <- data.frame(
        gene = character(0), module = character(0),
        deg = integer(0), wgcna = integer(0),
        viper = integer(0), decoupler = integer(0),
        score = integer(0), stat = numeric(0),
        stringsAsFactors = FALSE
    )
    if (!length(genes)) return(empty)

    wgcna_genes <- unique(as.character(wgcna_genes))
    viper_tfs <- unique(as.character(viper_tfs))
    decoupler_tfs <- unique(as.character(decoupler_tfs))

    score <- data.frame(
        gene = genes,
        module = module,
        deg = 1L,
        wgcna = as.integer(genes %in% wgcna_genes),
        viper = as.integer(genes %in% viper_tfs),
        decoupler = as.integer(genes %in% decoupler_tfs),
        stringsAsFactors = FALSE
    )
    if (!is.null(lrt_stats) && length(lrt_stats)) {
        score$stat <- as.numeric(unname(lrt_stats[genes]))
    } else {
        score$stat <- NA_real_
    }
    score$score <- score$deg + score$wgcna + score$viper + score$decoupler
    # Prefer genes with hits beyond LRT membership alone when available.
    keep <- score$score > 1L
    if (!any(keep)) keep <- rep(TRUE, nrow(score))
    ranked <- score[keep, , drop = FALSE]
    ranked <- ranked[order(
        -ranked$score,
        -replace(ranked$stat, is.na(ranked$stat), -Inf),
        ranked$gene
    ), , drop = FALSE]
    ranked <- ranked[seq_len(min(top_n, nrow(ranked))), , drop = FALSE]
    rownames(ranked) <- NULL
    ranked
}

#' Multi-hit ranking for LRT temporal clusters (tables, no barplots)
#'
#' For each degPatterns cluster, ranks the top `top_n` genes by co-occurrence
#' across LRT DEG membership, WGCNA modules, decoupleR TFs, and VIPER TFs.
#'
#' @param degpatterns_out Output of [run_degpatterns_pipeline()].
#' @param diffexp_out Output of [run_diffexp_lrt_pipeline()].
#' @param wgcna_out Output of [run_wgcna_pipeline()], or `NULL`.
#' @param viper_out Output of [run_viper_pipeline()], or `NULL`.
#' @param decoupler_out Output of [run_decoupler_pipeline()], or `NULL`.
#' @param alpha Significance threshold for TF layers. Default `0.05`.
#' @param top_n Genes per cluster table. Default `10`.
#' @param org Organism; default auto-detect.
#' @param save_dir Output directory. Default `getwd()`.
#' @return A list with `by_module`, `layers_used`, `mode = "lrt_module"`,
#'   `org`, `alpha`, `top_n`.
#' @export
run_multihit_lrt_pipeline <- function(
    degpatterns_out,
    diffexp_out,
    wgcna_out = NULL,
    viper_out = NULL,
    decoupler_out = NULL,
    alpha = 0.05,
    top_n = 10,
    org = NULL,
    save_dir = getwd()) {

    alpha <- validate_alpha(alpha)
    top_n <- max(1L, as.integer(top_n)[1])
    save_dir <- validate_save_dir(save_dir)

    clusters <- tryCatch(
        lrt_clusters_as_wgcna_like(degpatterns_out),
        error = function(e) NULL
    )
    if (is.null(clusters) || is.null(clusters$module_df) ||
            !nrow(clusters$module_df)) {
        stop("No degPatterns clusters available for LRT multi-hit ranking.",
             call. = FALSE)
    }
    if (is.null(org) || is.na(org)) {
        org <- .multihit_org(NULL, viper_out, decoupler_out, diffexp_out)
    }

    wgcna_genes <- character(0)
    if (!is.null(wgcna_out) && !is.null(wgcna_out$module_df) &&
            nrow(wgcna_out$module_df)) {
        gcol <- if ("gene" %in% colnames(wgcna_out$module_df)) {
            "gene"
        } else {
            colnames(wgcna_out$module_df)[[1]]
        }
        # Exclude grey / unassigned when present
        md <- wgcna_out$module_df
        if ("module" %in% colnames(md)) {
            md <- md[!tolower(as.character(md$module)) %in% c("grey", "gray"), ,
                     drop = FALSE]
        }
        wgcna_genes <- unique(as.character(md[[gcol]]))
    }

    layers_used <- c("LRT DEG")
    if (length(wgcna_genes)) layers_used <- c(layers_used, "WGCNA")
    if (!is.null(viper_out)) layers_used <- c(layers_used, "VIPER")
    if (!is.null(decoupler_out)) layers_used <- c(layers_used, "decoupleR")

    mods <- unique(as.character(clusters$module_df$module))
    by_module <- list()
    for (mod in mods) {
        message("Multi-hit gene ranking (LRT module): ", mod, " ...")
        genes <- as.character(
            clusters$module_df$gene[clusters$module_df$module == mod]
        )
        # LRT stats for this parent LRT
        lrt_label <- sub(" cluster .*", "", mod)
        lrt_stats <- NULL
        if (!is.null(diffexp_out$by_comparison[[lrt_label]])) {
            res <- diffexp_out$by_comparison[[lrt_label]]
            gcol <- if ("gene" %in% colnames(res)) "gene" else colnames(res)[[1]]
            if ("stat" %in% colnames(res)) {
                lrt_stats <- stats::setNames(
                    as.numeric(res$stat), as.character(res[[gcol]])
                )
            }
        }

        viper_tfs <- character(0)
        if (!is.null(viper_out) && !is.null(viper_out$contrast_acts[[mod]])) {
            viper_tfs <- .sig_tf_symbols(
                viper_out$contrast_acts[[mod]], alpha, org, "viper"
            )
        } else if (!is.null(viper_out) && length(viper_out$contrast_acts)) {
            # Fall back: any significant TF across cluster-level runs
            viper_tfs <- unique(unlist(lapply(
                viper_out$contrast_acts,
                function(a) .sig_tf_symbols(a, alpha, org, "viper")
            ), use.names = FALSE))
        }

        decoupler_tfs <- character(0)
        if (!is.null(decoupler_out) &&
                !is.null(decoupler_out$contrast_acts[[mod]])) {
            decoupler_tfs <- .sig_tf_symbols(
                decoupler_out$contrast_acts[[mod]], alpha, org, "decoupler"
            )
        } else if (!is.null(decoupler_out) &&
                       length(decoupler_out$contrast_acts)) {
            decoupler_tfs <- unique(unlist(lapply(
                decoupler_out$contrast_acts,
                function(a) .sig_tf_symbols(a, alpha, org, "decoupler")
            ), use.names = FALSE))
        }

        ranked <- rank_multihit_genes_module(
            module = mod,
            module_genes = genes,
            lrt_stats = lrt_stats,
            wgcna_genes = wgcna_genes,
            viper_tfs = viper_tfs,
            decoupler_tfs = decoupler_tfs,
            top_n = top_n
        )
        by_module[[mod]] <- ranked

        safe <- gsub("[^A-Za-z0-9]+", "_", mod)
        mod_dir <- file.path(save_dir, safe)
        dir.create(mod_dir, recursive = TRUE, showWarnings = FALSE)
        write_output(ranked, file.path(mod_dir, paste0(safe, "_top_genes.csv")),
                     format = "csv")
    }

    list(
        by_module = by_module,
        by_comparison = NULL,
        plots = NULL,
        layers_used = layers_used,
        mode = "lrt_module",
        org = org,
        alpha = alpha,
        top_n = top_n,
        save_dir = save_dir
    )
}
