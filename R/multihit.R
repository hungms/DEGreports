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
.multihit_org <- function(gsea_out = NULL, viper_out = NULL, decoupler_tf_out = NULL,
                          diffexp_out = NULL) {
    org <- tryCatch(gsea_out$org, error = function(e) NULL)
    if (is.null(org) || is.na(org)) {
        org <- tryCatch(gsea_out$results$org[[1]], error = function(e) NULL)
    }
    if (is.null(org) || is.na(org)) {
        org <- tryCatch(viper_out$org, error = function(e) NULL)
    }
    if (is.null(org) || is.na(org)) {
        org <- tryCatch(decoupler_tf_out$org, error = function(e) NULL)
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

#' Analysis layers scored by the Wald / DEG-table multi-hit ranking
#' @noRd
.multihit_layers <- c("deg", "gsea", "viper", "decoupler-tf", "decoupler-ccc", "cspa")

#' Genes that are a hit in each analysis, for one contrast
#'
#' The single definition of a multi-hit "hit", shared by the ranking
#' ([rank_multihit_genes_contrast()]) and the method-overlap UpSet plot
#' ([plot_multihit_upset()]). Thresholds are documented in
#' [rank_multihit_genes_contrast()].
#'
#' @inheritParams rank_multihit_genes_contrast
#' @return Named list of character vectors (`deg`, `gsea`, `viper`,
#'   `decoupler-tf`, `decoupler-ccc`, `cspa`); a method that was not run
#'   gives `character(0)`.
#' @export
multihit_hit_sets <- function(contrast, diffexp_out, gsea_out = NULL,
                              viper_out = NULL, decoupler_tf_out = NULL,
                              decoupler_ccc_out = NULL, cspa_out = NULL,
                              alpha = 0.05, lfc_threshold = 0.5,
                              tf_min_score = 1, ccc_min_score = 1,
                              viper_min_nes = 2, org = NULL) {
    sets <- stats::setNames(rep(list(character(0)), length(.multihit_layers)),
                            .multihit_layers)
    de <- diffexp_out$by_comparison[[contrast]]
    if (is.null(de) || !nrow(de)) return(sets)
    if (is.null(org) || is.na(org)) {
        org <- .multihit_org(gsea_out, viper_out, decoupler_tf_out, diffexp_out)
    }
    gene_col <- if ("gene" %in% colnames(de)) "gene" else colnames(de)[[1]]
    genes <- as.character(de[[gene_col]])

    is_deg <- !is.na(de$padj) & de$padj < alpha &
        !is.na(de$log2FoldChange) & abs(de$log2FoldChange) > lfc_threshold
    sets$deg <- unique(genes[is_deg])

    if (!is.null(cspa_out) && length(cspa_out$cspa_symbols)) {
        sets$cspa <- .match_cspa_genes(sets$deg, cspa_out$cspa_symbols)
    }
    if (!is.null(decoupler_tf_out)) {
        sets[["decoupler-tf"]] <- .sig_tf_symbols(
            decoupler_tf_out$contrast_acts[[contrast]], alpha, org, "decoupler",
            min_score = tf_min_score)
    }
    if (!is.null(decoupler_ccc_out)) {
        sets[["decoupler-ccc"]] <- .ccc_sig_genes(
            decoupler_ccc_out, contrast, alpha = alpha, min_score = ccc_min_score)
    }
    if (!is.null(viper_out)) {
        sets$viper <- .sig_tf_symbols(viper_out$contrast_acts[[contrast]], alpha,
                                      org, "viper", min_score = viper_min_nes)
    }
    if (!is.null(gsea_out) && !is.null(gsea_out$by_comparison[[contrast]])) {
        sig <- pool_sig_gsea_pathways(gsea_out$by_comparison[[contrast]], alpha = alpha)
        if (nrow(sig) > 0L) {
            sets$gsea <- unique(unlist(lapply(sig$leadingEdge, .parse_leading_edge),
                                       use.names = FALSE))
        }
    }
    sets
}

#' UpSet plot of how many genes each combination of methods shares
#'
#' Only genes shared by two or more methods are shown: single-method
#' combinations are dropped. Columns run left to right from the most methods
#' intersected to the fewest, with ties broken by gene count.
#'
#' @param sets Named list from [multihit_hit_sets()]. Empty methods are dropped.
#' @param title Optional plot title (the contrast label).
#' @return A `ggplot` wrapping the drawn ComplexHeatmap UpSet, with
#'   `attr(, "plot_size")` (inches) -- or `NULL` when no gene is a hit in
#'   more than one method.
#' @export
plot_multihit_upset <- function(sets, title = NULL) {
    sets <- sets[lengths(sets) > 0L]
    if (length(sets) < 2L) return(NULL)
    names(sets) <- unname(.multihit_layer_labels[names(sets)])
    m <- ComplexHeatmap::make_comb_mat(sets)
    m <- m[ComplexHeatmap::comb_degree(m) >= 2L]
    n_comb <- length(ComplexHeatmap::comb_size(m))
    if (!n_comb) return(NULL)
    width <- max(4.5, min(11, 2.5 + 0.35 * n_comb))
    height <- 1.6 + 0.3 * length(sets) + 1.4
    ht <- ComplexHeatmap::UpSet(
        m,
        set_order = names(sets),
        comb_order = order(-ComplexHeatmap::comb_degree(m),
                           -ComplexHeatmap::comb_size(m)),
        top_annotation = ComplexHeatmap::upset_top_annotation(m, add_numbers = TRUE),
        right_annotation = ComplexHeatmap::upset_right_annotation(m, add_numbers = TRUE),
        column_title = title
    )
    grob <- grid::grid.grabExpr(ComplexHeatmap::draw(ht), width = width, height = height)
    p <- ggplot2::ggplot() + ggplot2::annotation_custom(grob) + ggplot2::theme_void()
    attr(p, "plot_size") <- list(width = width, height = height)
    p
}

#' Display names for the multi-hit layers
#' @noRd
.multihit_layer_labels <- c(
    deg = "DEG", gsea = "GSEA LE", viper = "VIPER",
    `decoupler-tf` = "decoupleR TF", `decoupler-ccc` = "decoupleR CCC",
    cspa = "CSPA"
)

#' Rank genes by cross-analysis hit frequency for one contrast
#'
#' Scores each gene +1 per method it is a hit in (see [multihit_hit_sets()]):
#' \itemize{
#'   \item DEG: `padj < alpha` and `|log2FoldChange| > lfc_threshold`
#'   \item CSPA: a DEG (as above) that is a CSPA surface protein
#'   \item decoupleR-TF: `padj < alpha` and `|score| > tf_min_score`
#'   \item decoupleR-CCC: ligand or receptor of a pair with `padj < alpha`
#'     and `|score| > ccc_min_score`
#'   \item VIPER: `padj < alpha` and `|NES| > viper_min_nes` (human TFs
#'     mapped to organism symbols when needed)
#'   \item GSEA: in **any** significant pathway leading edge (`padj < alpha`;
#'     all collections pooled)
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
#' @param decoupler_tf_out Output of [run_decoupler_tf_pipeline()], or `NULL`.
#' @param decoupler_ccc_out Output of [run_decoupler_ccc_pipeline()], or `NULL`.
#' @param cspa_out Output of [run_cspa_marker_pipeline()], or `NULL`.
#' @param alpha Adjusted-p threshold for every method. Default `0.05`.
#' @param lfc_threshold Absolute log2 fold change a DEG must exceed. Default `0.5`.
#' @param tf_min_score Absolute decoupleR-TF ULM score cutoff. Default `1`.
#' @param ccc_min_score Absolute decoupleR-CCC ULM score cutoff. Default `1`.
#' @param viper_min_nes Absolute VIPER NES cutoff. Default `2`.
#' @param top_n Number of genes to return **per direction** (up and down).
#'   Default `10` (up to `2 * top_n` genes total).
#' @param org Organism (`"human"` / `"mouse"`). Default `NULL` (auto-detect).
#' @return A data frame with columns `gene`, `comparison`, `direction`,
#'   `deg`, `gsea`, `viper`, `decoupler-tf`, `decoupler-ccc`, `cspa`, `score`, plus DEG
#'   tie-break columns when present. Rows are ordered up-then-down, each
#'   block ranked by score.
#' @export
rank_multihit_genes_contrast <- function(
    contrast,
    diffexp_out,
    gsea_out = NULL,
    viper_out = NULL,
    decoupler_tf_out = NULL,
    decoupler_ccc_out = NULL,
    cspa_out = NULL,
    alpha = 0.05,
    lfc_threshold = 0.5,
    tf_min_score = 1,
    ccc_min_score = 1,
    viper_min_nes = 2,
    top_n = 10,
    org = NULL) {

    alpha <- validate_alpha(alpha)
    top_n <- max(1L, as.integer(top_n)[1])
    if (is.null(org) || is.na(org)) {
        org <- .multihit_org(gsea_out, viper_out, decoupler_tf_out, diffexp_out)
    }
    sets <- multihit_hit_sets(
        contrast, diffexp_out, gsea_out = gsea_out, viper_out = viper_out,
        decoupler_tf_out = decoupler_tf_out, decoupler_ccc_out = decoupler_ccc_out,
        cspa_out = cspa_out, alpha = alpha, lfc_threshold = lfc_threshold,
        tf_min_score = tf_min_score, ccc_min_score = ccc_min_score,
        viper_min_nes = viper_min_nes, org = org
    )

    empty <- data.frame(
        gene = character(0), comparison = character(0),
        direction = character(0),
        deg = integer(0), gsea = integer(0), viper = integer(0),
        `decoupler-tf` = integer(0), `decoupler-ccc` = integer(0),
        cspa = integer(0),
        score = integer(0),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    de <- diffexp_out$by_comparison[[contrast]]
    if (is.null(de) || !nrow(de)) {
        return(empty)
    }
    gene_col <- if ("gene" %in% colnames(de)) "gene" else colnames(de)[[1]]
    genes <- as.character(de[[gene_col]])
    score <- data.frame(gene = genes, comparison = contrast,
                        stringsAsFactors = FALSE, check.names = FALSE)
    for (ly in .multihit_layers) {
        score[[ly]] <- as.integer(genes %in% sets[[ly]])
    }
    score$stat <- if ("stat" %in% colnames(de)) as.numeric(de$stat) else NA_real_
    score$log2FoldChange <- if ("log2FoldChange" %in% colnames(de)) {
        as.numeric(de$log2FoldChange)
    } else {
        NA_real_
    }

    score$score <- score$deg + score$gsea + score$viper +
        score[["decoupler-tf"]] + score[["decoupler-ccc"]] + score$cspa
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
        "deg", "gsea", "viper", "decoupler-tf", "decoupler-ccc", "cspa", "score",
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
    layers <- c("deg", "gsea", "viper", "decoupler-tf", "decoupler-ccc", "cspa")
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
        viper = "#C44E52", `decoupler-tf` = "#8172B2",
        `decoupler-ccc` = "#64B5CD", cspa = "#CCB974"
    )
    lab <- .multihit_layer_labels

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
#'   `plots` (named list of barplot ggplots), `hit_sets` (per contrast, from
#'   [multihit_hit_sets()]), `upsets` (named list of [plot_multihit_upset()]
#'   plots), `layers_used`, `org`, `alpha`, `thresholds`.
#' @export
run_multihit_pipeline <- function(
    diffexp_out,
    gsea_out = NULL,
    viper_out = NULL,
    decoupler_tf_out = NULL,
    decoupler_ccc_out = NULL,
    cspa_out = NULL,
    contrasts = NULL,
    alpha = 0.05,
    lfc_threshold = 0.5,
    tf_min_score = 1,
    ccc_min_score = 1,
    viper_min_nes = 2,
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
        org <- .multihit_org(gsea_out, viper_out, decoupler_tf_out, diffexp_out)
    }
    if (is.null(contrasts)) {
        contrasts <- names(diffexp_out$by_comparison)
    }

    layers_used <- c("DEG")
    if (!is.null(gsea_out)) layers_used <- c(layers_used, "GSEA")
    if (!is.null(viper_out)) layers_used <- c(layers_used, "VIPER")
    if (!is.null(decoupler_tf_out)) layers_used <- c(layers_used, "decoupler-tf")
    if (!is.null(decoupler_ccc_out)) layers_used <- c(layers_used, "decoupler-ccc")
    if (!is.null(cspa_out)) layers_used <- c(layers_used, "CSPA")

    thresholds <- list(alpha = alpha, lfc_threshold = lfc_threshold,
                       tf_min_score = tf_min_score, ccc_min_score = ccc_min_score,
                       viper_min_nes = viper_min_nes)
    by_comparison <- list()
    plots <- list()
    hit_sets <- list()
    upsets <- list()
    for (comp in contrasts) {
        message("Multi-hit gene ranking: ", comp, " ...")
        args <- c(list(contrast = comp, diffexp_out = diffexp_out,
                       gsea_out = gsea_out, viper_out = viper_out,
                       decoupler_tf_out = decoupler_tf_out,
                       decoupler_ccc_out = decoupler_ccc_out,
                       cspa_out = cspa_out, org = org), thresholds)
        ranked <- do.call(rank_multihit_genes_contrast, c(args, list(top_n = top_n)))
        sets <- do.call(multihit_hit_sets, args)
        p <- plot_multihit_genes(ranked, title = comp)
        up <- tryCatch(plot_multihit_upset(sets, title = comp), error = function(e) {
            warning("UpSet plot failed for \"", comp, "\": ", conditionMessage(e),
                    call. = FALSE)
            NULL
        })
        by_comparison[[comp]] <- ranked
        plots[[comp]] <- p
        hit_sets[[comp]] <- sets
        upsets[comp] <- list(up)

        comp_dir <- file.path(save_dir, comp)
        dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
        stem <- contrast_file_stem(comp)
        write_output(ranked, file.path(comp_dir, paste0(stem, "_top_genes.csv")),
                     format = "csv")
        long <- data.frame(
            gene = unlist(sets, use.names = FALSE),
            method = rep(unname(.multihit_layer_labels[names(sets)]), lengths(sets)),
            stringsAsFactors = FALSE
        )
        write_output(long, file.path(comp_dir, paste0(stem, "_hit_sets.csv")),
                     format = "csv")
        if (save_plot && !is.null(up)) {
            sz <- attr(up, "plot_size")
            save_ggplot(up, file.path(comp_dir, paste0(stem, "_upset.pdf")),
                        width = sz$width, height = sz$height)
        }
        if (save_plot && !is.null(p)) {
            sz <- attr(p, "plot_size")
            w <- if (!is.null(sz)) sz$width else 4.9
            h <- if (!is.null(sz)) sz$height else max(2.8, 0.2 * nrow(ranked) + 1.5)
            save_ggplot(p, file.path(comp_dir, paste0(stem, "_multihit_barplot.pdf")),
                        width = w, height = h)
        }
    }

    list(
        by_comparison = by_comparison,
        plots = plots,
        hit_sets = hit_sets,
        upsets = upsets,
        layers_used = layers_used,
        org = org,
        alpha = alpha,
        thresholds = thresholds,
        top_n = top_n,
        save_dir = save_dir
    )
}

#' Collect significant TF symbols from a contrast/cluster acts table
#' @noRd
#'
#' VIPER: `padj < alpha` and `|nes| > min_score`; decoupleR: `padj < alpha`
#' and `|score| > min_score` (`padj` is BH over `p_value` when absent).
#' @noRd
.sig_tf_symbols <- function(acts, alpha, org, kind = c("viper", "decoupler"),
                            min_score = if (identical(kind[1], "viper")) 2 else 1) {
    kind <- match.arg(kind)
    if (is.null(acts) || !nrow(acts)) return(character(0))
    if (identical(kind, "viper")) {
        if (!all(c("tf", "padj", "nes") %in% colnames(acts))) return(character(0))
        keep <- !is.na(acts$padj) & acts$padj < alpha & abs(acts$nes) > min_score
        return(.map_human_symbols_to_org(as.character(acts$tf[keep]), org))
    }
    if (!all(c("source", "p_value", "score") %in% colnames(acts))) return(character(0))
    padj <- if ("padj" %in% colnames(acts)) acts$padj else
        stats::p.adjust(acts$p_value, method = "BH")
    keep <- !is.na(padj) & padj < alpha & abs(acts$score) > min_score
    unique(as.character(acts$source[keep]))
}

#' Rank genes within one LRT temporal cluster by multi-analysis hits
#'
#' Scores each cluster gene as:
#' \itemize{
#'   \item LRT DEG: +1 (always, for genes in the cluster of significant LRT genes)
#'   \item WGCNA: +1 if the gene is assigned to any WGCNA module
#'   \item VIPER: +1 if the gene is a significant TF for this cluster
#'   \item decoupleR-TF: +1 if the gene is a significant TF for this cluster
#'   \item decoupleR-CCC: +1 if the gene is in a significant L–R pair for
#'     the parent LRT
#' }
#' Returns the top `top_n` genes by score (then by LRT `stat` when available).
#'
#' @param module Module / cluster label (`"<lrt> cluster <n>"`).
#' @param module_genes Character vector of genes in the cluster.
#' @param lrt_stats Optional named numeric vector of LRT statistics.
#' @param wgcna_genes Character vector of genes in any WGCNA module.
#' @param viper_tfs Character vector of significant VIPER TFs (organism symbols).
#' @param decoupler_tfs Character vector of significant decoupleR TFs.
#' @param ccc_genes Character vector of genes in significant CCC pairs.
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
    ccc_genes = character(0),
    top_n = 10) {

    top_n <- max(1L, as.integer(top_n)[1])
    genes <- unique(as.character(module_genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    empty <- data.frame(
        gene = character(0), module = character(0),
        deg = integer(0), wgcna = integer(0),
        viper = integer(0), `decoupler-tf` = integer(0),
        `decoupler-ccc` = integer(0),
        score = integer(0), stat = numeric(0),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    if (!length(genes)) return(empty)

    wgcna_genes <- unique(as.character(wgcna_genes))
    viper_tfs <- unique(as.character(viper_tfs))
    decoupler_tfs <- unique(as.character(decoupler_tfs))
    ccc_genes <- unique(as.character(ccc_genes))

    score <- data.frame(
        gene = genes,
        module = module,
        deg = 1L,
        wgcna = as.integer(genes %in% wgcna_genes),
        viper = as.integer(genes %in% viper_tfs),
        `decoupler-tf` = as.integer(genes %in% decoupler_tfs),
        `decoupler-ccc` = as.integer(genes %in% ccc_genes),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    if (!is.null(lrt_stats) && length(lrt_stats)) {
        score$stat <- as.numeric(unname(lrt_stats[genes]))
    } else {
        score$stat <- NA_real_
    }
    score$score <- score$deg + score$wgcna + score$viper +
        score[["decoupler-tf"]] + score[["decoupler-ccc"]]
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
#' across LRT DEG membership, WGCNA modules, decoupleR-TF, decoupleR-CCC,
#' and VIPER TFs.
#'
#' @param degpatterns_out Output of [run_degpatterns_pipeline()].
#' @param diffexp_out Output of [run_diffexp_lrt_pipeline()].
#' @param wgcna_out Output of [run_wgcna_pipeline()], or `NULL`.
#' @param viper_out Output of [run_viper_pipeline()], or `NULL`.
#' @param decoupler_tf_out Output of [run_decoupler_tf_pipeline()], or `NULL`.
#' @param decoupler_ccc_out Output of [run_decoupler_ccc_pipeline()], or `NULL`.
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
    decoupler_tf_out = NULL,
    decoupler_ccc_out = NULL,
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
        org <- .multihit_org(NULL, viper_out, decoupler_tf_out, diffexp_out)
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
    if (!is.null(decoupler_tf_out)) layers_used <- c(layers_used, "decoupler-tf")
    if (!is.null(decoupler_ccc_out)) layers_used <- c(layers_used, "decoupler-ccc")

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
        if (!is.null(decoupler_tf_out) &&
                !is.null(decoupler_tf_out$contrast_acts[[mod]])) {
            decoupler_tfs <- .sig_tf_symbols(
                decoupler_tf_out$contrast_acts[[mod]], alpha, org, "decoupler"
            )
        } else if (!is.null(decoupler_tf_out) &&
                       length(decoupler_tf_out$contrast_acts)) {
            decoupler_tfs <- unique(unlist(lapply(
                decoupler_tf_out$contrast_acts,
                function(a) .sig_tf_symbols(a, alpha, org, "decoupler")
            ), use.names = FALSE))
        }

        ccc_genes <- character(0)
        if (!is.null(decoupler_ccc_out)) {
            ccc_genes <- .ccc_sig_genes(decoupler_ccc_out, mod, alpha = alpha)
        }

        ranked <- rank_multihit_genes_module(
            module = mod,
            module_genes = genes,
            lrt_stats = lrt_stats,
            wgcna_genes = wgcna_genes,
            viper_tfs = viper_tfs,
            decoupler_tfs = decoupler_tfs,
            ccc_genes = ccc_genes,
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
