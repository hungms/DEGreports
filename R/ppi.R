#' Initialize a STRINGdb object for human or mouse
#'
#' @param organism `"human"` or `"mouse"`.
#' @param score_threshold Minimum combined STRING score (0--1000). Default
#'   `400` (medium confidence).
#' @param cache_dir Directory for STRINGdb downloads (aliases / links).
#'   Created if missing.
#' @param version STRING database version. Default `"12.0"`.
#' @return A `STRINGdb` reference-class object.
#' @noRd
.init_string_db <- function(organism = c("human", "mouse"),
                            score_threshold = 400L,
                            cache_dir = tempdir(),
                            version = "12.0") {
    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
        stop("Package \"STRINGdb\" is required for PPI connectivity analysis. ",
             "Install with BiocManager::install(\"STRINGdb\").", call. = FALSE)
    }
    organism <- match.arg(organism)
    tax <- if (identical(organism, "mouse")) 10090L else 9606L
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    STRINGdb::STRINGdb$new(
        version = version,
        species = tax,
        score_threshold = as.integer(score_threshold)[1],
        input_directory = normalizePath(cache_dir, mustWork = FALSE)
    )
}

#' Map gene symbols to STRING IDs and return undirected edges among them
#'
#' @param string_db Object from [.init_string_db()].
#' @param genes Character vector of gene symbols.
#' @return A list with `edges` (`from`, `to`, `combined_score` as gene
#'   symbols), `mapping` (gene / STRING_id), and `n_mapped`. `edges` may
#'   have zero rows when no STRING links exist among the mapped genes.
#' @noRd
.string_edges_for_genes <- function(string_db, genes) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    empty <- list(
        edges = data.frame(from = character(), to = character(),
                           combined_score = integer(), stringsAsFactors = FALSE),
        mapping = data.frame(gene = character(), STRING_id = character(),
                             stringsAsFactors = FALSE),
        n_mapped = 0L
    )
    if (length(genes) < 2L) return(empty)

    mapped <- quiet(string_db$map(
        data.frame(gene = genes, stringsAsFactors = FALSE),
        "gene",
        removeUnmappedRows = TRUE
    ))
    if (is.null(mapped) || !nrow(mapped) || !"STRING_id" %in% colnames(mapped)) {
        return(empty)
    }
    mapped <- unique(mapped[, c("gene", "STRING_id"), drop = FALSE])
    mapped$gene <- as.character(mapped$gene)
    mapped$STRING_id <- as.character(mapped$STRING_id)
    # Prefer one STRING ID per gene symbol (first hit).
    mapped <- mapped[!duplicated(mapped$gene), , drop = FALSE]
    if (nrow(mapped) < 2L) {
        return(list(edges = empty$edges, mapping = mapped, n_mapped = nrow(mapped)))
    }

    ints <- quiet(string_db$get_interactions(mapped$STRING_id))
    if (is.null(ints) || !nrow(ints)) {
        return(list(edges = empty$edges, mapping = mapped, n_mapped = nrow(mapped)))
    }
    id2gene <- stats::setNames(mapped$gene, mapped$STRING_id)
    from <- unname(id2gene[as.character(ints$from)])
    to <- unname(id2gene[as.character(ints$to)])
    score <- as.integer(ints$combined_score)
    ok <- !is.na(from) & !is.na(to) & from != to
    from <- from[ok]
    to <- to[ok]
    score <- score[ok]
    if (!length(from)) {
        return(list(edges = empty$edges, mapping = mapped, n_mapped = nrow(mapped)))
    }

    # Canonical undirected edge (alphabetical endpoints); keep max score.
    a <- ifelse(from < to, from, to)
    b <- ifelse(from < to, to, from)
    edges <- data.frame(from = a, to = b, combined_score = score,
                        stringsAsFactors = FALSE)
    edges <- edges[order(-edges$combined_score), , drop = FALSE]
    edges <- edges[!duplicated(paste(edges$from, edges$to, sep = "\t")), ,
                   drop = FALSE]
    rownames(edges) <- NULL
    list(edges = edges, mapping = mapped, n_mapped = nrow(mapped))
}

#' Choose seed genes for a DE contrast PPI plot
#'
#' Prefers significant DEGs (`sig`); if fewer than 2, falls back to the top
#' `max_genes` by `rank_score`. Caps at `max_genes` (ranked by
#' `rank_score` when available).
#'
#' @noRd
.ppi_seed_genes_de <- function(res_df, gene_col = "gene", max_genes = Inf) {
    genes <- as.character(res_df[[gene_col]])
    ok <- !is.na(genes) & nzchar(genes)
    res_df <- res_df[ok, , drop = FALSE]
    if (!nrow(res_df)) return(character(0))

    if ("sig" %in% colnames(res_df)) {
        seed_idx <- which(!is.na(res_df$sig) & as.logical(res_df$sig))
    } else {
        seed_idx <- integer(0)
    }
    if (length(seed_idx) < 2L) {
        seed_idx <- seq_len(nrow(res_df))
    }
    seed <- res_df[seed_idx, , drop = FALSE]
    if ("rank_score" %in% colnames(seed)) {
        seed <- seed[order(-abs(seed$rank_score)), , drop = FALSE]
    }
    out <- unique(as.character(seed[[gene_col]]))
    if (!is.null(max_genes) && is.finite(max_genes)) {
        out <- utils::head(out, as.integer(max_genes))
    }
    out
}

#' Choose seed genes for a WGCNA module PPI plot
#'
#' Module members sorted by absolute kME, capped at `max_genes`.
#'
#' @noRd
.ppi_seed_genes_module <- function(module_df, module, gene_col = "gene",
                                   max_genes = Inf) {
    sub <- module_df[as.character(module_df$module) == as.character(module), ,
                     drop = FALSE]
    if (!nrow(sub)) return(character(0))
    if (!gene_col %in% colnames(sub)) {
        gene_col <- setdiff(colnames(sub),
                            c("module", "module_label", "module_color", "kME"))[1]
    }
    if ("kME" %in% colnames(sub)) {
        sub <- sub[order(-abs(sub$kME)), , drop = FALSE]
    }
    genes <- as.character(sub[[gene_col]])
    genes <- unique(genes[!is.na(genes) & nzchar(genes)])
    if (!is.null(max_genes) && is.finite(max_genes)) {
        genes <- utils::head(genes, as.integer(max_genes))
    }
    genes
}

#' Build an induced PPI subgraph and rank hubs by degree
#'
#' @param genes Character vector of seed gene symbols.
#' @param net Undirected edge list (`from`, `to`) -- typically from
#'   [.string_edges_for_genes()].
#' @param top_n_hubs Number of highest-degree genes to flag as hubs.
#' @param max_plot_nodes If the connected subgraph exceeds this size, keep
#'   only the highest-degree nodes (and edges among them) for plotting.
#' @return A list with `nodes` (full degree table for every gene with ≥1
#'   edge), `plot_nodes` / `edges` (plot-capped subgraph), `hubs`,
#'   `n_seed`, `n_with_edge`, `n_plotted`, or `NULL` if fewer than 2 genes
#'   share an edge.
#' @noRd
.build_ppi_subgraph <- function(genes, net, top_n_hubs = 10L,
                                max_plot_nodes = 50L) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    n_seed <- length(genes)
    if (n_seed < 2L || is.null(net) || !nrow(net)) return(NULL)

    edges_all <- net[net$from %in% genes & net$to %in% genes, , drop = FALSE]
    if (!nrow(edges_all)) return(NULL)

    verts_all <- sort(unique(c(edges_all$from, edges_all$to)))
    n_with_edge <- length(verts_all)
    deg_all <- table(factor(c(edges_all$from, edges_all$to), levels = verts_all))
    deg_all <- as.integer(deg_all)
    names(deg_all) <- verts_all

    hub_n <- min(as.integer(top_n_hubs), length(verts_all))
    hubs <- names(sort(deg_all, decreasing = TRUE))[seq_len(hub_n)]

    nodes_all <- data.frame(
        gene = verts_all,
        degree = as.integer(deg_all[verts_all]),
        hub = verts_all %in% hubs,
        stringsAsFactors = FALSE
    )
    nodes_all <- nodes_all[order(-nodes_all$degree, nodes_all$gene), , drop = FALSE]
    rownames(nodes_all) <- NULL

    edges <- edges_all
    verts <- verts_all
    deg <- deg_all
    if (length(verts) > max_plot_nodes) {
        keep <- names(sort(deg, decreasing = TRUE))[seq_len(max_plot_nodes)]
        edges <- edges[edges$from %in% keep & edges$to %in% keep, , drop = FALSE]
        verts <- sort(unique(c(edges$from, edges$to)))
        deg <- table(factor(c(edges$from, edges$to), levels = verts))
        deg <- as.integer(deg)
        names(deg) <- verts
    }
    if (length(verts) < 2L || !nrow(edges)) return(NULL)

    plot_nodes <- data.frame(
        gene = verts,
        degree = as.integer(deg[verts]),
        # Keep hub flags from the full network ranking (not re-ranked after plot cap).
        hub = verts %in% hubs,
        stringsAsFactors = FALSE
    )
    plot_nodes <- plot_nodes[order(-plot_nodes$degree, plot_nodes$gene), , drop = FALSE]
    rownames(plot_nodes) <- NULL

    list(
        nodes = nodes_all,
        plot_nodes = plot_nodes,
        edges = edges,
        hubs = hubs,
        n_seed = n_seed,
        n_with_edge = n_with_edge,
        n_plotted = nrow(plot_nodes)
    )
}

#' Add Louvain community membership (BioNERO-style node colouring)
#' @noRd
.ppi_add_communities <- function(nodes, edges) {
    if (!requireNamespace("igraph", quietly = TRUE)) return(nodes)
    g <- igraph::graph_from_data_frame(
        edges[, c("from", "to"), drop = FALSE],
        directed = FALSE,
        vertices = nodes["gene"]
    )
    g <- igraph::simplify(g)
    set.seed(1L)
    com <- tryCatch(
        igraph::cluster_louvain(g),
        error = function(e) igraph::cluster_fast_greedy(g)
    )
    mem <- igraph::membership(com)
    nodes$community <- factor(as.integer(mem[nodes$gene]))
    nodes
}

#' Qualitative palette for community colouring (BioNERO-like)
#' @noRd
.ppi_community_palette <- function(n) {
    base <- c(
        "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
        "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
        "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
        "#D55E00", "#CC79A7", "#999999", "#66C2A5", "#FC8D62"
    )
    if (n <= length(base)) return(base[seq_len(n)])
    grDevices::colorRampPalette(base)(n)
}

#' Layout an igraph for ggplot
#' @noRd
.ppi_layout_df <- function(nodes, edges) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package \"igraph\" is required for PPI connectivity plots. ",
             "Install with install.packages(\"igraph\").", call. = FALSE)
    }
    g <- igraph::graph_from_data_frame(
        edges[, c("from", "to"), drop = FALSE],
        directed = FALSE,
        vertices = nodes["gene"]
    )
    set.seed(1L)
    lay <- igraph::layout_with_fr(g)
    colnames(lay) <- c("x", "y")
    lay_df <- data.frame(gene = igraph::V(g)$name, lay, stringsAsFactors = FALSE)
    nodes <- merge(nodes, lay_df, by = "gene", sort = FALSE)
    edge_df <- data.frame(
        from = edges$from,
        to = edges$to,
        x = nodes$x[match(edges$from, nodes$gene)],
        y = nodes$y[match(edges$from, nodes$gene)],
        xend = nodes$x[match(edges$to, nodes$gene)],
        yend = nodes$y[match(edges$to, nodes$gene)],
        stringsAsFactors = FALSE
    )
    if ("combined_score" %in% colnames(edges)) {
        edge_df$combined_score <- edges$combined_score
    }
    list(nodes = nodes, edges = edge_df)
}

#' Plot a STRING PPI network in BioNERO style
#'
#' Renders an undirected protein-interaction graph with Fruchterman--Reingold
#' layout. Node size scales with degree; nodes are coloured by Louvain
#' community (BioNERO `plot_ppi` style); the top hubs by degree are labelled.
#' A `text` aesthetic is attached for [make_interactive_ppi()] / plotly hover.
#'
#' @param nodes Node table from [.build_ppi_subgraph()] (`gene`, `degree`,
#'   `hub`), optionally with a `community` column.
#' @param edges Edge table (`from`, `to`).
#' @param title Optional plot title.
#' @param top_n_hubs Number of hubs to label (uses `nodes$hub` when present,
#'   otherwise the top degrees). Default `10`.
#' @param add_color_legend Show community colour legend. Default `TRUE`.
#' @param show_hub_labels Draw static hub labels (`ggrepel`). Set `FALSE`
#'   when building a hover-only interactive plot. Default `TRUE`.
#' @return A `ggplot` object, or `NULL` if the graph is empty.
#' @export
plot_ppi_connectivity <- function(nodes, edges, title = NULL,
                                  top_n_hubs = 10L,
                                  add_color_legend = TRUE,
                                  show_hub_labels = TRUE) {
    if (is.null(nodes) || !nrow(nodes) || is.null(edges) || !nrow(edges)) {
        return(NULL)
    }
    if (!"community" %in% colnames(nodes)) {
        nodes <- .ppi_add_communities(nodes, edges)
    }
    lay <- .ppi_layout_df(nodes, edges)
    nodes <- lay$nodes
    edge_df <- lay$edges

    if (!"hub" %in% colnames(nodes)) {
        ord <- order(-nodes$degree, nodes$gene)
        nodes$hub <- FALSE
        nodes$hub[ord[seq_len(min(as.integer(top_n_hubs), nrow(nodes)))]] <- TRUE
    }
    hub_df <- nodes[nodes$hub, , drop = FALSE]
    hub_df <- hub_df[order(-hub_df$degree, hub_df$gene), , drop = FALSE]
    hub_df <- utils::head(hub_df, as.integer(top_n_hubs))

    nodes$community <- factor(nodes$community)
    n_com <- nlevels(nodes$community)
    pal <- .ppi_community_palette(n_com)
    names(pal) <- levels(nodes$community)

    if ("combined_score" %in% colnames(edge_df) &&
        any(!is.na(edge_df$combined_score))) {
        sc <- edge_df$combined_score
        edge_df$edge_width <- 0.2 + 1.2 * (sc - min(sc, na.rm = TRUE)) /
            max(1, max(sc, na.rm = TRUE) - min(sc, na.rm = TRUE))
    } else {
        edge_df$edge_width <- 0.45
    }

    tip <- paste0(
        "<b>", nodes$gene, "</b>",
        "<br>degree: ", nodes$degree,
        "<br>community: ", as.character(nodes$community),
        "<br>hub: ", ifelse(nodes$hub, "yes", "no")
    )
    if ("direction" %in% colnames(nodes)) {
        tip <- paste0(tip, "<br>DE: ", nodes$direction)
    }
    if ("kME" %in% colnames(nodes)) {
        tip <- paste0(tip, "<br>kME: ",
                      formatC(nodes$kME, format = "f", digits = 3))
    }
    nodes$.tooltip <- tip

    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = edge_df,
            ggplot2::aes(x = .data$x, y = .data$y,
                         xend = .data$xend, yend = .data$yend,
                         linewidth = .data$edge_width),
            colour = "grey75", alpha = 0.55, show.legend = FALSE
        ) +
        ggplot2::geom_point(
            data = nodes,
            ggplot2::aes(x = .data$x, y = .data$y,
                         size = .data$degree, colour = .data$community),
            alpha = 0.9
        ) +
        ggplot2::scale_colour_manual(
            values = pal, name = "Community",
            guide = if (isTRUE(add_color_legend)) "legend" else "none"
        ) +
        ggplot2::scale_size_continuous(range = c(2.5, 10), name = "Degree") +
        ggplot2::scale_linewidth_identity()

    if (isTRUE(show_hub_labels) && nrow(hub_df)) {
        p <- p + ggrepel::geom_label_repel(
            data = hub_df,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$gene),
            size = 3, colour = "grey25", fill = "white",
            label.size = 0.15, max.overlaps = Inf, seed = 1L,
            box.padding = 0.4, point.padding = 0.25,
            segment.colour = "grey50", min.segment.length = 0,
            show.legend = FALSE
        )
    }

    p <- p +
        ggplot2::coord_equal() +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, colour = "grey30"),
            legend.position = "right",
            plot.margin = ggplot2::margin(6, 6, 6, 6)
        ) +
        ggplot2::labs(subtitle = paste0(
            "STRING PPI (BioNERO-style); top ", nrow(hub_df),
            " hubs labelled (", nrow(nodes), " genes, ",
            nrow(edge_df), " edges)"
        ))

    if (!is.null(title) && nzchar(title)) {
        p <- p + ggplot2::labs(title = title)
    }
    # Stash layout tables so make_interactive_ppi() can reuse coordinates.
    attr(p, "ppi_nodes") <- nodes
    attr(p, "ppi_edges") <- edge_df
    attr(p, "ppi_palette") <- pal
    p
}

#' Interactive plotly PPI network with per-node gene hover labels
#'
#' Builds a plotly scatter network from the same layout used by
#' [plot_ppi_connectivity()]. Hover any node to see its gene symbol,
#' degree, community, and (when present) DE direction / kME. Top hubs
#' keep a permanent text label.
#'
#' @param nodes,edges As for [plot_ppi_connectivity()], or a ggplot from
#'   that function (layout is read from its attributes).
#' @param title Optional title.
#' @param top_n_hubs Number of hubs to keep permanently labelled. Default
#'   `10`.
#' @param width,height Widget size in pixels. Defaults `720`, `620`.
#' @return A `plotly` htmlwidget, or `NULL` if the graph is empty.
#' @export
make_interactive_ppi <- function(nodes, edges = NULL, title = NULL,
                                 top_n_hubs = 10L,
                                 width = 720, height = 620) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Package \"plotly\" is required for interactive PPI plots.",
             call. = FALSE)
    }

    if (inherits(nodes, "ggplot")) {
        p_gg <- nodes
        nodes <- attr(p_gg, "ppi_nodes")
        edges <- attr(p_gg, "ppi_edges")
        pal <- attr(p_gg, "ppi_palette")
        if (is.null(title) || !nzchar(title)) {
            title <- p_gg$labels$title
        }
    } else {
        pal <- NULL
        p_gg <- plot_ppi_connectivity(
            nodes, edges, title = title, top_n_hubs = top_n_hubs,
            show_hub_labels = FALSE
        )
        if (is.null(p_gg)) return(NULL)
        nodes <- attr(p_gg, "ppi_nodes")
        edges <- attr(p_gg, "ppi_edges")
        pal <- attr(p_gg, "ppi_palette")
    }
    if (is.null(nodes) || !nrow(nodes) || is.null(edges) || !nrow(edges)) {
        return(NULL)
    }

    nodes$community <- factor(nodes$community)
    if (is.null(pal) || !length(pal)) {
        pal <- .ppi_community_palette(nlevels(nodes$community))
        names(pal) <- levels(nodes$community)
    }
    nodes$colour <- unname(pal[as.character(nodes$community)])
    deg <- pmax(nodes$degree, 1)
    nodes$marker_size <- 8 + 18 * sqrt(deg) / max(sqrt(deg))

    hub_n <- min(as.integer(top_n_hubs), nrow(nodes))
    hub_genes <- utils::head(nodes$gene[order(-nodes$degree, nodes$gene)], hub_n)
    nodes$label <- ifelse(nodes$gene %in% hub_genes, nodes$gene, "")

    # Edge path as a single trace with NA breaks (plotly idiom).
    edge_x <- as.numeric(rbind(edges$x, edges$xend, NA_real_))
    edge_y <- as.numeric(rbind(edges$y, edges$yend, NA_real_))

    p <- plotly::plot_ly(width = width, height = height)
    p <- plotly::add_trace(
        p,
        x = edge_x, y = edge_y,
        type = "scatter", mode = "lines",
        line = list(color = "rgba(160,160,160,0.55)", width = 0.7),
        hoverinfo = "skip", showlegend = FALSE,
        name = "edges"
    )
    # One trace per community so the legend stays readable.
    for (com in levels(nodes$community)) {
        sub <- nodes[as.character(nodes$community) == com, , drop = FALSE]
        if (!nrow(sub)) next
        p <- plotly::add_trace(
            p,
            x = sub$x, y = sub$y,
            type = "scatter", mode = "markers+text",
            text = sub$label,
            textposition = "top middle",
            textfont = list(size = 11, color = "#333333"),
            marker = list(
                size = sub$marker_size,
                color = sub$colour[1],
                line = list(color = "rgba(40,40,40,0.55)", width = 0.6)
            ),
            hovertext = if (".tooltip" %in% colnames(sub)) sub$.tooltip else paste0("<b>", sub$gene, "</b>"),
            hoverinfo = "text",
            name = paste("Community", com),
            showlegend = TRUE
        )
    }
    p <- plotly::layout(
        p,
        title = list(text = title, x = 0.5, xanchor = "center"),
        xaxis = list(visible = FALSE, scaleanchor = "y", scaleratio = 1),
        yaxis = list(visible = FALSE),
        legend = list(orientation = "v", x = 1.02, y = 1, xanchor = "left"),
        margin = list(l = 20, r = 120, t = 50, b = 20),
        hovermode = "closest",
        plot_bgcolor = "white",
        paper_bgcolor = "white"
    )
    plotly::config(p, displaylogo = FALSE, responsive = FALSE)
}

#' Annotate PPI nodes with DE direction (CSV metadata; plot uses communities)
#' @noRd
.ppi_add_de_colours <- function(nodes, res_df, gene_col = "gene",
                                alpha = 0.05, lfc_threshold = 0.5) {
    dir <- classify_de_direction(res_df, alpha = alpha, lfc_threshold = lfc_threshold)
    map <- stats::setNames(as.character(dir), as.character(res_df[[gene_col]]))
    d <- unname(map[nodes$gene])
    d[is.na(d)] <- "ns"
    nodes$direction <- d
    nodes
}

#' Annotate PPI nodes with kME for CSV export
#' @noRd
.ppi_add_kme <- function(nodes, module_df, module, gene_col = "gene") {
    sub <- module_df[as.character(module_df$module) == as.character(module), ,
                     drop = FALSE]
    if (!gene_col %in% colnames(sub)) {
        gene_col <- setdiff(colnames(sub),
                            c("module", "module_label", "module_color", "kME"))[1]
    }
    if (!"kME" %in% colnames(sub)) return(nodes)
    map <- stats::setNames(sub$kME, as.character(sub[[gene_col]]))
    nodes$kME <- unname(map[nodes$gene])
    nodes
}

#' Build one STRING PPI result (plot + tables) for a gene set
#' @noRd
.ppi_result_for_genes <- function(genes, string_db, top_n_hubs, max_plot_nodes,
                                  title, annotate = NULL, n_sig = NA_integer_) {
    str <- .string_edges_for_genes(string_db, genes)
    sub <- .build_ppi_subgraph(
        genes, str$edges, top_n_hubs = top_n_hubs, max_plot_nodes = max_plot_nodes
    )
    if (is.null(sub)) {
        return(list(
            plot = NULL,
            nodes = NULL, edges = NULL, hubs = character(0),
            n_sig = n_sig,
            n_seed = length(unique(genes)),
            n_mapped = str$n_mapped,
            n_with_edge = 0L, n_plotted = 0L,
            max_plot_nodes = as.integer(max_plot_nodes),
            mapping = str$mapping
        ))
    }
    # Plot uses the degree-capped subgraph; the HTML/CSV degree table keeps
    # every gene with ≥1 STRING edge among seeds.
    plot_nodes <- .ppi_add_communities(sub$plot_nodes, sub$edges)
    nodes <- sub$nodes
    if (is.function(annotate)) {
        plot_nodes <- annotate(plot_nodes)
        nodes <- annotate(nodes)
    }
    # Carry community labels onto the full table where the gene was plotted.
    if ("community" %in% colnames(plot_nodes)) {
        com_map <- stats::setNames(
            as.character(plot_nodes$community), plot_nodes$gene
        )
        nodes$community <- unname(com_map[nodes$gene])
    }
    plot <- plot_ppi_connectivity(
        plot_nodes, sub$edges, title = title, top_n_hubs = top_n_hubs,
        show_hub_labels = TRUE
    )
    # Native STRING image should use only plotted genes, not all seeds.
    # Interactive plotly is built at report-render time from `plot` attrs
    # via make_interactive_ppi() (avoids serializing large htmlwidgets
    # through the pipeline return value).
    mapping_plot <- str$mapping[str$mapping$gene %in% plot_nodes$gene, , drop = FALSE]
    list(
        plot = plot,
        nodes = nodes,
        edges = sub$edges,
        hubs = sub$hubs,
        n_sig = n_sig,
        n_seed = sub$n_seed,
        n_mapped = str$n_mapped,
        n_with_edge = sub$n_with_edge,
        n_plotted = sub$n_plotted,
        max_plot_nodes = as.integer(max_plot_nodes),
        mapping = mapping_plot
    )
}

#' Save STRINGdb's native network image for a gene set (optional companion)
#' @noRd
.save_string_native_plot <- function(string_db, mapping, path) {
    if (is.null(mapping) || !nrow(mapping) || nrow(mapping) < 2L) return(invisible(NULL))
    ids <- unique(as.character(mapping$STRING_id))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (length(ids) < 2L) return(invisible(NULL))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tryCatch({
        grDevices::png(path, width = 900, height = 900, res = 120)
        string_db$plot_network(ids)
        grDevices::dev.off()
    }, error = function(e) {
        if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
        message("STRINGdb native plot failed (", conditionMessage(e), ").")
    })
    invisible(path)
}

#' Run STRING PPI connectivity analysis for DE contrasts and WGCNA modules
#'
#' For every differential-expression contrast and (optionally) every WGCNA
#' module, maps seed genes to [STRINGdb], retrieves interactions among them,
#' and draws a BioNERO-style connectivity plot (Louvain communities, degree-
#' sized nodes, top hubs labelled). Also writes STRINGdb's native network
#' PNG alongside each ggplot PDF. Imports packages `STRINGdb` and
#' `igraph` (and network access on first run to download STRING files).
#'
#' @param diffexp_out Output of [run_diffexp_pipeline()].
#' @param wgcna_out Optional output of [run_wgcna_pipeline()].
#' @param org Organism (`"human"` / `"mouse"`), or `NULL` to detect from
#'   gene symbols in `diffexp_out`.
#' @param top_n_hubs Number of highest-degree genes to label per plot.
#'   Default `10`.
#' @param max_genes Maximum seed genes per contrast/module. Default `Inf`
#'   keeps every significant DEG (or every module member). Set a finite
#'   integer to cap by `rank_score` / `|kME|`.
#' @param max_plot_nodes Maximum nodes drawn per plot. When more genes have
#'   STRING edges among the seed set, the highest-degree genes are kept.
#'   Default `50`.
#' @param score_threshold Minimum STRING combined score (0--1000). Default
#'   `400`.
#' @param gene_col Gene-id column name. Default `"gene"`.
#' @param alpha,lfc_threshold Thresholds used only to annotate DE direction
#'   on node tables. Defaults `0.05`, `0.5`.
#' @param save_dir Directory for **contrast** PPI folders
#'   (`<save_dir>/<comparison>/ppi/`). Pipeline uses `deseq2/`. Default `getwd()`.
#' @param module_save_dir Directory for **LRT-cluster** (non-WGCNA) module PPI
#'   folders (`<module_save_dir>/<module>/ppi/`). Default `NULL` uses `save_dir`
#'   (pipeline uses `ppi/`).
#' @param wgcna_module_save_dir Base directory for **WGCNA** module PPI
#'   (`<wgcna_module_save_dir>/<Module_N>/ppi/`). Default `NULL` writes WGCNA
#'   modules under `module_save_dir` (pipeline uses `wgcna/`).
#' @param save_plot Logical; write PDFs/PNGs alongside CSVs. Default `TRUE`.
#' @return A list with `org`, `by_comparison`, `by_module`, and STRING
#'   parameters. The full organism-wide edge list is not retained.
#' @export
run_ppi_pipeline <- function(
    diffexp_out,
    wgcna_out = NULL,
    org = NULL,
    top_n_hubs = 10L,
    max_genes = Inf,
    max_plot_nodes = 50L,
    score_threshold = 400L,
    gene_col = "gene",
    alpha = 0.05,
    lfc_threshold = 0.5,
    save_dir = getwd(),
    module_save_dir = NULL,
    wgcna_module_save_dir = NULL,
    save_plot = TRUE) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package \"igraph\" is required for PPI connectivity analysis. ",
             "Install with install.packages(\"igraph\").", call. = FALSE)
    }
    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
        stop("Package \"STRINGdb\" is required for PPI connectivity analysis. ",
             "Install with BiocManager::install(\"STRINGdb\").", call. = FALSE)
    }

    save_dir <- validate_save_dir(save_dir)
    if (is.null(module_save_dir)) {
        module_save_dir <- save_dir
    } else {
        module_save_dir <- validate_save_dir(module_save_dir)
    }
    if (!is.null(wgcna_module_save_dir)) {
        wgcna_module_save_dir <- validate_save_dir(wgcna_module_save_dir)
    }
    top_n_hubs <- as.integer(top_n_hubs)[1]
    max_plot_nodes <- as.integer(max_plot_nodes)[1]
    if (!is.null(max_genes) && is.finite(max_genes)) {
        max_genes <- as.integer(max_genes)[1]
        if (is.na(max_genes) || max_genes < 2L) {
            stop("`max_genes` must be >= 2, Inf, or NULL.", call. = FALSE)
        }
    } else {
        max_genes <- Inf
    }
    if (is.na(top_n_hubs) || top_n_hubs < 1L) {
        stop("`top_n_hubs` must be a positive integer.", call. = FALSE)
    }
    if (is.na(max_plot_nodes) || max_plot_nodes < 2L) {
        stop("`max_plot_nodes` must be >= 2.", call. = FALSE)
    }

    by_comp <- diffexp_out$by_comparison
    if (is.null(by_comp) || !length(by_comp)) {
        stop("`diffexp_out$by_comparison` is empty; cannot run PPI analysis.",
             call. = FALSE)
    }

    if (is.null(org) || !nzchar(as.character(org)[1]) ||
        is.na(as.character(org)[1])) {
        sample_genes <- unlist(lapply(by_comp, function(df) {
            as.character(df[[gene_col]])
        }), use.names = FALSE)
        org <- detect_organism(sample_genes)
    } else {
        org <- match.arg(as.character(org)[1], c("human", "mouse"))
    }

    cache_dir <- file.path(module_save_dir, "stringdb_cache")
    message("Initializing STRINGdb for organism=\"", org,
            "\" (score >= ", score_threshold, ") ...")
    string_db <- .init_string_db(
        organism = org,
        score_threshold = score_threshold,
        cache_dir = cache_dir
    )

    # Contrast folders use the raw deseq2 contrast label; module folders
    # sanitize names that contain spaces (e.g. "Module 1" -> "Module_1").
    safe_mod_dir <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
    .write_ppi_outputs <- function(out, out_dir, save_plot) {
        write_output(out$nodes, file.path(out_dir, "ppi_nodes.csv"))
        write_output(out$edges, file.path(out_dir, "ppi_edges.csv"))
        # Full degree-ranked node table (same as the HTML hub table / CSV button).
        hubs_out <- out$nodes[order(-out$nodes$degree, out$nodes$gene), , drop = FALSE]
        write_output(hubs_out, file.path(out_dir, "ppi_hubs.csv"))
        if (save_plot) {
            save_ggplot(out$plot, file.path(out_dir, "ppi_connectivity.pdf"),
                        width = 6.5, height = 5.5)
            .save_string_native_plot(
                string_db, out$mapping,
                file.path(out_dir, "stringdb_network.png")
            )
        }
    }

    by_comparison <- list()
    message("Building STRING PPI plots for ", length(by_comp),
            " DE contrast(s) ...")
    for (comp in names(by_comp)) {
        res <- by_comp[[comp]]
        n_sig <- if ("sig" %in% colnames(res)) {
            sum(!is.na(res$sig) & as.logical(res$sig))
        } else {
            NA_integer_
        }
        genes <- .ppi_seed_genes_de(res, gene_col = gene_col, max_genes = max_genes)
        out <- .ppi_result_for_genes(
            genes, string_db, top_n_hubs, max_plot_nodes, title = comp,
            n_sig = n_sig,
            annotate = function(nodes) {
                .ppi_add_de_colours(nodes, res, gene_col = gene_col,
                                    alpha = alpha, lfc_threshold = lfc_threshold)
            }
        )
        if (is.null(out$plot)) {
            message("  Skipping \"", comp, "\": fewer than 2 mapped genes share ",
                    "a STRING edge (mapped ", out$n_mapped, "/", out$n_seed, ").")
            next
        }
        message("  ", comp, ": sig=", out$n_sig,
                " -> seed=", out$n_seed,
                " -> STRING-mapped=", out$n_mapped,
                " -> with edges=", out$n_with_edge,
                " -> plotted=", out$n_plotted)
        by_comparison[[comp]] <- out
        .write_ppi_outputs(out, file.path(save_dir, comp, "ppi"), save_plot)
    }

    by_module <- list()
    if (!is.null(wgcna_out) && !is.null(wgcna_out$module_df)) {
        mod_df <- wgcna_out$module_df
        modules <- setdiff(unique(as.character(mod_df$module)), "Unassigned")
        sizes <- sort(table(mod_df$module[mod_df$module != "Unassigned"]),
                      decreasing = TRUE)
        modules <- intersect(names(sizes), modules)
        message("Building STRING PPI plots for ", length(modules),
                " WGCNA module(s) ...")
        for (mod in modules) {
            genes <- .ppi_seed_genes_module(
                mod_df, mod, gene_col = gene_col, max_genes = max_genes
            )
            out <- .ppi_result_for_genes(
                genes, string_db, top_n_hubs, max_plot_nodes, title = mod,
                annotate = function(nodes) {
                    .ppi_add_kme(nodes, mod_df, mod, gene_col = gene_col)
                }
            )
            if (is.null(out$plot)) {
                message("  Skipping \"", mod, "\": fewer than 2 mapped genes share ",
                        "a STRING edge (mapped ", out$n_mapped, "/", out$n_seed, ").")
                next
            }
            by_module[[mod]] <- out
            mod_base <- if (!is.null(wgcna_module_save_dir) &&
                            grepl("^Module ", as.character(mod), perl = TRUE)) {
                wgcna_module_save_dir
            } else {
                module_save_dir
            }
            .write_ppi_outputs(out, file.path(mod_base, safe_mod_dir(mod), "ppi"), save_plot)
        }
    }

    list(
        org = org,
        by_comparison = by_comparison,
        by_module = by_module,
        top_n_hubs = top_n_hubs,
        max_genes = max_genes,
        score_threshold = as.integer(score_threshold)[1],
        source = "STRINGdb"
    )
}

#' Remove the top-level `ppi/` folder (LRT-cluster PPI + STRINGdb cache)
#'
#' Called automatically at the end of [run_deseq2()] and
#' [`run_deseq2(mode = "lrt")`][run_deseq2] after HTML rendering. Contrast PPI under
#' `deseq2/<contrast>/ppi/` and WGCNA module PPI under
#' `wgcna/<Module_N>/ppi/` are kept.
#'
#' @param save_dir Pipeline root directory.
#' @return Invisibly `NULL`.
#' @export
cleanup_ppi_dir <- function(save_dir) {
    save_dir <- validate_save_dir(save_dir)
    ppi_dir <- file.path(save_dir, "ppi")
    if (dir.exists(ppi_dir)) {
        unlink(ppi_dir, recursive = TRUE, force = TRUE)
    }
    invisible(NULL)
}
