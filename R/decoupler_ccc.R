#' Run ligand–receptor (cell–cell communication) ULM with decoupleR
#'
#' Follows the bulk RNA-seq
#' [ligand–receptor section](https://decoupler.readthedocs.io/en/v1.9.2/notebooks/bulk.html)
#' of the decoupler tutorial: LIANA-style ligand–receptor pairs are formatted
#' as a prior-knowledge network (`source` = interaction, `target` = ligand or
#' receptor gene, `mor` = 1) and scored with Univariate Linear Models
#' (`decoupleR::run_ulm`, `minsize = 2`) on each contrast's DESeq2 Wald
#' `stat` vector. This is independent of CollecTRI TF activity
#' ([run_decoupler_tf_pipeline()]).
#'
#' Pass `cluster_diffexp` to additionally score each **degPatterns cluster**
#' of an LRT separately, the way [run_decoupler_tf_pipeline()] and
#' [run_viper_pipeline()] already treat temporal clusters.
#'
#' Requires Import `decoupleR`. The L–R resource is loaded from (in order)
#' an optional local CSV, [liana](https://saezlab.github.io/liana/) Consensus
#' if installed, or OmniPath ligand–receptor extra interactions.
#'
#' @param diffexp_out Output of [run_diffexp_pipeline()] or
#'   [prepare_diffexp_from_deg()] (needs `$by_comparison` with a `stat`
#'   column per gene).
#' @param cluster_diffexp Optional output of [lrt_clusters_diffexp_for_tf()],
#'   i.e. one signature per degPatterns cluster (`"<lrt> cluster <n>"`). Each
#'   cluster is then scored on its own, so a pair is only reported for a
#'   cluster when **both** partners are members of it — the whole-test
#'   ULM above cannot say which temporal pattern a pair belongs to. Written
#'   to `<lrt>/cluster_<n>/`. Default `NULL` (no cluster-level scoring).
#' @param n_lr Number of top / bottom ligand–receptor pairs per direction
#'   in the GSEA-style barplot. Default `10`.
#' @param minsize Minimum unique genes per interaction (passed to
#'   `decoupleR::run_ulm`). Default `2` (a ligand and a receptor).
#' @param gene_col Gene-id column name in diffexp tables. Default `"gene"`.
#' @param lr_file Optional path to a local CSV with columns `source`,
#'   `target`, `mor` (interaction, gene, weight). Default `NULL` fetches
#'   LIANA / OmniPath.
#' @param lr_net Optional pre-built network data frame (`source`, `target`,
#'   `mor`). When supplied, `lr_file` and remote fetches are skipped
#'   (used in tests).
#' @param save_dir Output directory. Default `getwd()`.
#' @param save_plot Logical; write PDFs alongside CSVs. Default `TRUE`.
#' @return A list with `org`, `net`, `contrast_acts`, `contrast_barplots`,
#'   `n_lr`, `minsize`, and — when `cluster_diffexp` is supplied —
#'   `cluster_acts`, `cluster_barplots` (both keyed by cluster label) and
#'   `clusters_by_test` (test label -> its cluster labels, including
#'   clusters that scored no pair).
#' @export
run_decoupler_ccc_pipeline <- function(
    diffexp_out,
    cluster_diffexp = NULL,
    n_lr = 10,
    minsize = 2,
    gene_col = "gene",
    lr_file = NULL,
    lr_net = NULL,
    save_dir = getwd(),
    save_plot = TRUE) {

    if (!requireNamespace("decoupleR", quietly = TRUE)) {
        stop("Package \"decoupleR\" is required for ligand–receptor analysis. ",
             "Install with BiocManager::install(\"decoupleR\").", call. = FALSE)
    }

    save_dir <- validate_save_dir(save_dir)
    n_lr <- as.integer(n_lr)[1]
    minsize <- as.integer(minsize)[1]
    if (is.na(n_lr) || n_lr < 1L) {
        stop("`n_lr` must be a positive integer.", call. = FALSE)
    }
    if (is.na(minsize) || minsize < 1L) {
        stop("`minsize` must be a positive integer.", call. = FALSE)
    }
    validate_logical(save_plot, "save_plot")

    by_comp <- diffexp_out$by_comparison
    if (is.null(by_comp) || !length(by_comp)) {
        stop("`diffexp_out$by_comparison` is empty; cannot run ligand–receptor ULM.",
             call. = FALSE)
    }

    sample_genes <- unlist(lapply(by_comp, function(df) {
        as.character(df[[gene_col]])
    }), use.names = FALSE)
    org <- detect_organism(sample_genes)

    if (!is.null(lr_net)) {
        net <- .validate_lr_net(lr_net)
    } else {
        message("Fetching ligand–receptor network for organism=\"", org, "\" ...")
        net <- .get_liana_lr_network(org, file = lr_file)
    }
    write_output(net, file.path(save_dir, "liana_lr_network.csv"))

    message("Running contrast-level ULM on ligand–receptor pairs (",
            length(by_comp), " contrast(s)) ...")
    contrast_acts <- list()
    contrast_barplots <- list()
    for (comp in names(by_comp)) {
        acts <- .ccc_ulm_signature(by_comp[[comp]], net, minsize = minsize,
                                   gene_col = gene_col, label = comp)
        if (is.null(acts)) next
        contrast_acts[[comp]] <- acts
        comp_dir <- file.path(save_dir, comp)
        write_output(acts, file.path(comp_dir, "ulm_lr.csv"))
        bp <- plot_lr_contrast_barplot(acts, n_lr = n_lr, title = comp)
        contrast_barplots[[comp]] <- bp
        .ccc_save_barplot(bp, comp_dir, save_plot)
    }

    cluster_acts <- list()
    cluster_barplots <- list()
    clusters_by_test <- list()
    by_cluster <- cluster_diffexp$by_comparison
    if (!is.null(cluster_diffexp) && !length(by_cluster)) {
        warning("`cluster_diffexp` has no `by_comparison` signatures; ",
                "skipping cluster-level ligand–receptor ULM.", call. = FALSE)
    }
    if (length(by_cluster)) {
        message("Running cluster-level ULM on ligand–receptor pairs (",
                length(by_cluster), " cluster(s)) ...")
        for (label in names(by_cluster)) {
            loc <- .ccc_cluster_location(label)
            clusters_by_test[[loc$test]] <- c(clusters_by_test[[loc$test]], label)
            acts <- .ccc_ulm_signature(by_cluster[[label]], net,
                                       minsize = minsize, gene_col = gene_col,
                                       label = label)
            if (is.null(acts)) {
                message("  ", label, ": no pair has both partners in this cluster.")
                next
            }
            cluster_acts[[label]] <- acts
            cl_dir <- file.path(save_dir, loc$test, loc$dir)
            write_output(acts, file.path(cl_dir, "ulm_lr.csv"))
            bp <- plot_lr_contrast_barplot(acts, n_lr = n_lr, title = label)
            cluster_barplots[[label]] <- bp
            .ccc_save_barplot(bp, cl_dir, save_plot)
        }
    }

    list(
        org = org,
        net = net,
        contrast_acts = contrast_acts,
        contrast_barplots = contrast_barplots,
        cluster_acts = cluster_acts,
        cluster_barplots = cluster_barplots,
        clusters_by_test = clusters_by_test,
        n_lr = n_lr,
        minsize = minsize
    )
}

#' Score one gene-level signature against the ligand–receptor network
#'
#' Returns `NULL` — with a warning for a malformed signature or a failed
#' fit — rather than stopping, so one unscorable cluster does not take the
#' whole module down. A cluster signature holds only that cluster's genes,
#' so `minsize = 2` keeps just the pairs whose ligand *and* receptor are
#' cluster members.
#' @noRd
.ccc_ulm_signature <- function(res, net, minsize, gene_col, label) {
    if (is.null(res) || !nrow(res)) return(NULL)
    missing_cols <- setdiff(c("stat", gene_col), colnames(res))
    if (length(missing_cols)) {
        warning("Skipping \"", label, "\": no `",
                paste(missing_cols, collapse = "`/`"), "` column.",
                call. = FALSE)
        return(NULL)
    }
    genes <- as.character(res[[gene_col]])
    ok <- !is.na(res$stat) & !is.na(genes) & nzchar(genes)
    if (!any(ok)) return(NULL)
    deg_mat <- matrix(res$stat[ok], ncol = 1L,
                      dimnames = list(genes[ok], "stat"))
    # `run_ulm()` errors when no source keeps `minsize` targets, which for a
    # cluster signature is the ordinary case rather than a fault -- most
    # pairs are split across clusters. Check first so only real failures warn.
    per_source <- table(net$source[net$target %in% rownames(deg_mat)])
    if (!any(per_source >= minsize)) return(NULL)
    acts <- tryCatch(
        as.data.frame(decoupleR::run_ulm(
            mat = deg_mat, net = net,
            .source = "source", .target = "target", .mor = "mor",
            minsize = minsize
        )),
        error = function(e) {
            warning("Ligand–receptor ULM failed for \"", label, "\": ",
                    conditionMessage(e), call. = FALSE)
            NULL
        }
    )
    if (is.null(acts) || !nrow(acts)) return(NULL)
    acts$comparison <- label
    if ("p_value" %in% colnames(acts)) {
        acts$padj <- stats::p.adjust(acts$p_value, method = "BH")
    }
    acts
}

#' Write a ligand–receptor barplot at its suggested size
#' @noRd
.ccc_save_barplot <- function(bp, out_dir, save_plot) {
    if (!isTRUE(save_plot) || is.null(bp)) return(invisible(NULL))
    sz <- attr(bp, "plot_size")
    save_ggplot(bp, file.path(out_dir, "ulm_lr_barplot.pdf"),
                width = if (!is.null(sz)) sz$width else 4.2,
                height = if (!is.null(sz)) sz$height else 3.2)
    invisible(NULL)
}

#' Split a `"<lrt> cluster <n>"` label into its test and output folder
#'
#' Cluster outputs nest inside the parent test's folder (`<lrt>/cluster_2/`),
#' as the cluster ORA and CSPA outputs do, so the test's report page can link
#' them with a relative path.
#' @noRd
.ccc_cluster_location <- function(label) {
    label <- as.character(label)[1]
    if (!grepl(" cluster ", label, fixed = TRUE)) {
        return(list(test = label, dir = "cluster_all"))
    }
    id <- sub("^.* cluster ", "", label)
    list(
        test = sub(" cluster .*", "", label),
        dir = paste0("cluster_", gsub("[^A-Za-z0-9_]+", "_", id))
    )
}

#' Validate source/target/mor ligand–receptor network
#' @noRd
.validate_lr_net <- function(net) {
    net <- as.data.frame(net, stringsAsFactors = FALSE)
    need <- c("source", "target", "mor")
    if (!all(need %in% colnames(net))) {
        stop("Ligand–receptor network must have columns source, target, mor.",
             call. = FALSE)
    }
    net <- net[, need, drop = FALSE]
    net$source <- as.character(net$source)
    net$target <- as.character(net$target)
    net$mor <- as.numeric(net$mor)
    net <- net[!is.na(net$mor) & nzchar(net$source) & nzchar(net$target), ,
               drop = FALSE]
    net <- unique(net)
    if (!nrow(net)) {
        stop("Ligand–receptor network is empty after filtering.", call. = FALSE)
    }
    net
}

#' Split complex members (LIANA / OmniPath style)
#' @noRd
.explode_complex_genes <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(character(0))
    parts <- unlist(strsplit(x, "[:_]", perl = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    unique(parts[nzchar(parts)])
}

#' Build decoupleR net from ligand / receptor pair columns
#' @noRd
.lr_pairs_to_net <- function(ligand, receptor) {
    ligand <- as.character(ligand)
    receptor <- as.character(receptor)
    ok <- !is.na(ligand) & !is.na(receptor) & nzchar(ligand) & nzchar(receptor)
    ligand <- ligand[ok]
    receptor <- receptor[ok]
    if (!length(ligand)) {
        stop("No ligand–receptor pairs to format.", call. = FALSE)
    }
    rows <- vector("list", length(ligand))
    for (i in seq_along(ligand)) {
        lig <- .explode_complex_genes(ligand[[i]])
        rec <- .explode_complex_genes(receptor[[i]])
        genes <- unique(c(lig, rec))
        if (length(genes) < 2L) next
        interaction <- paste(ligand[[i]], receptor[[i]], sep = "&")
        rows[[i]] <- data.frame(
            source = interaction,
            target = genes,
            mor = 1,
            stringsAsFactors = FALSE
        )
    }
    net <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
    if (is.null(net) || !nrow(net)) {
        stop("Failed to build a ligand–receptor network (need ≥2 genes per pair).",
             call. = FALSE)
    }
    unique(net)
}

#' Map human LR network genes to mouse symbols
#' @noRd
.map_lr_net_to_org <- function(net, org) {
    net <- .validate_lr_net(net)
    if (is.null(org) || is.na(org) || identical(org, "human")) {
        return(net)
    }
    mapped <- .map_human_symbols_to_org(unique(net$target), org)
    if (!length(mapped)) {
        warning("Could not map ligand–receptor genes to organism=\"", org,
                "\"; using human symbols.", call. = FALSE)
        return(net)
    }
    # babelgene returns mouse symbols only; rebuild by mapping each target.
    if (!requireNamespace("babelgene", quietly = TRUE) ||
            !identical(org, "mouse")) {
        return(net)
    }
    ortho <- tryCatch(
        babelgene::orthologs(unique(net$target), species = "mouse", human = TRUE),
        error = function(e) NULL
    )
    if (is.null(ortho) || !nrow(ortho)) return(net)
    human_col <- intersect(c("human_symbol", "human_ensembl"), colnames(ortho))
    mouse_col <- intersect(c("symbol", "ensembl"), colnames(ortho))
    if (!length(human_col) || !length(mouse_col)) return(net)
    lut <- stats::setNames(
        as.character(ortho[[mouse_col[[1]]]]),
        as.character(ortho[[human_col[[1]]]])
    )
    net$target <- unname(lut[net$target])
    net <- net[!is.na(net$target) & nzchar(net$target), , drop = FALSE]
    # Drop interactions that lost all but one gene after mapping.
    n_per <- table(net$source)
    keep <- names(n_per)[as.integer(n_per) >= 2L]
    net[net$source %in% keep, , drop = FALSE]
}

#' Fetch LIANA Consensus or OmniPath ligand–receptor extra as ULM net
#' @noRd
.get_liana_lr_network <- function(organism = c("human", "mouse"),
                                  file = NULL) {
    organism <- match.arg(organism)

    if (!is.null(file)) {
        message("Loading ligand–receptor network from ", file, " ...")
        return(.validate_lr_net(utils::read.csv(
            file, stringsAsFactors = FALSE, check.names = FALSE
        )))
    }

    net <- NULL
    if (requireNamespace("liana", quietly = TRUE)) {
        net <- tryCatch({
            res <- liana::select_resource("Consensus")
            df <- if (is.list(res) && !is.data.frame(res)) res[[1]] else res
            df <- as.data.frame(df)
            lig_col <- intersect(
                c("ligand", "ligand_complex", "source_genesymbol"),
                colnames(df)
            )
            rec_col <- intersect(
                c("receptor", "receptor_complex", "target_genesymbol"),
                colnames(df)
            )
            if (!length(lig_col) || !length(rec_col)) {
                stop("liana Consensus resource lacks ligand/receptor columns.")
            }
            .lr_pairs_to_net(df[[lig_col[[1]]]], df[[rec_col[[1]]]])
        }, error = function(e) {
            message("liana::select_resource() failed (",
                    conditionMessage(e), "); trying OmniPath ...")
            NULL
        })
    }

    if (is.null(net) && requireNamespace("OmnipathR", quietly = TRUE)) {
        net <- tryCatch({
            raw <- as.data.frame(OmnipathR::import_ligrecextra_interactions())
            lig_col <- intersect(
                c("source_genesymbol", "source"), colnames(raw)
            )
            rec_col <- intersect(
                c("target_genesymbol", "target"), colnames(raw)
            )
            if (!length(lig_col) || !length(rec_col)) {
                stop("OmniPath ligrecextra table lacks gene-symbol columns.")
            }
            .lr_pairs_to_net(raw[[lig_col[[1]]]], raw[[rec_col[[1]]]])
        }, error = function(e) {
            message("OmnipathR::import_ligrecextra_interactions() failed (",
                    conditionMessage(e), ").")
            NULL
        })
    }

    if (is.null(net) || !nrow(net)) {
        stop("Could not fetch a ligand–receptor resource. Install liana ",
             "(https://saezlab.github.io/liana/) or OmnipathR, or pass ",
             "`lr_file` / `lr_net`.", call. = FALSE)
    }

    net <- .validate_lr_net(net)
    .map_lr_net_to_org(net, organism)
}

#' Ligand–receptor ULM barplot (GSEA NES style)
#'
#' Top / bottom `n_lr` interactions by ULM `score`, with BH-adjusted
#' `padj` significance stars, matching [plot_gsea_barplot()].
#'
#' @param contrast_acts Long ULM result for one contrast (`source`, `score`,
#'   and `p_value` or `padj`).
#' @param n_lr Pairs per direction. Default `10`.
#' @param title Optional plot title (contrast label).
#' @return A `ggplot` object with `attr(,"plot_size")`, or `NULL` if empty.
#' @export
plot_lr_contrast_barplot <- function(contrast_acts, n_lr = 10, title = NULL) {
    if (is.null(contrast_acts) || !nrow(contrast_acts)) return(NULL)
    df <- contrast_acts
    df$pathway <- as.character(df$source)
    df$NES <- as.numeric(df$score)
    if (!"padj" %in% colnames(df)) {
        pv <- if ("p_value" %in% colnames(df)) df$p_value else NA_real_
        df$padj <- stats::p.adjust(pv, method = "BH")
    }
    p <- plot_gsea_barplot(df, top_n = n_lr, title = title)
    sz <- attr(p, "plot_size")
    p <- p +
        ggplot2::ylab("ULM score") +
        ggplot2::labs(fill = "score")
    attr(p, "plot_size") <- sz
    p
}

#' Genes belonging to significant ligand–receptor pairs
#' @noRd
.ccc_sig_genes <- function(ccc_out, contrast, alpha = 0.05, min_score = 1) {
    if (is.null(ccc_out) || is.null(ccc_out$contrast_acts)) {
        return(character(0))
    }
    acts <- ccc_out$contrast_acts[[contrast]]
    if (is.null(acts) && grepl(" cluster ", contrast, fixed = TRUE)) {
        parent <- sub(" cluster .*", "", contrast)
        acts <- ccc_out$contrast_acts[[parent]]
    }
    if (is.null(acts) || !nrow(acts)) return(character(0))
    pcol <- if ("padj" %in% colnames(acts)) "padj" else
        if ("p_value" %in% colnames(acts)) "p_value" else NULL
    if (is.null(pcol) || !"source" %in% colnames(acts)) return(character(0))
    keep <- !is.na(acts[[pcol]]) & acts[[pcol]] < alpha
    if ("score" %in% colnames(acts)) keep <- keep & abs(acts$score) > min_score
    src <- as.character(acts$source[keep])
    if (!length(src)) return(character(0))
    net <- ccc_out$net
    if (is.null(net) || !nrow(net)) return(character(0))
    unique(as.character(net$target[net$source %in% src]))
}
