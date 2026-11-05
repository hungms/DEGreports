#' Run CollecTRI + ULM transcription-factor activity analysis (decoupleR)
#'
#' Follows the bulk RNA-seq TF activity vignette for
#' [`decoupleR`](https://bioconductor.org/packages/decoupleR/): retrieves the
#' CollecTRI network from OmniPath, then runs Univariate Linear Models
#' (`run_ulm`) on (1) sample-level VST (or batch-corrected VST when `batch`
#' is set) and (2) each contrast's DESeq2 Wald `stat` vector.
#' DESeq2 Wald `stat` vector. Produces a sample×TF activity heatmap (top
#' variable TFs), per-contrast TF activity barplots, and writes CSVs under
#' `save_dir`.
#'
#' Requires Imports packages `decoupleR` and `OmnipathR`.
#'
#' @param dds A QC-filtered `DESeqDataSet`, or `NULL` to skip sample-level
#'   ULM and run **contrast-level** scoring only (needs
#'   `diffexp_out$by_comparison` with a `stat` column).
#' @param var Grouping column in `colData(dds)` used only for heatmap
#'   context when `dds` is supplied. Ignored when `dds = NULL`. Required
#'   as the protected group when `batch` is set.
#' @param batch Optional batch column(s), including the `pair` column of a
#'   paired design. When set (and `dds` is supplied),
#'   sample-level ULM uses batch-corrected expression. Contrast-level ULM
#'   still uses DESeq2 `stat` from the design (batch as a covariate, not
#'   corrected counts).
#' @param diffexp_out Output of [run_diffexp_pipeline()] or
#'   [prepare_diffexp_from_deg()] (needs `$by_comparison` with a `stat`
#'   column per gene).
#' @param n_tfs Number of top TFs to show in each contrast barplot.
#'   Default `25`.
#' @param n_tfs_heatmap Number of most variable TFs in the sample-level
#'   activity heatmap. Default `50`. Ignored when `dds = NULL`.
#' @param minsize Minimum number of targets per TF (passed to
#'   `decoupleR::run_ulm`). Default `5`.
#' @param gene_col Gene-id column name in diffexp tables. Default `"gene"`.
#' @param collectri_file Optional path to a local CollecTRI CSV with columns
#'   `source`, `target`, `mor`. Default `NULL` loads the date-stamped network
#'   shipped under `inst/extdata/` (or downloads it from the GitHub repo),
#'   then falls back to OmniPath only if that fails.
#' @param save_dir Output directory. Default `getwd()`.
#' @param save_plot Logical; write PDFs alongside CSVs. Default `TRUE`.
#' @return A list with `net`, `sample_acts`, `contrast_acts` (named list),
#'   `sample_heatmap`, `contrast_barplots` (named list), and `org`.
#' @export
run_decoupler_tf_pipeline <- function(
    dds = NULL,
    var = NULL,
    diffexp_out,
    n_tfs = 25,
    n_tfs_heatmap = 50,
    minsize = 5,
    gene_col = "gene",
    collectri_file = NULL,
    batch = NULL,
    save_dir = getwd(),
    save_plot = TRUE) {

    if (!requireNamespace("decoupleR", quietly = TRUE)) {
        stop("Package \"decoupleR\" is required for TF activity analysis. ",
             "Install with BiocManager::install(\"decoupleR\").", call. = FALSE)
    }

    if (!is.null(dds)) dds <- validate_dds(dds)
    save_dir <- validate_save_dir(save_dir)
    n_tfs <- as.integer(n_tfs)[1]
    n_tfs_heatmap <- as.integer(n_tfs_heatmap)[1]
    if (is.na(n_tfs) || n_tfs < 1) {
        stop("`n_tfs` must be a positive integer.", call. = FALSE)
    }
    if (is.na(n_tfs_heatmap) || n_tfs_heatmap < 1) {
        stop("`n_tfs_heatmap` must be a positive integer.", call. = FALSE)
    }

    by_comp <- diffexp_out$by_comparison
    if (is.null(by_comp) || !length(by_comp)) {
        stop("`diffexp_out$by_comparison` is empty; cannot run contrast-level ULM.",
             call. = FALSE)
    }

    if (!is.null(dds)) {
        org <- detect_organism(rownames(dds))
    } else {
        sample_genes <- unlist(lapply(by_comp, function(df) {
            as.character(df[[gene_col]])
        }), use.names = FALSE)
        org <- detect_organism(sample_genes)
    }
    message("Fetching CollecTRI network for organism=\"", org, "\" ...")
    net <- .get_collectri_network(org, file = collectri_file)
    write_output(net, file.path(save_dir, "collectri_network.csv"))

    sample_acts <- NULL
    sample_heatmap <- NULL
    expr_correction <- NA_character_
    if (!is.null(dds)) {
        if (!is.null(batch) && is.null(var)) {
            stop("`var` must be supplied when `batch` is set so the batch correction can protect the biological group.",
                 call. = FALSE)
        }
        message("Running sample-level ULM on sample expression ...")
        mat <- sample_level_expr_matrix(dds, var = var, batch = batch)
        expr_correction <- attr(mat, "expr_correction")
        storage.mode(mat) <- "double"
        sample_acts <- as.data.frame(decoupleR::run_ulm(
            mat = mat, net = net,
            .source = "source", .target = "target", .mor = "mor",
            minsize = minsize
        ))
        write_output(sample_acts, file.path(save_dir, "ulm_sample_activities.csv"))

        sample_heatmap <- plot_tf_sample_heatmap(sample_acts, n_tfs = n_tfs_heatmap)
        if (save_plot && !is.null(sample_heatmap)) {
            gt <- sample_heatmap$gtable
            w <- 0.5 * min(grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE), 10)
            h <- 0.5 * min(grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE), 10)
            grDevices::pdf(file.path(save_dir, "ulm_sample_tf_heatmap.pdf"), width = w, height = h)
            grid::grid.newpage()
            grid::grid.draw(gt)
            grDevices::dev.off()
        }
    } else {
        message("Skipping sample-level ULM (dds = NULL; contrast-level only).")
    }

    message("Running contrast-level ULM on DESeq2 Wald statistics (",
            length(by_comp), " contrast(s)) ...")
    contrast_acts <- list()
    contrast_barplots <- list()
    for (comp in names(by_comp)) {
        res <- by_comp[[comp]]
        if (!"stat" %in% colnames(res)) {
            warning("Skipping \"", comp, "\": no `stat` column.", call. = FALSE)
            next
        }
        genes <- as.character(res[[gene_col]])
        ok <- !is.na(res$stat) & !is.na(genes) & nzchar(genes)
        deg_mat <- matrix(res$stat[ok], ncol = 1L,
                          dimnames = list(genes[ok], "stat"))
        acts <- as.data.frame(decoupleR::run_ulm(
            mat = deg_mat, net = net,
            .source = "source", .target = "target", .mor = "mor",
            minsize = minsize
        ))
        acts$padj <- stats::p.adjust(acts$p_value, method = "BH")
        acts$comparison <- comp
        contrast_acts[[comp]] <- acts
        # Per-comparison subfolder keyed by the raw contrast label (same
        # names as deseq2/<label>/).
        comp_dir <- file.path(save_dir, comp)
        write_output(acts, file.path(comp_dir, "ulm_contrast.csv"))
        bp <- plot_tf_contrast_barplot(acts, n_tfs = n_tfs, title = comp)
        contrast_barplots[[comp]] <- bp
        if (save_plot && !is.null(bp)) {
            save_ggplot(bp, file.path(comp_dir, "ulm_contrast_barplot.pdf"),
                        width = 4, height = 2.25)
        }
    }

    list(
        org = org,
        net = net,
        sample_acts = sample_acts,
        contrast_acts = contrast_acts,
        sample_heatmap = sample_heatmap,
        contrast_barplots = contrast_barplots,
        n_tfs = n_tfs,
        n_tfs_heatmap = n_tfs_heatmap,
        expr_correction = expr_correction
    )
}

#' Read a local CollecTRI CSV (source/target/mor)
#' @noRd
.read_collectri_file <- function(file) {
    file <- normalizePath(file, mustWork = TRUE)
    net <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    need <- c("source", "target", "mor")
    if (!all(need %in% colnames(net))) {
        stop("CollecTRI file must have columns source, target, mor: ", file,
             call. = FALSE)
    }
    net <- net[, need, drop = FALSE]
    net$mor <- as.numeric(net$mor)
    net <- net[!is.na(net$mor) & nzchar(net$source) & nzchar(net$target), ,
               drop = FALSE]
    unique(net)
}

#' Fetch CollecTRI as source/target/mor, with local + OmniPath fallbacks
#'
#' Order: explicit `file` → bundled/GitHub `inst/extdata` CSV →
#' [decoupleR::get_collectri()] → OmniPath REST API.
#' @noRd
.get_collectri_network <- function(organism = c("human", "mouse"),
                                   file = NULL) {
    organism <- match.arg(organism)

    if (!is.null(file)) {
        message("Loading CollecTRI from ", file, " ...")
        return(.read_collectri_file(file))
    }
    local <- tryCatch(
        .collectri_default_path(organism, download = TRUE),
        error = function(e) {
            message("Bundled CollecTRI unavailable (", conditionMessage(e), ").")
            NULL
        }
    )
    if (!is.null(local)) {
        message("Loading CollecTRI from ", local, " ...")
        return(.read_collectri_file(local))
    }

    if (!requireNamespace("OmnipathR", quietly = TRUE)) {
        message("Package \"OmnipathR\" not installed; skipping ",
                "decoupleR::get_collectri() and using OmniPath REST ...")
    } else {
        net <- tryCatch(
            as.data.frame(decoupleR::get_collectri(
                organism = organism, split_complexes = FALSE
            )),
            error = function(e) {
                message("decoupleR::get_collectri() failed (",
                        conditionMessage(e),
                        "); falling back to OmniPath REST API ...")
                NULL
            }
        )
        if (!is.null(net) && all(c("source", "target", "mor") %in% colnames(net))) {
            net <- net[, c("source", "target", "mor"), drop = FALSE]
            return(unique(net))
        }
    }

    tax <- if (identical(organism, "mouse")) 10090L else 9606L
    url <- paste0(
        "https://omnipathdb.org/interactions?",
        "genesymbols=1&resources=CollecTRI&organisms=", tax
    )
    tmp <- tempfile(fileext = ".tsv")
    utils::download.file(url, tmp, quiet = TRUE)
    raw <- utils::read.delim(tmp, check.names = FALSE, stringsAsFactors = FALSE)
    unlink(tmp)
    if (!nrow(raw)) {
        stop("CollecTRI download returned 0 rows from OmniPath (", url, ").",
             call. = FALSE)
    }
    stim <- tolower(as.character(raw$is_stimulation)) %in% c("true", "1", "t")
    inh <- tolower(as.character(raw$is_inhibition)) %in% c("true", "1", "t")
    mor <- ifelse(stim & !inh, 1L, ifelse(inh & !stim, -1L, NA_integer_))
    # Ambiguous edges: prefer consensus flags when present.
    if ("consensus_stimulation" %in% colnames(raw)) {
        cstim <- tolower(as.character(raw$consensus_stimulation)) %in% c("true", "1", "t")
        cinh <- tolower(as.character(raw$consensus_inhibition)) %in% c("true", "1", "t")
        amb <- is.na(mor)
        mor[amb & cstim & !cinh] <- 1L
        mor[amb & cinh & !cstim] <- -1L
    }
    net <- data.frame(
        source = raw$source_genesymbol,
        target = raw$target_genesymbol,
        mor = as.numeric(mor),
        stringsAsFactors = FALSE
    )
    net <- net[!is.na(net$mor) & nzchar(net$source) & nzchar(net$target), , drop = FALSE]
    # Collapse duplicate source→target edges (ULM rejects repeated edges).
    if (nrow(net)) {
        net <- stats::aggregate(mor ~ source + target, data = net, FUN = function(x) {
            # Prefer a non-zero signed consensus; otherwise mean.
            sx <- unique(sign(x[x != 0]))
            if (length(sx) == 1) as.numeric(sx) else mean(x)
        })
    }
    if (!nrow(net)) {
        stop("Failed to build a non-empty CollecTRI network.", call. = FALSE)
    }
    net
}

#' Sample×TF activity heatmap (top variable TFs)
#'
#' @param sample_acts Long data frame from `decoupleR::run_ulm` (columns
#'   `source`, `condition`, `score`).
#' @param n_tfs Number of most variable TFs to keep. Default `50`.
#' @return A `pheatmap` object, or `NULL` if fewer than 2 TFs/samples.
#' @export
plot_tf_sample_heatmap <- function(sample_acts, n_tfs = 50) {
    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("Package \"tidyr\" is required for TF heatmaps.", call. = FALSE)
    }
    if (is.null(sample_acts) || nrow(sample_acts) == 0) return(NULL)

    tf_sd <- stats::aggregate(score ~ source, data = sample_acts, FUN = stats::sd)
    tf_sd <- tf_sd[order(-abs(tf_sd$score)), , drop = FALSE]
    tfs <- utils::head(tf_sd$source, n_tfs)
    sub <- sample_acts[sample_acts$source %in% tfs, , drop = FALSE]
    wide <- tidyr::pivot_wider(sub, id_cols = "condition",
                               names_from = "source", values_from = "score")
    mat <- as.matrix(wide[, -1, drop = FALSE])
    rownames(mat) <- wide$condition
    if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)

    mat <- scale(mat)
    mat[is.na(mat)] <- 0
    palette_length <- 100L
    my_color <- grDevices::colorRampPalette(c("Darkblue", "white", "red"))(palette_length)
    my_breaks <- c(
        seq(-3, 0, length.out = ceiling(palette_length / 2) + 1),
        seq(0.05, 3, length.out = floor(palette_length / 2))
    )
    pheatmap::pheatmap(
        mat,
        color = my_color,
        breaks = my_breaks,
        border_color = NA,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        silent = TRUE
    )
}

#' Top TF activity barplot for one contrast (ULM on DE stats)
#'
#' @param contrast_acts Long ULM result for a single contrast.
#' @param n_tfs Number of top TFs (both signs) to show. Default `25`.
#' @param title Optional plot title.
#' @return A `ggplot` object, or `NULL` if empty.
#' @export
plot_tf_contrast_barplot <- function(contrast_acts, n_tfs = 25, title = NULL) {
    if (is.null(contrast_acts) || nrow(contrast_acts) == 0) return(NULL)
    df <- contrast_acts
    df$rnk <- NA_real_
    pos <- df$score > 0
    if (any(pos)) df$rnk[pos] <- rank(-df$score[pos])
    if (any(!pos)) df$rnk[!pos] <- rank(-abs(df$score[!pos]))
    df <- df[order(df$rnk), , drop = FALSE]
    df <- utils::head(df, n_tfs)
    df$source <- factor(df$source, levels = df$source[order(df$score)])

    ggplot2::ggplot(df, ggplot2::aes(x = .data$source, y = .data$score, fill = .data$score)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_gradient2(
            low = "darkblue", mid = "whitesmoke", high = "indianred", midpoint = 0
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(face = "bold"),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
            axis.text.y = ggplot2::element_text(face = "bold"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
        ) +
        ggplot2::xlab("TFs") +
        ggplot2::ylab("ULM score") +
        ggplot2::ggtitle(title)
}
