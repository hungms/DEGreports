#' Load the human B-cell regulon from the bcellViper data package
#' @noRd
.load_bcell_regulon <- function() {
    if (!requireNamespace("bcellViper", quietly = TRUE)) {
        stop("Package \"bcellViper\" is required for VIPER analysis. ",
             "Install with BiocManager::install(\"bcellViper\").", call. = FALSE)
    }
    path <- system.file("data", "bcellViper.rda", package = "bcellViper")
    if (!nzchar(path) || !file.exists(path)) {
        stop("Could not locate bcellViper.rda inside the bcellViper package.",
             call. = FALSE)
    }
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (!exists("regulon", envir = env, inherits = FALSE)) {
        stop("bcellViper.rda did not contain a `regulon` object.", call. = FALSE)
    }
    env$regulon
}

#' Map a gene × sample matrix (or named numeric vector) to human symbols
#'
#' Human input is returned unchanged. Mouse symbols are mapped via
#' `babelgene::orthologs(..., species = "mouse", human = FALSE)`. When
#' several mouse genes map to one human gene, values are averaged.
#'
#' @param x Numeric matrix with gene rownames, or a named numeric vector.
#' @param org `"human"` or `"mouse"` (from [detect_organism()]).
#' @return Matrix or named vector in human gene space.
#' @noRd
.map_expression_to_human <- function(x, org) {
    is_vec <- is.null(dim(x))
    if (is_vec) {
        genes <- names(x)
        mat <- matrix(as.numeric(x), ncol = 1L, dimnames = list(genes, "value"))
    } else {
        mat <- as.matrix(x)
        genes <- rownames(mat)
    }
    if (is.null(genes) || !length(genes)) {
        stop("Expression input must have gene names.", call. = FALSE)
    }

    if (identical(org, "human")) {
        if (is_vec) {
            return(stats::setNames(as.numeric(mat[, 1]), rownames(mat)))
        }
        return(mat)
    }
    if (!identical(org, "mouse")) {
        stop("VIPER gene mapping supports organism \"human\" or \"mouse\" only; got \"",
             org, "\".", call. = FALSE)
    }
    if (!requireNamespace("babelgene", quietly = TRUE)) {
        stop("Package \"babelgene\" is required to map mouse genes to the human ",
             "bcellViper regulon. Install with install.packages(\"babelgene\").",
             call. = FALSE)
    }

    # babelgene rejects mixed symbol/Ensembl inputs; keep MGI-style symbols only.
    mouse_re <- "^[A-Z][A-Za-z0-9]+([._-][A-Za-z0-9]+)*$"
    keep_sym <- grepl(mouse_re, genes) & !grepl("^ENS", genes, ignore.case = TRUE)
    if (!any(keep_sym)) {
        stop("No mouse-style gene symbols available for babelgene ortholog mapping.",
             call. = FALSE)
    }
    if (!all(keep_sym)) {
        message("Dropping ", sum(!keep_sym),
                " non-symbol gene id(s) before mouse→human ortholog mapping.")
    }
    genes_query <- unique(genes[keep_sym])

    ortho <- babelgene::orthologs(genes_query, species = "mouse", human = FALSE)
    if (is.null(ortho) || !nrow(ortho)) {
        stop("babelgene returned no mouse→human orthologs for the input genes.",
             call. = FALSE)
    }
    need <- c("symbol", "human_symbol")
    if (!all(need %in% colnames(ortho))) {
        stop("Unexpected babelgene::orthologs() columns; expected symbol/human_symbol.",
             call. = FALSE)
    }
    ortho <- ortho[
        !is.na(ortho$human_symbol) & nzchar(as.character(ortho$human_symbol)) &
            !is.na(ortho$symbol) & nzchar(as.character(ortho$symbol)),
        , drop = FALSE
    ]
    if ("support_n" %in% colnames(ortho)) {
        ortho <- ortho[order(-ortho$support_n), , drop = FALSE]
    }
    # One human symbol per mouse gene (highest support), then collapse
    # multiple mouse genes that share a human ortholog by mean.
    ortho <- ortho[!duplicated(ortho$symbol), , drop = FALSE]
    ortho <- ortho[ortho$symbol %in% genes, , drop = FALSE]
    if (!nrow(ortho)) {
        stop("No mouse genes overlapped babelgene ortholog mappings.", call. = FALSE)
    }

    mapped_rows <- lapply(unique(ortho$human_symbol), function(hs) {
        mice <- ortho$symbol[ortho$human_symbol == hs]
        sub <- mat[mice, , drop = FALSE]
        if (nrow(sub) == 1L) {
            as.numeric(sub[1, ])
        } else {
            as.numeric(colMeans(sub, na.rm = TRUE))
        }
    })
    out <- do.call(rbind, mapped_rows)
    rownames(out) <- unique(ortho$human_symbol)
    colnames(out) <- colnames(mat)
    storage.mode(out) <- "double"

    if (is_vec) {
        return(stats::setNames(as.numeric(out[, 1]), rownames(out)))
    }
    out
}

#' Convert an msviper result to a tidy data frame
#' @noRd
.msviper_to_df <- function(msv, comparison = NULL) {
    nes <- msv$es$nes
    pval <- msv$es$p.value
    df <- data.frame(
        tf = names(nes),
        nes = as.numeric(nes),
        p.value = as.numeric(pval[names(nes)]),
        stringsAsFactors = FALSE
    )
    df$padj <- stats::p.adjust(df$p.value, method = "BH")
    df <- df[order(df$p.value, -abs(df$nes)), , drop = FALSE]
    if (!is.null(comparison)) df$comparison <- comparison
    rownames(df) <- NULL
    df
}

#' Run VIPER / msVIPER with the bcellViper B-cell regulon
#'
#' Uses the human B-cell transcriptional regulatory network from the
#' Bioconductor [bcellViper](https://bioconductor.org/packages/bcellViper/)
#' data package with [viper](https://califano.c2b2.columbia.edu/viper)
#' ([Alvarez et al., Nat Genet 2016](https://doi.org/10.1038/ng.3593)):
#' \itemize{
#'   \item Sample-level protein activity via `viper()` on scaled VST expression
#'   \item Contrast-level master-regulator analysis via `msviper()` on each
#'     comparison's DESeq2 Wald `stat` vector
#' }
#' Mouse datasets are mapped to human symbols with `babelgene` before
#' scoring against the human regulon.
#'
#' Requires Suggests packages `viper`, `bcellViper`, and (for mouse)
#' `babelgene`.
#'
#' @param dds A QC-filtered `DESeqDataSet`.
#' @param var Grouping column (unused for scoring; kept for API symmetry
#'   with other modules).
#' @param diffexp_out Output of [run_diffexp_pipeline()] (needs
#'   `$by_comparison` with a `stat` column per gene).
#' @param n_tfs Number of top master regulators shown in each contrast
#'   barcode plot (`plot.msviper`). Default `10`.
#' @param n_tfs_heatmap Number of most variable TFs in the sample-level
#'   activity heatmap. Default `50`.
#' @param minsize Minimum regulon size passed to `viper` / `msviper`.
#'   Default `25` (VIPER default); mouse ortholog mapping may leave fewer
#'   targets -- lower if many regulons are filtered out.
#' @param gene_col Gene-id column name in diffexp tables. Default `"gene"`.
#' @param save_dir Output directory. Default `getwd()`.
#' @param save_plot Logical; write PDFs alongside CSVs. Default `TRUE`.
#' @return A list with `org`, `regulon_n`, `sample_acts` (TF × sample
#'   matrix), `sample_heatmap`, `contrast_acts`, `contrast_msviper`,
#'   `n_tfs`, and `n_tfs_heatmap`.
#' @export
run_viper_pipeline <- function(
    dds,
    var,
    diffexp_out,
    n_tfs = 10,
    n_tfs_heatmap = 50,
    minsize = 25,
    gene_col = "gene",
    save_dir = getwd(),
    save_plot = TRUE) {

    if (!requireNamespace("viper", quietly = TRUE)) {
        stop("Package \"viper\" is required for VIPER analysis. ",
             "Install with BiocManager::install(\"viper\").", call. = FALSE)
    }

    dds <- validate_dds(dds)
    save_dir <- validate_save_dir(save_dir)
    minsize_missing <- missing(minsize)
    n_tfs <- as.integer(n_tfs)[1]
    n_tfs_heatmap <- as.integer(n_tfs_heatmap)[1]
    minsize <- as.integer(minsize)[1]
    if (is.na(n_tfs) || n_tfs < 1) {
        stop("`n_tfs` must be a positive integer.", call. = FALSE)
    }
    if (is.na(n_tfs_heatmap) || n_tfs_heatmap < 1) {
        stop("`n_tfs_heatmap` must be a positive integer.", call. = FALSE)
    }
    if (is.na(minsize) || minsize < 1) {
        stop("`minsize` must be a positive integer.", call. = FALSE)
    }

    org <- detect_organism(rownames(dds))
    # Ortholog mapping thins regulon target overlap; use a slightly lower
    # default minsize for mouse unless the caller overrode minsize.
    if (identical(org, "mouse") && minsize_missing) {
        minsize <- 15L
    }
    message("Loading bcellViper human B-cell regulon (organism=\"", org,
            "\", minsize=", minsize, ") ...")
    regulon <- .load_bcell_regulon()

    message("Running sample-level VIPER on scaled VST expression ...")
    vsd <- quiet(safe_vst(dds, blind = TRUE))
    mat <- as.matrix(SummarizedExperiment::assay(vsd))
    storage.mode(mat) <- "double"
    mat_h <- .map_expression_to_human(mat, org)
    # Scale genes across samples (VIPER vignette-style z-scores).
    mat_z <- t(scale(t(mat_h)))
    mat_z[!is.finite(mat_z)] <- 0

    sample_acts <- viper::viper(
        eset = mat_z,
        regulon = regulon,
        method = "scale",
        minsize = minsize,
        eset.filter = TRUE,
        verbose = FALSE
    )
    # viper returns TF × sample matrix of NES.
    write_output(
        data.frame(tf = rownames(sample_acts), sample_acts, check.names = FALSE),
        file.path(save_dir, "viper_sample_activities.csv")
    )

    sample_heatmap <- plot_viper_sample_heatmap(sample_acts, n_tfs = n_tfs_heatmap)
    if (save_plot && !is.null(sample_heatmap)) {
        gt <- sample_heatmap$gtable
        w <- 0.5 * min(grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE), 10)
        h <- 0.5 * min(grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE), 10)
        grDevices::pdf(file.path(save_dir, "viper_sample_tf_heatmap.pdf"),
                       width = max(4, w), height = max(3, h))
        grid::grid.newpage()
        grid::grid.draw(gt)
        grDevices::dev.off()
    }

    by_comp <- diffexp_out$by_comparison
    if (is.null(by_comp) || !length(by_comp)) {
        stop("`diffexp_out$by_comparison` is empty; cannot run contrast-level msVIPER.",
             call. = FALSE)
    }

    message("Running contrast-level msVIPER on DESeq2 Wald statistics (",
            length(by_comp), " contrast(s)) ...")
    contrast_acts <- list()
    contrast_msviper <- list()
    for (comp in names(by_comp)) {
        res <- by_comp[[comp]]
        if (!"stat" %in% colnames(res)) {
            warning("Skipping \"", comp, "\": no `stat` column.", call. = FALSE)
            next
        }
        genes <- as.character(res[[gene_col]])
        ok <- !is.na(res$stat) & !is.na(genes) & nzchar(genes)
        ges <- stats::setNames(as.numeric(res$stat[ok]), genes[ok])
        # Collapse duplicate gene ids (keep first).
        ges <- ges[!duplicated(names(ges))]
        ges_h <- .map_expression_to_human(ges, org)

        msv <- tryCatch(
            viper::msviper(
                ges = ges_h,
                regulon = regulon,
                minsize = minsize,
                ges.filter = TRUE,
                verbose = FALSE
            ),
            error = function(e) {
                warning("msVIPER failed for \"", comp, "\": ",
                        conditionMessage(e), call. = FALSE)
                NULL
            }
        )
        if (is.null(msv)) next

        acts <- .msviper_to_df(msv, comparison = comp)
        contrast_acts[[comp]] <- acts
        contrast_msviper[[comp]] <- msv
        safe <- gsub("[^A-Za-z0-9]+", "_", comp)
        write_output(acts, file.path(save_dir, paste0("msviper_contrast_", safe, ".csv")))

        if (save_plot) {
            pdf_path <- file.path(save_dir, paste0("msviper_contrast_", safe, "_barcode.pdf"))
            # Half the previous barcode footprint (was 8 x max(4, 0.35*n+1.5)).
            grDevices::pdf(pdf_path, width = 4, height = max(2, 0.175 * n_tfs + 0.75))
            tryCatch(
                plot(msv, mrs = n_tfs, cex = 0.55),
                error = function(e) {
                    warning("msVIPER plot failed for \"", comp, "\": ",
                            conditionMessage(e), call. = FALSE)
                }
            )
            grDevices::dev.off()
        }
    }

    list(
        org = org,
        regulon_n = length(regulon),
        sample_acts = sample_acts,
        sample_heatmap = sample_heatmap,
        contrast_acts = contrast_acts,
        contrast_msviper = contrast_msviper,
        n_tfs = n_tfs,
        n_tfs_heatmap = n_tfs_heatmap,
        minsize = minsize,
        network = "bcellViper"
    )
}

#' Sample×TF VIPER activity heatmap (top variable TFs)
#'
#' @param sample_acts TF × sample numeric matrix from `viper::viper()`.
#' @param n_tfs Number of most variable TFs to keep. Default `50`.
#' @return A `pheatmap` object, or `NULL` if fewer than 2 TFs/samples.
#' @export
plot_viper_sample_heatmap <- function(sample_acts, n_tfs = 50) {
    if (is.null(sample_acts) || !length(sample_acts)) return(NULL)
    mat <- as.matrix(sample_acts)
    if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)

    tf_var <- apply(mat, 1L, stats::var, na.rm = TRUE)
    tf_var[!is.finite(tf_var)] <- 0
    keep <- utils::head(names(sort(tf_var, decreasing = TRUE)), n_tfs)
    mat <- mat[keep, , drop = FALSE]
    # Rows = samples, columns = TFs (match TF-activity page convention).
    mat <- t(mat)
    mat <- scale(mat)
    mat[!is.finite(mat)] <- 0

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

#' Emit a VIPER / msVIPER base-graphics barcode plot with a PDF download
#'
#' @param msv An `msviper` object.
#' @param filename Output file name *without* extension.
#' @param mrs Number of master regulators to show. Default `10`.
#' @param width,height Device size in inches.
#' @param dpi Rasterization dpi for the on-page figure.
#' @return Invisibly `NULL`; `cat()`s HTML via [emit_widget()].
#' @export
emit_msviper_plot <- function(msv, filename, mrs = 10,
                              width = 4, height = 2.5, dpi = 150) {
    if (!requireNamespace("htmltools", quietly = TRUE)) {
        stop("Package \"htmltools\" is required to emit VIPER plots.", call. = FALSE)
    }
    if (!requireNamespace("downloadthis", quietly = TRUE)) {
        stop("Package \"downloadthis\" is required for download buttons.", call. = FALSE)
    }

    draw <- function() {
        old_par <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(old_par), add = TRUE)
        plot(msv, mrs = mrs, cex = 0.55)
    }

    tmp_pdf <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp_pdf, width = width, height = height)
    tryCatch(draw(), error = function(e) {
        plot.new()
        graphics::title(main = paste("VIPER plot failed:", conditionMessage(e)))
    })
    grDevices::dev.off()

    button <- downloadthis::download_file(
        path = tmp_pdf,
        output_name = filename,
        button_label = "Download PDF",
        button_type = "default",
        has_icon = TRUE,
        self_contained = TRUE
    )

    px_width <- as.integer(round(width * dpi))
    px_height <- as.integer(round(height * dpi))
    content <- htmltools::plotTag(
        draw(),
        alt = filename,
        width = px_width,
        height = px_height,
        pixelratio = dpi / 72,
        suppressSize = "xy",
        attribs = list(
            width = px_width,
            height = px_height,
            style = sprintf(
                "max-width: none !important; width: %dpx; height: %dpx;",
                px_width, px_height
            )
        )
    )
    emit_widget(wrap_with_corner_button(content, button))
}
