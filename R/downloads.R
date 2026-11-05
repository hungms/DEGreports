#' Link to an on-disk output file from a folder-bundled HTML report
#'
#' Emits a plain relative `<a href="...">` styled like the download
#' buttons, for CSVs/PDFs the pipeline already wrote beside the HTML
#' (e.g. `deseq2/<label>/<label>_diffexp.csv`). Prefer this over
#' [download_table_button()] for genome-scale tables so the HTML does
#' not embed a second full copy of the data.
#'
#' @param path Relative path from the HTML page directory to the file
#'   (e.g. `"WT__CC_vs_CB/WT__CC_vs_CB_diffexp.csv"`).
#' @param label Link text. Default `"Download CSV"`.
#' @return An `htmltools` tag suitable for [emit_widget()] /
#'   [emit_widget_with_download()].
#' @examples
#' \dontrun{
#' link_output_file("IgA_vs_IgM/IgA_vs_IgM_diffexp.csv")
#' }
#' @export
link_output_file <- function(path, label = "Download CSV") {
    if (!requireNamespace("htmltools", quietly = TRUE)) {
        stop("Package \"htmltools\" is required for file links.", call. = FALSE)
    }
    path <- as.character(path)[1]
    if (is.na(path) || !nzchar(path)) {
        stop("`path` must be a non-empty character string.", call. = FALSE)
    }
    htmltools::tags$a(
        href = path,
        class = "btn btn-default btn-sm",
        download = basename(path),
        htmltools::tags$i(class = "fa fa-download"),
        label
    )
}

#' Wrap a download button for consistent flex-row alignment
#' @noRd
.normalize_dl_button <- function(button) {
    htmltools::tags$span(class = "degreports-dl-btn-wrap", button)
}

#' Create a download button for a data frame (CSV)
#'
#' Thin wrapper around `downloadthis::download_this()` that embeds a
#' data frame as a self-contained "Download CSV" button. Prefer
#' [link_output_file()] for large tables that already exist on disk in
#' the results folder; use this helper for small / ad-hoc frames that
#' are not written as separate CSVs.
#'
#' @param df A data frame (or an object coercible to one) to offer as a
#'   CSV download.
#' @param filename Output file name *without* extension, e.g.
#'   `"condition_A_vs_B_deg_results"`. Default `"data"`.
#' @param button_label Button text. Default `"Download CSV"`.
#' @return An `htmltools` tag (invisibly printable in a knitr `results =
#'   'asis'` chunk via [emit_widget()], or returned directly from a
#'   normal chunk).
#' @examples
#' \dontrun{
#' download_table_button(res_df, filename = "A_vs_B_deg_results")
#' }
#' @export
download_table_button <- function(df, filename = "data", button_label = "Download CSV") {
    if (!requireNamespace("downloadthis", quietly = TRUE)) {
        stop("Package \"downloadthis\" is required for download buttons. Install it with ",
             "`install.packages(\"downloadthis\")`.", call. = FALSE)
    }
    downloadthis::download_this(
        as.data.frame(df),
        output_name = filename,
        output_extension = ".csv",
        button_label = button_label,
        button_type = "default",
        has_icon = TRUE,
        self_contained = TRUE,
        class = "btn-sm"
    )
}

#' Create a download button for a ggplot figure (PDF)
#'
#' Renders a "Download PDF" button next to any static `ggplot` figure shown in
#' the bundled HTML report. Inside a [render_report()] render the PDF is written
#' to the page's asset directory and linked; standalone, it is embedded via
#' `downloadthis::download_file()` (see `.pdf_download_tag()`). The PDF page
#' matches the on-page HTML size (see `dpi` / `widget`). Not applicable to the
#' `plotly`/`DT` *interactive* widgets' built-in camera/CSV controls --
#' this is for static `ggplot` figures (GSEA barplots, WGCNA figures,
#' PCA scree and loadings, etc). For interactive PCA / MA / volcano,
#' pass `widget` so the PDF matches the plotly canvas.
#'
#' @param plot A `ggplot` object.
#' @param filename Output file name *without* extension. Default `"plot"`.
#' @param button_label Button text. Default `"Download PDF"`.
#' @param width,height On-page figure size in inches (the same values
#'   passed to [capture_plot_tag()] / [emit_plot_with_download()]).
#'   Defaults `7`, `5`.
#' @param dpi Raster dpi used for the HTML `<img>` (default `150`, same
#'   as [capture_plot_tag()]). The PDF page is sized so that at 96 CSS
#'   pixels per inch it matches that on-page box, and the grob is scaled
#'   to fill it (same aspect ratio and relative layout as the HTML).
#' @param widget Optional `htmlwidget` (e.g. a `plotly` PCA / MA /
#'   volcano). When supplied, PDF size is taken from the widget's pixel
#'   width/height so it matches the interactive figure on the page.
#' @param pdf_width,pdf_height Optional explicit PDF page size in inches,
#'   overriding the size derived from the HTML box. Use for figures whose
#'   on-page shape suits a scrollable web page but not a printable page
#'   (see [pca_loadings_figure_size()]). Default `NULL` (match the HTML).
#' @return An `htmltools` tag (see [download_table_button()]).
#' @examples
#' \dontrun{
#' download_plot_button(gsea_barplot, filename = "A_vs_B_HALLMARK_gsea_barplot")
#' }
#' @export
download_plot_button <- function(plot, filename = "plot", button_label = "Download PDF",
                                  width = 7, height = 5, dpi = 150, widget = NULL,
                                  pdf_width = NULL, pdf_height = NULL) {
    if (!methods::is(plot, "ggplot")) {
        stop("`plot` must be a ggplot object.", call. = FALSE)
    }
    if (!is.null(widget)) {
        sz <- .htmlwidget_pdf_inches(widget)
        width <- sz$width
        height <- sz$height
        dpi <- 96
    }
    pdf_sz <- .override_pdf_in(.html_display_to_pdf_in(width, height, dpi), pdf_width, pdf_height)
    tmp_pdf <- tempfile(fileext = ".pdf")
    .save_ggplot_html_matched(plot, tmp_pdf, width, height, pdf_sz$width, pdf_sz$height)
    .pdf_download_tag(tmp_pdf, filename, button_label)
}

#' Create a download button for any grid-drawable figure (PDF)
#'
#' A more general counterpart to [download_plot_button()] for figures
#' that are not `ggplot` objects and so cannot go through
#' `ggplot2::ggsave()` -- notably the `pheatmap` grob returned by
#' [plot_sample_similarity()]. Draws the object to a PDF via the standard
#' `pdf()` graphics device and offers that file the same way
#' [download_plot_button()] does.
#'
#' @param drawable Any object with a `print()`/`plot()` method that
#'   renders to the active graphics device when drawn inside a `pdf()`
#'   device (e.g. a `pheatmap` grob, a base or `grid` plot object).
#' @param filename Output file name *without* extension. Default `"plot"`.
#' @param button_label Button text. Default `"Download PDF"`.
#' @param width,height On-page figure size in inches (same values used
#'   for the HTML display box). Defaults `7`, `7`.
#' @param dpi Raster dpi used for the HTML `<img>`. Default `150`. The
#'   PDF page matches that CSS box at 96 px/in; the grob is scaled to
#'   fill the page.
#' @param pdf_width,pdf_height Optional explicit PDF page size in inches,
#'   overriding the size derived from the HTML box (see
#'   [download_plot_button()]). Default `NULL` (match the HTML).
#' @return An `htmltools` tag (see [download_table_button()]).
#' @examples
#' \dontrun{
#' download_grid_plot_button(similarity_plot, filename = "sample_similarity_heatmap")
#' }
#' @export
download_grid_plot_button <- function(drawable, filename = "plot", button_label = "Download PDF",
                                       width = 7, height = 7, dpi = 150,
                                       pdf_width = NULL, pdf_height = NULL) {
    pdf_sz <- .override_pdf_in(.html_display_to_pdf_in(width, height, dpi), pdf_width, pdf_height)
    tmp_pdf <- tempfile(fileext = ".pdf")
    if (is.function(drawable)) {
        grDevices::pdf(tmp_pdf, width = pdf_sz$width, height = pdf_sz$height)
        drawable()
        grDevices::dev.off()
    } else {
        .save_gtable_html_matched(drawable, tmp_pdf, pdf_sz$width, pdf_sz$height)
    }
    .pdf_download_tag(tmp_pdf, filename, button_label)
}

#' Offer a just-written PDF for download, from disk inside a report
#'
#' The counterpart of [figure_img_tag()] for the PDF behind a download button.
#' While [render_report()] is rendering a page the file is moved into that
#' page's `<page>_figures/` directory and linked by relative path; outside a
#' render it falls back to `downloadthis::download_file()`, which inlines it as
#' a base64 `data:` URI so a standalone call still yields a self-contained tag.
#'
#' Inlining cost roughly a third of every contrast page: the MA and volcano
#' PDFs went in as 2.6 MB of base64 inside two `onclick` handlers, next to the
#' identical files the pipeline had already written in the same folder. That is
#' the same weight that made pandoc's post-processing of figures unusable, and
#' four workers converting such documents at once is what made pandoc fail
#' outright on individual pages.
#' @noRd
.pdf_download_tag <- function(pdf_path, filename, button_label = "Download PDF") {
    href <- .move_report_asset(pdf_path, filename, ".pdf")
    if (!is.null(href)) return(link_output_file(href, label = button_label))
    if (!requireNamespace("downloadthis", quietly = TRUE)) {
        stop("Package \"downloadthis\" is required for download buttons outside ",
             "a report render. Install it with `install.packages(\"downloadthis\")`.",
             call. = FALSE)
    }
    downloadthis::download_file(
        path = pdf_path,
        output_name = filename,
        button_label = button_label,
        button_type = "default",
        has_icon = TRUE,
        self_contained = TRUE,
        class = "btn-sm"
    )
}

#' Move a file into the page asset directory of the report being rendered
#'
#' @param path Existing file to move.
#' @param stem Desired base name, without extension; sanitised, and suffixed if
#'   another asset already claimed it (two figures on one page can share a
#'   `filename`).
#' @return The `href` to use, relative to the page's HTML file, or `NULL` when
#'   no asset directory is configured (caller should inline instead).
#' @noRd
.move_report_asset <- function(path, stem, ext) {
    dir <- getOption("degreports.figure_dir")
    href_base <- getOption("degreports.figure_href")
    if (is.null(dir) || is.null(href_base) || !file.exists(path)) return(NULL)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    stem <- gsub("[^A-Za-z0-9._+-]+", "_", basename(as.character(stem)[1]))
    if (!nzchar(stem)) stem <- "download"
    name <- paste0(stem, ext)
    i <- 1L
    while (file.exists(file.path(dir, name))) {
        i <- i + 1L
        name <- paste0(stem, "-", i, ext)
    }
    dest <- file.path(dir, name)
    if (!isTRUE(file.rename(path, dest))) {
        # file.rename() cannot cross devices, and the source is in tempdir().
        if (!isTRUE(file.copy(path, dest, overwrite = TRUE))) return(NULL)
        unlink(path)
    }
    paste0(href_base, "/", name)
}

#' PDF inches matching an HTML figure shown at `width`×`height` inches and `dpi`
#'
#' HTML `<img>` CSS size is `width * dpi` by `height * dpi` pixels.
#' PDF viewers treat 1 inch as 96 CSS pixels, so the downloaded page is
#' scaled by `dpi / 96`.
#' @noRd
.html_display_to_pdf_in <- function(width, height, dpi = 150) {
    dpi <- as.numeric(dpi)[1]
    if (!is.finite(dpi) || dpi <= 0) dpi <- 96
    list(
        width = as.numeric(width)[1] * dpi / 96,
        height = as.numeric(height)[1] * dpi / 96
    )
}

#' Replace either dimension of an HTML-matched PDF size with an explicit one
#'
#' Kept separate so a caller can override only the width or only the height
#' and leave the other tracking the on-page figure.
#' @noRd
.override_pdf_in <- function(pdf_sz, pdf_width = NULL, pdf_height = NULL) {
    if (!is.null(pdf_width)) pdf_sz$width <- as.numeric(pdf_width)[1]
    if (!is.null(pdf_height)) pdf_sz$height <- as.numeric(pdf_height)[1]
    pdf_sz
}

#' PDF inches matching an htmlwidget's on-page pixel size
#' @noRd
.htmlwidget_pdf_inches <- function(widget, fallback_width = 7, fallback_height = 5) {
    w <- if (!is.null(widget$width)) as.numeric(widget$width)[1] else NA_real_
    h <- if (!is.null(widget$height)) as.numeric(widget$height)[1] else NA_real_
    if ((!is.finite(w) || w <= 0) && !is.null(widget$x$layout$width)) {
        w <- as.numeric(widget$x$layout$width)[1]
    }
    if ((!is.finite(h) || h <= 0) && !is.null(widget$x$layout$height)) {
        h <- as.numeric(widget$x$layout$height)[1]
    }
    if (!is.finite(w) || w <= 0) w <- fallback_width * 96
    if (!is.finite(h) || h <= 0) h <- fallback_height * 96
    list(width = w / 96, height = h / 96)
}

#' Convert a gtable's row/column sizes to npc fractions of itself
#'
#' Absolute units (pt, mm, grobwidth) stay in proportion; `null` units
#' are resolved against the *current* device, so call this on a device
#' whose inches match the intended layout size.
#' @noRd
.gtable_to_npc <- function(gt) {
    w <- as.numeric(grid::convertWidth(gt$widths, "in", valueOnly = TRUE))
    h <- as.numeric(grid::convertHeight(gt$heights, "in", valueOnly = TRUE))
    w[!is.finite(w) | w < 0] <- 0
    h[!is.finite(h) | h < 0] <- 0
    tw <- sum(w)
    th <- sum(h)
    if (!is.finite(tw) || tw <= 0) tw <- 1
    if (!is.finite(th) || th <= 0) th <- 1
    # unit(0, "npc") viewports error ("Viewport has zero dimension(s)");
    # keep true zeros in points so gtable skips them.
    w_npc <- w / tw
    h_npc <- h / th
    gt$widths <- grid::unit(ifelse(w <= 0, 0, w_npc), ifelse(w <= 0, "pt", "npc"))
    gt$heights <- grid::unit(ifelse(h <= 0, 0, h_npc), ifelse(h <= 0, "pt", "npc"))
    gt
}

#' Multiply fontsize / cex on a grob tree so text scales with the layout
#' @noRd
.scale_grob_fontsize <- function(g, scale) {
    if (is.null(g) || !is.finite(scale) || scale == 1) return(g)
    if (inherits(g, "gtable")) {
        g$grobs <- lapply(g$grobs, .scale_grob_fontsize, scale = scale)
        return(g)
    }
    if (inherits(g, "gTree") && !is.null(g$children) && length(g$children)) {
        g$children <- lapply(g$children, .scale_grob_fontsize, scale = scale)
    }
    if (!is.null(g$gp)) {
        if (!is.null(g$gp$fontsize) && is.numeric(g$gp$fontsize)) {
            g$gp$fontsize <- g$gp$fontsize * scale
        }
        if (!is.null(g$gp$cex) && is.numeric(g$gp$cex)) {
            g$gp$cex <- g$gp$cex * scale
        }
        if (!is.null(g$gp$lwd) && is.numeric(g$gp$lwd)) {
            g$gp$lwd <- g$gp$lwd * scale
        }
    }
    g
}

#' Draw an npc-layout gtable onto a PDF page, uniformly scaled to fit
#' @noRd
.draw_npc_gtable_pdf <- function(gt_npc, path, layout_w, layout_h, pdf_w, pdf_h) {
    scale <- min(pdf_w / layout_w, pdf_h / layout_h)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    gt_npc <- .scale_grob_fontsize(gt_npc, scale)
    grDevices::pdf(path, width = pdf_w, height = pdf_h)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
        width = grid::unit(layout_w * scale, "in"),
        height = grid::unit(layout_h * scale, "in")
    ))
    grid::grid.draw(gt_npc)
    grid::popViewport()
    grDevices::dev.off()
}

#' Measure a gtable on a scratch device of a given size
#' @noRd
.gtable_measure <- function(gt, device_in) {
    scratch <- tempfile(fileext = ".pdf")
    grDevices::pdf(scratch, width = max(device_in, 1), height = max(device_in, 1))
    on.exit({
        grDevices::dev.off()
        unlink(scratch)
    }, add = TRUE)
    tryCatch(
        list(
            w = as.numeric(grid::convertWidth(gt$widths, "in", valueOnly = TRUE)),
            h = as.numeric(grid::convertHeight(gt$heights, "in", valueOnly = TRUE))
        ),
        error = function(e) NULL
    )
}

#' Split a gtable's rows/columns into fixed inches and device-elastic ones
#'
#' A `pheatmap` gtable has no natural size: its matrix row and column are
#' `npc`/`null`, so they absorb whatever page they are handed while the
#' surrounding furniture (dendrograms, labels, legends) stays a fixed number of
#' inches. Measuring at two device sizes separates the two -- a slope of zero
#' means fixed inches, anything else means elastic.
#' @noRd
.gtable_extent_parts <- function(gt) {
    small <- .gtable_measure(gt, 6)
    large <- .gtable_measure(gt, 18)
    part <- function(a, b) {
        if (is.null(a) || is.null(b) || !length(a)) {
            return(list(fixed = 0, elastic = TRUE))
        }
        slope <- (b - a) / 12
        slope[!is.finite(slope)] <- 0
        fixed <- a - slope * 6
        fixed[!is.finite(fixed) | fixed < 0] <- 0
        list(fixed = sum(fixed[slope <= 1e-6]), elastic = any(slope > 1e-6))
    }
    list(width = part(small$w, large$w), height = part(small$h, large$h))
}

#' Cell counts and label-driven minimum sizes for a pheatmap gtable
#'
#' The title and the colour legend are not sized by the gtable layout -- the
#' title is a text grob centred over the matrix column and the legend is drawn
#' inside the (elastic) matrix cell -- so neither shows up in `sum(gt$widths)`
#' or `sum(gt$heights)`. Both are returned here as floors on the page.
#' @noRd
.gtable_cell_metrics <- function(gt) {
    grob_named <- function(name) {
        i <- which(gt$layout$name == name)
        if (!length(i)) NULL else gt$grobs[[i[[1]]]]
    }
    text_lines <- function(g) {
        if (is.null(g) || is.null(g$label)) return(1L)
        max(1L, max(lengths(strsplit(as.character(g$label), "\n", fixed = TRUE))))
    }
    fontsize <- function(g, default = 8) {
        fs <- if (is.null(g)) NULL else g$gp$fontsize
        if (!length(fs) || !is.finite(fs[[1]])) default else as.numeric(fs[[1]])
    }
    rows <- grob_named("row_names")
    cols <- grob_named("col_names")
    n_row <- if (!is.null(rows$label)) length(rows$label) else NA_integer_
    n_col <- if (!is.null(cols$label)) length(cols$label) else NA_integer_
    if (is.na(n_col)) {
        anno <- grob_named("col_annotation")
        if (!is.null(anno$x)) n_col <- length(unique(as.numeric(anno$x)))
    }

    scratch <- tempfile(fileext = ".pdf")
    grDevices::pdf(scratch, width = 12, height = 12)
    title_width <- tryCatch({
        main <- grob_named("main")
        if (is.null(main)) 0 else
            grid::convertWidth(grid::grobWidth(main), "in", valueOnly = TRUE) + 0.2
    }, error = function(e) 0)
    grDevices::dev.off()
    unlink(scratch)
    if (!is.finite(title_width) || title_width < 0) title_width <- 0

    # A legend's tick labels are stacked, so it needs one line-height each.
    legend_height <- vapply(
        which(gt$layout$name %in% c("legend", "annotation_legend")),
        function(i) {
            kids <- gt$grobs[[i]]$children
            texts <- if (length(kids)) Filter(function(k) !is.null(k$label), kids) else list()
            if (!length(texts)) return(0)
            max(vapply(texts, function(k) {
                length(k$label) * fontsize(k, 10) * 1.6 / 72
            }, numeric(1)))
        }, numeric(1)
    )

    # A cell must be at least as tall as its (possibly wrapped) row label, and
    # at least as wide as one line of the rotated column label.
    list(
        n_row = n_row,
        n_col = n_col,
        min_height = fontsize(rows) * text_lines(rows) * 1.35 / 72,
        min_width = fontsize(cols) * 1.35 / 72,
        title_width = title_width,
        legend_height = if (length(legend_height)) max(legend_height) else 0
    )
}

#' Page size in inches for a `pheatmap` / `gtable` figure
#'
#' Sizes the page from the figure's content -- fixed furniture plus one cell per
#' row and column -- rather than from `sum(gt$heights)`, which for a `pheatmap`
#' just reports back whatever device it was measured on (see
#' `.gtable_extent_parts()`). Cells are never smaller than the labels they carry.
#'
#' @param drawable A `pheatmap` object or a bare `gtable`.
#' @param cell_width,cell_height Target inches per matrix column / row. Raised
#'   to fit the column / row labels when those need more room.
#' @param max_width,max_height Page caps in inches.
#' @return A list with the page `width`/`height` in inches, and `fit`: the
#'   uniform scale needed to bring the figure's layout onto that page, which is
#'   `1` unless an all-absolute gtable overruns the caps.
#' @export
gtable_page_size <- function(drawable, cell_width = 0.28, cell_height = 0.16,
                             max_width = 12, max_height = 24) {
    gt <- if (is.list(drawable) && !is.null(drawable$gtable)) drawable$gtable else drawable
    parts <- .gtable_extent_parts(gt)
    cells <- .gtable_cell_metrics(gt)

    side <- function(part, n, cell, min_cell, floor_in, cap) {
        if (!part$elastic) return(max(part$fixed, floor_in, 0.5))
        if (is.na(n) || n < 1L) return(cap)
        part$fixed + max(n * max(cell, min_cell), floor_in)
    }
    w <- side(parts$width, cells$n_col, cell_width, cells$min_width, 0, max_width)
    h <- side(parts$height, cells$n_row, cell_height, cells$min_height,
              cells$legend_height, max_height)
    w <- min(max(w, 1), max_width)
    h <- min(max(h, 1), max_height)
    w <- min(w + .gtable_title_deficit(gt, cells$title_width, w, h), max_width)

    # An all-absolute gtable can still overrun that page; scale rather than let
    # grid clip it.
    fit <- 1
    layout <- .gtable_layout_extent(gt, w, h)
    if (!is.null(layout)) {
        fit <- min(1, w / layout$width, h / layout$height)
        if (!is.finite(fit) || fit <= 0) fit <- 1
    }
    list(width = w, height = h, fit = fit,
         layout_width = if (is.null(layout)) w else layout$width,
         layout_height = if (is.null(layout)) h else layout$height)
}

#' Extra page width a `pheatmap` title needs to stop overhanging
#'
#' `pheatmap` centres the title over the matrix column, not over the figure, and
#' the matrix column is left of centre (row labels and legend sit to its right).
#' A title wider than that column therefore runs off the left edge, and simply
#' widening the page does not help: the surplus goes to the matrix column and
#' carries the title's centre right with it. Returning the shortfall lets the
#' caller add it once, which works whether the surplus widens the matrix column
#' (elastic gtable) or just pads both margins (absolute one).
#' @noRd
.gtable_title_deficit <- function(gt, title_width, page_w, page_h) {
    if (!is.finite(title_width) || title_width <= 0) return(0)
    mat <- gt$layout[gt$layout$name == "matrix", , drop = FALSE]
    if (!nrow(mat)) return(0)
    scratch <- tempfile(fileext = ".pdf")
    grDevices::pdf(scratch, width = max(page_w, 1), height = max(page_h, 1))
    on.exit({
        grDevices::dev.off()
        unlink(scratch)
    }, add = TRUE)
    widths <- tryCatch(
        as.numeric(grid::convertWidth(gt$widths, "in", valueOnly = TRUE)),
        error = function(e) NULL
    )
    if (is.null(widths)) return(0)
    lead <- if (mat$l[[1]] > 1L) sum(widths[seq_len(mat$l[[1]] - 1L)]) else 0
    span <- sum(widths[seq(mat$l[[1]], mat$r[[1]])])
    max(0, title_width - 2 * lead - span)
}

#' Total layout extent of a gtable on a page of the given inches
#' @noRd
.gtable_layout_extent <- function(gt, page_w, page_h) {
    scratch <- tempfile(fileext = ".pdf")
    grDevices::pdf(scratch, width = max(page_w, 1), height = max(page_h, 1))
    out <- tryCatch(
        list(
            width = grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE),
            height = grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
        ),
        error = function(e) NULL
    )
    grDevices::dev.off()
    unlink(scratch)
    if (is.null(out) || !is.finite(out$width) || !is.finite(out$height) ||
        out$width <= 0 || out$height <= 0) {
        return(NULL)
    }
    out
}

#' Zero-argument drawer for a gtable, uniformly scaled onto its page
#'
#' Drawn as-is when the layout fits the page. When it does not, the layout is
#' converted to npc fractions and drawn inside a viewport scaled to fit, with
#' font sizes scaled to match. The alternative -- `grid.draw()` onto a page
#' smaller than the layout -- silently clips, and because gtable centres its
#' layout it eats the figure from *both* ends: the column dendrogram off the top
#' and the column labels off the bottom.
#' @noRd
.gtable_fitted_drawer <- function(gt, size) {
    if (size$fit >= 1 - 1e-6) {
        return(function() {
            grid::grid.newpage()
            grid::grid.draw(gt)
        })
    }
    function() {
        gt_npc <- tryCatch(.gtable_to_npc(gt), error = function(e) NULL)
        grid::grid.newpage()
        if (is.null(gt_npc)) {
            grid::grid.draw(gt)
            return(invisible(NULL))
        }
        grid::pushViewport(grid::viewport(
            width = grid::unit(size$layout_width * size$fit, "in"),
            height = grid::unit(size$layout_height * size$fit, "in")
        ))
        grid::grid.draw(.scale_grob_fontsize(gt_npc, size$fit))
        grid::popViewport()
    }
}

#' Emit a `pheatmap` heatmap as an on-page figure with its download buttons
#'
#' The bundled report templates all showed heatmaps by measuring the gtable,
#' capping that size, and then `grid.draw()`-ing it -- which clips whenever the
#' measurement is off or the cap bites (see `.gtable_fitted_drawer()`). This
#' measures on a matching device via [gtable_page_size()], scales to fit, and
#' sizes the downloaded PDF to the same page.
#'
#' @param heatmap A `pheatmap` object (or a bare `gtable`).
#' @param filename Base name for the alt text and the downloaded PDF.
#' @param buttons Extra buttons (a list of tags) to show beside the
#'   "Download PDF" button, e.g. a [link_output_file()] CSV link.
#' @param scale On-page size relative to the figure's true physical size.
#'   `1` (default) shows it at the size it would print.
#' @param cell_width,cell_height,max_width,max_height Passed to
#'   [gtable_page_size()].
#' @param responsive If `TRUE`, let the figure shrink to fit a narrow report
#'   panel (`max-width: 100%`) instead of holding a fixed pixel size.
#' @return Invisibly `NULL`; called for its side effect of `cat()`-ing HTML.
#' @export
emit_heatmap_with_download <- function(heatmap, filename, buttons = list(),
                                        scale = 1, cell_width = 0.28,
                                        cell_height = 0.16, max_width = 12,
                                        max_height = 24, responsive = FALSE) {
    if (is.null(heatmap)) return(invisible(NULL))
    gt <- if (is.list(heatmap) && !is.null(heatmap$gtable)) heatmap$gtable else heatmap
    size <- gtable_page_size(heatmap, cell_width = cell_width,
                             cell_height = cell_height,
                             max_width = max_width, max_height = max_height)

    css_w <- as.integer(round(size$width * 96 * scale))
    css_h <- as.integer(round(size$height * 96 * scale))
    style <- if (isTRUE(responsive)) {
        sprintf("max-width: 100%% !important; height: auto !important; width: %dpx;", css_w)
    } else {
        sprintf("max-width: none !important; width: %dpx; height: %dpx;", css_w, css_h)
    }
    img <- figure_img_tag(
        .gtable_fitted_drawer(gt, size),
        px_width = size$width * 72,
        px_height = size$height * 72,
        pixelratio = 2,
        alt = filename,
        attribs = list(width = css_w, height = css_h, style = style)
    )
    pdf_button <- download_grid_plot_button(
        heatmap, filename = filename,
        pdf_width = size$width, pdf_height = size$height
    )
    emit_widget(wrap_with_corner_button(img, c(list(pdf_button), as.list(buttons))))
}

#' Print a ggplot onto a PDF page matching the HTML CSS box
#'
#' `ggplotGrob` + npc conversion fails for `ggrepel` / `coord_equal`
#' plots ("Viewport has zero dimension(s)"). Printing onto the target
#' device lets ggplot2 lay the figure out at the HTML-matched page size.
#' @noRd
.save_ggplot_html_matched <- function(plot, path, layout_w, layout_h, pdf_w, pdf_h) {
    pdf_w <- max(as.numeric(pdf_w)[1], 0.5)
    pdf_h <- max(as.numeric(pdf_h)[1], 0.5)
    ggplot2::ggsave(
        filename = path,
        plot = plot,
        device = grDevices::pdf,
        width = pdf_w,
        height = pdf_h,
        units = "in",
        limitsize = FALSE
    )
}

#' Scale a pheatmap/gtable grob to fit a PDF page (contain, no clip)
#' @noRd
.save_gtable_html_matched <- function(drawable, path, pdf_w, pdf_h) {
    gt <- if (is.list(drawable) && !is.null(drawable$gtable)) drawable$gtable else drawable
    pdf_w <- max(as.numeric(pdf_w)[1], 0.5)
    pdf_h <- max(as.numeric(pdf_h)[1], 0.5)
    tmp_layout <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp_layout, width = max(pdf_w, 1), height = max(pdf_h, 1))
    native_w <- tryCatch(
        grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE),
        error = function(e) pdf_w
    )
    native_h <- tryCatch(
        grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE),
        error = function(e) pdf_h
    )
    if (!is.finite(native_w) || native_w <= 0) native_w <- pdf_w
    if (!is.finite(native_h) || native_h <= 0) native_h <- pdf_h
    gt_npc <- tryCatch(.gtable_to_npc(gt), error = function(e) NULL)
    grDevices::dev.off()
    unlink(tmp_layout)
    ok <- FALSE
    if (!is.null(gt_npc)) {
        ok <- isTRUE(tryCatch({
            .draw_npc_gtable_pdf(gt_npc, path, native_w, native_h, pdf_w, pdf_h)
            TRUE
        }, error = function(e) FALSE))
    }
    if (!ok) {
        w <- max(pdf_w, native_w)
        h <- max(pdf_h, native_h)
        grDevices::pdf(path, width = w, height = h)
        grid::grid.newpage()
        grid::grid.draw(gt)
        grDevices::dev.off()
    }
}

#' Render any grid-drawable object (ggplot or grid/pheatmap grob) to a
#' self-contained base64 `<img>` tag
#'
#' Used to give a figure a stable `htmltools` tag *before* rendering, so
#' it can be composed into a single flex/positioned container alongside
#' its download button via [wrap_with_corner_button()] -- knitr's usual
#' implicit `print(plot)`-to-device figure capture happens in a separate
#' output stream from any `results = 'asis'` HTML emitted in the same or
#' another chunk, so the two cannot otherwise be nested inside one `<div>`.
#'
#' @param drawable A `ggplot` object, or any other object with a
#'   `print()`/`plot()` method that renders to the active graphics
#'   device (e.g. a `pheatmap` grob).
#' @param width,height Figure size in inches. Defaults `7`, `5`.
#' @param dpi Resolution (dots per inch) used to rasterize the figure at
#'   its logical (CSS) display size. Default `150`.
#' @param pixelratio Physical-to-logical pixel ratio used when actually
#'   drawing the bitmap (`htmltools::plotTag()`'s own default is `2`,
#'   i.e. draw at twice the CSS size for crisp rendering on high-DPI/Retina
#'   displays -- the previous hardcoded `pixelratio = 1` here produced
#'   visibly softer, blurrier figures than every interactive `plotly`
#'   widget in the same report). Default `2`.
#' @param alt Alt text for the `<img>` tag. Default `""`.
#' @return An `htmltools` tag wrapping a base64-encoded PNG `<img>`.
#' @details
#' The `<img>` is sized via HTML `width`/`height` *attributes* (not an
#' inline CSS `style` size). This matters because the bundled report's
#' global stylesheet forces `img { max-width: 100% !important; }` so
#' figures shrink to fit inside narrower report panels/tabs -- with a
#' `width`/`height` `style` (as `htmltools::plotTag()` sets by default),
#' `max-width` overrides only the width half of that fixed style while
#' the height stays pinned at its original literal pixel value, visibly
#' stretching/squishing the image out of its true aspect ratio wherever
#' the panel is narrower than the figure's native size (e.g. multi-panel
#' WGCNA plots in a narrower tab). `width`/`height` HTML attributes,
#' rather than a `style`, let the browser derive and preserve the
#' intrinsic aspect ratio automatically as it shrinks the image to fit.
#' @export
capture_plot_tag <- function(drawable, width = 7, height = 5, dpi = 150, pixelratio = 2, alt = "") {
    if (!requireNamespace("htmltools", quietly = TRUE)) {
        stop("Package \"htmltools\" is required to capture plot tags.", call. = FALSE)
    }
    px_width <- as.integer(round(width * dpi))
    px_height <- as.integer(round(height * dpi))
    figure_img_tag(function() print(drawable), px_width, px_height,
                   pixelratio = pixelratio, alt = alt)
}

#' Rasterize a figure into an HTML `<img>`, written to disk inside a report
#'
#' Drop-in replacement for `htmltools::plotTag()` in the bundled report
#' templates. Geometry is identical to what `plotTag()` would produce -- the
#' only difference is where the PNG lives. While [render_report()] is
#' rendering a page the bitmap is written to a `<page>_figures/` directory
#' beside that page's HTML and referenced by relative path; outside a render
#' it falls back to `plotTag()`'s inline base64 `data:` URI, so a standalone
#' call still yields a self-contained tag.
#'
#' @details
#' Inlining is what a single-file HTML page wants, but it makes
#' `rmarkdown::render()` unusable for figure-heavy pages. For a
#' non-self-contained document, `rmarkdown:::html_document_base()`'s
#' post-processor rewrites every `<img src>` through `utils::URLdecode()`,
#' which walks the string one byte at a time and regrows its output vector on
#' each step. That is quadratic, so it is free for a short file path and
#' catastrophic for a base64 payload: measured on the WGCNA page's own
#' figures, 40 KB takes 0.6 s and its largest (1.5 MB) figure takes about 14
#' minutes. Knitting and pandoc both completed and the HTML was already on
#' disk; the pipeline then sat in that post-processor for over half an hour
#' with no output, which is indistinguishable from a hang.
#'
#' @param draw Zero-argument function that draws the figure to the active
#'   graphics device (e.g. `function() grid::grid.draw(gt)`).
#' @param px_width,px_height Logical (CSS) pixel size of the bitmap.
#' @param pixelratio Physical-to-logical pixel ratio used when drawing, as in
#'   [capture_plot_tag()]. Default `2`.
#' @param alt Alt text for the `<img>`. Default `""`.
#' @param attribs Extra `<img>` attributes, or `NULL` (default) to pin the
#'   display size to `px_width`/`px_height` the way [capture_plot_tag()] does.
#' @return An `htmltools` `<img>` tag.
#' @export
figure_img_tag <- function(draw, px_width, px_height, pixelratio = 2, alt = "",
                           attribs = NULL) {
    if (!requireNamespace("htmltools", quietly = TRUE)) {
        stop("Package \"htmltools\" is required to capture plot tags.", call. = FALSE)
    }
    px_width <- as.integer(round(px_width))
    px_height <- as.integer(round(px_height))
    if (is.null(attribs)) {
        # Pin both the intrinsic bitmap size and the CSS display size, and
        # override the report theme's `img { max-width: 100% !important; }` so
        # figures stay at a fixed pixel size instead of scaling with the
        # browser window / tab width.
        attribs <- list(
            width = px_width,
            height = px_height,
            style = sprintf(
                "max-width: none !important; width: %dpx; height: %dpx;",
                px_width, px_height
            )
        )
    }
    src <- .write_report_figure(draw, px_width, px_height, pixelratio)
    if (is.null(src)) {
        return(htmltools::plotTag(
            draw(),
            alt = alt,
            width = px_width,
            height = px_height,
            pixelratio = pixelratio,
            suppressSize = "xy",
            attribs = attribs
        ))
    }
    htmltools::browsable(
        do.call(htmltools::tags$img, c(list(src = src, alt = alt), attribs))
    )
}

#' Write one report figure next to the page being rendered
#'
#' @return The `src` to use, relative to the page's HTML file, or `NULL` when
#'   no figure directory is configured (caller should inline instead).
#' @noRd
.write_report_figure <- function(draw, px_width, px_height, pixelratio) {
    fig_dir <- getOption("degreports.figure_dir")
    fig_href <- getOption("degreports.figure_href")
    if (is.null(fig_dir) || is.null(fig_href)) return(NULL)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    dest <- tempfile("fig-", tmpdir = fig_dir, fileext = ".png")
    htmltools::capturePlot(
        draw(),
        filename = dest,
        device = htmltools::defaultPngDevice(),
        width = px_width * pixelratio,
        height = px_height * pixelratio,
        res = 72 * pixelratio
    )
    paste0(fig_href, "/", basename(dest))
}

#' Compose a figure/table and its download button(s) into one container,
#' with the button(s) right-aligned directly above the figure/table
#'
#' Wraps `content` (a static plot `<img>` tag from [capture_plot_tag()],
#' a `plotly`/`DT` htmlwidget, or any other `htmltools`-taggable object)
#' and one or more download buttons in a single `<div>`: a right-aligned
#' button row in normal document flow, immediately above the figure/table
#' (pushing it down rather than overlaying on top of it), followed by the
#' content itself. Multiple buttons (e.g. a plot's PDF button alongside a
#' companion table's CSV button) sit side by side in that row.
#'
#' @param content The figure or table to display (an `htmltools` tag,
#'   `plotly` widget, `DT::datatable`, or similar).
#' @param button One button tag (from [download_table_button()],
#'   [download_plot_button()], [download_grid_plot_button()], or
#'   [link_output_file()]), or a list of several, to place in the row
#'   above the figure/table.
#' @param center If `TRUE` (the default), the figure/table content is
#'   centered on the page. Download buttons are always right-aligned to
#'   the full content width (not just the figure's shrink-wrapped edge).
#'   If `FALSE`, content spans the full available width (needed for
#'   responsive `DT` tables).
#' @return An `htmltools` tag (see [download_table_button()]).
#' @examples
#' \dontrun{
#' wrap_with_corner_button(
#'   capture_plot_tag(my_ggplot, width = 7, height = 5),
#'   download_plot_button(my_ggplot, filename = "my_plot")
#' )
#' }
#' @export
wrap_with_corner_button <- function(content, button, center = TRUE) {
    if (!requireNamespace("htmltools", quietly = TRUE)) {
        stop("Package \"htmltools\" is required to position download buttons.", call. = FALSE)
    }
    buttons <- if (inherits(button, "list") && !inherits(button, "shiny.tag")) button else list(button)
    buttons <- lapply(buttons, .normalize_dl_button)
    style_tag <- htmltools::tags$style(htmltools::HTML(paste(
        ".degreports-dl-row {",
        "  display: flex !important;",
        "  flex-direction: row !important;",
        "  justify-content: flex-end !important;",
        "  align-items: center !important;",
        "  align-content: center !important;",
        "  gap: 6px !important;",
        "  margin-bottom: 6px !important;",
        "  width: 100% !important;",
        "}",
        ".degreports-dl-btn-wrap {",
        "  display: inline-flex !important;",
        "  align-items: center !important;",
        "  align-self: center !important;",
        "  vertical-align: middle !important;",
        "  line-height: 0 !important;",
        "  margin: 0 !important;",
        "  padding: 0 !important;",
        "}",
        ".degreports-dl-row .btn,",
        ".degreports-dl-row a.btn,",
        ".degreports-dl-row button.btn,",
        ".degreports-dl-row .btn-group,",
        ".degreports-dl-row .btn-group > .btn {",
        "  padding: 4px 10px !important;",
        "  font-size: 12px !important;",
        "  line-height: 1.4 !important;",
        "  height: 28px !important;",
        "  min-height: 28px !important;",
        "  max-height: 28px !important;",
        "  box-sizing: border-box !important;",
        "  display: inline-flex !important;",
        "  align-items: center !important;",
        "  align-self: center !important;",
        "  vertical-align: middle !important;",
        "  margin: 0 !important;",
        "  white-space: nowrap !important;",
        "}",
        ".degreports-dl-row .btn i,",
        ".degreports-dl-row a.btn i {",
        "  font-size: 12px !important;",
        "  line-height: 1 !important;",
        "  margin: 0 0.35em 0 0 !important;",
        "  width: auto !important;",
        "  height: auto !important;",
        "}"
    )))
    button_row <- htmltools::tags$div(
        class = "degreports-dl-row",
        style_tag,
        buttons
    )
    content_div <- if (isTRUE(center)) {
        htmltools::tags$div(
            style = "display: flex; justify-content: center; width: 100%;",
            htmltools::tags$div(style = "display: inline-block; text-align: left;", content)
        )
    } else {
        htmltools::tags$div(style = "width: 100%;", content)
    }
    htmltools::tags$div(
        style = "width: 100%;",
        button_row,
        content_div
    )
}

#' Emit a static plot (ggplot or grid-drawable) together with its PDF
#' download button, positioned in the top-right corner
#'
#' Convenience wrapper combining [capture_plot_tag()],
#' [download_plot_button()]/[download_grid_plot_button()], and
#' [wrap_with_corner_button()] for the common case of a single static
#' figure with a single PDF download button. Intended for use in a
#' knitr `results = 'asis'` chunk.
#'
#' @param plot A `ggplot` object, or any other grid-drawable object
#'   (e.g. a `pheatmap` grob).
#' @param filename Output file name *without* extension.
#' @param width,height Figure size in inches. Defaults `7`, `5`.
#' @param button_label Button text. Default `"Download PDF"`.
#' @param dpi Resolution used to rasterize the on-page figure, and to
#'   size the downloaded PDF so it matches that CSS box. Default `150`.
#' @param pdf_width,pdf_height Optional explicit PDF page size in inches, for
#'   figures whose on-page shape does not suit a printable page (see
#'   [pca_loadings_figure_size()]). Default `NULL` (match the HTML box).
#' @return Invisibly `NULL`; called for its side effect of `cat()`-ing
#'   the rendered HTML (see [emit_widget()]).
#' @export
emit_plot_with_download <- function(plot, filename, width = 7, height = 5,
                                     button_label = "Download PDF", dpi = 150,
                                     pdf_width = NULL, pdf_height = NULL) {
    button <- if (methods::is(plot, "ggplot")) {
        download_plot_button(plot, filename = filename, button_label = button_label,
                             width = width, height = height, dpi = dpi,
                             pdf_width = pdf_width, pdf_height = pdf_height)
    } else {
        download_grid_plot_button(plot, filename = filename, button_label = button_label,
                                  width = width, height = height, dpi = dpi,
                                  pdf_width = pdf_width, pdf_height = pdf_height)
    }
    content <- capture_plot_tag(plot, width = width, height = height, dpi = dpi, alt = filename)
    emit_widget(wrap_with_corner_button(content, button))
}

#' Emit an interactive htmlwidget (`plotly`/`DT`) together with one or
#' more download buttons, positioned in the top-right corner
#'
#' Thin wrapper around [wrap_with_corner_button()] + [emit_widget()] for
#' the common case of an interactive widget (already a taggable object,
#' unlike a static `ggplot`) paired with a CSV/PDF download button.
#'
#' @param widget A `plotly` or `DT::datatable` htmlwidget.
#' @param button One button tag, or a list of several, to place in the
#'   row above the widget (see [wrap_with_corner_button()]).
#' @param center If `TRUE`, center the widget on the page under a
#'   full-width right-aligned button row. Use for [make_interactive()]
#'   plotly widgets. If `FALSE` (the default), content spans full width
#'   -- use for `DT` tables.
#' @return Invisibly `NULL`; called for its side effect of `cat()`-ing
#'   the rendered HTML (see [emit_widget()]).
#' @export
emit_widget_with_download <- function(widget, button, center = FALSE) {
    emit_widget(wrap_with_corner_button(widget, button, center = center))
}
