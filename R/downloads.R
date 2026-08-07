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
        style = "margin: 0 0 0.4em 0;",
        htmltools::tags$i(class = "fa fa-download", style = "margin-right: 0.35em;"),
        label
    )
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
        self_contained = TRUE
    )
}

#' Create a download button for a ggplot figure (PDF)
#'
#' Thin wrapper around `downloadthis::download_this()` that renders a
#' small self-contained "Download plot" button next to any static
#' `ggplot` figure shown in the bundled HTML report, exporting it as a
#' vector PDF at the given size on click. Not applicable to the
#' `plotly`/`DT` *interactive* widgets in the report (MA/volcano plots,
#' DEG/GSEA tables) -- those are already exportable via their own
#' built-in camera/CSV controls; this is for the static `ggplot` figures
#' (GSEA barplots, WGCNA diagnostic/heatmap/barplot figures, PCA scree
#' and loadings plots, etc).
#'
#' @param plot A `ggplot` object.
#' @param filename Output file name *without* extension. Default `"plot"`.
#' @param button_label Button text. Default `"Download PDF"`.
#' @param width,height Figure size in inches, passed to `ggplot2::ggsave()`.
#'   Defaults `7`, `5`.
#' @return An `htmltools` tag (see [download_table_button()]).
#' @examples
#' \dontrun{
#' download_plot_button(gsea_barplot, filename = "A_vs_B_HALLMARK_gsea_barplot")
#' }
#' @export
download_plot_button <- function(plot, filename = "plot", button_label = "Download PDF",
                                  width = 7, height = 5) {
    if (!requireNamespace("downloadthis", quietly = TRUE)) {
        stop("Package \"downloadthis\" is required for download buttons. Install it with ",
             "`install.packages(\"downloadthis\")`.", call. = FALSE)
    }
    if (!methods::is(plot, "ggplot")) {
        stop("`plot` must be a ggplot object.", call. = FALSE)
    }
    downloadthis::download_this(
        plot,
        output_name = filename,
        output_extension = ".pdf",
        button_label = button_label,
        button_type = "default",
        has_icon = TRUE,
        self_contained = TRUE,
        ggsave_args = list(width = width, height = height, units = "in")
    )
}

#' Create a download button for any grid-drawable figure (PDF)
#'
#' A more general counterpart to [download_plot_button()] for figures
#' that are not `ggplot` objects and so cannot go through
#' `ggplot2::ggsave()` -- notably the `pheatmap` grob returned by
#' [plot_sample_similarity()]. Draws the object to a temporary PDF via
#' the standard `pdf()` graphics device and wraps that file with
#' `downloadthis::download_file()`.
#'
#' @param drawable Any object with a `print()`/`plot()` method that
#'   renders to the active graphics device when drawn inside a `pdf()`
#'   device (e.g. a `pheatmap` grob, a base or `grid` plot object).
#' @param filename Output file name *without* extension. Default `"plot"`.
#' @param button_label Button text. Default `"Download PDF"`.
#' @param width,height Figure size in inches, passed to `grDevices::pdf()`.
#'   Defaults `7`, `7`.
#' @return An `htmltools` tag (see [download_table_button()]).
#' @examples
#' \dontrun{
#' download_grid_plot_button(similarity_plot, filename = "sample_similarity_heatmap")
#' }
#' @export
download_grid_plot_button <- function(drawable, filename = "plot", button_label = "Download PDF",
                                       width = 7, height = 7) {
    if (!requireNamespace("downloadthis", quietly = TRUE)) {
        stop("Package \"downloadthis\" is required for download buttons. Install it with ",
             "`install.packages(\"downloadthis\")`.", call. = FALSE)
    }
    tmp_pdf <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp_pdf, width = width, height = height)
    grDevices::dev.control("enable")
    print(drawable)
    grDevices::dev.off()

    downloadthis::download_file(
        path = tmp_pdf,
        output_name = filename,
        button_label = button_label,
        button_type = "default",
        has_icon = TRUE,
        self_contained = TRUE
    )
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
    # Pin both the intrinsic bitmap size and the CSS display size, and
    # override the report theme's `img { max-width: 100% !important; }` so
    # figures stay at a fixed pixel size instead of scaling with the
    # browser window / tab width.
    htmltools::plotTag(
        print(drawable),
        alt = alt,
        width = px_width,
        height = px_height,
        pixelratio = pixelratio,
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
    button_row <- htmltools::tags$div(
        style = paste(
            "display: flex; flex-direction: row; justify-content: flex-end;",
            "align-items: center; gap: 6px; margin-bottom: 6px; width: 100%;"
        ),
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
#' @param dpi Resolution used to rasterize the on-page figure. Default
#'   `150`, matching [capture_plot_tag()]'s own default (the previous
#'   `120` here, and `pixelratio = 1` in the older [capture_plot_tag()],
#'   produced visibly softer on-page figures than the same figure's own
#'   downloaded PDF).
#' @return Invisibly `NULL`; called for its side effect of `cat()`-ing
#'   the rendered HTML (see [emit_widget()]).
#' @export
emit_plot_with_download <- function(plot, filename, width = 7, height = 5,
                                     button_label = "Download PDF", dpi = 150) {
    button <- if (methods::is(plot, "ggplot")) {
        download_plot_button(plot, filename = filename, button_label = button_label, width = width, height = height)
    } else {
        download_grid_plot_button(plot, filename = filename, button_label = button_label, width = width, height = height)
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
