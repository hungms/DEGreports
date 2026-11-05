# Tests for the download / link helpers (download_table_button,
# download_plot_button, download_grid_plot_button, link_output_file).
# Embedded downloadthis buttons remain available for small ad-hoc tables;
# large pipeline outputs should prefer link_output_file() to on-disk CSVs.

skip_if_not_installed("downloadthis")
skip_if_not_installed("htmltools")

test_that("link_output_file emits a relative href download link", {
    tag <- link_output_file("IgA_vs_IgM/IgA_vs_IgM_diffexp.csv", label = "Download full CSV")
    expect_true(inherits(tag, "shiny.tag"))
    html_txt <- as.character(htmltools::as.tags(tag))
    expect_true(grepl("Download full CSV", html_txt))
    expect_true(grepl('href="IgA_vs_IgM/IgA_vs_IgM_diffexp.csv"', html_txt))
    expect_false(grepl("base64", html_txt))
})

test_that("link_output_file rejects empty paths", {
    expect_error(link_output_file(""), "non-empty")
})

test_that("download_table_button wraps a data frame as a self-contained CSV download button", {
    df <- data.frame(gene = c("TP53", "BRCA1"), log2FoldChange = c(1.5, -2.1))
    btn <- download_table_button(df, filename = "my_results")
    expect_true(inherits(btn, "shiny.tag") || inherits(btn, "shiny.tag.list") || inherits(btn, "html"))
    html_txt <- as.character(htmltools::as.tags(btn))
    expect_true(grepl("Download CSV", html_txt))
    # self-contained: the CSV payload is embedded (base64) rather than
    # referencing an external file/URL.
    expect_true(grepl("base64", html_txt))
})

test_that("download_table_button accepts a custom button label", {
    df <- data.frame(x = 1:3)
    btn <- download_table_button(df, filename = "data", button_label = "Get the data")
    html_txt <- as.character(htmltools::as.tags(btn))
    expect_true(grepl("Get the data", html_txt))
})

test_that("HTML display inches convert to PDF inches at 96 CSS pixels per inch", {
    expect_equal(
        DEGreports:::.html_display_to_pdf_in(4, 3, dpi = 96),
        list(width = 4, height = 3)
    )
    expect_equal(
        DEGreports:::.html_display_to_pdf_in(2, 1, dpi = 192),
        list(width = 4, height = 2)
    )
})

.pdf_mediabox_inches <- function(path) {
    raw <- readBin(path, what = "raw", n = file.info(path)$size)
    hit <- grepRaw("/MediaBox", raw, fixed = TRUE)
    if (!length(hit)) return(NULL)
    window <- raw[seq.int(hit[[1L]], min(length(raw), hit[[1L]] + 80L))]
    window <- window[window >= as.raw(32) & window <= as.raw(126)]
    txt <- rawToChar(window)
    m <- regexec("\\[\\s*([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)\\s*\\]", txt)
    g <- regmatches(txt, m)[[1]]
    if (length(g) < 5L) return(NULL)
    list(
        width = (as.numeric(g[4]) - as.numeric(g[2])) / 72,
        height = (as.numeric(g[5]) - as.numeric(g[3])) / 72
    )
}

test_that("ggplot download PDF page matches the HTML CSS box and is not a tiny inset", {
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    tmp <- tempfile(fileext = ".pdf")
    on.exit(unlink(tmp), add = TRUE)
    pdf_sz <- DEGreports:::.html_display_to_pdf_in(4, 3, dpi = 150)
    DEGreports:::.save_ggplot_html_matched(p, tmp, 4, 3, pdf_sz$width, pdf_sz$height)
    box <- .pdf_mediabox_inches(tmp)
    skip_if(is.null(box), "could not parse PDF MediaBox")
    expect_equal(box$width, pdf_sz$width, tolerance = 0.05)
    expect_equal(box$height, pdf_sz$height, tolerance = 0.05)
    expect_gt(file.info(tmp)$size, 1000)
})

test_that("ggplot download PDF succeeds for ggrepel / coord_equal plots", {
    skip_if_not_installed("ggrepel")
    set.seed(1)
    df <- data.frame(x = rnorm(30), y = rnorm(30), lab = paste0("g", 1:30))
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y, label = lab)) +
        ggplot2::geom_point() +
        ggrepel::geom_text_repel(max.overlaps = Inf, seed = 1L) +
        ggplot2::coord_equal() +
        ggplot2::theme_void()
    tmp <- tempfile(fileext = ".pdf")
    on.exit(unlink(tmp), add = TRUE)
    pdf_sz <- DEGreports:::.html_display_to_pdf_in(6.5, 5.5, dpi = 150)
    expect_error(
        DEGreports:::.save_ggplot_html_matched(p, tmp, 6.5, 5.5, pdf_sz$width, pdf_sz$height),
        NA
    )
    expect_true(file.exists(tmp) && file.info(tmp)$size > 1000)
})
test_that("pheatmap download PDF fits a smaller page without relying on native clip", {
    skip_if_not_installed("pheatmap")
    hm <- pheatmap::pheatmap(matrix(stats::rnorm(36), 6, 6), silent = TRUE,
                             cellwidth = 20, cellheight = 20)
    tmp <- tempfile(fileext = ".pdf")
    on.exit(unlink(tmp), add = TRUE)
    DEGreports:::.save_gtable_html_matched(hm, tmp, pdf_w = 3, pdf_h = 3)
    box <- .pdf_mediabox_inches(tmp)
    skip_if(is.null(box), "could not parse PDF MediaBox")
    expect_equal(box$width, 3, tolerance = 0.05)
    expect_equal(box$height, 3, tolerance = 0.05)
    expect_gt(file.info(tmp)$size, 1000)
})

test_that("download_plot_button accepts dpi and widget arguments", {
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    btn <- download_plot_button(p, filename = "my_plot", width = 4, height = 3, dpi = 150)
    html_txt <- as.character(htmltools::as.tags(btn))
    expect_true(grepl("Download PDF", html_txt))
    expect_true(grepl("base64", html_txt))
})

test_that("download_plot_button rejects non-ggplot input", {
    expect_error(download_plot_button(list(a = 1), filename = "x"), "must be a ggplot object")
})

test_that("download_grid_plot_button wraps a non-ggplot grid-drawable (e.g. a pheatmap grob) as a PDF download button", {
    skip_if_not_installed("pheatmap")
    p <- pheatmap::pheatmap(matrix(stats::rnorm(16), 4, 4), silent = TRUE)
    btn <- download_grid_plot_button(p, filename = "my_heatmap")
    html_txt <- as.character(htmltools::as.tags(btn))
    expect_true(grepl("Download PDF", html_txt))
    expect_true(grepl("base64", html_txt))
})

# --- Top-right corner button positioning helpers -----------------------

test_that("capture_plot_tag rasterizes a ggplot as a self-contained base64 <img> tag", {
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    tag <- capture_plot_tag(p, width = 4, height = 3, alt = "my_plot")
    expect_true(inherits(tag, "shiny.tag") || inherits(tag, "shiny.tag.list") || inherits(tag, "html"))
    html_txt <- as.character(htmltools::as.tags(tag))
    expect_true(grepl("<img", html_txt))
    expect_true(grepl("base64", html_txt))
    expect_true(grepl("my_plot", html_txt))
})

test_that("capture_plot_tag pins CSS width/height so the report theme cannot squash the figure", {
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    tag <- capture_plot_tag(p, width = 10, height = 5, dpi = 150)
    html_txt <- as.character(htmltools::as.tags(tag))
    expect_true(grepl('width="1500"', html_txt, fixed = TRUE))
    expect_true(grepl('height="750"', html_txt, fixed = TRUE))
    expect_true(grepl("width: 1500px", html_txt, fixed = TRUE))
    expect_true(grepl("height: 750px", html_txt, fixed = TRUE))
    expect_true(grepl("max-width: none", html_txt, fixed = TRUE))
})

test_that("capture_plot_tag defaults to a high pixelratio for crisp on-page rendering", {
    # Regression: the previous hardcoded pixelratio = 1 produced visibly
    # blurrier static figures than every interactive plotly widget in the
    # same report (which render at native screen resolution).
    expect_equal(formals(capture_plot_tag)$pixelratio, 2)
    expect_equal(formals(capture_plot_tag)$dpi, 150)
})

test_that("emit_plot_with_download defaults to the same higher dpi as capture_plot_tag", {
    expect_equal(formals(emit_plot_with_download)$dpi, 150)
})

test_that("wrap_with_corner_button places the button in a right-aligned row directly above the content (normal document flow, not an overlay)", {
    df <- data.frame(x = 1:3)
    btn <- download_table_button(df, filename = "data")
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    content <- capture_plot_tag(p, width = 4, height = 3, alt = "plot")
    wrapped <- wrap_with_corner_button(content, btn)
    html_txt <- as.character(htmltools::as.tags(wrapped))
    # button row is a normal in-flow flex row, right-aligned -- NOT an
    # absolutely-positioned overlay on top of the figure (regression:
    # buttons used to sit on top of/inside the plot).
    expect_true(grepl("justify-content: flex-end", html_txt))
    expect_false(grepl("position: absolute", html_txt))
    # the button row's markup must precede the figure's <img> markup, so
    # it renders above the figure rather than below or on top of it.
    expect_true(regexpr("Download CSV", html_txt)[1] < regexpr("<img", html_txt)[1])
    expect_true(grepl("<img", html_txt))
    expect_true(grepl("Download CSV", html_txt))
})

test_that("wrap_with_corner_button aligns PDF and CSV buttons on the same baseline", {
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    btn_pdf <- download_plot_button(p, filename = "plot", width = 4, height = 3)
    btn_csv <- link_output_file("scores.csv", label = "Download scores CSV")
    wrapped <- wrap_with_corner_button(htmltools::tags$div("content"), list(btn_pdf, btn_csv))
    html_txt <- as.character(htmltools::as.tags(wrapped))
    expect_true(grepl("degreports-dl-btn-wrap", html_txt))
    expect_true(grepl("align-items: center", html_txt))
    expect_true(grepl("Download PDF", html_txt))
    expect_true(grepl("Download scores CSV", html_txt))
})

test_that("wrap_with_corner_button places multiple buttons side by side in the row above the content", {
    df <- data.frame(x = 1:3)
    btn1 <- download_table_button(df, filename = "data1", button_label = "First button")
    btn2 <- download_table_button(df, filename = "data2", button_label = "Second button")
    wrapped <- wrap_with_corner_button(htmltools::tags$div("content"), list(btn1, btn2))
    html_txt <- as.character(htmltools::as.tags(wrapped))
    expect_true(grepl("First button", html_txt))
    expect_true(grepl("Second button", html_txt))
    expect_true(grepl("flex-direction: row", html_txt))
    expect_true(grepl("degreports-dl-row", html_txt))
    expect_true(grepl("height: 28px", html_txt))
})

test_that("wrap_with_corner_button centers static content by default but spans full width when center = FALSE", {
    content <- htmltools::tags$div("content")
    btn <- download_table_button(data.frame(x = 1), filename = "data")
    centered <- wrap_with_corner_button(content, btn, center = TRUE)
    full_width <- wrap_with_corner_button(content, btn, center = FALSE)
    centered_txt <- as.character(htmltools::as.tags(centered))
    full_txt <- as.character(htmltools::as.tags(full_width))
    expect_true(grepl("justify-content: center", centered_txt))
    expect_true(grepl("display: inline-block", centered_txt))
    expect_false(grepl("justify-content: center", full_txt))
    expect_true(grepl("width: 100%", full_txt))
})

test_that("emit_plot_with_download cats a self-contained figure + corner PDF button", {
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    html_txt <- paste(capture.output(emit_plot_with_download(p, filename = "my_plot", width = 4, height = 3)), collapse = "\n")
    expect_true(grepl("<img", html_txt))
    expect_true(grepl("Download PDF", html_txt))
    expect_true(grepl("justify-content: flex-end", html_txt))
})

test_that("emit_plot_with_download works with a non-ggplot grid-drawable (e.g. pheatmap)", {
    skip_if_not_installed("pheatmap")
    hm <- pheatmap::pheatmap(matrix(stats::rnorm(16), 4, 4), silent = TRUE)
    html_txt <- paste(capture.output(emit_plot_with_download(hm, filename = "my_heatmap", width = 4, height = 4)), collapse = "\n")
    expect_true(grepl("<img", html_txt))
    expect_true(grepl("Download PDF", html_txt))
})

test_that("emit_widget_with_download cats an interactive widget + corner button, full width (not centered)", {
    skip_if_not_installed("plotly")
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    widget <- plotly::ggplotly(p)
    btn <- download_plot_button(p, filename = "my_plot")
    html_txt <- paste(capture.output(emit_widget_with_download(widget, btn)), collapse = "\n")
    expect_true(grepl("justify-content: flex-end", html_txt))
    expect_true(grepl("width: 100%", html_txt))
    expect_false(grepl("text-align: center", html_txt))
})

test_that("emit_widget renders tag HTML without indentation (regression: nested divs must not trigger Pandoc's indented-code-block detection)", {
    # htmltools::as.character() pretty-prints nested tags with 2-space
    # indentation per nesting level. Once wrap_with_corner_button's nested
    # <div>s (outer wrapper > inner flex container > button row) are
    # serialized this way, lines reach >=4 leading spaces, which Pandoc's
    # Markdown parser treats as an indented code block -- silently
    # escaping the whole figure/button into literal <pre><code> text
    # instead of live HTML. emit_widget() must render with
    # indent = FALSE so no line is ever indented, regardless of nesting
    # depth.
    skip_if_not_installed("plotly")
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    widget <- plotly::ggplotly(p)
    btn <- download_plot_button(p, filename = "my_plot")
    wrapped <- wrap_with_corner_button(widget, btn, center = FALSE)
    html_txt <- paste(capture.output(emit_widget(wrapped)), collapse = "\n")
    lines <- strsplit(html_txt, "\n")[[1]]
    leading_spaces <- vapply(lines, function(l) attr(regexpr("^ *", l), "match.length"), integer(1))
    expect_true(all(leading_spaces == 0), info = paste("Indented line(s) found at positions:", paste(which(leading_spaces > 0), collapse = ", ")))
})

test_that("emit_widget output survives a real Pandoc Markdown round-trip without being escaped into a code block (end-to-end regression for the download-button-as-code-block bug)", {
    skip_if_not_installed("plotly")
    skip_if(!rmarkdown::pandoc_available(), "pandoc not available")
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    widget <- plotly::ggplotly(p)
    btn <- download_plot_button(p, filename = "my_plot")
    wrapped <- wrap_with_corner_button(widget, btn, center = FALSE)

    rmd_lines <- c(
        "---",
        "title: repro",
        "output: html_document",
        "---",
        "",
        "```{r setup, include=FALSE}",
        "knitr::opts_chunk$set(echo = FALSE)",
        "```",
        "",
        "## Section",
        "",
        "```{r plot1, results='asis'}",
        "DEGreports:::emit_widget(wrapped)",
        "```"
    )
    tmp_rmd <- tempfile(fileext = ".Rmd")
    writeLines(rmd_lines, tmp_rmd)
    out_dir <- tempfile()
    dir.create(out_dir)
    out_html <- rmarkdown::render(
        tmp_rmd, output_dir = out_dir, quiet = TRUE,
        envir = list2env(list(wrapped = wrapped))
    )
    html_txt <- paste(readLines(out_html), collapse = "\n")
    expect_false(grepl("<pre><code>", html_txt, fixed = TRUE))
    expect_true(grepl("justify-content: flex-end", html_txt, fixed = TRUE))
})

# Report figures must not be inlined as base64 while a page is rendering:
# rmarkdown's non-self-contained post-processor pushes every `<img src>`
# through utils::URLdecode(), which is quadratic in the src length, so a
# base64 payload turns a finished render into an apparent hang.

test_that("figure_img_tag inlines base64 when no report figure directory is configured", {
    skip_if_not_installed("ggplot2")
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    withr::local_options(list(degreports.figure_dir = NULL,
                              degreports.figure_href = NULL))
    html_txt <- as.character(htmltools::as.tags(capture_plot_tag(p, width = 2, height = 2)))
    expect_true(grepl("src=\"data:image/png;base64,", html_txt, fixed = TRUE))
})

test_that("figure_img_tag writes a PNG beside the page and references it by relative path", {
    skip_if_not_installed("ggplot2")
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    page_dir <- withr::local_tempdir()
    withr::local_options(list(
        degreports.figure_dir = file.path(page_dir, "qc_figures"),
        degreports.figure_href = "qc_figures"
    ))
    html_txt <- as.character(htmltools::as.tags(capture_plot_tag(p, width = 2, height = 2)))

    expect_false(grepl("base64", html_txt, fixed = TRUE))
    expect_match(html_txt, 'src="qc_figures/fig-[^"]+\\.png"')
    written <- list.files(file.path(page_dir, "qc_figures"), pattern = "\\.png$")
    expect_length(written, 1L)
    # The whole point: the src stays short enough for URLdecode() to be free.
    src <- sub('^.*src="', "", html_txt)
    expect_lt(nchar(sub('".*$', "", src)), 100L)
})

test_that("a PDF download button writes the PDF beside the page and links it", {
    skip_if_not_installed("ggplot2")
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    page_dir <- withr::local_tempdir()
    withr::local_options(list(
        degreports.figure_dir = file.path(page_dir, "qc_figures"),
        degreports.figure_href = "qc_figures"
    ))
    html_txt <- as.character(htmltools::as.tags(
        download_plot_button(p, filename = "pca_plot", width = 3, height = 2)
    ))
    # Regression: the PDF went in as base64 inside an onclick handler, which
    # was ~2.6 MB per DESeq2 contrast page on top of the identical file the
    # pipeline had already written in the same folder.
    expect_false(grepl("base64", html_txt, fixed = TRUE))
    expect_match(html_txt, 'href="qc_figures/pca_plot.pdf"', fixed = TRUE)
    expect_match(html_txt, "Download PDF")
    expect_true(file.exists(file.path(page_dir, "qc_figures", "pca_plot.pdf")))
})

test_that("two download buttons on one page cannot claim the same file name", {
    skip_if_not_installed("ggplot2")
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    page_dir <- withr::local_tempdir()
    withr::local_options(list(
        degreports.figure_dir = file.path(page_dir, "figs"),
        degreports.figure_href = "figs"
    ))
    href <- function() {
        txt <- as.character(htmltools::as.tags(download_plot_button(p, filename = "plot")))
        sub('".*$', "", sub('^.*href="', "", txt))
    }
    expect_equal(href(), "figs/plot.pdf")
    expect_equal(href(), "figs/plot-2.pdf")
    expect_length(list.files(file.path(page_dir, "figs"), pattern = "\\.pdf$"), 2L)
})

test_that("a PDF download button still embeds itself outside a report render", {
    skip_if_not_installed("ggplot2")
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    withr::local_options(list(degreports.figure_dir = NULL,
                              degreports.figure_href = NULL))
    html_txt <- as.character(htmltools::as.tags(
        download_plot_button(p, filename = "pca_plot")
    ))
    expect_true(grepl("base64", html_txt, fixed = TRUE))
})

test_that("figure_img_tag keeps the CSS box pinned to the requested display size", {
    skip_if_not_installed("ggplot2")
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    page_dir <- withr::local_tempdir()
    withr::local_options(list(
        degreports.figure_dir = file.path(page_dir, "figs"),
        degreports.figure_href = "figs"
    ))
    html_txt <- as.character(htmltools::as.tags(
        capture_plot_tag(p, width = 3, height = 2, dpi = 100)
    ))
    expect_match(html_txt, 'width="300"')
    expect_match(html_txt, 'height="200"')
    expect_match(html_txt, "width: 300px; height: 200px;", fixed = TRUE)
})

test_that("emit_plot_with_download can size the downloaded PDF independently of the HTML box", {
    skip_if_not_installed("ggplot2")
    html_only <- .html_display_to_pdf_in(16.4, 5.5, dpi = 150)
    overridden <- .override_pdf_in(html_only, pdf_width = 9, pdf_height = 12)
    expect_equal(overridden$width, 9)
    expect_equal(overridden$height, 12)
    # Overriding one dimension leaves the other tracking the HTML box.
    expect_equal(.override_pdf_in(html_only, pdf_height = 12)$width, html_only$width)
})

# A pheatmap gtable's matrix row/column are `npc`, so it has no natural size:
# sum(gt$heights) just reports back whatever device is current. Sizing a page
# from that measurement, or scaling the gtable's heights, leaves the layout
# larger than the page -- and gtable centres its layout, so grid clips the
# figure from both ends at once (column dendrogram off the top, column labels
# off the bottom).

pheatmap_fixture <- function(n_row = 40L, n_col = 8L) {
    set.seed(42)
    m <- matrix(stats::rnorm(n_row * n_col), nrow = n_row,
                dimnames = list(paste0("set_", seq_len(n_row)),
                                paste0("sample_", seq_len(n_col))))
    pheatmap::pheatmap(m, silent = TRUE)
}

test_that("gtable_page_size sizes an elastic pheatmap from content, not the current device", {
    skip_if_not_installed("pheatmap")
    hm <- pheatmap_fixture()

    measure_on <- function(w, h) {
        grDevices::pdf(NULL, width = w, height = h)
        on.exit(grDevices::dev.off())
        gtable_page_size(hm)
    }
    small <- measure_on(4, 3)
    large <- measure_on(20, 30)
    expect_equal(small$width, large$width)
    expect_equal(small$height, large$height)
    # Content-driven: one cell per row plus the fixed furniture.
    expect_gt(large$height, 40 * 0.16)
})

test_that("gtable_page_size gives a page the layout fits inside, so nothing is clipped", {
    skip_if_not_installed("pheatmap")
    hm <- pheatmap_fixture()
    size <- gtable_page_size(hm)
    expect_equal(size$fit, 1)
    grDevices::pdf(NULL, width = size$width, height = size$height)
    layout_h <- grid::convertHeight(sum(hm$gtable$heights), "in", valueOnly = TRUE)
    layout_w <- grid::convertWidth(sum(hm$gtable$widths), "in", valueOnly = TRUE)
    grDevices::dev.off()
    expect_lte(layout_h, size$height + 1e-6)
    expect_lte(layout_w, size$width + 1e-6)
})

test_that("gtable_page_size respects the caps and reports the scale needed to fit", {
    skip_if_not_installed("pheatmap")
    hm <- pheatmap_fixture(n_row = 300L)
    size <- gtable_page_size(hm, max_height = 6)
    expect_lte(size$height, 6)
    expect_lte(size$fit, 1)
})

test_that("plot_gsva_heatmap returns a gtable that fits its own page", {
    skip_if_not_installed("pheatmap")
    set.seed(7)
    scores <- matrix(stats::rnorm(30 * 6), nrow = 30,
                     dimnames = list(paste0("HALLMARK_SET_", seq_len(30)),
                                     paste0("s", seq_len(6))))
    hm <- plot_gsva_heatmap(scores, rep(c("a", "b"), each = 3),
                            collection = "HALLMARK")
    size <- gtable_page_size(hm)
    grDevices::pdf(NULL, width = size$width, height = size$height)
    layout_h <- grid::convertHeight(sum(hm$gtable$heights), "in", valueOnly = TRUE)
    grDevices::dev.off()
    expect_lte(layout_h, size$height + 1e-6)
    # The column dendrogram and column labels are on the page, not cropped off.
    expect_true(all(c("col_tree", "col_names") %in% hm$gtable$layout$name))
})

test_that("emit_heatmap_with_download emits the figure plus its download buttons", {
    skip_if_not_installed("pheatmap")
    hm <- pheatmap_fixture(n_row = 10L, n_col = 4L)
    page_dir <- withr::local_tempdir()
    withr::local_options(list(
        degreports.figure_dir = file.path(page_dir, "figs"),
        degreports.figure_href = "figs"
    ))
    html_txt <- paste(capture.output(
        emit_heatmap_with_download(hm, filename = "fixture_heatmap",
                                   buttons = list(htmltools::tags$span("extra")))
    ), collapse = "")

    expect_match(html_txt, 'src="figs/fig-[^"]+\\.png"')
    expect_match(html_txt, "Download PDF")
    expect_match(html_txt, "extra")
    expect_length(list.files(file.path(page_dir, "figs"), pattern = "\\.png$"), 1L)
})

test_that("gtable_page_size leaves room for a title wider than the matrix column", {
    skip_if_not_installed("pheatmap")
    set.seed(3)
    m <- matrix(stats::rnorm(20), nrow = 2,
                dimnames = list(c("condition", "batch"), paste0("PC", 1:10)))
    hm <- pheatmap::pheatmap(
        m, cluster_rows = FALSE, cluster_cols = FALSE, silent = TRUE,
        main = "-log10(p) PC\u2013covariate association (after ComBat-seq)"
    )
    bare <- pheatmap::pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE,
                               silent = TRUE)
    titled <- gtable_page_size(hm)
    untitled <- gtable_page_size(bare)
    expect_gt(titled$width, untitled$width)

    # The title is centred on the matrix column, so that column plus twice the
    # margin to its left has to cover the title.
    grDevices::pdf(NULL, width = titled$width, height = titled$height)
    widths <- as.numeric(grid::convertWidth(hm$gtable$widths, "in", valueOnly = TRUE))
    title_w <- grid::convertWidth(
        grid::grobWidth(hm$gtable$grobs[[which(hm$gtable$layout$name == "main")]]),
        "in", valueOnly = TRUE
    )
    grDevices::dev.off()
    mat <- hm$gtable$layout[hm$gtable$layout$name == "matrix", ]
    lead <- sum(widths[seq_len(mat$l[[1]] - 1L)])
    span <- sum(widths[seq(mat$l[[1]], mat$r[[1]])])
    expect_gte(span + 2 * lead, title_w)
})

test_that("gtable_page_size leaves room for the colour legend's stacked labels", {
    skip_if_not_installed("pheatmap")
    set.seed(4)
    m <- matrix(stats::rnorm(20), nrow = 2,
                dimnames = list(c("condition", "batch"), paste0("PC", 1:10)))
    hm <- pheatmap::pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE,
                             silent = TRUE)
    size <- gtable_page_size(hm)
    legend <- hm$gtable$grobs[[which(hm$gtable$layout$name == "legend")]]
    n_ticks <- max(vapply(
        Filter(function(k) !is.null(k$label), legend$children),
        function(k) length(k$label), integer(1)
    ))
    # Two heatmap rows alone would be well under an inch tall.
    expect_gt(size$height, n_ticks * 10 / 72)
})
