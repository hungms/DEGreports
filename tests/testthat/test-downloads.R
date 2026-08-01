# Tests for the download-button helpers (download_table_button,
# download_plot_button, download_grid_plot_button), all thin wrappers
# around downloadthis:: that must self-contain their payload so the
# bundled HTML report works fully offline as a single file.

skip_if_not_installed("downloadthis")

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

test_that("download_plot_button wraps a ggplot as a self-contained PDF download button", {
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    btn <- download_plot_button(p, filename = "my_plot")
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

test_that("capture_plot_tag emits real HTML width/height attributes (not a forced inline style) so the report's global img{max-width:100%} rule scales height in step with width", {
    # Regression: htmltools::plotTag()'s default inline
    # style="width:...px;height:...px" fights the report's global
    # `img { max-width: 100% !important; }` CSS rule -- width shrinks to
    # fit the container but the literal-pixel height in the style
    # attribute does not, squishing/stretching the image. HTML width/
    # height *attributes* (not a forced style) let the browser preserve
    # intrinsic aspect ratio when only width is constrained.
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    tag <- capture_plot_tag(p, width = 10, height = 5, dpi = 150)
    html_txt <- as.character(htmltools::as.tags(tag))
    expect_true(grepl('width="1500"', html_txt, fixed = TRUE))
    expect_true(grepl('height="750"', html_txt, fixed = TRUE))
    expect_false(grepl("style=\"width:", html_txt, fixed = TRUE))
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

test_that("wrap_with_corner_button places multiple buttons side by side in the row above the content", {
    df <- data.frame(x = 1:3)
    btn1 <- download_table_button(df, filename = "data1", button_label = "First button")
    btn2 <- download_table_button(df, filename = "data2", button_label = "Second button")
    wrapped <- wrap_with_corner_button(htmltools::tags$div("content"), list(btn1, btn2))
    html_txt <- as.character(htmltools::as.tags(wrapped))
    expect_true(grepl("First button", html_txt))
    expect_true(grepl("Second button", html_txt))
    expect_true(grepl("flex-direction: row", html_txt))
})

test_that("wrap_with_corner_button centers static content by default but spans full width when center = FALSE", {
    content <- htmltools::tags$div("content")
    btn <- download_table_button(data.frame(x = 1), filename = "data")
    centered <- wrap_with_corner_button(content, btn, center = TRUE)
    full_width <- wrap_with_corner_button(content, btn, center = FALSE)
    expect_true(grepl("text-align: center", as.character(htmltools::as.tags(centered))))
    expect_true(grepl("display: inline-block", as.character(htmltools::as.tags(centered))))
    expect_false(grepl("text-align: center", as.character(htmltools::as.tags(full_width))))
    expect_true(grepl("width: 100%", as.character(htmltools::as.tags(full_width))))
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
        "deseq2reports:::emit_widget(wrapped)",
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
