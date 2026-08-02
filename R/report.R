#' Registry of report pages: label, template, relative directory, and HTML file
#' @noRd
.report_pages <- function() {
    list(
        Index = list(
            key = "index", dir = ".", file = "index.html",
            template = "report_index.Rmd"
        ),
        QC = list(
            key = "qc", dir = "qc", file = "qc.html",
            template = "report_qc.Rmd"
        ),
        DESeq2 = list(
            key = "deseq2", dir = "deseq2", file = "deseq2.html",
            template = "report_deseq2.Rmd"
        ),
        clusterProfiler = list(
            key = "clusterProfiler", dir = "clusterprofiler", file = "clusterprofiler.html",
            template = "report_clusterProfiler.Rmd"
        ),
        WGCNA = list(
            key = "wgcna", dir = "wgcna", file = "wgcna.html",
            template = "report_wgcna.Rmd"
        ),
        decoupleR = list(
            key = "decoupleR", dir = "decoupler", file = "decoupler.html",
            template = "report_decoupleR.Rmd"
        ),
        VIPER = list(
            key = "viper", dir = "viper", file = "viper.html",
            template = "report_viper.Rmd"
        )
    )
}

#' Relative href from one report page directory to another
#' @noRd
.relative_href <- function(from_dir, to_dir, to_file) {
    if (identical(from_dir, to_dir)) {
        return(to_file)
    }
    if (identical(from_dir, ".")) {
        if (identical(to_dir, ".")) return(to_file)
        return(paste(to_dir, to_file, sep = "/"))
    }
    if (identical(to_dir, ".")) {
        return(paste("..", to_file, sep = "/"))
    }
    paste("..", to_dir, to_file, sep = "/")
}

#' Resolve the bundled Rmd template path for a report page
#' @noRd
.report_template_path <- function(template) {
    template_path <- system.file("rmarkdown", template, package = "deseq2reports")
    pkg_path <- tryCatch(find.package("deseq2reports"), error = function(e) "")
    src_candidate <- file.path(pkg_path, "rmarkdown", template)
    inst_candidate <- file.path(pkg_path, "inst", "rmarkdown", template)
    if (file.exists(inst_candidate)) {
        template_path <- inst_candidate
    } else if (file.exists(src_candidate)) {
        template_path <- src_candidate
    }
    if (!nzchar(template_path) || !file.exists(template_path)) {
        stop("Could not locate the bundled ", template, " in the installed package.",
             call. = FALSE)
    }
    template_path
}

#' Relative href map for the multi-page HTML report
#'
#' Returns a named character vector of relative paths from `from` to every
#' report page (including itself). Used by [report_nav_links()] and the
#' index page's section list so links stay correct when HTML files live in
#' per-module subfolders (`qc/qc.html`, `deseq2/deseq2.html`, ...).
#'
#' @param from One of `"Index"`, `"QC"`, `"DESeq2"`, `"clusterProfiler"`,
#'   `"WGCNA"`, `"decoupleR"`, `"VIPER"`.
#' @return Named character vector of hrefs, names = page labels.
#' @export
report_page_hrefs <- function(from = "Index") {
    pages <- .report_pages()
    if (!from %in% names(pages)) {
        stop("`from` must be one of: ", paste(names(pages), collapse = ", "),
             call. = FALSE)
    }
    from_dir <- pages[[from]]$dir
    stats::setNames(
        vapply(names(pages), function(label) {
            .relative_href(from_dir, pages[[label]]$dir, pages[[label]]$file)
        }, character(1)),
        names(pages)
    )
}

#' Build the shared top-of-page navigation bar for the multi-page HTML report
#'
#' Renders a single-line strip of cross-links between the report's HTML
#' pages (`index`, `qc`, `deseq2`, `clusterProfiler`, `wgcna`,
#' `decoupleR`, `viper`), used identically at the top of every page so a
#' reader can jump between sections. The current page is shown as plain
#' bold text (not a link to itself); any page named in `disabled` (e.g.
#' `"WGCNA"` when `wgcna_out` was `NULL` and that page was never rendered)
#' is shown greyed out instead of linked.
#'
#' @param current One of `"Index"`, `"QC"`, `"DESeq2"`, `"clusterProfiler"`,
#'   `"WGCNA"`, `"decoupleR"`, `"VIPER"` -- the page this bar is being
#'   placed on.
#' @param disabled Character vector of page labels (from the same set as
#'   `current`) that have no corresponding file and should be greyed out
#'   rather than linked. Default `character(0)`.
#' @param hrefs Optional named character vector of relative hrefs (as from
#'   [report_page_hrefs()]). Default `NULL` computes hrefs relative to
#'   `current`.
#' @return A length-1 character string of HTML (an already-closed `<div>`),
#'   suitable for `cat()`-ing into a `results='asis'` chunk.
#' @examples
#' cat(report_nav_links("QC"))
#' cat(report_nav_links("Index", disabled = "WGCNA"))
#' @export
report_nav_links <- function(current, disabled = character(0), hrefs = NULL) {
    pages <- .report_pages()
    if (!current %in% names(pages)) {
        stop("`current` must be one of: ", paste(names(pages), collapse = ", "),
             call. = FALSE)
    }
    if (is.null(hrefs)) {
        hrefs <- report_page_hrefs(current)
    } else {
        missing <- setdiff(names(pages), names(hrefs))
        if (length(missing)) {
            stop("`hrefs` is missing page label(s): ",
                 paste(missing, collapse = ", "), call. = FALSE)
        }
    }
    parts <- vapply(names(pages), function(label) {
        if (identical(label, current)) {
            sprintf("<strong>%s</strong>", label)
        } else if (label %in% disabled) {
            sprintf('<span style="color:#999;">%s</span>', label)
        } else {
            sprintf('<a href="%s">%s</a>', hrefs[[label]], label)
        }
    }, character(1))
    paste0(
        '<div style="margin-bottom: 1.5em; padding-bottom: 0.6em; ',
        'border-bottom: 1px solid #ddd; font-size: 1.05em;">',
        paste(parts, collapse = " &nbsp;|&nbsp; "),
        "</div>"
    )
}

#' Sanitize a string for use as an HTML id (vectorized)
#' @noRd
.html_id <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9]+", "-", x)
    x <- gsub("^-+|-+$", "", x)
    x[!nzchar(x) | is.na(x)] <- "tab"
    # HTML ids must start with a letter.
    needs_prefix <- grepl("^[0-9]", x)
    x[needs_prefix] <- paste0("id-", x[needs_prefix])
    tolower(x)
}

#' Open a Bootstrap 3 tabset (nav pills written into the HTML, not via `{.tabset}`)
#'
#' rmarkdown's `{.tabset}` conversion is client-side JS and often fails on
#' very large self-contained pages (e.g. WGCNA). Emitting the nav markup
#' directly makes tabs visible as soon as the HTML opens.
#'
#' @param tabset_id Unique HTML id for the tabset wrapper.
#' @param labels Character vector of tab labels (shown in the nav).
#' @return Character vector of pane ids (same length as `labels`), for use
#'   with [bs_tab_pane_start()].
#' @export
bs_tabset_start <- function(tabset_id, labels) {
    labels <- as.character(labels)
    if (!length(labels)) stop("`labels` must be non-empty.", call. = FALSE)
    tabset_id <- .html_id(tabset_id)
    pane_ids <- paste0(tabset_id, "-", .html_id(labels))
    # Ensure uniqueness if two labels sanitize identically.
    if (anyDuplicated(pane_ids)) {
        pane_ids <- paste0(pane_ids, "-", seq_along(pane_ids))
    }
    esc <- if (requireNamespace("htmltools", quietly = TRUE)) {
        htmltools::htmlEscape
    } else {
        function(x) gsub("&", "&amp;", gsub("<", "&lt;", gsub(">", "&gt;", x)))
    }
    cat('<div id="', tabset_id, '" class="deseq2reports-bs-tabset" style="margin: 0.8em 0 1.2em 0;">\n', sep = "")
    cat('<ul class="nav nav-tabs" role="tablist">\n')
    for (i in seq_along(labels)) {
        active <- if (i == 1L) ' class="active"' else ""
        cat(sprintf(
            '<li role="presentation"%s><a href="#%s" role="tab" data-toggle="tab">%s</a></li>\n',
            active, pane_ids[i], esc(labels[i])
        ))
    }
    cat('</ul>\n<div class="tab-content">\n')
    pane_ids
}

#' Open one Bootstrap 3 tab pane
#' @param pane_id Pane id returned by [bs_tabset_start()].
#' @param active Is this the initially visible pane? Default `FALSE`.
#' @export
bs_tab_pane_start <- function(pane_id, active = FALSE) {
    cls <- if (isTRUE(active)) "tab-pane fade in active" else "tab-pane fade"
    cat(sprintf('<div role="tabpanel" class="%s" id="%s" style="padding-top: 1em;">\n',
                cls, pane_id))
    invisible(NULL)
}

#' Close a Bootstrap 3 tab pane opened by [bs_tab_pane_start()]
#' @export
bs_tab_pane_end <- function() {
    cat("</div>\n")
    invisible(NULL)
}

#' Close a Bootstrap 3 tabset opened by [bs_tabset_start()]
#' @export
bs_tabset_end <- function() {
    cat("</div></div>\n")
    invisible(NULL)
}

#' Emit a fenced code block into a `results='asis'` report chunk
#'
#' Used to highlight a few well-known package functions under each report
#' section without changing the surrounding prose.
#'
#' @param ... Character strings (code lines). Vectors are flattened.
#' @param lang Fence language tag. Default `"r"`.
#' @param caption Optional one-line note printed above the fence
#'   (plain text, not markdown-bolded).
#' @return Invisibly `NULL`; called for its `cat()` side effect.
#' @export
emit_code_block <- function(..., lang = "r", caption = NULL) {
    lines <- unlist(list(...), use.names = FALSE)
    lines <- as.character(lines)
    if (!is.null(caption) && nzchar(caption)) {
        cat("\n\n", caption, "\n\n", sep = "")
    }
    cat("\n\n```", lang, "\n", sep = "")
    cat(paste(lines, collapse = "\n"))
    cat("\n```\n\n")
    invisible(NULL)
}

#' Shared params bag for every report page
#' @noRd
.report_shared_params <- function(
    dds,
    var,
    distance_out,
    diffexp_out,
    gsea_out,
    wgcna_out = NULL,
    decoupler_out = NULL,
    ppi_out = NULL,
    viper_out = NULL,
    alpha = 0.05,
    lfc_threshold = 0.5,
    run_args = NULL,
    disabled_labels = character(0)) {

    org_guess <- tryCatch(gsea_out$results$org[1], error = function(e) NA_character_)
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(decoupler_out$org, error = function(e) NA_character_)
    }
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(viper_out$org, error = function(e) NA_character_)
    }
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(ppi_out$org, error = function(e) NA_character_)
    }
    list(
        var = var,
        design_str = paste(deparse(DESeq2::design(dds)), collapse = " "),
        org = if (!is.null(org_guess) && !is.na(org_guess)) org_guess else NA_character_,
        n_samples = ncol(dds),
        n_genes = nrow(dds),
        qc_out = S4Vectors::metadata(dds)$deseq2reports_qc,
        distance_out = distance_out,
        diffexp_out = diffexp_out,
        gsea_out = gsea_out,
        wgcna_out = wgcna_out,
        decoupler_out = decoupler_out,
        ppi_out = ppi_out,
        viper_out = viper_out,
        alpha = alpha,
        lfc_threshold = lfc_threshold,
        run_args = run_args,
        package_version = as.character(utils::packageVersion("deseq2reports")),
        disabled_labels = disabled_labels
    )
}

#' Render one report HTML page into its module folder
#'
#' @param page Page label (`"Index"`, `"QC"`, ...) or key (`"index"`, `"qc"`, ...).
#' @param shared_params Params list from [.report_shared_params()].
#' @param save_dir Pipeline root directory (pages write into subfolders).
#' @return Absolute path to the rendered HTML file.
#' @noRd
.render_report_page <- function(page, shared_params, save_dir) {
    pages <- .report_pages()
    by_key <- stats::setNames(names(pages), vapply(pages, `[[`, character(1), "key"))
    if (page %in% names(by_key)) {
        label <- by_key[[page]]
    } else if (page %in% names(pages)) {
        label <- page
    } else {
        stop("`page` must be a report page label or key.", call. = FALSE)
    }
    spec <- pages[[label]]
    out_dir <- if (identical(spec$dir, ".")) {
        save_dir
    } else {
        validate_save_dir(file.path(save_dir, spec$dir))
    }
    page_params <- c(
        shared_params,
        list(
            current_page = label,
            page_hrefs = report_page_hrefs(label)
        )
    )
    rmarkdown::render(
        input = .report_template_path(spec$template),
        output_file = spec$file,
        output_dir = out_dir,
        params = page_params,
        envir = new.env(parent = globalenv()),
        quiet = TRUE
    )
    normalizePath(file.path(out_dir, spec$file), mustWork = TRUE)
}

#' Render the Index HTML page into `save_dir`
#' @inheritParams render_report
#' @param disabled_labels Page labels to grey out in the nav bar.
#' @return Path to `index.html`.
#' @export
render_index_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_out = NULL, ppi_out = NULL, viper_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_out, ppi_out, viper_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("Index", shared, save_dir)
}

#' Render the QC HTML page into `save_dir/qc/`
#' @inheritParams render_index_page
#' @return Path to `qc/qc.html`.
#' @export
render_qc_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_out = NULL, ppi_out = NULL, viper_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_out, ppi_out, viper_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("QC", shared, save_dir)
}

#' Render the DESeq2 HTML page into `save_dir/deseq2/`
#' @inheritParams render_index_page
#' @return Path to `deseq2/deseq2.html`.
#' @export
render_deseq2_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_out = NULL, ppi_out = NULL, viper_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_out, ppi_out, viper_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("DESeq2", shared, save_dir)
}

#' @rdname render_deseq2_page
#' @export
render_diffexp_page <- render_deseq2_page

#' Render the clusterProfiler HTML page into `save_dir/clusterprofiler/`
#' @inheritParams render_index_page
#' @return Path to `clusterprofiler/clusterprofiler.html`.
#' @export
render_clusterProfiler_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_out = NULL, ppi_out = NULL, viper_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_out, ppi_out, viper_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("clusterProfiler", shared, save_dir)
}

#' @rdname render_clusterProfiler_page
#' @export
render_pathway_page <- render_clusterProfiler_page

#' Render the WGCNA HTML page into `save_dir/wgcna/`
#' @inheritParams render_index_page
#' @return Path to `wgcna/wgcna.html`.
#' @export
render_wgcna_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_out = NULL, ppi_out = NULL, viper_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    if (is.null(wgcna_out)) {
        stop("`wgcna_out` must be supplied to render the WGCNA page.", call. = FALSE)
    }
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_out, ppi_out, viper_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("WGCNA", shared, save_dir)
}

#' Render the decoupleR HTML page into `save_dir/decoupler/`
#' @inheritParams render_index_page
#' @return Path to `decoupler/decoupler.html`.
#' @export
render_decoupleR_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_out = NULL, ppi_out = NULL, viper_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    if (is.null(decoupler_out)) {
        stop("`decoupler_out` must be supplied to render the decoupleR page.",
             call. = FALSE)
    }
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_out, ppi_out, viper_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("decoupleR", shared, save_dir)
}

#' @rdname render_decoupleR_page
#' @export
render_tf_activity_page <- render_decoupleR_page

#' Render the VIPER HTML page into `save_dir/viper/`
#' @inheritParams render_index_page
#' @return Path to `viper/viper.html`.
#' @export
render_viper_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_out = NULL, ppi_out = NULL, viper_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    if (is.null(viper_out)) {
        stop("`viper_out` must be supplied to render the VIPER page.",
             call. = FALSE)
    }
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_out, ppi_out, viper_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("VIPER", shared, save_dir)
}

#' Render the pipeline's multi-page HTML summary report
#'
#' Knits the package's bundled R Markdown templates against the outputs of
#' [run_qc()], [run_distance_pca()], [run_diffexp_pipeline()],
#' [run_gsea_pipeline()], [run_wgcna_pipeline()],
#' [run_decoupler_pipeline()], and [run_viper_pipeline()], writing each
#' section page into its module folder and `index.html` at the pipeline
#' root:
#' \itemize{
#'   \item `index.html`
#'   \item `qc/qc.html`
#'   \item `deseq2/deseq2.html`
#'   \item `clusterprofiler/clusterprofiler.html`
#'   \item `wgcna/wgcna.html` (omitted when `wgcna_out = NULL`)
#'   \item `decoupler/decoupler.html` (omitted when `decoupler_out = NULL`)
#'   \item `viper/viper.html` (omitted when `viper_out = NULL`)
#' }
#' Each file is self-contained (`self_contained = TRUE`); cross-links use
#' relative paths so the folder tree can be copied and opened offline.
#'
#' @param dds The (QC-filtered) `DESeqDataSet`, used to report sample/gene
#'   counts and the design formula in the report header. If it carries
#'   `S4Vectors::metadata(dds)$deseq2reports_qc` (attached automatically by
#'   [run_qc()]), a QC section is also included showing the number/percent
#'   of genes removed at each filtering step, the library-size plot, and
#'   the count-distribution diagnostic plot (with the expression cutoff
#'   marked) -- omitted if that metadata is absent, e.g. when `dds` was
#'   QC-filtered manually rather than via [run_qc()].
#' @param var Column name in `colData(dds)` used for grouping.
#' @param distance_out Output of [run_distance_pca()].
#' @param diffexp_out Output of [run_diffexp_pipeline()].
#' @param gsea_out Output of [run_gsea_pipeline()].
#' @param wgcna_out Output of [run_wgcna_pipeline()], or `NULL` (default)
#'   to omit the WGCNA page entirely (its link is greyed out, not removed,
#'   in the navigation bar of the other pages).
#' @param decoupler_out Output of [run_decoupler_pipeline()], or `NULL`
#'   (default) to omit the decoupleR page.
#' @param ppi_out Output of [run_ppi_pipeline()], or `NULL` (default) to
#'   omit PPI connectivity sections on the DE and WGCNA pages.
#' @param viper_out Output of [run_viper_pipeline()], or `NULL` (default)
#'   to omit the VIPER page.
#' @param alpha Significance threshold to display/use for report summaries.
#'   Default `0.05`.
#' @param lfc_threshold Log2 fold-change threshold used for up/down
#'   colouring on MA/volcano plots. Default `0.5`.
#' @param run_args Optional named list of pipeline arguments (as built by
#'   [run_deseq2_pipeline()]) to display on the index page. Default `NULL`.
#' @param save_dir Pipeline root directory. Section pages are written into
#'   per-module subfolders under this path. Default `getwd()`.
#' @return Invisibly, a named character vector of rendered file paths, with
#'   names `"index"`, `"qc"`, `"deseq2"`, `"clusterProfiler"`, and
#'   (when rendered) `"wgcna"` / `"decoupleR"` / `"viper"`.
#' @examples
#' \dontrun{
#' render_report(dds, var = "condition", distance_out = distance_out,
#'               diffexp_out = diffexp_out, gsea_out = gsea_out,
#'               save_dir = "results")
#' }
#' @export
render_report <- function(
    dds,
    var,
    distance_out,
    diffexp_out,
    gsea_out,
    wgcna_out = NULL,
    decoupler_out = NULL,
    ppi_out = NULL,
    viper_out = NULL,
    alpha = 0.05,
    lfc_threshold = 0.5,
    run_args = NULL,
    save_dir = getwd()) {

    # Every bundled Rmd template embeds interactive plotly/DT widgets via
    # htmltools::tagList() at render time, so htmltools must be installed
    # even though this R function itself never calls it directly.
    if (!requireNamespace("htmltools", quietly = TRUE)) {
        stop("Package 'htmltools' is required to render the HTML report.", call. = FALSE)
    }

    save_dir <- validate_save_dir(save_dir)

    disabled_labels <- character(0)
    if (is.null(wgcna_out)) disabled_labels <- c(disabled_labels, "WGCNA")
    if (is.null(decoupler_out)) disabled_labels <- c(disabled_labels, "decoupleR")
    if (is.null(viper_out)) disabled_labels <- c(disabled_labels, "VIPER")

    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_out, ppi_out, viper_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )

    # Module pages first, then index (so the landing page is written last).
    render_order <- c("qc", "deseq2", "clusterProfiler")
    if (!is.null(wgcna_out)) render_order <- c(render_order, "wgcna")
    if (!is.null(decoupler_out)) render_order <- c(render_order, "decoupleR")
    if (!is.null(viper_out)) render_order <- c(render_order, "viper")
    render_order <- c(render_order, "index")

    out_paths <- stats::setNames(character(length(render_order)), render_order)
    for (key in render_order) {
        out_paths[[key]] <- .render_report_page(key, shared, save_dir)
    }

    message("HTML report written to: ", out_paths[["index"]],
            " (+ ", length(out_paths) - 1,
            " linked section page(s) in per-module folders)")
    invisible(out_paths)
}
