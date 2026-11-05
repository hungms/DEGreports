#' Can `rmarkdown` see a pandoc new enough to knit the report?
#'
#' `rmarkdown` looks for pandoc on `PATH`, in `RSTUDIO_PANDOC`, and in the
#' handful of locations RStudio and Quarto install it. A batch job, an
#' `Rscript` invocation on a cluster node, or a container therefore often
#' has no pandoc even though the same package library renders fine from the
#' desktop.
#' @noRd
.rmarkdown_pandoc_ok <- function() {
    isTRUE(tryCatch(rmarkdown::pandoc_available("1.12.3"), error = function(e) FALSE))
}

#' Hand `rmarkdown` a pandoc installed by the **pandoc** package, if any
#'
#' pandoc is a system dependency (declared in `SystemRequirements`), and on
#' a cluster without root the practical way to satisfy it is the **pandoc**
#' package's `pandoc_install()`, which drops a binary in the user data
#' directory. That binary is invisible to `rmarkdown` until it is activated,
#' which would otherwise have to be repeated by hand in every session.
#' @return `TRUE` if a managed pandoc was activated.
#' @noRd
.activate_managed_pandoc <- function() {
    if (!requireNamespace("pandoc", quietly = TRUE)) return(FALSE)
    tryCatch({
        if (!length(pandoc::pandoc_installed_versions())) return(FALSE)
        pandoc::pandoc_activate(rmarkdown = TRUE, quiet = TRUE)
        TRUE
    }, error = function(e) FALSE)
}

#' Is a pandoc new enough to knit the report available to this session?
#' @noRd
.pandoc_ok <- function() {
    if (.rmarkdown_pandoc_ok()) return(TRUE)
    if (.activate_managed_pandoc()) return(.rmarkdown_pandoc_ok())
    FALSE
}

#' Explain a missing pandoc and the ways to supply one
#' @noRd
.pandoc_missing_msg <- function() {
    paste0(
        "pandoc (>= 1.12.3) was not found, so the HTML report cannot be knitted. ",
        "It is a system dependency (see SystemRequirements), and sessions outside ",
        "RStudio -- an HPC batch job, Rscript, a container -- do not inherit RStudio's ",
        "bundled copy. Install one system-wide (`module load pandoc`, ",
        "`conda install -c conda-forge pandoc`, or https://pandoc.org/installing.html); ",
        "without admin rights, ",
        "`install.packages(c(\"pandoc\", \"gh\")); pandoc::pandoc_install()` ",
        "puts a binary in your user directory that this package activates on its own ",
        "from then on. An existing binary elsewhere can be pointed at with ",
        "`rmarkdown::find_pandoc(dir = \"/path/to/pandoc/bin\")`, or by exporting ",
        "RSTUDIO_PANDOC before starting R."
    )
}

#' Warn up front when a run will not be able to render its report
#'
#' Called at the top of the pipeline entry points so the missing dependency
#' is reported in the first seconds of a run rather than after every
#' analysis step has finished.
#' @param render_html The caller's `render_html` argument.
#' @return Invisibly, `TRUE` when a warning was issued.
#' @noRd
.warn_if_no_pandoc <- function(render_html) {
    if (!isTRUE(render_html) || .pandoc_ok()) return(invisible(FALSE))
    warning(.pandoc_missing_msg(),
            " The pipeline will still run to completion and write every table and ",
            "figure; only the report at the end is skipped.", call. = FALSE)
    invisible(TRUE)
}

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
        GSVA = list(
            key = "gsva", dir = "gsva", file = "gsva.html",
            template = "report_gsva.Rmd"
        ),
        WGCNA = list(
            key = "wgcna", dir = "wgcna", file = "wgcna.html",
            template = "report_wgcna.Rmd"
        ),
        `decoupleR-TF` = list(
            key = "decoupleR-TF", dir = "decoupler-tf", file = "decoupleR-TF.html",
            template = "report_decoupleR_tf.Rmd"
        ),
        `decoupleR-CCC` = list(
            key = "decoupleR-CCC", dir = "decoupler-ccc", file = "decoupleR-CCC.html",
            template = "report_decoupleR_ccc.Rmd"
        ),
        VIPER = list(
            key = "viper", dir = "viper", file = "viper.html",
            template = "report_viper.Rmd"
        ),
        Multihit = list(
            key = "multihit", dir = "multihit", file = "multihit.html",
            template = "report_multihit.Rmd"
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

#' Relative href map from a per-contrast page (possibly nested)
#'
#' @param module_dir Module folder name (e.g. `"deseq2"`).
#' @param nest_depth Folder depth under the module (`1` for flat contrasts,
#'   `2` for three-factor `comparison/stratum`).
#' @noRd
.contrast_page_hrefs <- function(module_dir, nest_depth = 1L) {
    nest_depth <- max(1L, as.integer(nest_depth)[1L])
    up_to_module <- paste(rep("..", nest_depth), collapse = "/")
    up_to_root <- paste(rep("..", nest_depth + 1L), collapse = "/")
    pages <- .report_pages()
    stats::setNames(
        vapply(names(pages), function(lab) {
            to <- pages[[lab]]
            if (identical(to$dir, ".")) {
                return(paste(up_to_root, to$file, sep = "/"))
            }
            if (identical(to$dir, module_dir)) {
                return(paste(up_to_module, to$file, sep = "/"))
            }
            paste(up_to_root, to$dir, to$file, sep = "/")
        }, character(1)),
        names(pages)
    )
}

.render_contrast_report_page <- function(module_key, contrast, shared_params, save_dir) {
    pages <- .report_pages()
    by_key <- stats::setNames(names(pages), vapply(pages, `[[`, character(1), "key"))
    label <- by_key[[module_key]]
    spec <- pages[[label]]
    out_dir <- validate_save_dir(file.path(save_dir, spec$dir, contrast))
    nest_depth <- .contrast_nest_depth(contrast)
    stem <- contrast_file_stem(contrast)
    page_params <- utils::modifyList(
        shared_params,
        list(
            current_page = label,
            page_hrefs = .contrast_page_hrefs(spec$dir, nest_depth = nest_depth),
            contrast = contrast
        )
    )
    .render_page_rmd(spec$template, paste0(stem, ".html"), out_dir, page_params)
}

#' Resolve the bundled Rmd template path for a report page
#' @noRd
.report_template_path <- function(template) {
    template_path <- system.file("rmarkdown", template, package = "DEGreports")
    pkg_path <- tryCatch(find.package("DEGreports"), error = function(e) "")
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

#' Knit one report template to HTML with private knitr intermediates
#'
#' `rmarkdown::render()` defaults `intermediates_dir` to `dirname(input)`, so
#' every page knitting from the same bundled template writes the same
#' `<template>.knit.md`/`.utf8.md` next to it. Two consequences, both bugs:
#' rendering per-contrast pages in parallel makes the workers overwrite each
#' other's `.knit.md` (every contrast HTML ends up byte-identical), and the
#' intermediates -- 19 MB for WGCNA -- land in the package source tree, which
#' stalls badly when that tree sits in a synced folder. Give each render its
#' own scratch directory instead.
#'
#' Also points [capture_plot_tag()] at a `figures/` directory beside the page
#' so its bitmaps are written there rather than inlined as base64; see
#' `figure_img_tag()` for why inlining makes `rmarkdown::render()` appear to
#' hang on figure-heavy pages. The directory is emptied first, so a re-run
#' replaces its figures instead of accumulating orphans.
#' @noRd
.render_page_rmd <- function(template, output_file, output_dir, page_params) {
    scratch <- tempfile("degreports-knit-")
    dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)

    fig_href <- paste0(tools::file_path_sans_ext(basename(output_file)), "_figures")
    fig_dir <- file.path(output_dir, fig_href)
    unlink(fig_dir, recursive = TRUE, force = TRUE)
    old_opts <- options(degreports.figure_dir = fig_dir,
                        degreports.figure_href = fig_href)
    on.exit(options(old_opts), add = TRUE)

    rmarkdown::render(
        input = .report_template_path(template),
        output_file = output_file,
        output_dir = output_dir,
        intermediates_dir = scratch,
        output_options = list(
            self_contained = FALSE,
            lib_dir = file.path(output_dir, "lib")
        ),
        params = page_params,
        envir = new.env(parent = globalenv()),
        quiet = TRUE
    )
    normalizePath(file.path(output_dir, output_file), mustWork = TRUE)
}

#' Give this (forked) process its own rmarkdown temp-file prefix
#'
#' rmarkdown writes header/body includes to `tempdir()` as
#' `rmarkdown-str<hex>.html` and, when a render finishes, deletes **every**
#' file matching that pattern there. Forked render workers share one
#' `tempdir()`, so a worker finishing deleted includes another worker's pandoc
#' was about to read ("withBinaryFile: does not exist"). A per-process prefix
#' makes each worker's cleanup match only its own files. Call only inside a
#' forked worker: it changes rmarkdown's namespace for the rest of that process.
#' @noRd
.isolate_rmarkdown_tmpfiles <- function() {
    ns <- asNamespace("rmarkdown")
    if (!exists("tmpfile_pattern", envir = ns, inherits = FALSE)) return(invisible(FALSE))
    utils::assignInNamespace("tmpfile_pattern",
                             paste0("rmarkdown-str-p", Sys.getpid(), "-"),
                             ns = "rmarkdown")
    invisible(TRUE)
}

#' Did a per-contrast render return a written HTML file?
#'
#' A worker can come back with a condition object (its render errored), a
#' `try-error` (`mclapply()`'s wrapper, including a child killed outright), or
#' `NULL`. Only a path to a file that exists counts.
#' @noRd
.rendered_ok <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x) && file.exists(x)
}

#' One-line reason from whatever a failed render came back as
#' @noRd
.render_failure_reason <- function(x) {
    if (inherits(x, "condition")) return(conditionMessage(x))
    if (inherits(x, "try-error")) {
        return(sub("^Error[^:]*:\\s*", "", trimws(as.character(x))))
    }
    if (is.null(x)) return("worker returned nothing (child process died)")
    "no HTML file was written"
}

#' Relative href map for the multi-page HTML report
#'
#' Returns a named character vector of relative paths from `from` to every
#' report page (including itself). Used by [report_nav_links()] and the
#' index page's section list so links stay correct when HTML files live in
#' per-module subfolders (`qc/qc.html`, `deseq2/deseq2.html`, ...).
#'
#' @param from One of `"Index"`, `"QC"`, `"DESeq2"`, `"clusterProfiler"`,
#'   `"GSVA"`, `"WGCNA"`, `"decoupleR-TF"`, `"decoupleR-CCC"`, `"VIPER"`,
#'   `"Multihit"`.
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
#' pages (`index`, `qc`, `deseq2`, `clusterProfiler`, `gsva`, `wgcna`, `decoupleR-TF`,
#' `decoupleR-CCC`, `viper`, `multihit`), used identically at the top of every page so a reader can jump
#' between sections. The current page is shown as plain bold text (not a
#' link to itself); any page named in `disabled` (e.g. `"WGCNA"` when
#' `wgcna_out` was `NULL` and that page was never rendered) is shown
#' greyed out instead of linked.
#'
#' @param current One of `"Index"`, `"QC"`, `"DESeq2"`, `"clusterProfiler"`,
#'   `"GSVA"`, `"WGCNA"`, `"decoupleR-TF"`, `"decoupleR-CCC"`, `"VIPER"`,
#'   `"Multihit"` -- the page this bar is being
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
    cat('<div id="', tabset_id, '" class="degreports-bs-tabset" style="margin: 0.8em 0 1.2em 0;">\n', sep = "")
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

#' Bullet list of links to per-contrast pages
#'
#' A stem like `CB__noK` would otherwise leave its `**` markers literal:
#' pandoc reads the `__` as an unclosed strong emphasis of its own.
#' @noRd
.emit_contrast_links <- function(labels) {
    for (comp in labels) {
        stem <- contrast_file_stem(comp)
        href <- file.path(comp, paste0(stem, ".html"))
        cat("- [**", gsub("_", "\\_", stem, fixed = TRUE), "**](", href, ")\n",
            sep = "")
    }
    invisible(NULL)
}

#' Emit a hub page's navigation to its per-contrast pages
#'
#' Groups the links by contrast class -- main effects, then each main
#' effect's within-stratum children, then the interaction contrasts -- so a
#' multi-factor grid does not arrive as one long undifferentiated list. A
#' single-factor (or LRT) grid has nothing to group and is emitted flat.
#'
#' @param specs Contrast specs from [expand_contrast_specs()], subset to the
#'   contrasts that actually have a page.
#' @param var The pipeline's `var`, for the section titles.
#' @param heading_level Markdown heading level for the group titles, so a
#'   caller under `### Contrasts` passes `4` (the default).
#' @return Invisibly `NULL`; called for the markdown it prints.
#' @export
emit_contrast_nav <- function(specs, var, heading_level = 4L) {
    if (is.null(specs) || !nrow(specs)) return(invisible(NULL))
    labels <- as.character(specs$label)
    cls <- if ("class" %in% names(specs)) {
        as.character(specs$class)
    } else {
        rep(NA_character_, length(labels))
    }
    nested <- grepl("/", labels, fixed = TRUE)
    if (!any(nested) && !any(cls %in% c("main", "interaction"))) {
        .emit_contrast_links(labels)
        cat("\n")
        return(invisible(NULL))
    }

    hashes <- paste(rep("#", max(1L, as.integer(heading_level))), collapse = "")
    heading <- function(text) cat("\n\n", hashes, " ", text, "\n\n", sep = "")

    main <- labels[!nested & cls %in% "main"]
    if (length(main)) {
        heading(contrast_class_title("main", var))
        .emit_contrast_links(main)
    }
    for (parent in unique(dirname(labels[nested]))) {
        heading(paste0(parent, " \u2014 by stratum"))
        .emit_contrast_links(labels[nested & dirname(labels) == parent])
    }
    interaction <- labels[!nested & cls %in% "interaction"]
    if (length(interaction)) {
        heading(contrast_class_title("interaction", var))
        .emit_contrast_links(interaction)
    }
    # Anything left is a within-stratum contrast with a flat label (the
    # two-factor convention): title each group by its class rather than
    # dropping them all under one "Other" heading.
    others <- setdiff(labels[!nested], c(main, interaction))
    if (length(others)) {
        other_cls <- cls[match(others, labels)]
        for (cl in unique(other_cls)) {
            title <- contrast_class_title(cl, var)
            if (is.na(title) || identical(title, cl)) title <- "Other"
            heading(title)
            .emit_contrast_links(others[other_cls %in% cl])
        }
    }
    cat("\n")
    invisible(NULL)
}

#' A comparison's pooled side, read as prose
#'
#' `"CB+CC"` becomes `"the average of CB and CC"`, since a one-vs-rest
#' interaction weights the pooled levels equally rather than pooling their
#' samples. A single level is returned unchanged.
#' @noRd
.ix_side_phrase <- function(x) {
    parts <- strsplit(as.character(x), "+", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    if (length(parts) < 2L) return(paste(parts, collapse = ""))
    paste0("the average of ", paste(parts[-length(parts)], collapse = ", "),
           " and ", parts[[length(parts)]])
}

#' Split an `ix_terms` cell into its tested factors and their comparisons
#'
#' `"celltype:CA_vs_CB+CC|K_mutation:K_vs_noK"` gives one entry per tested
#' factor, each with the column name, the comparison label, and its two
#' sides.
#' @noRd
.ix_parse_terms <- function(x) {
    parts <- strsplit(as.character(x), "|", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    lapply(parts, function(p) {
        comp <- sub("^[^:]*:", "", p)
        sides <- strsplit(comp, "_vs_", fixed = TRUE)[[1]]
        list(
            v = sub(":.*$", "", p),
            comp = comp,
            hi = sides[[1L]],
            lo = if (length(sides) > 1L) paste(sides[-1L], collapse = "_vs_") else ""
        )
    })
}

#' What question each interaction contrast answers
#'
#' A difference of differences is easy to misread, so the hub page states each
#' interaction row as a question and says what a positive `log2FoldChange`
#' means for it. Rows with two tested factors read "does the K effect differ
#' between C and noC, within CB"; the three-way row asks whether that
#' interaction itself differs across the third factor. The reading is
#' symmetric in the tested pair, which is why the "positive" column gives
#' both directions.
#'
#' @param specs Contrast specs from [expand_contrast_specs()]. Rows that are
#'   not interactions (or predate the `ix_terms` column) are dropped.
#' @param link Logical. If `TRUE` (default), the contrast column is a
#'   markdown link to that contrast's page, relative to the hub.
#' @return A data frame with columns `contrast`, `question`, `positive`; zero
#'   rows when `specs` holds no interaction contrasts.
#' @export
interaction_question_table <- function(specs, link = TRUE) {
    empty <- data.frame(contrast = character(0), question = character(0),
                        positive = character(0), stringsAsFactors = FALSE)
    if (is.null(specs) || !nrow(specs)) return(empty)
    if (!all(c("label", "class", "ix_terms") %in% names(specs))) return(empty)
    ix <- specs[specs$class %in% "interaction" & !is.na(specs$ix_terms), ,
                drop = FALSE]
    if (!nrow(ix)) return(empty)

    held <- contrast_stratum_labels(ix)
    code <- function(x) paste0("`", x, "`")
    rows <- lapply(seq_len(nrow(ix)), function(i) {
        terms <- .ix_parse_terms(ix$ix_terms[[i]])
        label <- as.character(ix$label[[i]])
        stem <- contrast_file_stem(label)
        first <- terms[[1L]]
        between <- paste0(first$v, " ", first$hi, " and ", .ix_side_phrase(first$lo))
        if (length(terms) >= 3L) {
            question <- paste0(
                "Does the ", terms[[2L]]$v, " \u00d7 ", terms[[3L]]$v,
                " interaction (", code(terms[[2L]]$comp), " \u00d7 ",
                code(terms[[3L]]$comp), ") differ between ", between, "?"
            )
            positive <- paste0(
                "That interaction is larger in ", first$v, " ", first$hi,
                " than in ", .ix_side_phrase(first$lo), "."
            )
        } else if (length(terms) == 2L) {
            within <- if (is.na(held[[i]])) "" else paste0(", within ", held[[i]])
            question <- paste0(
                "Does the ", terms[[2L]]$v, " effect (", code(terms[[2L]]$comp),
                ") differ between ", between, within, "?"
            )
            positive <- paste0(
                code(terms[[2L]]$comp), " is more up in ", first$v, " ", first$hi,
                " than in ", .ix_side_phrase(first$lo), " \u2014 equivalently, ",
                code(first$comp), " is more up in ", terms[[2L]]$v, " ",
                terms[[2L]]$hi, " than in ", .ix_side_phrase(terms[[2L]]$lo), "."
            )
        } else {
            question <- NA_character_
            positive <- NA_character_
        }
        data.frame(
            contrast = if (isTRUE(link)) {
                paste0("[", code(stem), "](", file.path(label, paste0(stem, ".html")), ")")
            } else {
                label
            },
            question = question,
            positive = positive,
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out <- out[!is.na(out$question), , drop = FALSE]
    rownames(out) <- NULL
    out
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
    decoupler_tf_out = NULL,
    ppi_out = NULL,
    viper_out = NULL,
    cspa_out = NULL,
    alpha = 0.05,
    lfc_threshold = 0.5,
    run_args = NULL,
    disabled_labels = character(0),
    gsea_overlap_out = NULL,
    multihit_out = NULL,
    degpatterns_out = NULL,
    cluster_ora_out = NULL,
    gsva_out = NULL,
    decoupler_ccc_out = NULL) {

    org_guess <- tryCatch(gsea_out$results$org[1], error = function(e) NA_character_)
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(gsea_out$org, error = function(e) NA_character_)
    }
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(cluster_ora_out$org, error = function(e) NA_character_)
    }
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(decoupler_tf_out$org, error = function(e) NA_character_)
    }
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(viper_out$org, error = function(e) NA_character_)
    }
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(ppi_out$org, error = function(e) NA_character_)
    }
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(cspa_out$org, error = function(e) NA_character_)
    }
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(gsva_out$org, error = function(e) NA_character_)
    }
    if (is.null(org_guess) || is.na(org_guess)) {
        org_guess <- tryCatch(multihit_out$org, error = function(e) NA_character_)
    }
    if ((is.null(org_guess) || is.na(org_guess)) &&
            !is.null(diffexp_out$results) &&
            "gene" %in% colnames(diffexp_out$results)) {
        org_guess <- tryCatch(
            detect_organism(as.character(diffexp_out$results$gene)),
            error = function(e) NA_character_
        )
    }

    if (is.null(dds)) {
        design_str <- "n/a (DEG table input; no DESeqDataSet)"
        n_samples <- NA_integer_
        n_genes <- if (!is.null(diffexp_out$results) &&
                       "gene" %in% colnames(diffexp_out$results)) {
            length(unique(as.character(diffexp_out$results$gene)))
        } else {
            NA_integer_
        }
        qc_out <- NULL
    } else {
        design_str <- paste(deparse(DESeq2::design(dds)), collapse = " ")
        n_samples <- ncol(dds)
        n_genes <- nrow(dds)
        qc_out <- S4Vectors::metadata(dds)$degreports_qc
    }

    list(
        var = var,
        design_str = design_str,
        org = if (!is.null(org_guess) && !is.na(org_guess)) org_guess else NA_character_,
        n_samples = n_samples,
        n_genes = n_genes,
        qc_out = qc_out,
        distance_out = distance_out,
        diffexp_out = diffexp_out,
        gsea_out = gsea_out,
        gsea_overlap_out = gsea_overlap_out,
        multihit_out = multihit_out,
        degpatterns_out = degpatterns_out,
        cluster_ora_out = cluster_ora_out,
        wgcna_out = wgcna_out,
        decoupler_tf_out = decoupler_tf_out,
        ppi_out = ppi_out,
        viper_out = viper_out,
        cspa_out = cspa_out,
        gsva_out = gsva_out,
        decoupler_ccc_out = decoupler_ccc_out,
        alpha = alpha,
        lfc_threshold = lfc_threshold,
        run_args = run_args,
        package_version = as.character(utils::packageVersion("DEGreports")),
        disabled_labels = disabled_labels
    )
}

#' Drop interactive plotly handles from a distance/PCA result (keeps tables)
#' @noRd
.strip_distance_interactive <- function(distance_out) {
    if (is.null(distance_out)) return(NULL)
    drop <- c(
        "pca_plot_interactive", "pca_plot_uncorrected_interactive",
        "similarity_plot", "similarity_plot_corrected"
    )
    for (nm in intersect(drop, names(distance_out))) {
        distance_out[[nm]] <- NULL
    }
    distance_out
}

#' Drop ggplot / gtable objects from a WGCNA pipeline result
#' @noRd
.strip_wgcna_plots <- function(wgcna_out) {
    if (is.null(wgcna_out)) return(NULL)
    drop <- c(
        "soft_threshold_plot", "module_size_plot", "module_trait_heatmap",
        "module_eigengene_heatmap", "dendrogram_plot", "ora_plots_by_module"
    )
    for (nm in intersect(drop, names(wgcna_out))) {
        wgcna_out[[nm]] <- NULL
    }
    wgcna_out
}

#' Drop network ggplot objects from a PPI result
#' @noRd
.strip_ppi_plots <- function(ppi_out) {
    if (is.null(ppi_out)) return(NULL)
    if (length(ppi_out$by_module)) {
        ppi_out$by_module <- lapply(ppi_out$by_module, function(x) {
            x$plot <- NULL
            x
        })
    }
    if (length(ppi_out$by_comparison)) {
        ppi_out$by_comparison <- lapply(ppi_out$by_comparison, function(x) {
            x$plot <- NULL
            x
        })
    }
    ppi_out
}

#' Drop heatmap grobs from a GSVA result
#' @noRd
.strip_gsva_plots <- function(gsva_out) {
    if (is.null(gsva_out) || !length(gsva_out$heatmaps)) return(gsva_out)
    gsva_out$heatmaps <- lapply(gsva_out$heatmaps, function(h) {
        if (is.null(h)) return(h)
        h$gtable <- NULL
        h
    })
    gsva_out
}

#' Resolve a render-order key to its page label
#' @noRd
.report_page_label <- function(page) {
    pages <- .report_pages()
    by_key <- stats::setNames(names(pages), vapply(pages, `[[`, character(1), "key"))
    if (page %in% names(by_key)) {
        return(by_key[[page]])
    }
    if (page %in% names(pages)) {
        return(page)
    }
    page
}

#' Slim shared report params before rendering a page (avoids serializing unused plots)
#' @noRd
.report_params_for_page <- function(page, shared) {
    label <- .report_page_label(page)
    out <- shared
    if (!label %in% c("QC", "Index") && !is.null(out$distance_out)) {
        out$distance_out <- .strip_distance_interactive(out$distance_out)
    }
    if (!label %in% c("WGCNA", "Index", "Multihit") && !is.null(out$wgcna_out)) {
        out$wgcna_out <- .strip_wgcna_plots(out$wgcna_out)
    }
    if (!label %in% c("WGCNA", "Index", "DESeq2") &&
        !is.null(out$ppi_out) && length(out$ppi_out$by_module)) {
        out$ppi_out$by_module <- NULL
    }
    if (!label %in% c("DESeq2", "Index") &&
        !is.null(out$ppi_out) && length(out$ppi_out$by_comparison)) {
        out$ppi_out$by_comparison <- NULL
    }
    if (!label %in% c("GSVA", "Index") && !is.null(out$gsva_out)) {
        out$gsva_out <- .strip_gsva_plots(out$gsva_out)
    }
    out
}

#' Free heavy plot objects from the live shared-params bag after a page renders
#' @noRd
.release_heavy_report_cache <- function(shared, rendered_key) {
    label <- .report_page_label(rendered_key)
    if (identical(label, "WGCNA")) {
        shared$wgcna_out <- .strip_wgcna_plots(shared$wgcna_out)
        if (!is.null(shared$ppi_out) && length(shared$ppi_out$by_module)) {
            shared$ppi_out <- .strip_ppi_plots(shared$ppi_out)
        }
    } else if (identical(label, "GSVA")) {
        shared$gsva_out <- .strip_gsva_plots(shared$gsva_out)
    } else if (identical(label, "QC")) {
        shared$distance_out <- .strip_distance_interactive(shared$distance_out)
    }
    shared
}

#' Human-readable wall time since `started`, for render progress messages
#' @noRd
.format_elapsed <- function(started) {
    secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (secs < 90) {
        return(paste0(format(round(secs, 1), nsmall = 1), "s"))
    }
    paste0(format(round(secs / 60, 1), nsmall = 1), " min")
}

#' Rendered page size, for render progress messages
#'
#' Reported because page size is what drives render time: pandoc and
#' `rmarkdown`'s post-processing both scale with it, so an unexpectedly slow
#' page shows up here as an unexpectedly large one.
#' @noRd
.format_page_size <- function(path) {
    bytes <- tryCatch(file.size(path), error = function(e) NA_real_)
    if (!length(bytes) || is.na(bytes)) return("")
    if (bytes < 1024^2) {
        return(paste0(", ", round(bytes / 1024), " KB"))
    }
    paste0(", ", format(round(bytes / 1024^2, 1), nsmall = 1), " MB")
}

#' Render one hub HTML page from an in-memory shared-params bag
#' @noRd
.render_one_hub_page <- function(key, shared_bag, save_dir) {
    message("  Rendering HTML page: ", key, " ...")
    utils::flush.console()
    started <- Sys.time()
    page_params <- .report_params_for_page(key, shared_bag)
    path <- tryCatch(
        .render_report_page(key, page_params, save_dir),
        error = function(e) {
            stop("Failed to render report page \"", key, "\": ",
                 conditionMessage(e), call. = FALSE)
        }
    )
    message("  Finished HTML page: ", key, " (", .format_elapsed(started),
            .format_page_size(path), ")")
    utils::flush.console()
    list(
        path = path,
        shared = .release_heavy_report_cache(shared_bag, key)
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
    page_params <- utils::modifyList(
        shared_params,
        list(
            current_page = label,
            page_hrefs = report_page_hrefs(label),
            contrast = NULL
        )
    )
    if (identical(label, "Index")) page_params$repro <- .repro_for_report(save_dir)
    # Folder-bundled delivery: each HTML lives next to its own `lib/` so
    # rmarkdown 2.x `relativeTo(output_dir, lib_dir)` stays valid (a parent
    # shared lib, or a symlink into it, is normalized to a non-descendant
    # path and fails). Offline browsing still works via relative `lib/` hrefs.
    .render_page_rmd(spec$template, spec$file, out_dir, page_params)
}

#' Render the Index HTML page into `save_dir`
#' @inheritParams render_report
#' @param disabled_labels Page labels to grey out in the nav bar.
#' @return Path to `index.html`.
#' @export
render_index_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL, cspa_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
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
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL, cspa_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
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
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL, cspa_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
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
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL, cspa_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
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
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL, cspa_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    if (is.null(wgcna_out)) {
        stop("`wgcna_out` must be supplied to render the WGCNA page.", call. = FALSE)
    }
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("WGCNA", shared, save_dir)
}

#' Render the decoupleR-TF HTML page into `save_dir/decoupler-tf/`
#' @inheritParams render_index_page
#' @return Path to `decoupler-tf/decoupleR-TF.html`.
#' @export
render_decoupleR_tf_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL, cspa_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    if (is.null(decoupler_tf_out)) {
        stop("`decoupler_tf_out` must be supplied to render the decoupleR-TF page.",
             call. = FALSE)
    }
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("decoupleR-TF", shared, save_dir)
}

#' @rdname render_decoupleR_tf_page
#' @export
render_tf_activity_page <- render_decoupleR_tf_page

#' Render the GSVA HTML page into `save_dir/gsva/`
#' @inheritParams render_index_page
#' @param gsva_out Output of [run_gsva_pipeline()].
#' @return Path to `gsva/gsva.html`.
#' @export
render_gsva_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL,
    cspa_out = NULL, gsva_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    if (is.null(gsva_out)) {
        stop("`gsva_out` must be supplied to render the GSVA page.",
             call. = FALSE)
    }
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
        disabled_labels,
        gsva_out = gsva_out
    )
    .render_report_page("GSVA", shared, save_dir)
}

#' Render the decoupleR-CCC HTML hub into `save_dir/decoupler-ccc/`
#' @inheritParams render_index_page
#' @param decoupler_ccc_out Output of [run_decoupler_ccc_pipeline()].
#' @return Path to `decoupler-ccc/decoupleR-CCC.html`.
#' @export
render_decoupleR_ccc_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL,
    cspa_out = NULL, decoupler_ccc_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    if (is.null(decoupler_ccc_out)) {
        stop("`decoupler_ccc_out` must be supplied to render the decoupleR-CCC page.",
             call. = FALSE)
    }
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
        disabled_labels,
        decoupler_ccc_out = decoupler_ccc_out
    )
    .render_report_page("decoupleR-CCC", shared, save_dir)
}

#' Render the VIPER HTML page into `save_dir/viper/`
#' @inheritParams render_index_page
#' @return Path to `viper/viper.html`.
#' @export
render_viper_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL, cspa_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    if (is.null(viper_out)) {
        stop("`viper_out` must be supplied to render the VIPER page.",
             call. = FALSE)
    }
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
        disabled_labels
    )
    .render_report_page("VIPER", shared, save_dir)
}

#' Render the Multihit HTML page into `save_dir/multihit/`
#' @inheritParams render_index_page
#' @param multihit_out Output of [run_multihit_pipeline()] /
#'   [run_multihit_lrt_pipeline()].
#' @return Path to `multihit/multihit.html`.
#' @export
render_multihit_page <- function(
    dds, var, distance_out, diffexp_out, gsea_out,
    wgcna_out = NULL, decoupler_tf_out = NULL, ppi_out = NULL, viper_out = NULL, cspa_out = NULL,
    multihit_out = NULL,
    alpha = 0.05, lfc_threshold = 0.5, run_args = NULL,
    disabled_labels = character(0), save_dir = getwd()) {
    if (is.null(multihit_out)) {
        stop("`multihit_out` must be supplied to render the Multihit page.",
             call. = FALSE)
    }
    save_dir <- validate_save_dir(save_dir)
    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
        disabled_labels,
        multihit_out = multihit_out
    )
    .render_report_page("Multihit", shared, save_dir)
}

#' Render the pipeline's multi-page HTML summary report
#'
#' Knits the package's bundled R Markdown templates against the outputs of
#' [run_qc()], [run_distance_pca()], [run_diffexp_pipeline()],
#' [run_gsea_pipeline()], [run_wgcna_pipeline()],
#' [run_decoupler_tf_pipeline()], and [run_viper_pipeline()], writing each
#' section page into its module folder and `index.html` at the pipeline
#' root:
#' \itemize{
#'   \item `index.html`
#'   \item `qc/qc.html`
#'   \item `deseq2/deseq2.html`
#'   \item `clusterprofiler/clusterprofiler.html`
#'   \item `wgcna/wgcna.html` (omitted when `wgcna_out = NULL`)
#'   \item `decoupler-tf/decoupleR-TF.html` (omitted when `decoupler_tf_out = NULL`)
#'   \item `decoupler-ccc/decoupleR-CCC.html` plus per-contrast pages
#'     (omitted when `decoupler_ccc_out = NULL`)
#'   \item `viper/viper.html` (omitted when `viper_out = NULL`)
#'   \item `multihit/multihit.html` (omitted when `multihit_out = NULL`)
#' }
#' Pages use folder-bundled delivery (`self_contained = FALSE`) with a
#' shared `lib/` asset directory under `save_dir`; cross-links and CSV
#' downloads use relative paths so the whole results tree can be copied
#' and opened offline (do not share a lone HTML file without its siblings).
#'
#' @param dds The (QC-filtered) `DESeqDataSet`, used to report sample/gene
#'   counts and the design formula in the report header. If it carries
#'   `S4Vectors::metadata(dds)$degreports_qc` (attached automatically by
#'   [run_qc()]), a QC section is also included showing the number/percent
#'   of genes removed at each filtering step, the library-size plot, and
#'   the count-distribution diagnostic plot (with the expression cutoff
#'   marked) -- omitted if that metadata is absent, e.g. when `dds` was
#'   QC-filtered manually rather than via [run_qc()]. Pass `NULL` for
#'   DEG-table-only reports ([run_deg()]): the QC page is skipped
#'   and greyed out in the nav.
#' @param var Column name in `colData(dds)` used for grouping (or a label
#'   such as `"comparison"` when `dds = NULL`).
#' @param distance_out Output of [run_distance_pca()], or `NULL` when
#'   starting from a DEG table.
#' @param diffexp_out Output of [run_diffexp_pipeline()],
#'   [run_diffexp_lrt_pipeline()], or [prepare_diffexp_from_deg()].
#' @param gsea_out Output of [run_gsea_pipeline()] (may be `NULL`).
#' @param gsea_overlap_out Output of [run_gsea_overlap_pipeline()], or
#'   `NULL` (default). When present, leading-edge overlap plots are shown
#'   at the bottom of the clusterProfiler page.
#' @param multihit_out Output of [run_multihit_pipeline()], or `NULL`
#'   (default) to omit the Multihit page.
#' @param wgcna_out Output of [run_wgcna_pipeline()], or `NULL` (default)
#'   to omit the WGCNA page entirely (its link is greyed out, not removed,
#'   in the navigation bar of the other pages).
#' @param decoupler_tf_out Output of [run_decoupler_tf_pipeline()], or `NULL`
#'   (default) to omit the decoupleR-TF page.
#' @param decoupler_ccc_out Output of [run_decoupler_ccc_pipeline()], or
#'   `NULL` (default) to omit the decoupleR-CCC page.
#' @param ppi_out Output of [run_ppi_pipeline()], or `NULL` (default) to
#'   omit PPI connectivity sections on the DE and WGCNA pages.
#' @param viper_out Output of [run_viper_pipeline()], or `NULL` (default)
#'   to omit the VIPER page.
#' @param cspa_out Output of [run_cspa_marker_pipeline()], or `NULL`.
#' @param degpatterns_out Optional degPatterns result list (LRT).
#' @param cluster_ora_out Optional cluster ORA result list (LRT).
#' @param alpha Significance threshold to display/use for report summaries.
#'   Default `0.05`.
#' @param lfc_threshold Log2 fold-change threshold used for up/down
#'   colouring on MA/volcano plots. Default `0.5`.
#' @param run_args Optional named list of pipeline arguments (as built by
#'   [run_deseq2()] / [run_deg()]) to display on the index page.
#'   Default `NULL`.
#' @param save_dir Pipeline root directory. Section pages are written into
#'   per-module subfolders under this path. Default `getwd()`.
#' @param parallel_render If `TRUE`, knit per-contrast HTML pages
#'   (`deseq2/<contrast>/`, `clusterprofiler/<contrast>/`,
#'   `decoupler-ccc/<contrast>/`) in parallel on Unix via
#'   `parallel::mclapply()`. Hub pages (`qc.html`, `wgcna.html`, …) always
#'   render in the current R session. Default `FALSE`.
#' @section When a contrast page fails:
#' A page that errors in a worker is retried once in the current session,
#' where pandoc's own output is visible and nothing else is competing for the
#' filesystem. If the retry fails as well, the page is reported as a warning
#' and the rest of the report (including the index) still renders, rather than
#' one page taking down a run whose analyses are already finished.
#' @param render_workers Number of parallel workers for `parallel_render`.
#'   Default `NULL` uses `min(4, parallel::detectCores() - 1)`; the cap keeps
#'   the forked knits from exhausting RAM and stalling later hub pages.
#' @section When pandoc is missing:
#' Every page is knitted through pandoc, so a session without one (an HPC
#' batch job or `Rscript` invocation, which do not see RStudio's bundled
#' copy) cannot render any of them. That is reported once, as a warning
#' naming the ways to supply a pandoc, and the function returns `NULL`
#' rather than erroring -- a pipeline that has already written its tables
#' and figures then still completes.
#' @return Invisibly, a named character vector of rendered file paths, with
#'   names `"index"`, `"qc"` (when `dds` is supplied), `"deseq2"`,
#'   `"clusterProfiler"`, and (when rendered) `"wgcna"` / `"decoupleR-TF"` /
#'   `"decoupleR-CCC"` / `"viper"` / `"multihit"`; or `NULL` when pandoc is
#'   unavailable.
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
    gsea_overlap_out = NULL,
    multihit_out = NULL,
    wgcna_out = NULL,
    decoupler_tf_out = NULL,
    ppi_out = NULL,
    viper_out = NULL,
    cspa_out = NULL,
    gsva_out = NULL,
    decoupler_ccc_out = NULL,
    degpatterns_out = NULL,
    cluster_ora_out = NULL,
    alpha = 0.05,
    lfc_threshold = 0.5,
    run_args = NULL,
    save_dir = getwd(),
    parallel_render = FALSE,
    render_workers = NULL,
    contrasts_only = FALSE) {

    # Every bundled Rmd template embeds interactive plotly/DT widgets via
    # htmltools::tagList() at render time, so htmltools must be installed
    # even though this R function itself never calls it directly.
    if (!requireNamespace("htmltools", quietly = TRUE)) {
        stop("Package 'htmltools' is required to render the HTML report.", call. = FALSE)
    }

    # Every page goes through pandoc, so a missing one fails all of them
    # identically. Reported once here, rather than as a page-level error
    # from whichever page happened to be knitted first.
    if (!.pandoc_ok()) {
        warning(.pandoc_missing_msg(), " Skipping the report; the analysis outputs ",
                "already written to disk are unaffected.", call. = FALSE)
        return(invisible(NULL))
    }

    save_dir <- validate_save_dir(save_dir)

    disabled_labels <- character(0)
    if (is.null(dds)) disabled_labels <- c(disabled_labels, "QC")
    if (is.null(gsva_out)) disabled_labels <- c(disabled_labels, "GSVA")
    if (is.null(wgcna_out)) disabled_labels <- c(disabled_labels, "WGCNA")
    if (is.null(decoupler_tf_out)) disabled_labels <- c(disabled_labels, "decoupleR-TF")
    if (is.null(decoupler_ccc_out)) disabled_labels <- c(disabled_labels, "decoupleR-CCC")
    if (is.null(viper_out)) disabled_labels <- c(disabled_labels, "VIPER")
    if (is.null(multihit_out)) disabled_labels <- c(disabled_labels, "Multihit")

    shared <- .report_shared_params(
        dds, var, distance_out, diffexp_out, gsea_out,
        wgcna_out, decoupler_tf_out, ppi_out, viper_out, cspa_out, alpha, lfc_threshold, run_args,
        disabled_labels,
        gsea_overlap_out = gsea_overlap_out,
        multihit_out = multihit_out,
        degpatterns_out = degpatterns_out,
        cluster_ora_out = cluster_ora_out,
        gsva_out = gsva_out,
        decoupler_ccc_out = decoupler_ccc_out
    )

    # Module hubs first, then per-contrast pages for heavy modules, then
    # remaining modules, then index last. Skip QC when there is no dds
    # (DEG-table pipeline).
    render_order <- character(0)
    if (!is.null(dds)) render_order <- c(render_order, "qc")
    render_order <- c(render_order, "deseq2", "clusterProfiler")
    if (!is.null(gsva_out)) render_order <- c(render_order, "gsva")
    if (!is.null(wgcna_out)) render_order <- c(render_order, "wgcna")
    if (!is.null(decoupler_tf_out)) render_order <- c(render_order, "decoupleR-TF")
    if (!is.null(decoupler_ccc_out)) render_order <- c(render_order, "decoupleR-CCC")
    if (!is.null(viper_out)) render_order <- c(render_order, "viper")
    if (!is.null(multihit_out)) render_order <- c(render_order, "multihit")
    render_order <- c(render_order, "index")

    out_paths <- stats::setNames(rep(NA_character_, length(render_order)), render_order)

    render_contrast_pages <- function(key, shared_bag, save_dir) {
        if (!key %in% c("deseq2", "clusterProfiler", "decoupleR-CCC")) {
            return(invisible(NULL))
        }
        comps <- if (identical(key, "deseq2")) {
            names(diffexp_out$by_comparison)
        } else if (identical(key, "decoupleR-CCC")) {
            names(decoupler_ccc_out$contrast_acts)
        } else if (!is.null(cluster_ora_out) &&
                   !is.null(cluster_ora_out$ora_by_lrt) &&
                   length(cluster_ora_out$ora_by_lrt)) {
            names(cluster_ora_out$ora_by_lrt)
        } else if (!is.null(gsea_out) && !is.null(gsea_out$by_comparison)) {
            names(gsea_out$by_comparison)
        } else {
            character(0)
        }
        if (!length(comps)) return(invisible(NULL))
        message("  Rendering ", length(comps), " per-contrast HTML page(s) for ",
                key, " ...")
        utils::flush.console()

        # Slim the params bag once, in the parent. Doing it per child made
        # every fork copy-on-write its way through the full shared bag.
        page_params <- .report_params_for_page(key, shared_bag)

        render_one <- function(comp) {
            started <- Sys.time()
            message("  Rendering HTML page: ", key, " / ", comp, " ...")
            utils::flush.console()
            path <- tryCatch(
                .render_contrast_report_page(key, comp, page_params, save_dir),
                error = function(e) {
                    stop("Failed to render ", key, " page \"", comp, "\": ",
                         conditionMessage(e), call. = FALSE)
                }
            )
            message("  Finished HTML page: ", key, " / ", comp,
                    " (", .format_elapsed(started), .format_page_size(path), ")")
            utils::flush.console()
            path
        }

        # Forks inherit this session's heap, and each one knits a multi-MB
        # self-contained page. Too many at once pushes the machine into swap
        # and the *next* hub render (WGCNA) looks hung, so cap the fan-out.
        n_workers <- if (is.null(render_workers)) {
            min(4L, max(1L, parallel::detectCores(logical = FALSE) - 1L, na.rm = TRUE))
        } else {
            as.integer(render_workers)[1]
        }
        if (is.na(n_workers) || n_workers < 1L) n_workers <- 1L
        used_parallel <- isTRUE(parallel_render) &&
            length(comps) > 1L &&
            n_workers > 1L &&
            identical(.Platform$OS.type, "unix")
        try_render <- function(comp) tryCatch(render_one(comp), error = identity)
        if (used_parallel) {
            invisible(gc(full = TRUE))
            res <- parallel::mclapply(
                comps, function(comp) {
                    .isolate_rmarkdown_tmpfiles()
                    try_render(comp)
                },
                mc.cores = min(n_workers, length(comps)),
                mc.preschedule = FALSE
            )
            failed <- !vapply(res, .rendered_ok, logical(1))
            if (any(failed)) {
                # A page can die in a worker for reasons that have nothing to
                # do with its own content -- several pandoc processes chewing
                # through multi-MB documents at once, or a fork that lost its
                # memory. Retry those alone in this session, where pandoc's own
                # stderr reaches the console instead of a dead child's.
                message("  ", sum(failed), " of ", length(comps), " ", key,
                        " page(s) failed under parallel render (",
                        .render_failure_reason(res[[which(failed)[1]]]),
                        "); retrying sequentially: ",
                        paste(comps[failed], collapse = ", "), " ...")
                utils::flush.console()
                for (i in which(failed)) {
                    invisible(gc(full = TRUE))
                    res[[i]] <- try_render(comps[i])
                }
            }
        } else {
            res <- lapply(comps, try_render)
        }
        failed <- !vapply(res, .rendered_ok, logical(1))
        if (any(failed)) {
            # Warn rather than abort: the remaining pages, the index, and
            # whatever the caller does after rendering are all worth more than
            # these, and the hub still lists every contrast that was run.
            warning(sum(failed), " of ", length(comps), " ", key,
                    " contrast page(s) could not be rendered and are missing ",
                    "from the report: ", paste(comps[failed], collapse = ", "),
                    ". First error: ",
                    .render_failure_reason(res[[which(failed)[1]]]),
                    call. = FALSE, immediate. = TRUE)
        }
        rm(page_params)
        invisible(gc(full = TRUE))
    }

    if (isTRUE(contrasts_only)) {
        for (key in intersect(render_order, c("deseq2", "clusterProfiler", "decoupleR-CCC"))) {
            render_contrast_pages(key, shared, save_dir)
        }
        message("HTML contrast pages updated under: ", save_dir)
        return(invisible(out_paths))
    }

    # Hubs always knit in this R session. Shipping each hub to a future
    # worker meant contrast pages waited until *all* hubs finished, and
    # WGCNA/GSVA workers could sit silent for many minutes (load_all +
    # huge RDS). `parallel_render` only parallelizes per-contrast pages.
    for (key in render_order) {
        if (key %in% c("gsva", "wgcna")) {
            message("  (the ", key, " page can take several minutes) ...")
        }
        res <- .render_one_hub_page(key, shared, save_dir)
        out_paths[[key]] <- res$path
        shared <- res$shared
        invisible(gc(full = TRUE))
        render_contrast_pages(key, shared, save_dir)
    }

    message("HTML report written to: ", out_paths[["index"]],
            " (+ ", max(0L, length(out_paths) - 1L),
            " linked section page(s) in per-module folders)")
    invisible(out_paths)
}
