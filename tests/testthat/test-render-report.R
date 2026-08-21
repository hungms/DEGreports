# Tests for report_nav_links() (the shared cross-page navigation bar) and
# render_report() (the multi-page HTML renderer: index + per-module pages).

# contrast_file_stem() takes one label at a time.
contrast_page_paths <- function(save_dir, comps, module = "deseq2") {
    vapply(comps, function(comp) {
        file.path(save_dir, module, comp,
                  paste0(contrast_file_stem(comp), ".html"))
    }, character(1), USE.NAMES = FALSE)
}

test_that("report_page_hrefs() returns nested relative paths from Index and sibling folders", {
    from_index <- report_page_hrefs("Index")
    expect_equal(from_index[["Index"]], "index.html")
    expect_equal(from_index[["QC"]], "qc/qc.html")
    expect_equal(from_index[["DESeq2"]], "deseq2/deseq2.html")
    expect_equal(from_index[["clusterProfiler"]], "clusterprofiler/clusterprofiler.html")
    expect_false("GSEA overlap" %in% names(from_index))
    expect_equal(from_index[["GSVA"]], "gsva/gsva.html")
    expect_equal(from_index[["WGCNA"]], "wgcna/wgcna.html")
    expect_equal(from_index[["decoupleR-TF"]], "decoupler-tf/decoupleR-TF.html")
    expect_equal(from_index[["decoupleR-CCC"]], "decoupler-ccc/decoupleR-CCC.html")
    expect_equal(from_index[["VIPER"]], "viper/viper.html")
    expect_equal(from_index[["Multihit"]], "multihit/multihit.html")

    from_qc <- report_page_hrefs("QC")
    expect_equal(from_qc[["Index"]], "../index.html")
    expect_equal(from_qc[["QC"]], "qc.html")
    expect_equal(from_qc[["WGCNA"]], "../wgcna/wgcna.html")
    expect_equal(from_qc[["VIPER"]], "../viper/viper.html")
    expect_equal(from_qc[["Multihit"]], "../multihit/multihit.html")
})

test_that("report_nav_links() marks the current page as plain bold text, not a link", {
    html <- report_nav_links("QC")
    expect_true(grepl("<strong>QC</strong>", html))
    expect_false(grepl('<a href="qc.html">QC</a>', html))
    # other pages use paths relative to the QC folder
    expect_true(grepl('<a href="../index.html">Index</a>', html))
    expect_true(grepl('<a href="../deseq2/deseq2.html">DESeq2</a>', html))
    expect_true(grepl('<a href="../clusterprofiler/clusterprofiler.html">clusterProfiler</a>', html))
    expect_true(grepl('<a href="../wgcna/wgcna.html">WGCNA</a>', html))
    expect_true(grepl('<a href="../multihit/multihit.html">Multihit</a>', html))
})

test_that("report_nav_links() greys out disabled pages instead of linking them", {
    html <- report_nav_links("Index", disabled = "WGCNA")
    expect_true(grepl("<strong>Index</strong>", html))
    expect_true(grepl('color:#999;">WGCNA</span>', html))
    expect_false(grepl('<a href="wgcna/wgcna.html">', html))
})

test_that("report_nav_links() rejects an unrecognised current-page label", {
    expect_error(report_nav_links("Not A Page"), "must be one of")
})

test_that("render_report() writes cross-linked, self-contained HTML files in per-module folders when wgcna_out is supplied", {
    skip_if_not_installed("rmarkdown")
    skip_if_not_installed("htmltools")
    skip_if_not_installed("plotly")
    skip_if_not_installed("DT")
    skip_if_not_installed("msigdbr")
    skip_if_not_installed("WGCNA")
    skip_if(!rmarkdown::pandoc_available(), "pandoc not available")

    dds <- simulate_dds_real_genes(n_genes = 120, n_per_group = 4, seed = 11)
    dds <- run_qc(dds, var = "condition", min_count = 1)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)

    distance_out <- run_distance_pca(dds, var = "condition")
    diffexp_out <- run_diffexp_pipeline(dds, var = "condition")
    gsea_out <- run_gsea_pipeline(diffexp_out$results, collections = c("HALLMARK"), run_clusterprofiler = FALSE)
    wgcna_out <- run_wgcna_pipeline(dds, var = "condition", save_dir = tempdir(),
                                     min_module_size = 10, ora_collection = c("GOBP"),
                                     save_plot = FALSE)

    out_dir <- withr::local_tempdir()
    paths <- render_report(
        dds, var = "condition",
        distance_out = distance_out, diffexp_out = diffexp_out,
        gsea_out = gsea_out, wgcna_out = wgcna_out,
        save_dir = out_dir
    )

    expect_setequal(names(paths), c("index", "qc", "deseq2", "clusterProfiler", "wgcna"))
    expect_true(all(file.exists(paths)))
    expect_true(all(file.info(paths)$size > 1000))
    expect_equal(basename(dirname(paths[["qc"]])), "qc")
    expect_equal(basename(dirname(paths[["deseq2"]])), "deseq2")
    expect_equal(basename(dirname(paths[["clusterProfiler"]])), "clusterprofiler")
    expect_equal(basename(dirname(paths[["wgcna"]])), "wgcna")

    index_html <- paste(readLines(paths[["index"]]), collapse = "\n")
    expect_true(grepl('href="qc/qc.html"', index_html))
    expect_true(grepl('href="wgcna/wgcna.html"', index_html))
    expect_false(grepl('color:#999;">WGCNA', index_html))  # not disabled when wgcna_out is present
    expect_true(grepl('color:#999;">Multihit', index_html))  # greyed when multihit_out is NULL

    qc_html <- paste(readLines(paths[["qc"]]), collapse = "\n")
    expect_true(grepl('href="../index.html"', qc_html))
    expect_true(grepl("<strong>QC</strong>", qc_html))

    # Regression: emit_widget()'s indented tag serialization used to make
    # Pandoc treat deeply-nested corner-button <div>s as an indented code
    # block, escaping the figure/button markup into literal <pre><code>
    # text instead of rendering it live.
    expect_false(grepl("<pre><code>", qc_html, fixed = TRUE))
    expect_true(grepl("justify-content: flex-end", qc_html, fixed = TRUE))
    diffexp_html <- paste(readLines(paths[["deseq2"]]), collapse = "\n")
    expect_false(grepl("<pre><code>", diffexp_html, fixed = TRUE))
    wgcna_html <- paste(readLines(paths[["wgcna"]]), collapse = "\n")
    expect_false(grepl("<pre><code>", wgcna_html, fixed = TRUE))

    # Regression: even with folder-bundled (`self_contained = FALSE`) HTML,
    # pandoc can still inject a remote MathJax
    # <script src="https://mathjax.rstudio.com/..."> into every page; a
    # network-restricted viewer can then block navigation/rendering. None
    # of these reports use LaTeX math, so `mathjax: null` in each
    # template's YAML header must keep this out.
    for (html_txt in list(index_html, qc_html, diffexp_html, wgcna_html)) {
        expect_false(grepl("mathjax.rstudio.com", html_txt, fixed = TRUE))
    }

    # Per-contrast GSEA pages tab collections (HALLMARK / GOBP / …) via
    # Bootstrap tabs, not stacked "Collection: …" headings.
    gsea_comp <- names(gsea_out$by_comparison)[[1]]
    gsea_page <- file.path(out_dir, "clusterprofiler", gsea_comp,
                           paste0(contrast_file_stem(gsea_comp), ".html"))
    expect_true(file.exists(gsea_page))
    gsea_html <- paste(readLines(gsea_page, warn = FALSE), collapse = "\n")
    expect_true(grepl("deseq2reports-bs-tabset", gsea_html, fixed = TRUE))
    expect_true(grepl(">HALLMARK<", gsea_html, fixed = TRUE))
    expect_false(grepl("Collection: HALLMARK", gsea_html, fixed = TRUE))
})

test_that("render_report()'s wgcna page has a top-level 'Module genes' section tabbed by module, and per-module pathway enrichment tabbed module-first then collection", {
    skip_if_not_installed("rmarkdown")
    skip_if_not_installed("htmltools")
    skip_if_not_installed("plotly")
    skip_if_not_installed("DT")
    skip_if_not_installed("msigdbr")
    skip_if_not_installed("WGCNA")
    skip_if(!rmarkdown::pandoc_available(), "pandoc not available")

    dds <- simulate_dds_real_genes(n_genes = 120, n_per_group = 4, seed = 11)
    dds <- run_qc(dds, var = "condition", min_count = 1)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)

    distance_out <- run_distance_pca(dds, var = "condition")
    diffexp_out <- run_diffexp_pipeline(dds, var = "condition")
    gsea_out <- run_gsea_pipeline(diffexp_out$results, collections = c("HALLMARK"), run_clusterprofiler = FALSE)
    wgcna_out <- run_wgcna_pipeline(dds, var = "condition", save_dir = tempdir(),
                                     min_module_size = 10, ora_collection = c("GOBP", "HALLMARK"),
                                     save_plot = FALSE)

    out_dir <- withr::local_tempdir()
    paths <- render_report(
        dds, var = "condition",
        distance_out = distance_out, diffexp_out = diffexp_out,
        gsea_out = gsea_out, wgcna_out = wgcna_out,
        save_dir = out_dir
    )

    # Pandoc may soft-wrap tab labels ("Module 1" -> "Module\n1"); collapse
    # whitespace before substring searches so assertions stay stable.
    wgcna_html <- paste(readLines(paths[["wgcna"]], warn = FALSE), collapse = "\n")
    wgcna_flat <- gsub("[[:space:]]+", " ", wgcna_html)

    modules <- unique(wgcna_out$module_df$module)
    collections <- names(wgcna_out$ora_by_module)
    # Module genes tabs omit WGCNA's grey "Unassigned" catch-all.
    gene_modules <- setdiff(modules, "Unassigned")

    expect_true(grepl(">Module genes<", wgcna_flat, fixed = TRUE))
    html_len <- nchar(wgcna_flat)
    module_genes_pos <- regexpr(">Module genes<", wgcna_flat, fixed = TRUE)
    for (mod in gene_modules) {
        mod_pos <- regexpr(paste0(">", mod, "<"),
                           substring(wgcna_flat, module_genes_pos, html_len),
                           fixed = TRUE)
        expect_true(mod_pos > 0, info = paste("module heading not found after Module genes:", mod))
    }

    # per-module pathway enrichment: module tabs first, collections nested
    pme_pos <- regexpr(">Per-module pathway enrichment<", wgcna_flat, fixed = TRUE)
    expect_true(pme_pos > 0)
    pme_section <- substring(wgcna_flat, pme_pos, html_len)
    pme_len <- nchar(pme_section)
    for (mod in gene_modules) {
        mod_tab_pos <- regexpr(paste0(">", mod, "<"), pme_section, fixed = TRUE)
        expect_true(mod_tab_pos > 0, info = paste("module tab not found in pathway enrichment section:", mod))
        after_mod <- substring(pme_section, mod_tab_pos, pme_len)
        for (coll in collections) {
            # inner tabs use the bare collection name (not "Collection: ...")
            coll_pos <- regexpr(paste0(">", coll, "<"), after_mod, fixed = TRUE)
            expect_true(coll_pos > 0, info = paste("collection tab not nested under module:", mod, coll))
        }
    }
})

test_that("a contrast page that only fails inside a worker is retried in-session", {
    skip_if_not_installed("rmarkdown")
    skip_if_not_installed("DT")
    skip_if_not_installed("plotly")
    skip_if(!rmarkdown::pandoc_available(), "pandoc not available")
    skip_on_os("windows")

    dds <- simulate_dds(n_genes = 60, group_levels = c("groupA", "groupB", "groupC"),
                        n_per_group = 3, seed = 13)
    dds <- run_qc(dds, var = "condition", min_count = 1)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    out_dir <- withr::local_tempdir()
    distance_out <- run_distance_pca(dds, var = "condition")
    diffexp_out <- run_diffexp_pipeline(dds, var = "condition",
                                        save_dir = file.path(out_dir, "de"))
    comps <- names(diffexp_out$by_comparison)

    # Stand in for whatever kills pandoc when several forks run it at once:
    # fail in any child process, succeed when the parent retries alone.
    parent_pid <- Sys.getpid()
    real_page <- deseq2reports:::.render_contrast_report_page
    local_mocked_bindings(
        .render_contrast_report_page = function(...) {
            if (!identical(Sys.getpid(), parent_pid)) {
                stop("pandoc document conversion failed with error 1", call. = FALSE)
            }
            real_page(...)
        }
    )

    expect_no_warning(render_report(
        dds, var = "condition", distance_out = distance_out,
        diffexp_out = diffexp_out, gsea_out = NULL, save_dir = out_dir,
        parallel_render = TRUE, render_workers = 2
    ))
    expect_true(all(file.exists(contrast_page_paths(out_dir, comps))))
})

test_that("a contrast page that fails everywhere warns and leaves the rest of the report standing", {
    skip_if_not_installed("rmarkdown")
    skip_if_not_installed("DT")
    skip_if_not_installed("plotly")
    skip_if(!rmarkdown::pandoc_available(), "pandoc not available")

    dds <- simulate_dds(n_genes = 60, group_levels = c("groupA", "groupB", "groupC"),
                        n_per_group = 3, seed = 13)
    dds <- run_qc(dds, var = "condition", min_count = 1)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    out_dir <- withr::local_tempdir()
    distance_out <- run_distance_pca(dds, var = "condition")
    diffexp_out <- run_diffexp_pipeline(dds, var = "condition",
                                        save_dir = file.path(out_dir, "de"))
    comps <- names(diffexp_out$by_comparison)
    doomed <- comps[[1]]

    real_page <- deseq2reports:::.render_contrast_report_page
    local_mocked_bindings(
        .render_contrast_report_page = function(module_key, contrast, ...) {
            if (identical(contrast, doomed)) {
                stop("pandoc document conversion failed with error 1", call. = FALSE)
            }
            real_page(module_key, contrast, ...)
        }
    )

    paths <- NULL
    expect_warning(
        paths <- render_report(dds, var = "condition", distance_out = distance_out,
                               diffexp_out = diffexp_out, gsea_out = NULL,
                               save_dir = out_dir, parallel_render = FALSE),
        "could not be rendered"
    )
    # The run continues: index, hubs, and every other contrast page survive.
    expect_true(file.exists(paths[["index"]]))
    expect_true(all(file.exists(
        contrast_page_paths(out_dir, setdiff(comps, doomed))
    )))
    expect_false(file.exists(contrast_page_paths(out_dir, doomed)))
})

test_that("render_report(parallel_render = TRUE) writes one distinct page per contrast and no intermediates next to the template", {
    skip_if_not_installed("rmarkdown")
    skip_if_not_installed("DT")
    skip_if_not_installed("plotly")
    skip_if(!rmarkdown::pandoc_available(), "pandoc not available")
    skip_on_os("windows")

    dds <- simulate_dds(n_genes = 100,
                        group_levels = c("groupA", "groupB", "groupC"),
                        n_per_group = 3, seed = 13)
    dds <- run_qc(dds, var = "condition", min_count = 1)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    distance_out <- run_distance_pca(dds, var = "condition")
    diffexp_out <- run_diffexp_pipeline(dds, var = "condition")
    comps <- names(diffexp_out$by_comparison)
    expect_gt(length(comps), 1L)

    out_dir <- withr::local_tempdir()
    render_report(
        dds, var = "condition",
        distance_out = distance_out,
        diffexp_out = diffexp_out,
        gsea_out = NULL,
        save_dir = out_dir,
        parallel_render = TRUE,
        render_workers = 2
    )

    pages <- vapply(comps, function(comp) {
        file.path(out_dir, "deseq2", comp,
                  paste0(contrast_file_stem(comp), ".html"))
    }, character(1))
    expect_true(all(file.exists(pages)))

    # Regression: knitr's intermediates defaulted to the template's own
    # directory, so parallel workers clobbered a single shared
    # report_deseq2.knit.md and every contrast page came out byte-identical.
    digests <- vapply(pages, function(p) {
        paste(tools::md5sum(p), collapse = "")
    }, character(1))
    expect_equal(length(unique(digests)), length(comps))

    template_dir <- dirname(deseq2reports:::.report_template_path("report_deseq2.Rmd"))
    expect_length(list.files(template_dir, pattern = "\\.(knit|utf8)\\.md$"), 0L)
})

test_that("render_report() groups a three-factor grid into main effects and interactions", {
    skip_if_not_installed("rmarkdown")
    skip_if_not_installed("DT")
    skip_if_not_installed("plotly")
    skip_if(!rmarkdown::pandoc_available(), "pandoc not available")

    dds <- simulate_three_factor_dds(n_genes = 60, n_per_cell = 3)
    var <- c("celltype", "C_mutation", "K_mutation")
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    out_dir <- withr::local_tempdir()
    distance_out <- run_distance_pca(dds, var = var[1])
    diffexp_out <- run_diffexp_pipeline(
        dds, var = var, include_interactions = TRUE,
        save_dir = file.path(out_dir, "de"), save_plot = FALSE
    )

    paths <- render_report(
        dds, var = var, distance_out = distance_out,
        diffexp_out = diffexp_out, gsea_out = NULL, save_dir = out_dir
    )

    hub <- gsub("[[:space:]]+", " ",
                paste(readLines(paths[["deseq2"]], warn = FALSE), collapse = "\n"))
    expect_true(grepl("Main effects (across all other factors)", hub, fixed = TRUE))
    expect_true(grepl("Interactions (difference of differences)", hub, fixed = TRUE))
    expect_false(grepl("by stratum", hub, fixed = TRUE))
    # Each interaction is also spelled out as the question it answers.
    expect_true(grepl("What each interaction asks", hub, fixed = TRUE))
    expect_true(grepl(
        "Does the K_mutation effect (<code>K_vs_noK</code>) differ between C_mutation C and noC, within celltype = CB?",
        hub, fixed = TRUE
    ))

    specs <- diffexp_out$contrast_specs
    expect_true(all(file.exists(contrast_page_paths(out_dir, specs$label))))
    ix_page <- contrast_page_paths(out_dir, "C_vs_noC_x_K_vs_noK_in_CB")
    ix_html <- paste(readLines(ix_page, warn = FALSE), collapse = "\n")
    # The label is not a pair of sample groups, so the MA / volcano
    # annotations must not claim genes are "up in C".
    expect_false(grepl("up in C<", ix_html, fixed = TRUE))
})

test_that("render_report() omits the wgcna page and greys out its nav link when wgcna_out = NULL", {
    skip_if_not_installed("rmarkdown")
    skip_if_not_installed("htmltools")
    skip_if_not_installed("plotly")
    skip_if_not_installed("DT")
    skip_if_not_installed("msigdbr")
    skip_if(!rmarkdown::pandoc_available(), "pandoc not available")

    dds <- simulate_dds_real_genes(n_genes = 100, n_per_group = 4, seed = 12)
    dds <- run_qc(dds, var = "condition", min_count = 1)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)

    distance_out <- run_distance_pca(dds, var = "condition")
    diffexp_out <- run_diffexp_pipeline(dds, var = "condition")
    gsea_out <- run_gsea_pipeline(diffexp_out$results, collections = c("HALLMARK"), run_clusterprofiler = FALSE)

    out_dir <- withr::local_tempdir()
    paths <- render_report(
        dds, var = "condition",
        distance_out = distance_out, diffexp_out = diffexp_out,
        gsea_out = gsea_out, wgcna_out = NULL,
        save_dir = out_dir
    )

    expect_setequal(names(paths), c("index", "qc", "deseq2", "clusterProfiler"))
    expect_false(file.exists(file.path(out_dir, "wgcna.html")))
    expect_false(file.exists(file.path(out_dir, "wgcna", "wgcna.html")))

    index_html <- paste(readLines(paths[["index"]]), collapse = "\n")
    expect_true(grepl('color:#999;">WGCNA', index_html))
    expect_false(grepl('href="wgcna/wgcna.html">WGCNA', index_html))
    expect_true(grepl('color:#999;">Multihit', index_html))
    expect_true(grepl('color:#999;">VIPER', index_html))
})
