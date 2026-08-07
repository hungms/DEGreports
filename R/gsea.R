#' Map a package organism string to the species name `msigdbr` expects
#' @param org `"human"` or `"mouse"`.
#' @return `"Homo sapiens"` or `"Mus musculus"`.
#' @noRd
org_to_species <- function(org) {
    switch(org, human = "Homo sapiens", mouse = "Mus musculus",
           stop("Unsupported organism: ", org, call. = FALSE))
}

#' Friendly MSigDB collection aliases supported by this package
#'
#' Maps a short, memorable alias to the `msigdbr` collection/subcollection
#' pair it corresponds to, and to the gene-set-name prefix `msigdbr` uses
#' for that collection (needed to cosmetically strip the prefix in
#' [plot_gsea_barplot()]). Keeping this as a single named table means the
#' alias -> `msigdbr` mapping and the alias -> display-prefix mapping can
#' never drift apart.
#' @noRd
.msigdb_collection_map <- list(
    HALLMARK     = list(cat = "H",  subcat = NULL,             prefix = "HALLMARK_"),
    GOBP         = list(cat = "C5", subcat = "GO:BP",          prefix = "GOBP_"),
    GOCC         = list(cat = "C5", subcat = "GO:CC",          prefix = "GOCC_"),
    GOMF         = list(cat = "C5", subcat = "GO:MF",          prefix = "GOMF_"),
    KEGG_MEDICUS = list(cat = "C2", subcat = "CP:KEGG_MEDICUS", prefix = "KEGG_MEDICUS_"),
    # NOTE: msigdbr's `CP:KEGG_LEGACY` subcollection still names its gene
    # sets with the older bare "KEGG_" prefix (e.g. "KEGG_ABC_TRANSPORTERS"),
    # not "KEGG_LEGACY_" -- verified against the installed msigdbr db.
    KEGG_LEGACY  = list(cat = "C2", subcat = "CP:KEGG_LEGACY",  prefix = "KEGG_"),
    REACTOME     = list(cat = "C2", subcat = "CP:REACTOME",    prefix = "REACTOME_"),
    PID          = list(cat = "C2", subcat = "CP:PID",         prefix = "PID_")
)

#' Resolve a friendly collection alias to its `msigdbr` collection/subcollection
#' @param collection One of `names(.msigdb_collection_map)`, e.g. `"HALLMARK"`.
#' @return A list with elements `cat` and `subcat` (`subcat` may be `NULL`).
#' @noRd
resolve_msigdb_alias <- function(collection) {
    spec <- .msigdb_collection_map[[collection]]
    if (is.null(spec)) {
        stop("Unsupported MSigDB collection alias: \"", collection, "\". Supported aliases are: ",
             paste(names(.msigdb_collection_map), collapse = ", "), ".", call. = FALSE)
    }
    spec
}

#' Fetch an MSigDB gene set collection as a raw `msigdbr` data frame
#'
#' Wraps `msigdbr::msigdbr()`, handling both the legacy (`category`/
#' `subcategory`) and current (`collection`/`subcollection`) argument names
#' across `msigdbr` package versions. Collection identifiers are given as
#' friendly aliases (`"HALLMARK"`, `"GOBP"`, `"GOCC"`, `"GOMF"`,
#' `"KEGG_MEDICUS"`, `"REACTOME"`), resolved via `resolve_msigdb_alias()`. Shared by
#' [import_msigdbr()] (fgsea-style named-list format) and
#' [import_msigdbr_term2gene()] (clusterProfiler's `TERM2GENE` data-frame
#' format) so the collection is only fetched once per representation.
#'
#' @param org `"human"` or `"mouse"`.
#' @param collection A single collection alias, e.g. `"HALLMARK"` or
#'   `"GOBP"`. See `resolve_msigdb_alias()` for the supported aliases.
#' @return The raw `msigdbr` data frame, with a `gene_symbol` column
#'   normalized (renamed from `db_gene_symbol` on `msigdbr` versions that
#'   use that name).
#' @noRd
fetch_msigdbr_df <- function(org, collection) {
    org <- validate_org(org)
    species <- org_to_species(org)
    spec <- resolve_msigdb_alias(collection)
    cat <- spec$cat
    subcat <- spec$subcat

    msigdbr_args <- formals(msigdbr::msigdbr)
    df <- if ("collection" %in% names(msigdbr_args)) {
        if (is.null(subcat)) {
            msigdbr::msigdbr(species = species, collection = cat)
        } else {
            msigdbr::msigdbr(species = species, collection = cat, subcollection = subcat)
        }
    } else {
        if (is.null(subcat)) {
            msigdbr::msigdbr(species = species, category = cat)
        } else {
            msigdbr::msigdbr(species = species, category = cat, subcategory = subcat)
        }
    }

    if (!"gene_symbol" %in% colnames(df) && "db_gene_symbol" %in% colnames(df)) {
        df$gene_symbol <- df$db_gene_symbol
    }
    df
}

#' Fetch an MSigDB gene set collection as a named list of gene vectors
#'
#' @param org `"human"` or `"mouse"`.
#' @param collection A single collection alias, e.g. `"HALLMARK"` or
#'   `"GOBP"`. See `resolve_msigdb_alias()` for the supported aliases.
#' @return A named list; each element is a character vector of gene symbols
#'   for one gene set, named by `gs_name`.
#' @export
import_msigdbr <- function(org, collection) {
    df <- fetch_msigdbr_df(org, collection)
    split(df$gene_symbol, df$gs_name)
}

#' Fetch an MSigDB gene set collection as a clusterProfiler `TERM2GENE` table
#'
#' @param org `"human"` or `"mouse"`.
#' @param collection A single collection alias, e.g. `"HALLMARK"` or
#'   `"GOBP"`. See `resolve_msigdb_alias()` for the supported aliases.
#' @return A data frame with columns `term` (gene set name) and `gene`
#'   (gene symbol), the format `clusterProfiler::GSEA()` expects for its
#'   `TERM2GENE` argument.
#' @export
import_msigdbr_term2gene <- function(org, collection) {
    df <- fetch_msigdbr_df(org, collection)
    data.frame(term = df$gs_name, gene = df$gene_symbol, stringsAsFactors = FALSE)
}

#' Run GSEA (via `clusterProfiler::GSEA()`) for one ranked gene list against
#' one collection
#'
#' A second, independent GSEA implementation run alongside the package's
#' primary [run_gsea()]/`fgsea` pipeline, added specifically so the raw
#' `clusterProfiler` `gseaResult` S4 object -- not just a results data
#' frame -- is available for users who want it (e.g. for
#' `clusterProfiler`/`enrichplot` downstream visualization functions such
#' as `enrichplot::gseaplot2()`). Uses the same MSigDB gene sets (fetched
#' via [import_msigdbr_term2gene()]) and the same rank vector as the
#' `fgsea`-based pipeline, so results are directly comparable.
#'
#' @param ranks A named, sorted numeric vector as returned by
#'   [build_gsea_ranks()].
#' @param term2gene A `TERM2GENE` data frame as returned by
#'   [import_msigdbr_term2gene()] (columns `term`, `gene`).
#' @param min_size,max_size Minimum/maximum gene set size to test. Defaults
#'   `15`/`500`.
#' @param eps Boundary for estimating GSEA p-values; passed to
#'   `clusterProfiler::GSEA()`. Default `1e-10`.
#' @param seed Logical or integer; passed to `clusterProfiler::GSEA()`'s
#'   `seed` argument for reproducibility. Default `TRUE`.
#' @return The raw `gseaResult` S4 object from `clusterProfiler::GSEA()`
#'   (`pvalueCutoff = 1`, i.e. **unfiltered** -- every tested gene set is
#'   returned regardless of significance), or `NULL` if no gene set had
#'   enough member genes overlapping `ranks` to clear `min_size` (a
#'   warning is issued in that case rather than propagating
#'   `clusterProfiler`'s error).
#' @export
run_clusterprofiler_gsea <- function(ranks, term2gene, min_size = 15, max_size = 500,
                                      eps = 1e-10, seed = TRUE) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package \"clusterProfiler\" is required for run_clusterprofiler_gsea() ",
             "but is not installed.", call. = FALSE)
    }
    if (!is.numeric(ranks) || is.null(names(ranks))) {
        stop("`ranks` must be a named numeric vector (see build_gsea_ranks()).", call. = FALSE)
    }
    tryCatch(
        clusterProfiler::GSEA(
            geneList = ranks,
            TERM2GENE = term2gene,
            minGSSize = min_size,
            maxGSSize = max_size,
            eps = eps,
            pvalueCutoff = 1,
            seed = seed,
            verbose = FALSE
        ),
        error = function(e) {
            warning("clusterProfiler::GSEA() found no gene sets with enough member genes ",
                    "overlapping the ranked list (after min_size/max_size filtering): ",
                    conditionMessage(e), call. = FALSE)
            NULL
        }
    )
}

#' Run GSEA (via `fgsea`) for one ranked gene list against one collection
#'
#' @param ranks A named, sorted numeric vector as returned by
#'   [build_gsea_ranks()].
#' @param gene_sets A named list of character vectors (gene sets), as
#'   returned by [import_msigdbr()].
#' @param min_size,max_size Minimum/maximum gene set size to test. Defaults
#'   `15`/`500`.
#' @param eps Boundary for estimating GSEA p-values; passed to
#'   `fgsea::fgsea()`. Default `1e-10`.
#' @return A data frame of `fgsea` results, sorted by ascending `pval`, with
#'   list-columns coerced to comma-separated strings for CSV export.
#' @export
run_gsea <- function(ranks, gene_sets, min_size = 15, max_size = 500, eps = 1e-10) {
    if (!is.numeric(ranks) || is.null(names(ranks))) {
        stop("`ranks` must be a named numeric vector (see build_gsea_ranks()).", call. = FALSE)
    }
    res <- fgsea::fgsea(
        pathways = gene_sets,
        stats = ranks,
        minSize = min_size,
        maxSize = max_size,
        eps = eps
    )
    res <- as.data.frame(res)
    res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ",")
    # Keep fgsea's native `pval` / `padj` (BH); do not recompute FDR separately.
    res[order(res$pval), ]
}

#' Figure size (inches) for a GSEA NES barplot given the number of pathways shown
#'
#' Height scales with pathway count but stays tighter than a fixed 7" panel
#' (which vertically stretches bar labels when ~20 pathways are shown).
#' Width also grows mildly with pathway count.
#' @noRd
.gsea_barplot_size <- function(n_pathways) {
    n <- max(1L, as.integer(n_pathways))
    # 40% smaller than the previous footprint.
    list(
        width = 0.6 * max(5.5, min(8, 5 + 0.05 * n)),
        height = 0.6 * max(2.0, min(4.2, 0.7 * (0.22 * n + 1.2)))
    )
}

#' Run GSEA for every comparison and every requested MSigDB collection
#'
#' @param diffexp_results Either the `results` data frame from
#'   [run_diffexp_pipeline()] (with a `comparison` column) or its
#'   `by_comparison` list.
#' @param collections Character vector of MSigDB collection aliases:
#'   `"HALLMARK"`, `"GOBP"`, `"GOCC"`, `"GOMF"`, `"KEGG_MEDICUS"`,
#'   `"KEGG_LEGACY"`, `"REACTOME"`. Default `"HALLMARK"` only when calling
#'   this function directly ([run_deseq2_pipeline()] passes its own
#'   6-collection default -- see its `collections` argument). Organism
#'   (human/mouse) is auto-detected from the gene symbols in
#'   `diffexp_results[[gene_col]]` via [detect_organism()] -- there is no
#'   `org` argument to set.
#' @param ranking Ranking metric, passed to [build_gsea_ranks()]. Default
#'   `"pxfc"`.
#' @param gene_col Name of the gene-identifier column in `diffexp_results`.
#'   Default `"gene"`.
#' @param min_size,max_size Gene set size bounds, passed to [run_gsea()].
#' @param top_n_plot Number of top pathways per direction to show in the
#'   GSEA barplot. Default `10`.
#' @param save_dir Directory to write per-comparison GSEA outputs to.
#'   Default `getwd()`.
#' @param save_plot Logical. If `TRUE`, saves a GSEA barplot per comparison
#'   per collection. Default `TRUE`.
#' @param run_clusterprofiler Logical. If `TRUE` (default `FALSE`), also
#'   runs a second, independent GSEA pass via
#'   [run_clusterprofiler_gsea()] (`clusterProfiler::GSEA()`) for every
#'   comparison/collection, using the same rank vector and gene sets. The
#'   raw `gseaResult` S4 objects are returned in `by_comparison_cp` and
#'   also serialized to `<comparison>/gsea/<collection>_clusterprofiler_gsea.rds`;
#'   an unfiltered `as.data.frame()` export of each is additionally
#'   written to `<comparison>/gsea/<collection>_clusterprofiler_gsea.csv`.
#'   This is a second implementation run *alongside* (not instead of) the
#'   package's primary `fgsea`-based results in `results`/`by_comparison`.
#' @return A list with elements `results` (combined data frame across all
#'   comparisons and collections, from the `fgsea`-based pipeline),
#'   `by_comparison` (nested list, indexed first by comparison then by
#'   collection, from the `fgsea`-based pipeline), and, when
#'   `run_clusterprofiler = TRUE`, `by_comparison_cp` (nested list of raw
#'   `gseaResult` S4 objects, indexed first by comparison then by
#'   collection; `NULL` for any comparison/collection where
#'   `clusterProfiler::GSEA()` found no testable gene sets).
#' @examples
#' \dontrun{
#' gsea_out <- run_gsea_pipeline(diffexp_out$results, collections = "HALLMARK")
#' }
#' @export
run_gsea_pipeline <- function(
    diffexp_results,
    collections = "HALLMARK",
    ranking = "log2FoldChange",
    gene_col = "gene",
    min_size = 15,
    max_size = 500,
    top_n_plot = 10,
    save_dir = getwd(),
    save_plot = TRUE,
    run_clusterprofiler = FALSE) {

    collections <- validate_msigdb_collections(collections)
    ranking <- validate_ranking(ranking)
    validate_logical(save_plot, "save_plot")
    validate_logical(run_clusterprofiler, "run_clusterprofiler")
    save_dir <- validate_save_dir(save_dir)

    if (run_clusterprofiler && !requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package \"clusterProfiler\" is required when run_clusterprofiler = TRUE ",
             "but is not installed.", call. = FALSE)
    }

    if (is.list(diffexp_results) && !is.data.frame(diffexp_results)) {
        diffexp_results <- do.call(rbind, diffexp_results)
    }
    if (!"comparison" %in% colnames(diffexp_results)) {
        stop("`diffexp_results` must contain a `comparison` column (as produced by run_diffexp_pipeline()).", call. = FALSE)
    }
    if (!gene_col %in% colnames(diffexp_results)) {
        stop("`diffexp_results` must contain a \"", gene_col, "\" gene-identifier column.", call. = FALSE)
    }
    org <- detect_organism(diffexp_results[[gene_col]])
    message("Detected organism: ", org)

    gene_set_cache <- stats::setNames(
        lapply(collections, function(coll) import_msigdbr(org, coll)),
        collections
    )
    term2gene_cache <- if (run_clusterprofiler) {
        stats::setNames(
            lapply(collections, function(coll) import_msigdbr_term2gene(org, coll)),
            collections
        )
    } else {
        NULL
    }

    comparisons <- unique(diffexp_results$comparison)
    by_comparison <- list()
    ranks_by_comparison <- list()
    by_comparison_cp <- if (run_clusterprofiler) list() else NULL

    for (comparison in comparisons) {
        message("Running GSEA: ", comparison, " ...")
        comp_df <- diffexp_results[diffexp_results$comparison == comparison, ]
        ranks <- build_gsea_ranks(comp_df, gene_col = gene_col, ranking = ranking)
        ranks_by_comparison[[comparison]] <- ranks

        comp_dir <- file.path(save_dir, comparison, "gsea")
        by_collection <- list()
        by_collection_cp <- if (run_clusterprofiler) list() else NULL

        for (coll in collections) {
            gsea_res <- run_gsea(ranks, gene_set_cache[[coll]], min_size = min_size, max_size = max_size)
            # `df$col <- scalar` errors on a 0-row data.frame ("replacement
            # has 1 row, data has 0"), which fgsea can legitimately return
            # when no gene set in this collection has enough member genes
            # overlapping `ranks` to clear `min_size`. Guard explicitly so
            # a collection/comparison with no testable gene sets doesn't
            # crash the whole pipeline -- an empty-but-correctly-shaped
            # result is written instead.
            if (nrow(gsea_res) > 0) {
                gsea_res$collection <- coll
                gsea_res$comparison <- comparison
                gsea_res$org <- org
            } else {
                gsea_res$collection <- character(0)
                gsea_res$comparison <- character(0)
                gsea_res$org <- character(0)
            }

            coll_label <- gsub("[^A-Za-z0-9]+", "_", coll)
            write_output(gsea_res, file.path(comp_dir, paste0(coll_label, "_gsea.csv")), format = "csv")

            if (nrow(gsea_res) == 0) {
                warning("No gene sets in collection '", coll, "' had enough member genes overlapping ",
                        "the ranked list for comparison '", comparison, "' (after min_size/max_size ",
                        "filtering); an empty GSEA result was written for this collection/comparison.",
                        call. = FALSE)
            }

            if (save_plot && nrow(gsea_res) > 0) {
                p <- plot_gsea_barplot(gsea_res, top_n = top_n_plot, title = comparison)
                sz <- attr(p, "plot_size")
                if (is.null(sz)) sz <- list(width = 6.5, height = 5)
                save_ggplot(p, file.path(comp_dir, paste0(coll_label, "_gsea_barplot.pdf")),
                            width = sz$width, height = sz$height)
                p_bc <- tryCatch(
                    plot_gsea_barcode_panel(
                        gsea_res, ranks = ranks,
                        gene_sets = gene_set_cache[[coll]],
                        top_n = top_n_plot, title = comparison
                    ),
                    error = function(e) {
                        warning("GSEA barcode panel failed for ", comparison, "/", coll,
                                ": ", conditionMessage(e), call. = FALSE)
                        NULL
                    }
                )
                if (!is.null(p_bc)) {
                    sz_bc <- attr(p_bc, "plot_size")
                    if (is.null(sz_bc)) sz_bc <- list(width = 4.25, height = 2.25)
                    save_ggplot(p_bc, file.path(comp_dir, paste0(coll_label, "_gsea_barcode.pdf")),
                                width = sz_bc$width, height = sz_bc$height)
                }
            }

            by_collection[[coll]] <- gsea_res

            if (run_clusterprofiler) {
                cp_res <- run_clusterprofiler_gsea(ranks, term2gene_cache[[coll]],
                                                    min_size = min_size, max_size = max_size)
                if (!is.null(cp_res)) {
                    dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
                    saveRDS(cp_res, file.path(comp_dir, paste0(coll_label, "_clusterprofiler_gsea.rds")))
                    write_output(as.data.frame(cp_res),
                                 file.path(comp_dir, paste0(coll_label, "_clusterprofiler_gsea.csv")),
                                 format = "csv")
                }
                by_collection_cp[[coll]] <- cp_res
            }
        }
        by_comparison[[comparison]] <- by_collection
        if (run_clusterprofiler) by_comparison_cp[[comparison]] <- by_collection_cp
    }

    combined <- do.call(rbind, lapply(by_comparison, function(x) do.call(rbind, x)))
    rownames(combined) <- NULL

    out <- list(
        results = combined,
        by_comparison = by_comparison,
        ranks_by_comparison = ranks_by_comparison,
        gene_sets = gene_set_cache,
        ranking = ranking,
        org = org
    )
    if (run_clusterprofiler) out$by_comparison_cp <- by_comparison_cp
    out
}

#' Convert BH-adjusted p-values to conventional significance stars
#' @param x Numeric vector of `padj` values (may contain `NA`).
#' @return Character vector: `"***"` (< 0.001), `"**"` (< 0.01), `"*"`
#'   (< 0.05), or `""` otherwise (including `NA`).
#' @noRd
padj_to_stars <- function(x) {
    stars <- ifelse(is.na(x), "",
              ifelse(x < 0.001, "***",
              ifelse(x < 0.01, "**",
              ifelse(x < 0.05, "*", ""))))
    stars
}

#' Select the same pathways shown by [plot_gsea_barplot()]
#'
#' Top `top_n` pathways per NES direction by `|NES|`, with cosmetic
#' `pathway_label` (collection prefix stripped).
#' @noRd
.gsea_barplot_pathways <- function(gsea_res, top_n = 10, signif_only = FALSE) {
    df <- gsea_res
    if (!"padj" %in% colnames(df)) {
        stop("`gsea_res` must contain a `padj` column (as returned by fgsea / run_gsea()).",
             call. = FALSE)
    }
    if (signif_only) df <- df[!is.na(df$padj) & df$padj < 0.05, ]
    if (nrow(df) == 0) {
        stop("No pathways to plot (after signif_only filtering, if applied).", call. = FALSE)
    }
    df$direction <- ifelse(df$NES > 0, "up", "down")
    df <- df[order(-abs(df$NES)), ]
    top_df <- do.call(rbind, lapply(split(df, df$direction), function(d) {
        d[seq_len(min(top_n, nrow(d))), ]
    }))
    top_df$pathway_label <- as.character(top_df$pathway)
    if ("collection" %in% colnames(top_df)) {
        for (coll in unique(top_df$collection)) {
            spec <- .msigdb_collection_map[[coll]]
            if (!is.null(spec)) {
                rows <- top_df$collection == coll
                top_df$pathway_label[rows] <- sub(
                    paste0("^", spec$prefix), "", top_df$pathway_label[rows]
                )
            }
        }
    } else {
        for (spec in .msigdb_collection_map) {
            top_df$pathway_label <- sub(
                paste0("^", spec$prefix), "", top_df$pathway_label
            )
        }
    }
    top_df
}

#' GSEA barplot: top and bottom enriched pathways by normalized enrichment score
#'
#' Shows the `top_n` most positive (up) and `top_n` most negative (down)
#' pathways by NES, **regardless of significance** -- i.e. not restricted
#' to `padj < 0.05` by default -- with conventional significance stars
#' (`*` padj < 0.05, `**` < 0.01, `***` < 0.001) from fgsea's BH-adjusted
#' `padj` annotated just beyond each bar's tip. Pass `signif_only = TRUE`
#' to restrict to significant pathways first, as before. Bars are outlined
#' in black and coloured on a continuous red-blue diverging scale by NES
#' (red = positive/up, blue = negative/down). Pathway names are drawn as
#' in-plot text labels hugging the zero line on the side opposite their
#' bar (rather than conventional axis tick labels), matching common
#' published GSEA-barplot conventions.
#'
#' @param gsea_res A data frame as returned by [run_gsea()] (must contain
#'   `pathway`, `NES`, and `padj`).
#' @param top_n Number of pathways *per direction* (up/down by NES) to
#'   plot. Default `10`.
#' @param signif_only Logical. If `TRUE`, restricts to `padj < 0.05`
#'   before selecting the top/bottom pathways. Default `FALSE` (shows
#'   top/bottom by NES regardless of significance).
#' @param title Optional plot title (e.g. the comparison name), shown bold
#'   and centred. Default `NULL` (no title).
#' @return A `ggplot` object: horizontal barplot of NES, ordered by NES,
#'   coloured by a continuous red-blue diverging NES scale, with pathway
#'   names and significance stars labelled directly on the plot. The
#'   returned object carries `attr(,"plot_size")` with suggested
#'   `width`/`height` in inches scaled to the number of pathways shown.
#' @export
plot_gsea_barplot <- function(gsea_res, top_n = 10, signif_only = FALSE, title = NULL) {
    top_df <- .gsea_barplot_pathways(gsea_res, top_n = top_n, signif_only = signif_only)
    top_df$pathway <- factor(top_df$pathway, levels = top_df$pathway[order(top_df$NES)])
    top_df$stars <- padj_to_stars(top_df$padj)
    nes_limit <- max(abs(top_df$NES), na.rm = TRUE)

    # Nudge the star label beyond the end of each bar, on the side the bar
    # extends towards, so it doesn't overlap the column itself; `hjust`
    # (rather than just a bigger nudge) is what actually pushes the star
    # glyphs further away from the bar tip, since text width is not in
    # data units.
    label_nudge <- nes_limit * 0.06
    top_df$label_y <- top_df$NES + ifelse(top_df$NES >= 0, label_nudge, -label_nudge)
    top_df$star_hjust <- ifelse(top_df$NES >= 0, 0, 1)

    # Pathway name sits in the empty half of the panel (opposite its own
    # bar): right-aligned for up (positive-NES) pathways, left-aligned for
    # down (negative-NES) pathways. The anchor point is nudged away from
    # NES = 0 by the same margin as the star labels, so the text doesn't
    # crowd the zero-reference line -- hjust then extends it further still
    # into the empty half.
    top_df$name_hjust <- ifelse(top_df$NES >= 0, 1, 0)
    top_df$name_x <- ifelse(top_df$NES >= 0, -label_nudge, label_nudge)

    lvls <- .parse_contrast_levels(if (is.null(title) || !nzchar(title)) "up_vs_down" else title)
    up_label <- paste0(lvls$test, " \u2192")
    down_label <- paste0("\u2190 ", lvls$ref)

    p <- ggplot2::ggplot(top_df, ggplot2::aes(x = .data$pathway, y = .data$NES, fill = .data$NES)) +
        ggplot2::geom_col(width = 0.86, colour = "black", linewidth = 0.3) +
        ggplot2::geom_text(ggplot2::aes(y = .data$name_x, label = .data$pathway_label, hjust = .data$name_hjust),
                            size = 3.2) +
        ggplot2::geom_text(ggplot2::aes(y = .data$label_y, label = .data$stars, hjust = .data$star_hjust),
                            size = 4.5) +
        ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
        unclipped_coord_flip() +
        ggplot2::scale_y_continuous(limits = c(-nes_limit, nes_limit),
                                     expand = ggplot2::expansion(mult = c(0.25, 0.25))) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                                       midpoint = 0, limits = c(-nes_limit, nes_limit)) +
        # After coord_flip, NES is horizontal: left = down, right = up.
        ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
                           label = up_label, colour = "#B2182B",
                           size = 5, fontface = "bold") +
        ggplot2::annotate("text", x = Inf, y = -Inf, hjust = -0.05, vjust = 1.3,
                           label = down_label, colour = "#2166AC",
                           size = 5, fontface = "bold") +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
        ) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab("Normalized Enrichment Score (NES)") +
        ggplot2::labs(fill = "NES", title = title)

    attr(p, "plot_size") <- .gsea_barplot_size(nrow(top_df))
    p
}

#' GSEA gene-set barcode panel
#'
#' For the same pathways shown by [plot_gsea_barplot()], draws a multi-row
#' enrichment barcode: pathway names on the left, member-gene positions
#' along the ranked list in the middle (one track per pathway; red =
#' `NES > 0`, blue = `NES < 0`), and `padj` + `NES` on the right.
#'
#' @param gsea_res A data frame as returned by [run_gsea()].
#' @param ranks Named numeric gene ranking vector (e.g. from
#'   [build_gsea_ranks()]), decreasing order preferred but not required.
#' @param gene_sets Named list of character gene-symbol vectors (e.g. from
#'   [import_msigdbr()]).
#' @param top_n Pathways per NES direction, matching the barplot. Default
#'   `10`.
#' @param signif_only Passed through to pathway selection. Default `FALSE`.
#' @param title Optional title.
#' @return A `ggplot` object with `attr(,"plot_size")` width/height inches.
#' @export
plot_gsea_barcode_panel <- function(gsea_res, ranks, gene_sets,
                                    top_n = 10, signif_only = FALSE,
                                    title = NULL) {
    if (is.null(ranks) || !length(ranks) || is.null(names(ranks))) {
        stop("`ranks` must be a named numeric vector.", call. = FALSE)
    }
    if (is.null(gene_sets) || !length(gene_sets)) {
        stop("`gene_sets` must be a named list of gene sets.", call. = FALSE)
    }
    top_df <- .gsea_barplot_pathways(gsea_res, top_n = top_n, signif_only = signif_only)
    # High NES at top of panel.
    top_df <- top_df[order(-top_df$NES, top_df$pathway), , drop = FALSE]
    top_df$pathway_label <- factor(top_df$pathway_label, levels = rev(unique(top_df$pathway_label)))
    top_df$nes_dir <- ifelse(top_df$NES >= 0, "pos", "neg")

    ranks <- sort(ranks, decreasing = TRUE)
    n_rank <- length(ranks)
    rank_pos <- stats::setNames(seq_along(ranks), names(ranks))

    tick_list <- lapply(seq_len(nrow(top_df)), function(i) {
        pw <- as.character(top_df$pathway[i])
        lab <- as.character(top_df$pathway_label[i])
        members <- intersect(gene_sets[[pw]], names(ranks))
        if (!length(members)) return(NULL)
        data.frame(
            pathway_label = lab,
            gene = members,
            x = as.numeric(rank_pos[members]),
            nes_dir = top_df$nes_dir[i],
            stringsAsFactors = FALSE
        )
    })
    ticks <- do.call(rbind, tick_list)
    if (is.null(ticks) || !nrow(ticks)) {
        stop("No pathway member genes overlap the ranked list for barcode plotting.",
             call. = FALSE)
    }
    ticks$pathway_label <- factor(ticks$pathway_label, levels = levels(top_df$pathway_label))
    ticks$nes_dir <- factor(ticks$nes_dir, levels = c("pos", "neg"))
    # Single track per pathway row; colour by pathway NES sign.
    ticks$y <- as.numeric(ticks$pathway_label)
    ticks$ymin <- ticks$y - 0.35
    ticks$ymax <- ticks$y + 0.35

    meta <- top_df
    meta$y <- as.numeric(meta$pathway_label)
    meta$padj_lab <- ifelse(
        is.na(meta$padj), "NA",
        formatC(meta$padj, format = "e", digits = 2)
    )
    meta$nes_lab <- formatC(meta$NES, format = "f", digits = 2)

    x_padj <- n_rank * 1.18
    x_nes <- n_rank * 1.32
    x_max <- n_rank * 1.42

    p <- ggplot2::ggplot() +
        ggplot2::geom_rect(
            data = meta,
            ggplot2::aes(xmin = 0.5, xmax = n_rank + 0.5,
                         ymin = .data$y - 0.45, ymax = .data$y + 0.45),
            fill = NA, colour = "grey40", linewidth = 0.3
        ) +
        ggplot2::geom_segment(
            data = ticks,
            ggplot2::aes(x = .data$x, xend = .data$x,
                         y = .data$ymin, yend = .data$ymax,
                         colour = .data$nes_dir),
            linewidth = 0.35, alpha = 0.9
        ) +
        ggplot2::geom_text(
            data = meta,
            ggplot2::aes(x = 0, y = .data$y, label = .data$pathway_label),
            hjust = 1, size = 2.8, colour = "grey15"
        ) +
        ggplot2::geom_text(
            data = meta,
            ggplot2::aes(x = x_padj, y = .data$y, label = .data$padj_lab),
            hjust = 0, size = 2.6, colour = "grey20"
        ) +
        ggplot2::geom_text(
            data = meta,
            ggplot2::aes(x = x_nes, y = .data$y, label = .data$nes_lab),
            hjust = 0, size = 2.6, colour = "grey20"
        ) +
        ggplot2::annotate(
            "text", x = 0, y = max(meta$y) + 0.85, label = "Pathway",
            hjust = 1, size = 3.2, fontface = "bold"
        ) +
        ggplot2::annotate(
            "text", x = x_padj, y = max(meta$y) + 0.85, label = "padj",
            hjust = 0, size = 3.2, fontface = "bold"
        ) +
        ggplot2::annotate(
            "text", x = x_nes, y = max(meta$y) + 0.85, label = "NES",
            hjust = 0, size = 3.2, fontface = "bold"
        ) +
        ggplot2::scale_colour_manual(
            values = c(pos = "#B2182B", neg = "#2166AC"),
            labels = c(pos = "NES > 0", neg = "NES < 0"),
            name = NULL,
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            limits = c(-n_rank * 0.42, x_max),
            expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0.4, max(meta$y) + 1.2),
            expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(
            legend.position = "bottom",
            plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
            plot.margin = ggplot2::margin(8, 12, 8, 80)
        ) +
        ggplot2::labs(
            title = title,
            subtitle = "Gene enrichment barcode (red = NES > 0 / blue = NES < 0)"
        )

    n_pw <- nrow(meta)
    # Display / PDF size: half of the previous barcode footprint.
    attr(p, "plot_size") <- list(
        width = 4.25,
        height = max(1.31, 0.1575 * n_pw + 0.525)
    )
    p
}

#' Top and bottom pathways by NES, including leading-edge gene sets
#'
#' Selects the `top_n` pathways with the most positive NES and the `top_n`
#' with the most negative NES (regardless of significance), for reporting
#' alongside [plot_gsea_barplot()]. Unlike a plain top-N-by-p-value table,
#' this always surfaces both enrichment directions and keeps the
#' `leadingEdge` column (already comma-separated by [run_gsea()]).
#'
#' @param gsea_res A data frame as returned by [run_gsea()].
#' @param top_n Number of pathways *per direction* to return. Default `10`.
#' @return A data frame: up to `2 * top_n` rows, ordered by descending NES,
#'   with all columns of `gsea_res` (including `leadingEdge`).
#' @export
gsea_top_bottom_table <- function(gsea_res, top_n = 10) {
    df <- gsea_res[order(-gsea_res$NES), ]
    n <- nrow(df)
    top_n <- min(top_n, ceiling(n / 2))
    top_idx <- seq_len(min(top_n, n))
    bottom_idx <- seq(to = n, length.out = min(top_n, n))
    idx <- union(top_idx, bottom_idx)
    df[idx, ]
}

#' Interactive, sortable/searchable GSEA results table
#'
#' Wraps a GSEA result table (e.g. from [gsea_top_bottom_table()] or the
#' full `gsea_res` data frame) in a `DT::datatable()`, styled the same way
#' as [make_deg_datatable()]: sortable/searchable columns, a CSV export
#' button, and rounded numeric display (the underlying data is untouched).
#'
#' @param gsea_res A data frame as returned by [run_gsea()] or
#'   [gsea_top_bottom_table()].
#' @param digits Number of decimal places to round numeric columns to for
#'   display. Default `4`.
#' @param page_length Number of rows shown per page. Default `10`.
#' @return A `DT::datatable` htmlwidget.
#' @export
make_gsea_datatable <- function(gsea_res, digits = 4, page_length = 10) {
    if (!requireNamespace("DT", quietly = TRUE)) {
        stop("Package \"DT\" is required for interactive tables. Install it with ",
             "`install.packages(\"DT\")`.", call. = FALSE)
    }
    display_df <- gsea_res
    num_cols <- vapply(display_df, is.numeric, logical(1))
    display_df[num_cols] <- lapply(display_df[num_cols], round, digits = digits)

    DT::datatable(
        display_df,
        filter = "top",
        extensions = "Buttons",
        rownames = FALSE,
        options = list(
            pageLength = page_length,
            lengthMenu = c(10, 25, 50, 100, nrow(display_df)),
            scrollX = TRUE,
            dom = "Blfrtip",
            buttons = c("csv")
        )
    )
}
