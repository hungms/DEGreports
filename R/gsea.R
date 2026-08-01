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
    # Storey q-values alongside BH padj (used for barplot stars / filtering).
    if (nrow(res) > 0) {
        res$qvalue <- .compute_qvalues(res$pval)
    } else {
        res$qvalue <- numeric(0)
    }
    res[order(res$pval), ]
}

#' Compute Storey q-values when the `qvalue` package is available; otherwise
#' fall back to BH-adjusted p-values so a `qvalue` column is always present.
#' @noRd
.compute_qvalues <- function(pval) {
    if (requireNamespace("qvalue", quietly = TRUE)) {
        qv <- tryCatch(qvalue::qvalue(pval)$qvalues, error = function(e) NULL)
        if (!is.null(qv)) {
            return(as.numeric(qv))
        }
    }
    stats::p.adjust(pval, method = "BH")
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
    ranking = "pxfc",
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
    by_comparison_cp <- if (run_clusterprofiler) list() else NULL

    for (comparison in comparisons) {
        message("Running GSEA: ", comparison, " ...")
        comp_df <- diffexp_results[diffexp_results$comparison == comparison, ]
        ranks <- build_gsea_ranks(comp_df, gene_col = gene_col, ranking = ranking)

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

    out <- list(results = combined, by_comparison = by_comparison)
    if (run_clusterprofiler) out$by_comparison_cp <- by_comparison_cp
    out
}

#' Convert significance values (padj or qvalue) to conventional stars
#' @param x Numeric vector of adjusted p- / q-values (may contain `NA`).
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

#' GSEA barplot: top and bottom enriched pathways by normalized enrichment score
#'
#' Shows the `top_n` most positive (up) and `top_n` most negative (down)
#' pathways by NES, **regardless of significance** -- i.e. not restricted
#' to `qvalue < 0.05` by default -- with conventional significance stars
#' (`*` qvalue < 0.05, `**` < 0.01, `***` < 0.001) annotated just beyond
#' each bar's tip. Pass `signif_only = TRUE` to restrict to significant
#' pathways first, as before. Bars are outlined in black and coloured on
#' a continuous red-blue diverging scale by NES (red = positive/up, blue
#' = negative/down). Pathway names are drawn as in-plot text labels
#' hugging the zero line on the side opposite their bar (rather than
#' conventional axis tick labels), matching common published GSEA-barplot
#' conventions.
#'
#' @param gsea_res A data frame as returned by [run_gsea()] (must contain
#'   `pathway`, `NES`, and `qvalue` -- or `padj` as a fallback).
#' @param top_n Number of pathways *per direction* (up/down by NES) to
#'   plot. Default `10`.
#' @param signif_only Logical. If `TRUE`, restricts to `qvalue < 0.05`
#'   (falling back to `padj` if `qvalue` is absent) before selecting the
#'   top/bottom pathways. Default `FALSE` (shows top/bottom by NES
#'   regardless of significance).
#' @param title Optional plot title (e.g. the comparison name), shown bold
#'   and centred. Default `NULL` (no title).
#' @return A `ggplot` object: horizontal barplot of NES, ordered by NES,
#'   coloured by a continuous red-blue diverging NES scale, with pathway
#'   names and significance stars labelled directly on the plot. The
#'   returned object carries `attr(,"plot_size")` with suggested
#'   `width`/`height` in inches scaled to the number of pathways shown.
#' @export
plot_gsea_barplot <- function(gsea_res, top_n = 10, signif_only = FALSE, title = NULL) {
    df <- gsea_res
    signif_col <- if ("qvalue" %in% colnames(df)) "qvalue" else "padj"
    if (signif_only) df <- df[!is.na(df[[signif_col]]) & df[[signif_col]] < 0.05, ]
    if (nrow(df) == 0) {
        stop("No pathways to plot (after signif_only filtering, if applied).", call. = FALSE)
    }
    df$direction <- ifelse(df$NES > 0, "up", "down")
    df <- df[order(-abs(df$NES)), ]
    top_df <- do.call(rbind, lapply(split(df, df$direction), function(d) {
        d[seq_len(min(top_n, nrow(d))), ]
    }))
    top_df$pathway <- factor(top_df$pathway, levels = top_df$pathway[order(top_df$NES)])
    # Cosmetic-only display label: strip the collection's gene-set-name
    # prefix (e.g. "HALLMARK_", "GOBP_", "KEGG_MEDICUS_") -- the underlying
    # `pathway` column, and any downstream table, is untouched, this only
    # affects what's drawn on the bars. If `gsea_res$collection` is present
    # (as attached by run_gsea_pipeline()) the exact prefix for each row's
    # collection is looked up and stripped; otherwise every known prefix is
    # tried (a plain sub() on each is a no-op for non-matching rows).
    top_df$pathway_label <- as.character(top_df$pathway)
    if ("collection" %in% colnames(top_df)) {
        for (coll in unique(top_df$collection)) {
            spec <- .msigdb_collection_map[[coll]]
            if (!is.null(spec)) {
                rows <- top_df$collection == coll
                top_df$pathway_label[rows] <- sub(paste0("^", spec$prefix), "", top_df$pathway_label[rows])
            }
        }
    } else {
        for (spec in .msigdb_collection_map) {
            top_df$pathway_label <- sub(paste0("^", spec$prefix), "", top_df$pathway_label)
        }
    }
    top_df$stars <- padj_to_stars(top_df[[signif_col]])
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
    # Prefer column order ... padj, qvalue, ... when both are present.
    if (all(c("padj", "qvalue") %in% colnames(display_df))) {
        cols <- colnames(display_df)
        cols <- cols[cols != "qvalue"]
        padj_i <- match("padj", cols)
        cols <- append(cols, "qvalue", after = padj_i)
        display_df <- display_df[, cols, drop = FALSE]
    }
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
