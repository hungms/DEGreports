#' Validate a DESeq2 dataset object
#'
#' @param dds A `DESeqDataSet` object.
#' @return The validated `dds` object (unmodified).
#' @noRd
validate_dds <- function(dds) {
    if (!methods::is(dds, "DESeqDataSet")) {
        stop("`dds` must be a DESeqDataSet object (see DESeq2::DESeqDataSetFromMatrix()).", call. = FALSE)
    }
    if (ncol(dds) < 2) {
        stop("`dds` must contain at least 2 samples.", call. = FALSE)
    }
    dds
}

#' Validate grouping variable(s) and coerce each to a factor
#'
#' @param dds A `DESeqDataSet` object.
#' @param var Column name(s) in `colData(dds)` used for grouping. Length 1
#'   (single-factor), 2 (two-factor), or 3 (three-factor; typically with
#'   `include_interactions = TRUE`).
#' @return `dds` with each `colData(dds)[[var]]` coerced to a factor.
#' @noRd
validate_var <- function(dds, var) {
    if (!is.character(var) || length(var) < 1 || length(var) > 3 || anyNA(var)) {
        stop("`var` must be a character vector of length 1, 2, or 3 naming column(s) in colData(dds).",
             call. = FALSE)
    }
    if (anyDuplicated(var)) {
        stop("`var` must not contain duplicated column names.", call. = FALSE)
    }
    cd_names <- colnames(SummarizedExperiment::colData(dds))
    missing <- setdiff(var, cd_names)
    if (length(missing) > 0) {
        stop("`var` column(s) not in colData(dds): ", paste(missing, collapse = ", "),
             ". Available columns: ", paste(cd_names, collapse = ", "), call. = FALSE)
    }
    for (v in var) {
        grp <- SummarizedExperiment::colData(dds)[[v]]
        if (!is.factor(grp)) grp <- factor(grp)
        if (nlevels(grp) < 2) {
            stop("`var` (\"", v, "\") must have at least 2 levels; found ", nlevels(grp), ".",
                 call. = FALSE)
        }
        tab <- table(grp)
        if (any(tab < 2)) {
            warning("Level(s) with fewer than 2 samples in `", v, "`: ",
                    paste(names(tab)[tab < 2], collapse = ", "),
                    ". Comparisons involving these levels may be unreliable or fail.",
                    call. = FALSE)
        }
        SummarizedExperiment::colData(dds)[[v]] <- grp
    }
    dds
}

#' Validate and apply a design formula
#'
#' @param dds A `DESeqDataSet` object.
#' @param design A design formula, a character string coercible to one, or
#'   `NULL` (caller should have already resolved via [resolve_design()]).
#' @return `dds` with `design(dds)` set to the validated formula.
#' @noRd
validate_design <- function(dds, design) {
    if (is.null(design)) {
        stop("`design` resolved to `NULL`; pass an explicit design or let ",
             "`resolve_design(var, batch)` build one.", call. = FALSE)
    }
    if (is.character(design)) design <- stats::as.formula(design)
    if (!inherits(design, "formula")) {
        stop("`design` must be a formula or a character string coercible to one (e.g. \"~ condition\").", call. = FALSE)
    }
    design_vars <- all.vars(design)
    missing_vars <- setdiff(design_vars, colnames(SummarizedExperiment::colData(dds)))
    if (length(missing_vars) > 0) {
        stop("Design formula references column(s) not present in colData(dds): ",
             paste(missing_vars, collapse = ", "), call. = FALSE)
    }
    DESeq2::design(dds) <- design
    dds
}

#' Validate that a single string is one of a fixed set of allowed choices
#'
#' Shared by every "pick one of these options" argument (`org`, `shrinkage`,
#' `ranking`, `output_format`) so the same check and error
#' wording isn't repeated five times over.
#'
#' @param x The value to validate.
#' @param valid Character vector of allowed values.
#' @param name Argument name, used in the error message.
#' @param lower If `TRUE` (default), lower-case `x` before comparing/returning.
#' @return The validated (optionally lower-cased) string.
#' @noRd
validate_choice <- function(x, valid, name, lower = TRUE) {
    if (!is.character(x) || length(x) != 1 || is.na(x)) {
        stop("`", name, "` must be a single character string.", call. = FALSE)
    }
    if (lower) x <- tolower(x)
    if (!x %in% valid) {
        stop("`", name, "` must be one of: ", paste(valid, collapse = ", "), "; got \"", x, "\".", call. = FALSE)
    }
    x
}

#' Validate an organism identifier
#' @param org Character, "human" or "mouse".
#' @return The validated, lower-cased `org` string.
#' @noRd
validate_org <- function(org) {
    validate_choice(org, c("human", "mouse"), "org")
}

#' Validate a logical argument (or vector thereof)
#' @param x Logical value(s) to validate.
#' @param name Name to use in the error message.
#' @return Invisibly, `NULL`.
#' @noRd
validate_logical <- function(x, name = "argument") {
    if (!is.logical(x) || anyNA(x)) {
        stop("`", name, "` must be TRUE or FALSE (no NAs).", call. = FALSE)
    }
    invisible(NULL)
}

#' Validate the minimum count filtering threshold
#' @param min_count Numeric threshold, coerced to a non-negative integer.
#' @return The validated integer threshold.
#' @noRd
validate_min_count <- function(min_count) {
    if (!is.numeric(min_count) || length(min_count) != 1 || is.na(min_count) || min_count < 0) {
        stop("`min_count` must be a single non-negative number.", call. = FALSE)
    }
    as.integer(round(min_count))
}

#' Validate the significance threshold (adjusted p-value)
#' @param alpha Numeric in (0, 1).
#' @return The validated `alpha`.
#' @noRd
validate_alpha <- function(alpha) {
    if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
        stop("`alpha` must be a single number strictly between 0 and 1.", call. = FALSE)
    }
    alpha
}

#' Validate the log2 fold-change threshold used for significance calling
#' @param lfc_threshold Non-negative numeric.
#' @return The validated `lfc_threshold`.
#' @noRd
validate_lfc_threshold <- function(lfc_threshold) {
    if (!is.numeric(lfc_threshold) || length(lfc_threshold) != 1 || is.na(lfc_threshold) || lfc_threshold < 0) {
        stop("`lfc_threshold` must be a single non-negative number.", call. = FALSE)
    }
    lfc_threshold
}

#' Validate the LFC shrinkage method
#' @param shrinkage One of "ashr", "apeglm", "normal", "none".
#' @return The validated, lower-cased `shrinkage` string.
#' @noRd
validate_shrinkage <- function(shrinkage) {
    shrinkage <- validate_choice(shrinkage, c("ashr", "apeglm", "normal", "none"), "shrinkage")
    # DESeq2::lfcShrink() dispatches to these packages internally by name
    # (type = "ashr" / "apeglm") without them appearing in our own code as
    # `pkg::fun()` calls; check for their presence explicitly here so a
    # missing dependency fails with a clear message at argument-validation
    # time rather than deep inside lfcShrink().
    if (shrinkage == "ashr" && !requireNamespace("ashr", quietly = TRUE)) {
        stop("`shrinkage = \"ashr\"` requires the 'ashr' package to be installed.", call. = FALSE)
    }
    if (shrinkage == "apeglm" && !requireNamespace("apeglm", quietly = TRUE)) {
        stop("`shrinkage = \"apeglm\"` requires the 'apeglm' package to be installed.", call. = FALSE)
    }
    shrinkage
}

#' Validate the DEG/GSEA ranking metric
#' @param ranking One of "pxfc", "padj", "log2FoldChange", "stat".
#' @return The validated `ranking` string.
#' @noRd
validate_ranking <- function(ranking) {
    # not lower-cased: "log2FoldChange" is camelCase, so casing is significant
    validate_choice(ranking, c("pxfc", "padj", "log2FoldChange", "stat"), "ranking", lower = FALSE)
}

#' Validate the batch-correction column(s)
#' @param dds A `DESeqDataSet` object.
#' @param batch Column name(s) in `colData(dds)`, or `NULL`.
#' @return `dds` with each `colData(dds)[[batch]]` coerced to a factor (if
#'   not `NULL`).
#' @noRd
validate_batch <- function(dds, batch) {
    if (is.null(batch)) return(dds)
    if (!is.character(batch) || !length(batch) || anyNA(batch) ||
            any(!nzchar(batch)) || anyDuplicated(batch)) {
        stop("`batch` must be `NULL` or a character vector of distinct column names in colData(dds).",
             call. = FALSE)
    }
    cd_names <- colnames(SummarizedExperiment::colData(dds))
    missing <- setdiff(batch, cd_names)
    if (length(missing)) {
        stop("`batch` column(s) not in colData(dds): ", paste(missing, collapse = ", "),
             ". Available columns: ", paste(cd_names, collapse = ", "), call. = FALSE)
    }
    for (b in batch) {
        x <- SummarizedExperiment::colData(dds)[[b]]
        if (!is.factor(x)) x <- factor(x)
        if (nlevels(x) < 2) {
            stop("`batch` (\"", b, "\") must have at least 2 levels to be worth correcting for.",
                 call. = FALSE)
        }
        SummarizedExperiment::colData(dds)[[b]] <- x
    }
    dds
}

#' Validate the pairing column of a paired design
#'
#' `pair` names the unit measured under several conditions -- a donor,
#' patient or animal. It enters every model as a blocking term
#' (`~ pair + var`), so each comparison is made within units, and it is
#' removed from the sample-level matrix alongside `batch`.
#'
#' A paired design only estimates factors that vary *within* a unit, so
#' this errors when every unit sits in a single level of a tested factor
#' (the factor is between-unit and would be aliased with `pair`), warns
#' about units that contribute nothing to a comparison, and errors when a
#' `batch` column is constant within every unit (it is then already
#' absorbed by `pair`).
#'
#' @param dds A `DESeqDataSet` object.
#' @param pair Column name in `colData(dds)`, or `NULL`.
#' @param var,time Tested factor column(s) (`time` in LRT mode).
#' @param batch Batch column name(s), or `NULL`.
#' @return `dds` with `colData(dds)[[pair]]` coerced to a factor.
#' @noRd
validate_pair <- function(dds, pair, var = NULL, batch = NULL, time = NULL) {
    if (is.null(pair)) return(dds)
    if (!is.character(pair) || length(pair) != 1L || is.na(pair) || !nzchar(pair)) {
        stop("`pair` must be `NULL` or a single column name in colData(dds).", call. = FALSE)
    }
    cd <- SummarizedExperiment::colData(dds)
    if (!pair %in% colnames(cd)) {
        stop("`pair` (\"", pair, "\") is not a column of colData(dds). Available columns: ",
             paste(colnames(cd), collapse = ", "), call. = FALSE)
    }
    if (pair %in% c(var, batch, time)) {
        stop("`pair` (\"", pair, "\") must not also be given in `var`, `batch` or `time`.",
             call. = FALSE)
    }
    p <- cd[[pair]]
    if (anyNA(p)) {
        stop("`pair` (\"", pair, "\") has missing values; every sample needs a ", pair, ".",
             call. = FALSE)
    }
    p <- droplevels(factor(p))
    if (nlevels(p) < 2L) {
        stop("`pair` (\"", pair, "\") must have at least 2 levels.", call. = FALSE)
    }
    levels_per_unit <- function(col) {
        vapply(split(as.character(cd[[col]]), p), function(x) length(unique(x)), integer(1))
    }
    for (f in c(var, time)) {
        n_lv <- levels_per_unit(f)
        if (all(n_lv <= 1L)) {
            stop("`pair` (\"", pair, "\") is nested within `", f, "`: every ", pair,
                 " has samples in a single level of ", f, ", so a paired design cannot estimate the ",
                 f, " effect. Use `pair` only for factors measured within each ", pair,
                 " (e.g. before/after treatment); otherwise leave `pair` unset.", call. = FALSE)
        }
        single <- names(n_lv)[n_lv <= 1L]
        if (length(single)) {
            warning(length(single), " ", pair, " level(s) have samples in only one level of `", f,
                    "` and add nothing to its within-", pair, " comparison: ",
                    paste(utils::head(single, 10L), collapse = ", "),
                    if (length(single) > 10L) ", ..." else "", call. = FALSE)
        }
    }
    for (b in batch) {
        if (all(levels_per_unit(b) == 1L)) {
            stop("`batch` (\"", b, "\") is constant within every `pair` (\"", pair,
                 "\") level, so `pair` already accounts for it. Drop \"", b, "\" from `batch`.",
                 call. = FALSE)
        }
    }
    SummarizedExperiment::colData(dds)[[pair]] <- p
    dds
}

#' Validate a custom colour palette against the number of group levels
#' @param dds A `DESeqDataSet` object.
#' @param var Grouping column name.
#' @param pals Character vector of colours, or `NULL`.
#' @return The validated `pals` (or `NULL`).
#' @noRd
validate_pals <- function(dds, var, pals) {
    if (is.null(pals)) return(NULL)
    n_levels <- nlevels(SummarizedExperiment::colData(dds)[[var]])
    if (!is.character(pals) || length(pals) < n_levels) {
        stop("`pals` must be a character vector of at least ", n_levels,
             " colours (one per level of `", var, "`); got ", length(pals), ".", call. = FALSE)
    }
    pals
}

#' Validate/generate the comparison vector or contrast-spec table
#'
#' If `comparisons` is `NULL`, auto-expands the contrast grid via
#' [expand_contrast_specs()] (pairwise + one-vs-rest for one factor; every
#' factor's main effect for two or three, plus interaction contrasts when
#' `include_interactions = TRUE`). If supplied as a character vector, each
#' entry must be either a simple `"A_vs_B"` string (single-factor) or a
#' contrast `label` from the auto grid (two-factor), and the returned object
#' is the filtered specs table (or the character vector for single-factor
#' backward compatibility).
#'
#' @param dds A `DESeqDataSet` object.
#' @param var Grouping column name(s), length 1 or 2.
#' @param comparisons Character vector of comparison / label strings, a
#'   contrast-spec table (or list of [contrast_spec()] rows), or `NULL`.
#' @param include_interactions Logical passed to [expand_contrast_specs()].
#'   Default `FALSE`.
#' @param stratified Logical passed to [expand_contrast_specs()]. Default
#'   `FALSE`.
#' @return Always a validated `contrast_specs` table, whatever the arity of
#'   `var` and whatever form `comparisons` took. (It used to return a bare
#'   character vector for single-factor `var`, which every caller then had to
#'   branch on.)
#' @noRd
validate_comparisons <- function(dds, var, comparisons = NULL,
                                 include_interactions = FALSE,
                                 stratified = FALSE) {
    specs <- expand_contrast_specs(dds, var,
                                   include_interactions = include_interactions,
                                   stratified = stratified)

    if (is.null(comparisons)) return(specs)

    # A table (or a list of rows) is the supported way to run contrasts the
    # automatic grid does not generate: validate it against `dds` and use it.
    if (!is.character(comparisons)) {
        user <- as_contrast_specs(comparisons, var = var[[1L]])
        if (is.null(user)) {
            stop("`comparisons` is empty.", call. = FALSE)
        }
        return(validate_contrast_specs(user, dds = dds, var = var))
    }
    if (length(comparisons) == 0) {
        stop("`comparisons` must be `NULL`, a non-empty character vector, or a ",
             "contrast-spec table.", call. = FALSE)
    }
    comparisons <- unique(comparisons)

    # Labels of the automatic grid can always be named.
    known <- comparisons %in% specs$label
    if (all(known)) {
        return(specs[match(comparisons, specs$label), , drop = FALSE])
    }
    # For a single grouping factor, any well-formed comparison of its levels is
    # accepted even when the grid does not list it -- that is how a
    # pooled-vs-pooled comparison such as "d0+d1_vs_d7+d14" is requested.
    if (length(var) == 1L) {
        levels_var <- levels(factor(SummarizedExperiment::colData(dds)[[var]]))
        for (comp in comparisons) validate_comparison(comp, levels_var)
        out <- as_contrast_specs(comparisons, var = var)
        return(validate_contrast_specs(out, dds = dds, var = var))
    }
    stop("`comparisons` not in the auto contrast grid for `var = c(",
         paste(sprintf('"%s"', var), collapse = ", "), ")`: ",
         paste(comparisons[!known], collapse = ", "),
         ".\nValid labels include e.g.: ",
         paste(utils::head(specs$label, 8), collapse = ", "),
         "...\nTo run a contrast the grid does not generate, build it with ",
         "contrast_spec() and pass the table (see `add_contrasts()`).",
         call. = FALSE)
}

#' Validate a single comparison string against known levels
#' @param comparison A single `"A_vs_B"` or `"A_vs_B+C"` string.
#' @param levels_var Character vector of valid levels.
#' @return Invisibly, `NULL`.
#' @noRd
validate_comparison <- function(comparison, levels_var) {
    if (!grepl("_vs_", comparison, fixed = TRUE)) {
        stop("Comparison \"", comparison, "\" is not in \"A_vs_B\" format.", call. = FALSE)
    }
    parts <- strsplit(comparison, "_vs_", fixed = TRUE)[[1]]
    if (length(parts) != 2) {
        stop("Comparison \"", comparison, "\" must contain exactly one \"_vs_\" separator.", call. = FALSE)
    }
    groups <- unlist(strsplit(parts, "+", fixed = TRUE))
    missing <- setdiff(groups, levels_var)
    if (length(missing) > 0) {
        stop("Comparison \"", comparison, "\" references level(s) not found in the data: ",
             paste(missing, collapse = ", "), ". Available levels: ", paste(levels_var, collapse = ", "), call. = FALSE)
    }
    invisible(NULL)
}

#' Validate requested MSigDB collections
#' @param collections Character vector of MSigDB collection aliases:
#'   `"HALLMARK"`, `"GOBP"`, `"GOCC"`, `"GOMF"`, `"KEGG"`,
#'   `"REACTOME"`, `"PID"`. `"KEGG_LEGACY"` is accepted as a deprecated
#'   alias for `"KEGG"`.
#' @return The validated `collections` vector (with deprecated aliases
#'   remapped).
#' @noRd
validate_msigdb_collections <- function(collections) {
    supported <- c("HALLMARK", "GOBP", "GOCC", "GOMF", "KEGG",
                   "REACTOME", "PID")
    if (!is.character(collections) || length(collections) == 0) {
        stop("`collections` must be a non-empty character vector, one or more of: ",
             paste(supported, collapse = ", "), ".", call. = FALSE)
    }
    if ("KEGG_LEGACY" %in% collections) {
        message('Collection alias "KEGG_LEGACY" is now "KEGG"; remapping.')
        collections[collections == "KEGG_LEGACY"] <- "KEGG"
        collections <- unique(collections)
    }
    unsupported <- setdiff(collections, supported)
    if (length(unsupported) > 0) {
        stop("Unsupported MSigDB collection(s): ", paste(unsupported, collapse = ", "),
             ". Supported: ", paste(supported, collapse = ", "), ".", call. = FALSE)
    }
    collections
}

#' Validate / normalise the `remove_vdj` argument
#'
#' Accepts `NULL` / `FALSE` (off), `TRUE` / `"both"` (BCR + TCR), or one or
#' more of `"BCR"` / `"TCR"`. Returns a character vector of kinds to remove
#' (`character(0)` when disabled).
#' @noRd
validate_remove_vdj <- function(remove_vdj) {
    if (is.null(remove_vdj) || isFALSE(remove_vdj)) {
        return(character(0))
    }
    if (isTRUE(remove_vdj)) {
        return(c("BCR", "TCR"))
    }
    if (is.character(remove_vdj)) {
        kinds <- toupper(trimws(remove_vdj))
        kinds <- kinds[nzchar(kinds)]
        if (length(kinds) == 1 && identical(kinds, "BOTH")) {
            return(c("BCR", "TCR"))
        }
        unsupported <- setdiff(kinds, c("BCR", "TCR"))
        if (length(unsupported) > 0) {
            stop("`remove_vdj` must be NULL/FALSE/TRUE/\"both\", or one or more of ",
                 "\"BCR\", \"TCR\". Unsupported value(s): ",
                 paste(unsupported, collapse = ", "), ".", call. = FALSE)
        }
        return(unique(kinds))
    }
    stop("`remove_vdj` must be NULL/FALSE/TRUE/\"both\", or one or more of \"BCR\", \"TCR\".",
         call. = FALSE)
}

#' Validate the output file format
#' @param format One of "csv" (currently the only supported format).
#' @return The validated, lower-cased `format` string.
#' @noRd
validate_output_format <- function(format) {
    validate_choice(format, c("csv"), "output_format")
}

#' Validate/create the save directory
#' @param save_dir Path to a directory (created if it does not exist).
#' @return The validated, normalized `save_dir` path.
#' @noRd
validate_save_dir <- function(save_dir) {
    if (!is.character(save_dir) || length(save_dir) != 1 || is.na(save_dir) || !nzchar(save_dir)) {
        stop("`save_dir` must be a single non-empty character string.", call. = FALSE)
    }
    if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    }
    normalizePath(save_dir, mustWork = TRUE)
}
