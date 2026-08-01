#' Resolve the DESeq2 design formula from `var` / `batch` (or an override)
#'
#' Everyday use passes only `var` (length 1 or 2) and optional `batch`;
#' `design = NULL` then expands to `~ var`, `~ batch + var`, `~ A * B`, or
#' `~ batch + A * B` as appropriate. An explicit `design` is validated and
#' returned unchanged.
#'
#' @param var Character vector of length 1 or 2 naming grouping columns in
#'   `colData(dds)`.
#' @param batch Optional character vector of covariate column names, or `NULL`.
#' @param design Optional design formula / character string override, or
#'   `NULL` (default) to auto-build from `var` / `batch`.
#' @return A `formula`.
#' @export
resolve_design <- function(var, batch = NULL, design = NULL) {
    if (!is.null(design)) {
        if (is.character(design)) design <- stats::as.formula(design)
        if (!inherits(design, "formula")) {
            stop("`design` must be a formula, a character string coercible to one, or `NULL`.",
                 call. = FALSE)
        }
        return(design)
    }
    if (!is.character(var) || length(var) < 1 || length(var) > 2 || anyNA(var)) {
        stop("`var` must be a character vector of length 1 or 2.", call. = FALSE)
    }
    if (!is.null(batch)) {
        if (!is.character(batch) || length(batch) < 1 || anyNA(batch)) {
            stop("`batch` must be `NULL` or a non-empty character vector of column names.",
                 call. = FALSE)
        }
        if (length(intersect(batch, var)) > 0) {
            stop("`batch` columns must not overlap with `var`.", call. = FALSE)
        }
    }
    var_term <- if (length(var) == 1) var else paste(var, collapse = " * ")
    rhs <- if (is.null(batch)) {
        var_term
    } else {
        paste(c(batch, var_term), collapse = " + ")
    }
    stats::as.formula(paste("~", rhs), env = baseenv())
}

#' Expand the automatic contrast grid for one- or two-factor `var`
#'
#' Returns a data frame of contrast specifications consumed by
#' [run_diffexp_pipeline()] and the HTML report. Columns:
#' \describe{
#'   \item{label}{Disk / tab name for the contrast}
#'   \item{class}{`"single"`, `"within_A"`, `"within_B"`, or `"interaction"`}
#'   \item{type}{`"pairwise"`, `"one_vs_all"`, or `"interaction"`}
#'   \item{stratum_var}{Factor held fixed (two-factor within-* only), else `NA`}
#'   \item{stratum_level}{Level of `stratum_var`, else `NA`}
#'   \item{test_var}{Column used for the contrast (or first of `var` for interaction)}
#'   \item{comparison}{`"A_vs_B"` / `"A_vs_B+C"` string for within / single;
#'     interaction coefficient name for interaction rows}
#' }
#'
#' @param dds A `DESeqDataSet`.
#' @param var Character vector of length 1 or 2.
#' @return A data frame of contrast specs (0 rows only if somehow no levels).
#' @export
expand_contrast_specs <- function(dds, var) {
    dds <- validate_dds(dds)
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    if (!is.character(var) || length(var) < 1 || length(var) > 2 || anyNA(var)) {
        stop("`var` must be a character vector of length 1 or 2.", call. = FALSE)
    }
    missing <- setdiff(var, colnames(cd))
    if (length(missing) > 0) {
        stop("`var` column(s) not in colData(dds): ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    if (length(var) == 1) {
        levels_var <- levels(factor(cd[[var]]))
        comps <- generate_comparisons(levels_var)
        types <- classify_comparison_type(comps)
        return(data.frame(
            label = comps,
            class = "single",
            type = types,
            stratum_var = NA_character_,
            stratum_level = NA_character_,
            test_var = var,
            comparison = comps,
            stringsAsFactors = FALSE
        ))
    }

    a <- var[1]
    b <- var[2]
    a_levels <- levels(factor(cd[[a]]))
    b_levels <- levels(factor(cd[[b]]))
    rows <- list()

    # Within each A level: B pairwise (+ B one-vs-rest if |B| >= 3)
    for (al in a_levels) {
        b_comps <- generate_comparisons(b_levels)
        b_types <- classify_comparison_type(b_comps)
        for (i in seq_along(b_comps)) {
            rows[[length(rows) + 1]] <- data.frame(
                label = paste0(al, "__", b_comps[i]),
                class = "within_A",
                type = b_types[i],
                stratum_var = a,
                stratum_level = al,
                test_var = b,
                comparison = b_comps[i],
                stringsAsFactors = FALSE
            )
        }
    }

    # Within each B level: A pairwise (+ A one-vs-rest if |A| >= 3)
    for (bl in b_levels) {
        a_comps <- generate_comparisons(a_levels)
        a_types <- classify_comparison_type(a_comps)
        for (i in seq_along(a_comps)) {
            rows[[length(rows) + 1]] <- data.frame(
                label = paste0(bl, "__", a_comps[i]),
                class = "within_B",
                type = a_types[i],
                stratum_var = b,
                stratum_level = bl,
                test_var = a,
                comparison = a_comps[i],
                stringsAsFactors = FALSE
            )
        }
    }

    # Interaction coefficients from ~ A * B (ref = first factor level of each)
    # DESeq2 names: "AlevelB.BlevelB" for the interaction of non-ref levels.
    a_ref <- a_levels[1]
    b_ref <- b_levels[1]
    for (al in setdiff(a_levels, a_ref)) {
        for (bl in setdiff(b_levels, b_ref)) {
            coef <- paste0(a, al, ".", b, bl)
            rows[[length(rows) + 1]] <- data.frame(
                label = paste0("interaction__", a, al, "_x_", b, bl),
                class = "interaction",
                type = "interaction",
                stratum_var = NA_character_,
                stratum_level = NA_character_,
                test_var = a,
                comparison = coef,
                stringsAsFactors = FALSE
            )
        }
    }

    do.call(rbind, rows)
}

#' Human-readable section title for a contrast class
#'
#' @param class One of `"single"`, `"within_A"`, `"within_B"`, `"interaction"`.
#' @param var Character vector of length 1 or 2.
#' @return A character string.
#' @export
contrast_class_title <- function(class, var) {
    if (identical(class, "single") || identical(class, "pairwise") ||
        identical(class, "one_vs_all")) {
        return(NA_character_)
    }
    if (identical(class, "interaction")) return("Interaction")
    if (length(var) < 2) return(as.character(class))
    a <- var[1]
    b <- var[2]
    if (identical(class, "within_A")) {
        paste0("Within ", a, " \u2014 ", b, " effects")
    } else if (identical(class, "within_B")) {
        paste0("Within ", b, " \u2014 ", a, " effects")
    } else {
        as.character(class)
    }
}
