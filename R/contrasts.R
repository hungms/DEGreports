#' Resolve the DESeq2 design formula from `var` / `batch` (or an override)
#'
#' Everyday use passes only `var` (length 1 or 2) and optional `batch`;
#' `design = NULL` then expands to `~ var`, `~ batch + var`, and for two
#' factors either `~ A + B` / `~ batch + A + B` (additive; default) or
#' `~ A * B` / `~ batch + A * B` when `include_interactions = TRUE`. An
#' explicit `design` is validated and returned unchanged.
#'
#' @param var Character vector of length 1 or 2 naming grouping columns in
#'   `colData(dds)`.
#' @param batch Optional character vector of covariate column names, or `NULL`.
#' @param design Optional design formula / character string override, or
#'   `NULL` (default) to auto-build from `var` / `batch`.
#' @param include_interactions Logical. When `var` has length 2 and
#'   `design` is `NULL`, use `A * B` if `TRUE` and additive `A + B` if
#'   `FALSE` (default). Ignored for single-factor `var` or when `design`
#'   is supplied.
#' @return A `formula`.
#' @export
resolve_design <- function(var, batch = NULL, design = NULL,
                           include_interactions = FALSE) {
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
    if (length(var) == 1) {
        var_term <- var
    } else if (isTRUE(include_interactions)) {
        var_term <- paste(var, collapse = " * ")
    } else {
        var_term <- paste(var, collapse = " + ")
    }
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
#'     primary interaction coefficient name for interaction rows}
#'   \item{coef_subtract}{For pairwise interaction contrasts, the second
#'     coefficient subtracted from `comparison`; otherwise `NA`}
#' }
#'
#' @param dds A `DESeqDataSet`.
#' @param var Character vector of length 1 or 2.
#' @param include_interactions Logical. When `var` has length 2, include
#'   interaction contrasts if `TRUE`. Default `FALSE`. Interaction rows
#'   cover (i) each non-reference level of `var[1]` × each non-reference
#'   level of `var[2]` vs the reference interaction, and (ii) all pairwise
#'   differences among those non-reference interactions (e.g. CK vs C × CC).
#' @return A data frame of contrast specs (0 rows only if somehow no levels).
#' @export
expand_contrast_specs <- function(dds, var, include_interactions = FALSE) {
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
            coef_subtract = rep(NA_character_, length(comps)),
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
                coef_subtract = NA_character_,
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
                coef_subtract = NA_character_,
                stringsAsFactors = FALSE
            )
        }
    }

    if (isTRUE(include_interactions)) {
        # Interaction coefficients from ~ A * B (ref = first factor level of each)
        # DESeq2 names: "AlevelB.BlevelB" for the interaction of non-ref levels.
        a_ref <- a_levels[1]
        b_ref <- b_levels[1]
        a_nonref <- setdiff(a_levels, a_ref)
        b_nonref <- setdiff(b_levels, b_ref)
        .ix_coef <- function(al, bl) paste0(a, al, ".", b, bl)

        for (al in a_nonref) {
            for (bl in b_nonref) {
                coef <- .ix_coef(al, bl)
                rows[[length(rows) + 1]] <- data.frame(
                    label = paste0("interaction__", a, al, "_x_", b, bl),
                    class = "interaction",
                    type = "interaction",
                    stratum_var = NA_character_,
                    stratum_level = NA_character_,
                    test_var = a,
                    comparison = coef,
                    coef_subtract = NA_character_,
                    stringsAsFactors = FALSE
                )
            }
        }

        # Pairwise among non-ref A interactions at each non-ref B
        # e.g. (CK:CC) - (C:CC) = CK vs C differing by cell type
        # Order as later vs earlier non-ref (e.g. CK_vs_C), matching
        # combn column order with rows swapped for the contrast sign.
        if (length(a_nonref) >= 2L && length(b_nonref) >= 1L) {
            a_pairs <- utils::combn(a_nonref, 2L)
            for (bl in b_nonref) {
                for (j in seq_len(ncol(a_pairs))) {
                    a_hi <- a_pairs[2L, j]
                    a_lo <- a_pairs[1L, j]
                    rows[[length(rows) + 1]] <- data.frame(
                        label = paste0("interaction__", a, a_hi, "_vs_", a, a_lo, "_x_", b, bl),
                        class = "interaction",
                        type = "interaction",
                        stratum_var = NA_character_,
                        stratum_level = NA_character_,
                        test_var = a,
                        comparison = .ix_coef(a_hi, bl),
                        coef_subtract = .ix_coef(a_lo, bl),
                        stringsAsFactors = FALSE
                    )
                }
            }
        }

        # Pairwise among non-ref B interactions at each non-ref A
        if (length(b_nonref) >= 2L && length(a_nonref) >= 1L) {
            b_pairs <- utils::combn(b_nonref, 2L)
            for (al in a_nonref) {
                for (j in seq_len(ncol(b_pairs))) {
                    b_hi <- b_pairs[2L, j]
                    b_lo <- b_pairs[1L, j]
                    rows[[length(rows) + 1]] <- data.frame(
                        label = paste0("interaction__", a, al, "_x_", b, b_hi, "_vs_", b, b_lo),
                        class = "interaction",
                        type = "interaction",
                        stratum_var = NA_character_,
                        stratum_level = NA_character_,
                        test_var = a,
                        comparison = .ix_coef(al, b_hi),
                        coef_subtract = .ix_coef(al, b_lo),
                        stringsAsFactors = FALSE
                    )
                }
            }
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
