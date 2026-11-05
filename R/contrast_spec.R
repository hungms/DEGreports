# The single definition of what a contrast is.
#
# Every contrast in this package -- auto-generated or hand-written, Wald or
# LRT or read back from a precomputed DEG table -- is one row of a
# `contrast_specs` table. This file owns that schema: the field registry, the
# constructor, the validator, the string codecs the table's packed columns use,
# and the accessors that unpack them. Nothing else should parse a spec field.
#
# The table is a `data.frame` subclass, so everything that worked on the plain
# data frame these specs used to be -- `specs$label`, `specs[i, , drop =
# FALSE]`, `rbind()`, `knitr::kable()` -- still works unchanged.

# ---- field registry --------------------------------------------------------

#' Contrast classes and the fitting strategy each one implies
#'
#' `class` used to carry both "how is this fitted" and "which report section
#' does it belong to"; `fit_method` now names the former explicitly so the
#' dispatch in [run_diffexp_pipeline()] reads off a field instead of a chain
#' of `identical(class, ...)` tests.
#' @noRd
.contrast_fit_methods <- function() {
    c(
        single      = "levels",            # comparison of `var` levels, dds as given
        main        = "levels_additive",   # all samples, additive design rebuilt
        within_A    = "levels_subset",     # one stratum's samples only
        within_B    = "levels_subset",
        within_C    = "levels_subset",
        interaction = "coefs",             # combination of model coefficients
        lrt         = "lrt",               # full vs reduced likelihood ratio
        deg         = "precomputed"        # read from a supplied DEG table
    )
}

#' Fitting strategy for a contrast class
#'
#' Unknown classes (a `within_4` from a four-factor grid, say) fall back on
#' the `within_*` convention rather than erroring, so adding a class does not
#' require editing this map.
#' @param class Character vector of contrast classes.
#' @return Character vector of fit methods, same length as `class`.
#' @noRd
.fit_method_for_class <- function(class) {
    map <- .contrast_fit_methods()
    out <- unname(map[as.character(class)])
    unknown <- is.na(out)
    if (any(unknown)) {
        within <- grepl("^within(_|$)", as.character(class)[unknown])
        out[unknown][within] <- "levels_subset"
        out[unknown][!within] <- "levels"
    }
    out
}

# ---- string codecs ---------------------------------------------------------
#
# Four packed encodings live in a spec table. All four are read and written
# here and nowhere else:
#
#   coefficients   "a+b"            summed model coefficients
#   weights        "a=1;b=-0.5"     coefficients with fractional weights
#   terms          "B+C" / "b1+c1"  compound strata and multi-factor test_var
#   ix terms       "f:comp|f:comp"  which comparison each factor contributes
#
# `+` cannot occur in a DESeq2 coefficient name (they come from make.names()),
# and `validate_level_names()` keeps it out of factor levels, so `+` is safe
# as the separator for all of them.

#' Join values into a `+`-packed spec field
#' @noRd
.pack_terms <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(NA_character_)
    paste(x, collapse = "+")
}

#' Split a `+`-packed spec field back into a character vector
#'
#' `NA` and empty strings decode to `character(0)`, so a caller can treat
#' "no value" and "nothing left after unpacking" the same way.
#' @noRd
.unpack_terms <- function(x) {
    if (is.null(x)) return(character(0))
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(character(0))
    x <- unlist(strsplit(x, "+", fixed = TRUE), use.names = FALSE)
    x <- trimws(x)
    x[nzchar(x)]
}

#' Coefficient names from a `+`-packed spec field
#'
#' Kept as a separate name from [.unpack_terms()] because the call sites read
#' better, and because a future coefficient-name check belongs here.
#' @noRd
.clean_coef_names <- function(x) {
    .unpack_terms(x)
}

#' Encode a coefficient-weight map into one spec-table cell
#'
#' Fractional weights cannot be expressed by the `+`-joined
#' `comparison` / `coef_subtract` pair, which only ever sums coefficients, so
#' a pooled comparison stores `"name=weight;name=weight"` instead. Neither
#' `=` nor `;` can occur in a DESeq2 coefficient name, which comes from
#' [make.names()].
#' @noRd
.format_coef_weights <- function(w) {
    if (is.null(w) || !length(w)) return(NA_character_)
    num <- vapply(unname(w), format, character(1), digits = 12, trim = TRUE)
    paste0(names(w), "=", num, collapse = ";")
}

#' Decode [.format_coef_weights()] back into a named numeric vector
#' @return Named numeric vector, or `NULL` when there is nothing to decode.
#' @noRd
.parse_coef_weights <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(NULL)
    parts <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(NULL)
    kv <- strsplit(parts, "=", fixed = TRUE)
    bad <- vapply(kv, function(p) length(p) != 2L, logical(1))
    if (any(bad)) {
        stop("Malformed coefficient weight(s): ",
             paste(parts[bad], collapse = ", "),
             ". Expected \"name=weight\" pairs separated by \";\".", call. = FALSE)
    }
    val <- suppressWarnings(as.numeric(trimws(vapply(kv, `[[`, character(1), 2L))))
    nm <- trimws(vapply(kv, `[[`, character(1), 1L))
    if (anyNA(val)) {
        stop("Non-numeric coefficient weight(s) for: ",
             paste(nm[is.na(val)], collapse = ", "), ".", call. = FALSE)
    }
    stats::setNames(val, nm)
}

# ---- constructor / coercion ------------------------------------------------

#' Build one contrast specification
#'
#' The constructor for a single row of a `contrast_specs` table. Everyday
#' analysis never calls it -- [expand_contrast_specs()] builds the automatic
#' grid -- but it is the supported way to describe a contrast the automatic
#' grid does not generate, which can then be appended with
#' [add_contrasts()] and passed to any pipeline entry point as `comparisons`.
#'
#' Exactly one of the three ways of saying *what* is contrasted must be given:
#' \describe{
#'   \item{`comparison`}{An `"A_vs_B"` / `"A_vs_B+C"` string over the levels
#'     of `test_var`, fitted by [run_diffexp()] (`fit_method` `"levels"`,
#'     `"levels_additive"` or `"levels_subset"`).}
#'   \item{`coef` (+ optional `coef_subtract`)}{Model coefficient names as in
#'     `DESeq2::resultsNames()`, summed and optionally differenced.}
#'   \item{`coef_weights`}{A named numeric vector of coefficient weights, for
#'     a contrast that averages rather than sums -- a pooled one-vs-rest
#'     comparison, for instance.}
#' }
#'
#' @param label Identity of the contrast: the report tab name and the key in
#'   `by_comparison`. Must be unique within a spec table.
#' @param class Contrast class, which selects the report section:
#'   `"single"`, `"main"`, `"interaction"`, `"within_A"` / `"within_B"` /
#'   `"within_C"` (or `"within_<i>"`), `"lrt"`, `"deg"`. Default `"single"`.
#' @param type `"pairwise"`, `"one_vs_all"`, `"interaction"` or `"lrt"`.
#'   Default: inferred from `comparison` via [classify_comparison_type()], or
#'   `"interaction"` when coefficients are given.
#' @param test_var Column(s) in `colData(dds)` being contrasted.
#' @param comparison An `"A_vs_B"` / `"A_vs_B+C"` comparison string.
#' @param coef,coef_subtract Coefficient names summed in the numerator /
#'   denominator.
#' @param coef_weights Named numeric vector of coefficient weights (or its
#'   `"name=weight;..."` encoding).
#' @param stratum Named character vector of factors held fixed
#'   (`c(celltype = "CB")`), or `NULL`.
#' @param question One-line prose statement of what the contrast asks, shown
#'   in the contrast manifest and the report. Default: generated from the
#'   other fields by [contrast_question()].
#' @param path Output sub-directory for this contrast, relative to the
#'   module's `save_dir`. Default: a filesystem-safe form of `label`.
#' @param fit_method Fitting strategy override; normally derived from `class`.
#' @param ... Further fields (e.g. `full` / `reduced` / `test_terms` for an
#'   LRT row), recycled into the row as character columns.
#' @return A one-row `contrast_specs` table.
#' @seealso [add_contrasts()], [expand_contrast_specs()],
#'   [validate_contrast_specs()]
#' @examples
#' # A pooled-vs-pooled comparison the automatic grid never generates
#' contrast_spec(
#'   label = "earlyBoth_vs_lateBoth",
#'   test_var = "timepoint",
#'   comparison = "d0+d1_vs_d7+d14"
#' )
#' @export
contrast_spec <- function(label,
                          class = "single",
                          type = NULL,
                          test_var = NA_character_,
                          comparison = NA_character_,
                          coef = NULL,
                          coef_subtract = NULL,
                          coef_weights = NULL,
                          stratum = NULL,
                          question = NULL,
                          path = NULL,
                          fit_method = NULL,
                          ...) {

    if (!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label)) {
        stop("`label` must be a single non-empty character string.", call. = FALSE)
    }
    if (!is.character(class) || length(class) != 1L || is.na(class)) {
        stop("`class` must be a single character string.", call. = FALSE)
    }

    if (is.character(coef_weights)) coef_weights <- .parse_coef_weights(coef_weights)
    has_comparison <- !is.null(comparison) && length(comparison) == 1L &&
        !is.na(comparison) && nzchar(comparison)
    has_coef <- length(.clean_coef_names(coef)) > 0L
    has_weights <- length(coef_weights) > 0L
    n_given <- sum(has_comparison, has_coef, has_weights)
    if (n_given == 0L) {
        stop("Contrast \"", label, "\": give one of `comparison`, `coef` ",
             "(with optional `coef_subtract`), or `coef_weights`.", call. = FALSE)
    }
    if (n_given > 1L) {
        stop("Contrast \"", label, "\": give only one of `comparison`, `coef` ",
             "or `coef_weights`; got ", n_given, ".", call. = FALSE)
    }
    if (has_weights && (is.null(names(coef_weights)) || anyNA(coef_weights))) {
        stop("Contrast \"", label,
             "\": `coef_weights` must be a fully named, non-`NA` numeric vector.",
             call. = FALSE)
    }
    if (!is.null(stratum) &&
        (!is.character(stratum) || is.null(names(stratum)) || anyNA(stratum))) {
        stop("Contrast \"", label, "\": `stratum` must be a fully named, ",
             "non-`NA` character vector such as `c(celltype = \"CB\")`.",
             call. = FALSE)
    }

    if (is.null(type)) {
        type <- if (has_comparison) {
            classify_comparison_type(comparison)
        } else {
            "interaction"
        }
    }
    if (is.null(fit_method)) {
        fit_method <- if (has_comparison) .fit_method_for_class(class) else "coefs"
    }

    row <- data.frame(
        label = label,
        class = class,
        type = as.character(type)[[1L]],
        stratum_var = if (is.null(stratum)) NA_character_ else .pack_terms(names(stratum)),
        stratum_level = if (is.null(stratum)) NA_character_ else .pack_terms(unname(stratum)),
        test_var = .pack_terms(test_var),
        comparison = if (has_comparison) {
            as.character(comparison)[[1L]]
        } else if (has_coef) {
            .pack_terms(.clean_coef_names(coef))
        } else {
            NA_character_
        },
        coef_subtract = .pack_terms(.clean_coef_names(coef_subtract)),
        coef_weights = .format_coef_weights(coef_weights),
        stringsAsFactors = FALSE
    )
    extra <- list(...)
    for (nm in names(extra)) {
        row[[nm]] <- if (is.null(extra[[nm]])) NA_character_ else as.character(extra[[nm]])[[1L]]
    }
    row$path <- if (is.null(path)) contrast_path(label) else as.character(path)[[1L]]
    row$fit_method <- as.character(fit_method)[[1L]]
    row$question <- if (is.null(question)) NA_character_ else as.character(question)[[1L]]
    if (is.na(row$question)) row$question <- contrast_question(row)

    new_contrast_specs(row)
}

#' Low-level `contrast_specs` constructor
#'
#' Stamps the class onto a data frame of already-built rows without
#' re-validating, and fills in the derived fields (`path`, `fit_method`,
#' `question`) that the grid builders do not set themselves.
#' @param df A data frame with at least a `label` column.
#' @return A `contrast_specs` table.
#' @noRd
new_contrast_specs <- function(df) {
    if (is.null(df) || !nrow(df)) return(NULL)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!"label" %in% names(df)) {
        stop("A contrast-spec table needs a `label` column.", call. = FALSE)
    }
    if (!"class" %in% names(df)) df$class <- "single"
    # Fill the derived fields row by row, so a table that mixes hand-written
    # rows (which set them) with auto-generated ones (which do not) keeps the
    # explicit values and derives only the gaps.
    if (!"path" %in% names(df)) df$path <- NA_character_
    gap <- is.na(df$path) | !nzchar(df$path)
    if (any(gap)) df$path[gap] <- contrast_path(df$label[gap])
    if (!"fit_method" %in% names(df)) df$fit_method <- NA_character_
    gap <- is.na(df$fit_method) | !nzchar(df$fit_method)
    if (any(gap)) df$fit_method[gap] <- .fit_method_for_class(df$class[gap])
    if (!"question" %in% names(df)) df$question <- NA_character_
    gap <- is.na(df$question) | !nzchar(df$question)
    if (any(gap)) df$question[gap] <- contrast_question(df[gap, , drop = FALSE])
    rownames(df) <- NULL
    class(df) <- unique(c("contrast_specs", "data.frame"))
    df
}

#' Coerce comparisons to a contrast-spec table
#'
#' Accepts what a user can plausibly hand a pipeline entry point as
#' `comparisons`: an existing spec table (returned as-is after validation), a
#' character vector of comparison strings or labels, or a list of
#' [contrast_spec()] rows.
#'
#' @param x A `contrast_specs` table, data frame, character vector, or list
#'   of one-row spec tables.
#' @param var Grouping column name(s), used as `test_var` when `x` is a
#'   character vector. Default `NA`.
#' @param class Contrast class assigned to character input. Default
#'   `"single"`.
#' @return A `contrast_specs` table, or `NULL` for empty input.
#' @export
as_contrast_specs <- function(x, var = NA_character_, class = "single") {
    if (is.null(x) || (is.atomic(x) && !length(x))) return(NULL)
    if (inherits(x, "contrast_specs")) return(x)
    if (is.data.frame(x)) return(new_contrast_specs(x))
    if (is.list(x)) {
        parts <- lapply(x, as_contrast_specs, var = var, class = class)
        return(do.call(.rbind_specs, parts))
    }
    if (!is.character(x)) {
        stop("`comparisons` must be a character vector, a contrast-spec data ",
             "frame, or a list of contrast_spec() rows; got ",
             class(x)[[1L]], ".", call. = FALSE)
    }
    x <- unique(x[!is.na(x) & nzchar(x)])
    if (!length(x)) return(NULL)
    new_contrast_specs(data.frame(
        label = x,
        class = class,
        type = classify_comparison_type(x),
        stratum_var = NA_character_,
        stratum_level = NA_character_,
        test_var = .pack_terms(var),
        comparison = x,
        coef_subtract = NA_character_,
        stringsAsFactors = FALSE
    ))
}

#' Append contrasts to a grid
#'
#' The supported way to run something the automatic grid does not generate:
#' take the grid [expand_contrast_specs()] produced, add rows built with
#' [contrast_spec()] or [contrast_spec_from_coefs()], and hand the result to
#' any pipeline entry point as `comparisons`.
#'
#' @param specs A `contrast_specs` table (or `NULL` to start from nothing).
#' @param ... One-row spec tables, lists of them, or character comparison
#'   strings.
#' @param var Grouping column name(s) used for character input.
#' @return A `contrast_specs` table.
#' @examples
#' \dontrun{
#' grid <- expand_contrast_specs(dds, var = c("tissue", "genotype"))
#' grid <- add_contrasts(grid, contrast_spec(
#'   label = "BM_KO_vs_SP_WT",
#'   test_var = "group",
#'   comparison = "BM.KO_vs_SP.WT"
#' ))
#' run_deseq2(dds, var = c("tissue", "genotype"), comparisons = grid)
#' }
#' @seealso [contrast_spec()], [contrast_spec_from_coefs()]
#' @export
add_contrasts <- function(specs, ..., var = NA_character_) {
    extra <- lapply(list(...), as_contrast_specs, var = var)
    out <- do.call(.rbind_specs, c(list(as_contrast_specs(specs, var = var)), extra))
    validate_contrast_specs(out)
}

#' Describe a contrast as a weighted combination of model coefficients
#'
#' The most general way to ask a question of a fitted DESeq2 model: give the
#' coefficient names from `DESeq2::resultsNames(dds)` and the weight each
#' carries. Weights that are all `±1` are fitted as a list contrast, anything
#' else (an average over several levels, say) as a numeric contrast vector.
#'
#' @param label Contrast label.
#' @param weights Named numeric vector of coefficient weights, e.g.
#'   `c(tissueSP.genotypeKO = 1, tissueBM.genotypeKO = -1)`.
#' @param question Optional prose statement of what the contrast asks.
#' @param test_var Optional `colData` column(s) the contrast concerns, for the
#'   manifest.
#' @param ... Passed to [contrast_spec()].
#' @return A one-row `contrast_specs` table.
#' @examples
#' contrast_spec_from_coefs(
#'   label = "KO_effect_in_SP_vs_BM",
#'   weights = c(tissueSP.genotypeKO = 1),
#'   question = "Does the KO effect differ between SP and BM?"
#' )
#' @seealso [add_contrasts()], [run_diffexp_coef()]
#' @export
contrast_spec_from_coefs <- function(label, weights, question = NULL,
                                     test_var = NA_character_, ...) {
    if (is.character(weights)) weights <- .parse_coef_weights(weights)
    if (!is.numeric(weights) || is.null(names(weights)) || anyNA(weights)) {
        stop("`weights` must be a fully named, non-`NA` numeric vector of ",
             "coefficient weights.", call. = FALSE)
    }
    # ±1 throughout is expressible as a list contrast, which DESeq2 shrinks
    # more cheaply; anything else needs the weight vector.
    summable <- all(abs(abs(weights) - 1) < 1e-9)
    if (summable) {
        contrast_spec(
            label = label, class = "interaction", type = "interaction",
            test_var = test_var,
            coef = names(weights)[weights > 0],
            coef_subtract = names(weights)[weights < 0],
            question = question, ...
        )
    } else {
        contrast_spec(
            label = label, class = "interaction", type = "interaction",
            test_var = test_var, coef_weights = weights, question = question, ...
        )
    }
}

# ---- validation ------------------------------------------------------------

#' Validate a contrast-spec table
#'
#' Checks the structural invariants the fitting dispatch relies on: a unique,
#' non-empty `label` and `path` on every row, a known `fit_method`, exactly
#' one way of saying what is contrasted, and -- when `dds` is supplied -- that
#' every referenced `colData` column, factor level and model coefficient
#' actually exists. Errors name the offending row and field, since a spec
#' table is now something an analyst can write by hand.
#'
#' @param specs A `contrast_specs` table (or a data frame to coerce).
#' @param dds Optional `DESeqDataSet`, to check column / level / coefficient
#'   references. Default `NULL` (structure only).
#' @param var Optional grouping column name(s) the table should belong to.
#' @return The validated table, invisibly, as a `contrast_specs` object.
#' @export
validate_contrast_specs <- function(specs, dds = NULL, var = NULL) {
    specs <- as_contrast_specs(specs, var = if (is.null(var)) NA_character_ else var)
    if (is.null(specs)) {
        stop("No contrasts to run: the contrast-spec table is empty.", call. = FALSE)
    }

    bad_label <- is.na(specs$label) | !nzchar(specs$label)
    if (any(bad_label)) {
        stop("Contrast row(s) ", paste(which(bad_label), collapse = ", "),
             " have an empty `label`.", call. = FALSE)
    }
    dup <- duplicated(specs$label)
    if (any(dup)) {
        stop("Duplicated contrast label(s): ",
             paste(unique(specs$label[dup]), collapse = ", "),
             ". Every label keys one output folder and one report tab, so it ",
             "must be unique.", call. = FALSE)
    }
    dup_path <- duplicated(specs$path)
    if (any(dup_path)) {
        stop("Contrast label(s) ",
             paste(specs$label[specs$path %in% specs$path[dup_path]], collapse = ", "),
             " map to the same output path (\"",
             paste(unique(specs$path[dup_path]), collapse = "\", \""),
             "\"). Give them distinct `label`s or an explicit `path`.",
             call. = FALSE)
    }

    known <- c(unname(.contrast_fit_methods()), "levels_subset")
    unknown_fit <- !specs$fit_method %in% known
    if (any(unknown_fit)) {
        stop("Unknown `fit_method` on contrast(s) ",
             paste(specs$label[unknown_fit], collapse = ", "), ": ",
             paste(unique(specs$fit_method[unknown_fit]), collapse = ", "),
             ". Expected one of: ", paste(unique(known), collapse = ", "), ".",
             call. = FALSE)
    }

    for (i in seq_len(nrow(specs))) {
        row <- specs[i, , drop = FALSE]
        if (identical(row$fit_method[[1L]], "lrt")) next
        if (identical(row$fit_method[[1L]], "precomputed")) next
        coefs <- contrast_coefs(row)
        has_comparison <- identical(row$fit_method[[1L]], "levels") ||
            identical(row$fit_method[[1L]], "levels_additive") ||
            identical(row$fit_method[[1L]], "levels_subset")
        if (has_comparison) {
            if (is.na(row$comparison[[1L]]) || !nzchar(row$comparison[[1L]])) {
                stop("Contrast \"", row$label[[1L]], "\" (fit_method \"",
                     row$fit_method[[1L]], "\") has no `comparison` string.",
                     call. = FALSE)
            }
        } else if (!length(coefs$pos) && !length(coefs$weights)) {
            stop("Contrast \"", row$label[[1L]],
                 "\" has no numerator coefficient and no `coef_weights`.",
                 call. = FALSE)
        }
        if (length(coefs$weights) && (length(coefs$pos) || length(coefs$neg))) {
            stop("Contrast \"", row$label[[1L]], "\": give either `coef` / ",
                 "`coef_subtract` or `coef_weights`, not both.", call. = FALSE)
        }
        # Compare the packed fields directly: contrast_stratum() returns
        # nothing when they disagree, which would hide the inconsistency.
        n_sv <- length(.unpack_terms(row$stratum_var[[1L]]))
        n_sl <- length(.unpack_terms(row$stratum_level[[1L]]))
        if (n_sv != n_sl) {
            stop("Contrast \"", row$label[[1L]],
                 "\": `stratum_var` and `stratum_level` list different numbers ",
                 "of values (", n_sv, " vs ", n_sl, ").", call. = FALSE)
        }
    }

    if (!is.null(dds)) .validate_specs_against_dds(specs, dds)
    invisible(specs)
}

#' Check a spec table's column / level / coefficient references against `dds`
#' @noRd
.validate_specs_against_dds <- function(specs, dds) {
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    for (i in seq_len(nrow(specs))) {
        row <- specs[i, , drop = FALSE]
        lab <- row$label[[1L]]
        cols <- unique(c(.unpack_terms(row$test_var[[1L]]),
                         .unpack_terms(row$stratum_var[[1L]])))
        missing <- setdiff(cols, colnames(cd))
        if (length(missing)) {
            stop("Contrast \"", lab, "\" references colData column(s) not in ",
                 "`dds`: ", paste(missing, collapse = ", "), ".", call. = FALSE)
        }
        tv <- .unpack_terms(row$test_var[[1L]])
        if (length(tv) == 1L && !is.na(row$comparison[[1L]]) &&
            grepl("_vs_", row$comparison[[1L]], fixed = TRUE)) {
            validate_comparison(row$comparison[[1L]], levels(factor(cd[[tv]])))
        }
        stratum <- contrast_stratum(row)
        for (sv in names(stratum)) {
            lvls <- levels(factor(cd[[sv]]))
            if (!stratum[[sv]] %in% lvls) {
                stop("Contrast \"", lab, "\" holds ", sv, " = \"", stratum[[sv]],
                     "\", which is not a level of ", sv, ". Available: ",
                     paste(lvls, collapse = ", "), ".", call. = FALSE)
            }
        }
    }
    invisible(TRUE)
}

#' Reject factor levels that cannot survive label parsing or a file path
#'
#' A contrast label is built by pasting level names around `"_vs_"` and `"+"`,
#' and then used as an output directory, so a level containing either of those
#' -- or a path separator -- silently produces a mislabelled contrast or a
#' stray nested folder. Failing at grid-build time, naming the factor and the
#' level, is far more actionable.
#'
#' @param dds A `DESeqDataSet`.
#' @param var Character vector of `colData` column names to check.
#' @return Invisibly `TRUE` when every level is safe.
#' @export
validate_level_names <- function(dds, var) {
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    bad <- list()
    for (v in intersect(var, colnames(cd))) {
        lvls <- levels(factor(cd[[v]]))
        offenders <- lvls[grepl("[+/\\\\]|_vs_", lvls) | trimws(lvls) != lvls |
                              !nzchar(trimws(lvls))]
        if (length(offenders)) bad[[v]] <- offenders
    }
    if (!length(bad)) return(invisible(TRUE))
    detail <- paste(vapply(names(bad), function(v) {
        paste0("  ", v, ": ", paste0("\"", bad[[v]], "\"", collapse = ", "))
    }, character(1)), collapse = "\n")
    stop("Factor level name(s) cannot be used in a contrast label:\n", detail,
         "\n\nContrast labels are built by joining level names with \"_vs_\" and ",
         "\"+\", and are also used as output directory names, so a level ",
         "containing \"_vs_\", \"+\", \"/\", \"\\\\\" or surrounding whitespace ",
         "would produce a mislabelled contrast. Rename the level(s), e.g. ",
         "`colData(dds)$", names(bad)[[1L]], " <- factor(gsub(\"+\", \"_plus_\", ",
         "colData(dds)$", names(bad)[[1L]], ", fixed = TRUE))`.",
         call. = FALSE)
}

# ---- accessors -------------------------------------------------------------

#' Coefficients a contrast combines
#'
#' Unpacks the three coefficient fields of one spec row into the form the
#' fitting code wants, so no caller splits those strings itself.
#'
#' @param spec A one-row `contrast_specs` table.
#' @return A list with `pos` (summed in the numerator), `neg` (summed in the
#'   denominator) and `weights` (named numeric, for contrasts that average
#'   coefficients). A row that compares factor levels rather than coefficients
#'   gives three empty elements.
#' @export
contrast_coefs <- function(spec) {
    stopifnot(nrow(spec) == 1L)
    is_coefs <- identical(spec$fit_method[[1L]], "coefs") ||
        identical(spec$class[[1L]], "interaction")
    weights <- if ("coef_weights" %in% names(spec)) {
        .parse_coef_weights(spec$coef_weights[[1L]])
    } else {
        NULL
    }
    list(
        pos = if (is_coefs && is.null(weights)) .clean_coef_names(spec$comparison[[1L]]) else character(0),
        neg = if (is_coefs && is.null(weights)) .clean_coef_names(spec$coef_subtract[[1L]]) else character(0),
        weights = if (is.null(weights)) numeric(0) else weights
    )
}

#' Factors a contrast holds fixed
#'
#' @param spec A one-row `contrast_specs` table.
#' @return A named character vector (`c(celltype = "CB")`), empty when the
#'   contrast holds nothing fixed.
#' @export
contrast_stratum <- function(spec) {
    stopifnot(nrow(spec) == 1L)
    sv <- .unpack_terms(spec$stratum_var[[1L]])
    sl <- .unpack_terms(spec$stratum_level[[1L]])
    if (!length(sv) || length(sv) != length(sl)) {
        return(stats::setNames(character(0), character(0)))
    }
    stats::setNames(sl, sv)
}

#' Columns a contrast tests
#' @param spec A one-row `contrast_specs` table.
#' @return Character vector of `colData` column names.
#' @export
contrast_test_vars <- function(spec) {
    stopifnot(nrow(spec) == 1L)
    .unpack_terms(spec$test_var[[1L]])
}

#' Filesystem-safe output path for a contrast label
#'
#' `label` is the contrast's identity -- the report tab, the `by_comparison`
#' key -- and has historically doubled as its output directory, which couples
#' naming to what a filesystem will accept. `path` is that second job, split
#' out: nested labels (`"C_vs_noC/CB__noK"`) keep their nesting, everything
#' else is sanitised.
#'
#' @param label Character vector of contrast labels.
#' @return Character vector of relative paths.
#' @export
contrast_path <- function(label) {
    label <- as.character(label)
    parts <- strsplit(label, "/+", perl = TRUE)
    vapply(parts, function(p) {
        p <- p[nzchar(p)]
        p <- gsub("[[:space:]]+", "_", p)
        p <- gsub("[^A-Za-z0-9._+=-]", "_", p)
        if (!length(p)) return("contrast")
        paste(p, collapse = "/")
    }, character(1), USE.NAMES = FALSE)
}

# ---- prose -----------------------------------------------------------------

#' One-line statement of what a contrast asks
#'
#' The question an interaction contrast answers was already spelled out by
#' [interaction_question_table()]; this extends the same treatment to main,
#' within-stratum, single-factor and LRT rows, so every contrast in the
#' manifest and the report says what it tests in words.
#'
#' @param specs A `contrast_specs` table (or any data frame with the fields).
#' @return Character vector, one element per row.
#' @export
contrast_question <- function(specs) {
    if (is.null(specs) || !nrow(specs)) return(character(0))
    fld <- function(nm) {
        if (nm %in% names(specs)) as.character(specs[[nm]]) else rep(NA_character_, nrow(specs))
    }
    class_v <- fld("class")
    comp_v <- fld("comparison")
    tv_v <- fld("test_var")
    sv_v <- fld("stratum_var")
    sl_v <- fld("stratum_level")
    ix_v <- fld("ix_terms")
    terms_v <- fld("test_terms")

    phrase_sides <- function(comp) {
        if (is.na(comp) || !grepl("_vs_", comp, fixed = TRUE)) return(NA_character_)
        parts <- strsplit(comp, "_vs_", fixed = TRUE)[[1L]]
        test <- paste(.unpack_terms(parts[[1L]]), collapse = " + ")
        ref <- if (length(parts) > 1L) {
            paste(.unpack_terms(parts[[2L]]), collapse = " + ")
        } else {
            NA_character_
        }
        if (is.na(ref)) return(NA_character_)
        paste0(test, " versus ", ref)
    }

    vapply(seq_len(nrow(specs)), function(i) {
        cls <- class_v[[i]]
        stratum <- if (!is.na(sv_v[[i]]) && !is.na(sl_v[[i]])) {
            paste0(
                ", within ",
                paste(paste0(.unpack_terms(sv_v[[i]]), " = ",
                             .unpack_terms(sl_v[[i]])), collapse = " and ")
            )
        } else {
            ""
        }
        if (identical(cls, "lrt")) {
            return(paste0("Does dropping ",
                          if (is.na(terms_v[[i]])) "the tested term(s)" else terms_v[[i]],
                          " worsen the model fit?"))
        }
        if (identical(cls, "interaction")) {
            terms <- .unpack_ix_terms(ix_v[[i]])
            if (!length(terms)) {
                return(paste0("Does the interaction differ", stratum, "?"))
            }
            # "celltype:CC_vs_CB" -> factor "celltype", comparison "CC_vs_CB".
            fac <- names(terms)
            comp <- sub("^[^:]*:", "", unname(terms))
            # "CC vs CB" for the effect under test, "CC and CB" for the levels
            # it is compared across.
            as_vs <- function(x) gsub("_vs_", " vs ", gsub("+", " + ", x, fixed = TRUE),
                                      fixed = TRUE)
            as_and <- function(x) {
                sides <- strsplit(x, "_vs_", fixed = TRUE)[[1L]]
                paste(gsub("+", " + ", sides, fixed = TRUE), collapse = " and ")
            }
            effect <- paste0("the ", fac[[1L]], " ", as_vs(comp[[1L]]), " effect")
            if (length(terms) == 2L) {
                return(paste0("Does ", effect, " differ between ", fac[[2L]], " ",
                              as_and(comp[[2L]]), stratum, "?"))
            }
            # Three-way: the two-way interaction itself is what varies.
            inner <- paste0(effect, " across ", fac[[2L]], " ",
                            as_and(comp[[2L]]))
            return(paste0("Does the difference in ", inner, " itself differ between ",
                          fac[[3L]], " ", as_and(comp[[3L]]), stratum, "?"))
        }
        sides <- phrase_sides(comp_v[[i]])
        if (is.na(sides)) {
            return(paste0("Contrast ", comp_v[[i]], stratum, "."))
        }
        scope <- if (identical(cls, "main")) {
            paste0(", adjusted for the other factors")
        } else {
            stratum
        }
        paste0("Is ", tv_v[[i]], " ", sides, " differentially expressed", scope, "?")
    }, character(1))
}

#' Decode the `"factor:comparison|factor:comparison"` interaction-term field
#' @return Named character vector (names = factors), empty when absent.
#' @noRd
.unpack_ix_terms <- function(x) {
    if (is.null(x) || is.na(x) || !nzchar(x)) {
        return(stats::setNames(character(0), character(0)))
    }
    parts <- unlist(strsplit(as.character(x), "|", fixed = TRUE), use.names = FALSE)
    parts <- trimws(parts[nzchar(parts)])
    if (!length(parts)) return(stats::setNames(character(0), character(0)))
    nm <- sub(":.*$", "", parts)
    stats::setNames(parts, nm)
}

# ---- methods ---------------------------------------------------------------

#' @export
print.contrast_specs <- function(x, ...) {
    cat("<contrast_specs>: ", nrow(x), " contrast(s)\n", sep = "")
    by_class <- table(x$class)
    cat("  classes: ",
        paste(paste0(names(by_class), " (", as.integer(by_class), ")"),
              collapse = ", "), "\n", sep = "")
    show <- as.data.frame(x)[, intersect(c("label", "class", "fit_method",
                                           "test_var", "comparison"),
                                         names(x)), drop = FALSE]
    print(utils::head(show, 20L), row.names = FALSE)
    if (nrow(x) > 20L) cat("  ... ", nrow(x) - 20L, " more\n", sep = "")
    invisible(x)
}

#' @export
summary.contrast_specs <- function(object, ...) {
    data.frame(
        label = object$label,
        class = object$class,
        fit_method = object$fit_method,
        question = contrast_question(object),
        stringsAsFactors = FALSE
    )
}
