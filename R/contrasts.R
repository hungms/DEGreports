#' Validate a grouping-factor vector
#'
#' The one place the `var` contract is enforced. There is no upper limit on
#' the number of factors any more -- the grid builders are written for an
#' arbitrary number -- but every extra factor multiplies the contrast grid and
#' needs every level combination to be populated, so three is still as far as
#' most designs usefully go.
#' @noRd
.validate_var_vector <- function(var) {
    if (!is.character(var) || length(var) < 1L || anyNA(var) || !all(nzchar(var))) {
        stop("`var` must be a non-empty character vector of colData column names.",
             call. = FALSE)
    }
    if (anyDuplicated(var)) {
        stop("`var` must not contain duplicated column names: ",
             paste(unique(var[duplicated(var)]), collapse = ", "), ".", call. = FALSE)
    }
    var
}

#' Resolve the DESeq2 design formula from `var` / `batch` (or an override)
#'
#' Everyday use passes only `var` (length 1–3) and optional `batch`;
#' `design = NULL` then expands to `~ var`, `~ batch + var`, and for two
#' or three factors either an additive sum (default) or a full interaction
#' product (`A * B` / `A * B * C`) when `include_interactions = TRUE`. An
#' explicit `design` is validated and returned unchanged.
#'
#' @param var Character vector of length 1, 2, or 3 naming grouping columns
#'   in `colData(dds)`.
#' @param batch Optional character vector of covariate column names, or `NULL`.
#' @param design Optional design formula / character string override, or
#'   `NULL` (default) to auto-build from `var` / `batch`.
#' @param include_interactions Logical. When `var` has length ≥ 2 and
#'   `design` is `NULL`, use `A * B` / `A * B * C` if `TRUE` and additive
#'   `A + B` / `A + B + C` if `FALSE` (default). Ignored for single-factor
#'   `var` or when `design` is supplied.
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
    var <- .validate_var_vector(var)
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

#' Expand the automatic contrast grid for one-, two-, or three-factor `var`
#'
#' Returns a data frame of contrast specifications consumed by
#' [run_diffexp_pipeline()] and the HTML report. Columns:
#' \describe{
#'   \item{label}{Disk / tab name for the contrast}
#'   \item{class}{`"single"` (one-factor `var`), `"main"` (a factor's main
#'     effect across every level of the others), `"interaction"` (a
#'     difference of differences), or `"within_A"` / `"within_B"` /
#'     `"within_C"` (a factor inside one fixed level combination of the
#'     others, `stratified = TRUE` only)}
#'   \item{type}{`"pairwise"`, `"one_vs_all"`, or `"interaction"`}
#'   \item{stratum_var}{Factor(s) held fixed, else `NA`}
#'   \item{stratum_level}{Level(s) of `stratum_var`, else `NA`}
#'   \item{test_var}{Column(s) contrasted; `+`-joined for interaction rows}
#'   \item{comparison}{`"A_vs_B"` / `"A_vs_B+C"` string for main / within /
#'     single; `+`-joined coefficient names summed in the numerator for
#'     interaction rows, or `NA` for an interaction that needs
#'     `coef_weights`}
#'   \item{coef_subtract}{`+`-joined coefficient names summed in the
#'     denominator of an interaction contrast; otherwise `NA`}
#'   \item{coef_weights}{`"name=weight;..."` map for an interaction contrast
#'     that averages rather than sums coefficients (a pooled one-vs-rest
#'     comparison); `NA` otherwise. Interaction rows only}
#'   \item{ix_terms}{`"factor:comparison|factor:comparison"` for an
#'     interaction row, which [interaction_question_table()] turns into the
#'     question the contrast answers; `NA` otherwise}
#' }
#'
#' @section Which contrasts are generated:
#' A single-factor `var` gives every pairwise comparison plus, with three or
#' more levels, every one-vs-rest comparison.
#'
#' For two or three factors, **every factor's main effect is always
#' expanded**: its pairwise (and one-vs-rest) comparisons taken across all
#' levels of the other factors, e.g. `v1l1_vs_v1l2+v1l3` regardless of `v2`
#' and `v3`, and `v2l1_vs_v2l2` regardless of `v1` and `v3`. These are the
#' `"main"` rows, and [run_diffexp_pipeline()] fits them on all samples with
#' an additive design (`~ [batch +] other factors + tested factor`), so the
#' other factors are adjusted for rather than left in the residual. Labels
#' are flat (`"v2l1_vs_v2l2"`), as for a single-factor `var`.
#'
#' `include_interactions = TRUE` adds **difference-of-differences** rows,
#' each asking whether one factor's effect differs between two levels of
#' another: `C_vs_noC_x_K_vs_noK_in_CB` is "does the K effect differ between
#' C and noC, within CB" — which is the same question as "does the C effect
#' differ between K and noK, within CB", so a factor pair contributes one row
#' per stratum, not two. With three factors each pair is asked inside every
#' level of the third factor, plus one three-way row per level pairing
#' (`CC_vs_CB_x_C_vs_noC_x_K_vs_noK`: does the C × K interaction itself
#' differ between CC and CB). [run_diffexp_pipeline()] fits all of them from
#' a single fit of the full product design (`~ [batch +] A * B * C`) on all
#' samples, each contrast being a combination of that fit's coefficients.
#'
#' A factor with three or more levels is compared **one-vs-rest** in these
#' rows (`CA_vs_CB+CC_x_K_vs_noK_in_C`) and not pairwise, matching how its
#' main effects read and keeping the grid linear in the number of levels. The
#' pooled side is an *average*, each of its levels weighted `1/(k - 1)`
#' irrespective of how many samples it has, so the contrast is fitted from a
#' numeric weight vector over `resultsNames()` rather than a summed list
#' contrast.
#'
#' @param dds A `DESeqDataSet`.
#' @param var Character vector of length 1, 2, or 3.
#' @param include_interactions Logical. When `var` has length ≥ 2, add the
#'   difference-of-differences interaction contrasts to the main effects.
#'   Default `FALSE`, i.e. main effects only. Ignored for single-factor `var`.
#' @param stratified Logical. Add the within-stratum contrasts — each factor
#'   compared inside one fixed level combination of the others, fitted on
#'   that subset alone. Default `FALSE`; not exposed by the pipeline entry
#'   points, since an interaction contrast answers the same question with all
#'   samples behind it.
#' @return A data frame of contrast specs (0 rows only if somehow no levels).
#' @export
expand_contrast_specs <- function(dds, var, include_interactions = FALSE,
                                  stratified = FALSE) {
    dds <- validate_dds(dds)
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    var <- .validate_var_vector(var)
    missing <- setdiff(var, colnames(cd))
    if (length(missing) > 0) {
        stop("`var` column(s) not in colData(dds): ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    validate_level_names(dds, var)

    if (length(var) == 1) {
        levels_var <- levels(factor(cd[[var]]))
        comps <- generate_comparisons(levels_var)
        types <- classify_comparison_type(comps)
        return(new_contrast_specs(data.frame(
            label = comps,
            class = "single",
            type = types,
            stratum_var = NA_character_,
            stratum_level = NA_character_,
            test_var = var,
            comparison = comps,
            coef_subtract = rep(NA_character_, length(comps)),
            stringsAsFactors = FALSE
        )))
    }

    # Two factors name their comparisons alphabetically, three or more follow
    # factor level order. Kept as it was so existing labels do not move.
    main <- .main_effect_specs(
        cd, var,
        comps_fun = if (length(var) >= 3L) .comps_factor_order else generate_comparisons
    )
    parts <- list(main)
    if (isTRUE(stratified)) {
        parts <- c(parts, list(.within_stratum_specs(cd, var)))
    }
    if (isTRUE(include_interactions)) {
        parts <- c(parts, list(.interaction_specs(cd, var)))
    }
    do.call(.rbind_specs, parts)
}

#' Level comparisons of one factor, as weights on its levels
#'
#' Two levels give the single pairwise comparison, non-reference vs reference
#' ([.comps_factor_order()]'s direction convention). Three or more give one
#' *pooled* comparison per level -- that level against the mean of the rest,
#' each of the rest weighted `1/(k - 1)` so an unbalanced design cannot tilt
#' the comparison towards its larger groups. The pairwise combinations are
#' deliberately not enumerated: pooling asks the same question of a
#' many-levelled factor without multiplying the interaction grid.
#'
#' @return List of `list(levels, weights, label)`, `levels[[1]]` being the
#'   level under test and always the one weighted `+1`.
#' @noRd
.level_comparisons <- function(lvls) {
    lvls <- as.character(lvls)
    k <- length(lvls)
    if (k < 2L) return(list())
    if (k == 2L) {
        return(list(list(
            levels = c(lvls[2L], lvls[1L]),
            weights = c(1, -1),
            label = paste0(lvls[2L], "_vs_", lvls[1L])
        )))
    }
    lapply(seq_len(k), function(i) {
        rest <- lvls[-i]
        list(
            levels = c(lvls[i], rest),
            weights = c(1, rep(-1 / length(rest), length(rest))),
            label = paste0(lvls[i], "_vs_", paste(rest, collapse = "+"))
        )
    })
}

#' Name of one DESeq2 interaction coefficient, or `NULL` if it is zero
#'
#' DESeq2 names an interaction term by pasting `<column><level>` for each
#' participating factor, joined by `.`, in the order the factors appear in the
#' design formula -- e.g. `celltypeCC.C_mutationC.K_mutationK`. A term naming
#' any factor's *reference* level is not a coefficient at all: under treatment
#' coding it is structurally zero, so it drops out of the contrast.
#'
#' @param at Named character vector, one level per participating factor.
#' @param refs Named character vector of each factor's reference level.
#' @noRd
.ix_coef_name <- function(var, at, refs) {
    if (any(at == refs[names(at)])) return(NULL)
    ord <- order(match(names(at), var))
    paste(paste0(names(at)[ord], at[ord]), collapse = ".")
}

#' Every subset of the held-fixed factors, the empty one first
#'
#' A contrast fitted inside one fixed level needs the coefficient for that
#' level folded in as well; inside two, it needs both of them and their
#' product. Enumerating the subsets is what makes [.ix_weighted_terms()]
#' independent of how many factors are held fixed.
#' @param stratum Named character vector (factor = level), or `NULL`.
#' @return List of named character vectors, starting with `NULL`.
#' @noRd
.ix_stratum_subsets <- function(stratum) {
    if (is.null(stratum) || !length(stratum)) return(list(NULL))
    out <- list(NULL)
    for (k in seq_along(stratum)) {
        out <- c(out, utils::combn(seq_along(stratum), k,
                                   function(s) stratum[s], simplify = FALSE))
    }
    out
}

#' Weighted coefficient terms of one interaction contrast
#'
#' A difference of differences expands into one term per corner of the tested
#' levels, weighted by the product of that corner's level weights: `±1`
#' throughout for pairwise comparisons, and `∓1/(k - 1)` on each pooled level
#' of a [.level_comparisons()] one-vs-rest. With a factor held fixed at
#' `stratum`, each corner contributes both its own interaction term and the
#' higher-order term that includes the held-fixed level -- the latter
#' vanishing when that level is the reference, which is why the interaction in
#' the reference stratum is a single coefficient and the one in another
#' stratum is a sum.
#'
#' Weights are accumulated per coefficient rather than appended, because a
#' pooled comparison can reach the same coefficient from several corners: the
#' three-way `CA_vs_CB+CC_x_C_vs_noC_x_K_vs_noK`, for instance, picks up
#' `C_mutationC.K_mutationK` once at `+1` and twice at `-1/2`, so it cancels
#' exactly and must drop out.
#'
#' @param tested List of `list(v, levels, weights, label)`, length 2 (a
#'   two-way interaction) or 3 (a three-way one).
#' @param stratum Named length-1 character vector (factor = level), or `NULL`.
#' @return Named numeric vector of coefficient weights, in the order the
#'   corners were walked, with exact cancellations removed.
#' @noRd
.ix_weighted_terms <- function(var, refs, tested, stratum = NULL) {
    vs <- vapply(tested, `[[`, character(1), "v")
    corners <- expand.grid(lapply(tested, function(t) seq_along(t$levels)),
                           KEEP.OUT.ATTRS = FALSE)
    nm <- character(0)
    wt <- numeric(0)
    for (r in seq_len(nrow(corners))) {
        pick <- as.integer(corners[r, ])
        w <- prod(vapply(seq_along(tested), function(i) {
            tested[[i]]$weights[[pick[[i]]]]
        }, numeric(1)))
        corner <- stats::setNames(
            vapply(seq_along(tested), function(i) {
                tested[[i]]$levels[[pick[[i]]]]
            }, character(1)),
            vs
        )
        # With factors held fixed, each corner contributes its own interaction
        # term plus every higher-order term that adds some of the held-fixed
        # levels -- all subsets of the stratum, the empty one included. Terms
        # naming a reference level drop out inside .ix_coef_name().
        coefs <- unlist(
            lapply(.ix_stratum_subsets(stratum), function(s) {
                .ix_coef_name(var, c(corner, s), refs)
            }),
            use.names = FALSE
        )
        if (!length(coefs)) next
        nm <- c(nm, coefs)
        wt <- c(wt, rep(w, length(coefs)))
    }
    if (!length(nm)) return(stats::setNames(numeric(0), character(0)))
    out <- vapply(unique(nm), function(x) sum(wt[nm == x]), numeric(1))
    out[abs(out) > 1e-9]
}

# `.format_coef_weights()` / `.parse_coef_weights()` used to live here; they
# are now in R/contrast_spec.R with the other spec-field codecs.

#' Interaction ("difference of differences") contrast rows
#'
#' For every pair of factors, asks whether one factor's effect differs
#' between two levels of the other -- a question that is symmetric in the
#' pair, so each pair contributes one row per level combination rather than
#' two. With three factors each pair is additionally asked inside every level
#' of the third factor, and the three-way term ("does that interaction itself
#' differ between two levels of the third factor") is added on top.
#'
#' A factor with three or more levels is compared one-vs-rest rather than
#' pairwise ([.level_comparisons()]), so its part of the grid grows with the
#' number of levels instead of with their pairs.
#'
#' Every row is fitted from a single fit of the full product design by
#' [run_diffexp_pipeline()]. A comparison of two levels sums coefficients, so
#' the row stores them as `+`-joined names in `comparison` / `coef_subtract`;
#' a pooled one-vs-rest comparison has to *average* them instead, so the row
#' stores fractional weights in `coef_weights` and leaves the other two `NA`
#' rather than encoding a contrast it cannot express.
#' @noRd
.interaction_specs <- function(cd, var) {
    if (length(var) < 2L) return(NULL)
    lv <- lapply(var, function(v) levels(factor(cd[[v]])))
    names(lv) <- var
    refs <- vapply(lv, `[[`, character(1), 1L)
    comps_of <- lapply(lv, .level_comparisons)
    names(comps_of) <- var

    rows <- list()
    # "celltypeCA_vs_celltypeCB+celltypeCC", for the rare grid whose factors
    # share level names.
    comp_qual <- function(t) {
        paste0(t$v, t$levels[[1L]], "_vs_",
               paste0(t$v, t$levels[-1L], collapse = "+"))
    }
    add <- function(tested, stratum) {
        w <- .ix_weighted_terms(var, refs, tested, stratum)
        if (!length(w)) return(invisible(NULL))
        tested_vars <- vapply(tested, `[[`, character(1), "v")
        comps <- vapply(tested, `[[`, character(1), "label")
        label <- paste(comps, collapse = "_x_")
        label_qual <- paste(vapply(tested, comp_qual, character(1)), collapse = "_x_")
        if (!is.null(stratum)) {
            label <- paste0(label, "_in_", paste(unname(stratum), collapse = "_"))
            label_qual <- paste0(label_qual, "_in_",
                                 paste0(names(stratum), unname(stratum),
                                        collapse = "_"))
        }
        # Two factors sharing level names would make "no_vs_yes_x_no_vs_yes"
        # of it; spell the columns out instead.
        if (anyDuplicated(comps)) label <- label_qual
        summable <- all(abs(abs(w) - 1) < 1e-9)
        rows[[length(rows) + 1L]] <<- data.frame(
            label = label,
            label_qualified = label_qual,
            class = "interaction",
            type = "interaction",
            stratum_var = .pack_terms(names(stratum)),
            stratum_level = .pack_terms(unname(stratum)),
            test_var = .pack_terms(tested_vars),
            comparison = if (summable) {
                paste(names(w)[w > 0], collapse = "+")
            } else {
                NA_character_
            },
            coef_subtract = if (summable && any(w < 0)) {
                paste(names(w)[w < 0], collapse = "+")
            } else {
                NA_character_
            },
            coef_weights = if (summable) NA_character_ else .format_coef_weights(w),
            ix_terms = paste(paste0(tested_vars, ":", comps), collapse = "|"),
            stringsAsFactors = FALSE
        )
    }
    # Each tested factor's comparisons, tagged with the column they came from,
    # so both the coefficient names and the label follow `var` order.
    tested_of <- function(v, comp) c(list(v = v), comp)

    # Every subset of two or more factors interacts, asked inside every level
    # combination of the factors left over. For two factors that is the single
    # pair; for three, each pair inside every level of the third plus the
    # three-way term; and the same rule keeps going for more, which is why the
    # grid no longer needs an arity-specific branch.
    for (k in 2:length(var)) {
        combos <- utils::combn(seq_along(var), k)
        for (j in seq_len(ncol(combos))) {
            idx <- combos[, j]
            tested_vars <- var[idx]
            rest <- var[setdiff(seq_along(var), idx)]
            # Comparison combinations with the last tested factor varying
            # fastest, matching the nested loops this replaces.
            picks <- expand.grid(rev(lapply(comps_of[tested_vars], seq_along)),
                                 KEEP.OUT.ATTRS = FALSE)
            picks <- picks[, rev(seq_len(ncol(picks))), drop = FALSE]
            strata <- if (length(rest)) {
                expand.grid(lv[rest], stringsAsFactors = FALSE,
                            KEEP.OUT.ATTRS = FALSE)
            } else {
                NULL
            }
            for (p in seq_len(nrow(picks))) {
                tested <- lapply(seq_along(tested_vars), function(m) {
                    v <- tested_vars[[m]]
                    tested_of(v, comps_of[[v]][[picks[p, m]]])
                })
                if (is.null(strata)) {
                    add(tested = tested, stratum = NULL)
                    next
                }
                for (s in seq_len(nrow(strata))) {
                    add(
                        tested = tested,
                        stratum = stats::setNames(
                            as.character(strata[s, , drop = TRUE]), rest
                        )
                    )
                }
            }
        }
    }

    if (!length(rows)) return(NULL)
    out <- do.call(rbind, rows)
    # Factors that share level names ("yes"/"no" twice over) would otherwise
    # give two contrasts one label and one output folder.
    dup <- duplicated(out$label) | duplicated(out$label, fromLast = TRUE)
    out$label[dup] <- out$label_qualified[dup]
    out$label_qualified <- NULL
    rownames(out) <- NULL
    out
}

#' Main-effect contrast rows, one block per factor of a multi-factor `var`
#'
#' Each factor's own comparisons, taken across every level of the other
#' factors (the `"main"` class). Labels are the bare comparison strings, so
#' for a three-factor `var` a main effect and its per-stratum children share
#' a folder: `deseq2/C_vs_noC/C_vs_noC.html` alongside
#' `deseq2/C_vs_noC/CB__K/CB__K.html`.
#'
#' @param comps_fun Comparison generator, so each caller keeps its own level
#'   ordering convention ([generate_comparisons()] sorts alphabetically;
#'   [.comps_factor_order()] follows factor level order).
#' @noRd
.main_effect_specs <- function(cd, var, comps_fun) {
    rows <- lapply(var, function(v) {
        comps <- comps_fun(levels(factor(cd[[v]])))
        data.frame(
            label = comps,
            class = "main",
            type = classify_comparison_type(comps),
            stratum_var = NA_character_,
            stratum_level = NA_character_,
            test_var = v,
            comparison = comps,
            coef_subtract = NA_character_,
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    # Two factors can share level names ("yes"/"no" twice over), which would
    # otherwise give two different contrasts one label and one output folder.
    dup <- duplicated(out$label) | duplicated(out$label, fromLast = TRUE)
    out$label[dup] <- paste0(out$test_var[dup], "_", out$label[dup])
    rownames(out) <- NULL
    out
}

#' `rbind()` spec blocks that need not carry identical columns
#'
#' The three-factor grid adds `stratum_tag`, which the main-effect and
#' two-factor blocks have no use for.
#' @noRd
.rbind_specs <- function(...) {
    parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, list(...))
    if (!length(parts)) return(NULL)
    cols <- unique(unlist(lapply(parts, names)))
    parts <- lapply(parts, function(df) {
        df <- as.data.frame(df, stringsAsFactors = FALSE)
        for (col in setdiff(cols, names(df))) df[[col]] <- NA_character_
        df[, cols, drop = FALSE]
    })
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    new_contrast_specs(out)
}

#' Pairwise / one-vs-rest comparisons in factor level order
#'
#' Unlike [generate_comparisons()], which sorts level names alphabetically,
#' this names each pairwise comparison non-reference vs reference following
#' `levels()` order (first level = reference).
#' @noRd
.comps_factor_order <- function(lvls) {
    lvls <- as.character(lvls)
    if (length(lvls) < 2L) {
        stop("Need at least 2 levels to build a contrast.", call. = FALSE)
    }
    if (length(lvls) == 2L) {
        return(paste0(lvls[2L], "_vs_", lvls[1L]))
    }
    pairs <- utils::combn(seq_along(lvls), 2L)
    pairwise <- apply(pairs, 2L, function(ij) {
        paste0(lvls[ij[2L]], "_vs_", lvls[ij[1L]])
    })
    onetoall <- vapply(seq_along(lvls), function(i) {
        paste0(lvls[i], "_vs_", paste(lvls[-i], collapse = "+"))
    }, character(1))
    c(pairwise, onetoall)
}

#' Within-stratum contrast grid, for any number of grouping factors
#'
#' Each factor compared inside one fixed level combination of the others,
#' fitted on that subset alone. Only generated when
#' `expand_contrast_specs(stratified = TRUE)`; the interaction coefficients
#' themselves come from [.interaction_specs()].
#'
#' The two arities this package shipped with name their rows differently, and
#' both conventions are kept so existing output folders and report sections
#' do not move:
#' \describe{
#'   \item{two factors}{`class` names the factor held **fixed**
#'     (`"within_A"` = inside a level of `var[1]`, testing `var[2]`), the
#'     label is `"<level>__<comparison>"`, and comparisons are named
#'     alphabetically by [generate_comparisons()].}
#'   \item{three or more}{`class` names the factor being **tested**
#'     (`"within_A"` = testing `var[1]` inside the others), the label is the
#'     nested path `"<comparison>/<stratum_tag>"` with `stratum_tag` joining
#'     the held-fixed level values with `"__"`, and comparisons follow factor
#'     level order via [.comps_factor_order()].}
#' }
#' @noRd
.within_stratum_specs <- function(cd, var) {
    n <- length(var)
    if (n < 2L) return(NULL)
    lv <- lapply(var, function(v) levels(factor(cd[[v]])))
    names(lv) <- var
    # Class letters follow `var` order: A, B, C, ... then within_<i> beyond Z.
    letter <- function(i) if (i <= 26L) LETTERS[[i]] else as.character(i)

    rows <- list()
    # Block order is historical: with two factors the grid starts with the
    # comparisons of `var[2]` inside `var[1]`, with three or more it starts
    # with `var[1]` tested inside the others.
    for (i in if (n == 2L) rev(seq_along(var)) else seq_along(var)) {
        tested <- var[[i]]
        stratum_vars <- var[-i]
        comps <- if (n == 2L) {
            generate_comparisons(lv[[tested]])
        } else {
            .comps_factor_order(lv[[tested]])
        }
        types <- classify_comparison_type(comps)
        # Which factor the class names differs by arity (see above).
        class <- paste0("within_", letter(if (n == 2L) (3L - i) else i))
        grid <- expand.grid(lv[stratum_vars], stringsAsFactors = FALSE,
                            KEEP.OUT.ATTRS = FALSE)
        for (k in seq_len(nrow(grid))) {
            sl <- as.character(grid[k, , drop = TRUE])
            stratum_tag <- paste(sl, collapse = "__")
            for (j in seq_along(comps)) {
                row <- data.frame(
                    label = if (n == 2L) {
                        paste0(stratum_tag, "__", comps[j])
                    } else {
                        file.path(comps[j], stratum_tag)
                    },
                    class = class,
                    type = types[j],
                    stratum_var = .pack_terms(stratum_vars),
                    stratum_level = .pack_terms(sl),
                    test_var = tested,
                    comparison = comps[j],
                    coef_subtract = NA_character_,
                    stringsAsFactors = FALSE
                )
                if (n > 2L) row$stratum_tag <- stratum_tag
                rows[[length(rows) + 1L]] <- row
            }
        }
    }
    if (!length(rows)) return(NULL)
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

#' File-name stem for a contrast label
#'
#' Nested three-factor labels look like `"C_vs_noC/CB__noK"`; artifact
#' files use the basename (`"CB__noK_diffexp.csv"`). Flat labels
#' (`"IgG_vs_IgM"`) are returned unchanged.
#'
#' @param label Contrast label / folder key.
#' @return Character scalar, `basename(label)`.
#' @export
contrast_file_stem <- function(label) {
    basename(as.character(label)[[1L]])
}

#' Nesting depth of a contrast label (`1` = flat, `2` = nested)
#' @noRd
.contrast_nest_depth <- function(label) {
    parts <- strsplit(as.character(label)[[1L]], "/+", perl = TRUE)[[1L]]
    parts <- parts[nzchar(parts)]
    max(1L, length(parts))
}

#' `"factor = level"` label for each spec row's held-fixed factor
#'
#' Used by the report templates to group contrasts by the stratum they hold
#' fixed. Rows that hold nothing fixed (main effects, three-way
#' interactions) give `NA`.
#'
#' @param specs A contrast-specs data frame from [expand_contrast_specs()].
#' @return Character vector, one element per row of `specs`.
#' @export
contrast_stratum_labels <- function(specs) {
    if (is.null(specs) || !nrow(specs)) return(character(0))
    lvl <- as.character(specs$stratum_level)
    var <- as.character(specs$stratum_var)
    out <- paste0(var, " = ", lvl)
    out[is.na(lvl) | !nzchar(lvl) | is.na(var) | !nzchar(var)] <- NA_character_
    out
}

#' Human-readable section title for a contrast class
#'
#' @param class One of `"single"`, `"main"`, `"within_A"`, `"within_B"`,
#'   `"within_C"`, `"interaction"`.
#' @param var Character vector of length 1–3.
#' @return A character string.
#' @export
contrast_class_title <- function(class, var) {
    if (identical(class, "single") || identical(class, "pairwise") ||
        identical(class, "one_vs_all")) {
        return(NA_character_)
    }
    if (identical(class, "main")) return("Main effects (across all other factors)")
    if (identical(class, "interaction")) {
        return("Interactions (difference of differences)")
    }
    if (length(var) < 2) return(as.character(class))
    # "within_A" ... "within_Z", then "within_4" and up for wide designs.
    tag <- sub("^within_", "", as.character(class))
    i <- if (tag %in% LETTERS) match(tag, LETTERS) else suppressWarnings(as.integer(tag))
    if (is.na(i) || i > length(var)) return(as.character(class))

    if (length(var) == 2L) {
        # Two factors: the class names the factor held fixed.
        held <- var[[i]]
        tested <- var[[3L - i]]
        return(paste0("Within ", held, " \u2014 ", tested, " effects"))
    }
    # Three or more: the class names the factor being tested.
    paste0("Within ", paste(var[-i], collapse = " \u00d7 "), " \u2014 ", var[[i]],
           " effects")
}
