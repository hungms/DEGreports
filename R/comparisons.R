#' Generate all pairwise and one-vs-rest comparisons for a set of group levels
#'
#' Given the levels of a grouping variable, generates every pairwise contrast
#' (`"A_vs_B"`) plus, when there are more than two levels, every one-vs-rest
#' contrast (`"A_vs_B+C"`, meaning A vs. the pooled remaining levels). Pairwise
#' comparison names are ordered alphabetically (`"A_vs_B"`, never both
#' `"A_vs_B"` and `"B_vs_A"`) so that the same contrast is never generated
#' twice under two different names.
#'
#' @param var_levels Character vector of group levels (as they appear in
#'   `levels(colData(dds)[[var]])`). Duplicates are removed internally.
#' @return A character vector of comparison names in `"A_vs_B"` /
#'   `"A_vs_B+C"` format.
#' @examples
#' generate_comparisons(c("Control", "Treated"))
#' generate_comparisons(c("A", "B", "C"))
#' @export
generate_comparisons <- function(var_levels) {
    var_levels <- unique(var_levels)
    if (length(var_levels) < 2) {
        stop("`var_levels` must contain at least 2 unique levels; got ", length(var_levels), ".", call. = FALSE)
    }

    if (length(var_levels) == 2) {
        # canonical alphabetical order: second (alphabetically later) level is
        # the "test" group, first is the reference -- matches DESeq2 contrast
        # convention of contrast(var, level2, level1) but naming stays sorted.
        sorted_levels <- sort(as.character(var_levels))
        comparison_vec <- paste0(sorted_levels[2], "_vs_", sorted_levels[1])
        return(comparison_vec)
    }

    combos <- utils::combn(var_levels, 2)
    pairwise_vec <- apply(combos, 2, function(pair) {
        if (as.character(pair[2]) < as.character(pair[1])) {
            paste0(as.character(pair[2]), "_vs_", as.character(pair[1]))
        } else {
            paste0(as.character(pair[1]), "_vs_", as.character(pair[2]))
        }
    })

    onetoall_vec <- vapply(seq_along(var_levels), function(i) {
        paste0(var_levels[i], "_vs_", paste0(var_levels[-i], collapse = "+"))
    }, character(1))

    unique(c(pairwise_vec, onetoall_vec))
}

#' Classify comparison names as pairwise or one-vs-rest
#'
#' [generate_comparisons()] names a one-vs-rest contrast with a `"+"`
#' joining the pooled reference levels on the right-hand side of `"_vs_"`
#' (e.g. `"A_vs_B+C"`), while a pairwise contrast never contains a `"+"`
#' (e.g. `"A_vs_B"`). This is the sole structural marker distinguishing
#' the two classes, and is what the bundled HTML report uses to route
#' each comparison's differential-expression and GSEA results into its
#' "Pairwise" or "One vs. all" subsection.
#'
#' @param comparison_names Character vector of comparison names, in the
#'   `"A_vs_B"` / `"A_vs_B+C"` format produced by [generate_comparisons()].
#' @return A character vector the same length as `comparison_names`, each
#'   element `"pairwise"` or `"one_vs_all"`.
#' @examples
#' classify_comparison_type(c("B_vs_A", "A_vs_B+C", "C_vs_A+B"))
#' @export
classify_comparison_type <- function(comparison_names) {
    ifelse(grepl("+", comparison_names, fixed = TRUE), "one_vs_all", "pairwise")
}

#' Fail early with a clear message if the design matrix is rank-deficient
#'
#' Reconstructs the model matrix DESeq2 would use for `dds`'s current design
#' and checks whether it is full column rank *before* handing off to
#' `DESeq2::DESeq()` / `DESeq2::results()`. This mirrors the rank check
#' `DESeq2` performs internally, but raises a clearer, comparison-specific
#' error message pointing at the offending comparison and (for pooled
#' one-vs-rest comparisons) suggesting how to avoid it.
#'
#' @param dds_obj A `DESeqDataSet` object with `design(dds_obj)` already set.
#' @param comparison The comparison name currently being evaluated (used only
#'   for the error message).
#' @param group_var The grouping column name (used only for the error message).
#' @return Invisibly, `NULL`, if the model matrix is full rank.
#' @export
assert_full_rank_model_matrix <- function(dds_obj, comparison, group_var) {
    des <- DESeq2::design(dds_obj)
    if (!inherits(des, "formula")) {
        return(invisible(NULL))
    }
    mm <- stats::model.matrix(des, data = as.data.frame(SummarizedExperiment::colData(dds_obj)))
    rk <- qr(mm)$rank
    nc <- ncol(mm)
    if (rk >= nc) {
        return(invisible(NULL))
    }

    pooled <- grepl("+", comparison, fixed = TRUE)
    msg <- paste0(
        "Model matrix is not full rank for comparison \"", comparison, "\" ",
        "(rank ", rk, " < ", nc, " columns). DESeq2 cannot fit this contrast.\n\n"
    )
    if (pooled) {
        msg <- paste0(
            msg,
            "This is a pooled one-vs-rest comparison; pooling can make the design ",
            "matrix rank-deficient when other design terms are confounded with the ",
            "pooled grouping. Consider restricting to pairwise comparisons only ",
            "(e.g. `comps <- generate_comparisons(levels(colData(dds)[['",
            group_var, "']])); comps[!grepl(\"+\", comps, fixed = TRUE)]`)."
        )
    } else {
        msg <- paste0(
            msg,
            "Check that terms in the design are not linearly redundant with `",
            group_var, "` for this contrast."
        )
    }
    stop(msg, call. = FALSE)
}

#' Merge group levels for a pooled ("A_vs_B+C") comparison
#'
#' Creates a temporary grouping factor on a *copy* of `dds` in which the
#' pooled levels in `comparison` are collapsed into a single combined level
#' (named by joining the pooled level names with `"+"`), while leaving the
#' reference/original `dds` object passed in by the caller untouched. Levels
#' are matched exactly (not via regex/substring matching) so that a level
#' name which happens to be a substring of another level is never
#' mis-collapsed.
#'
#' @param dds A `DESeqDataSet` object.
#' @param var Grouping column name.
#' @param comparison A comparison string, e.g. `"A_vs_B+C"`.
#' @return A list with elements `dds` (a copy of `dds` with a new column
#'   `paste0(var, "__comparison")` holding the merged grouping factor) and
#'   `merged_var` (the name of that new column).
#' @noRd
merge_comparison_levels <- function(dds, var, comparison) {
    parts <- strsplit(comparison, "_vs_", fixed = TRUE)[[1]]
    test_groups <- unlist(strsplit(parts[1], "+", fixed = TRUE))
    ref_groups <- unlist(strsplit(parts[2], "+", fixed = TRUE))

    grp <- as.character(SummarizedExperiment::colData(dds)[[var]])
    merged <- rep(NA_character_, length(grp))
    merged[grp %in% test_groups] <- parts[1]
    merged[grp %in% ref_groups] <- parts[2]
    # any level not part of this comparison is dropped (NA); samples with NA
    # are excluded from this contrast's model fit.

    merged_var <- paste0(var, "__comparison")
    dds_sub <- dds[, !is.na(merged)]
    merged <- merged[!is.na(merged)]
    SummarizedExperiment::colData(dds_sub)[[merged_var]] <- factor(merged, levels = c(parts[2], parts[1]))

    list(dds = dds_sub, merged_var = merged_var, ref = parts[2], test = parts[1])
}
