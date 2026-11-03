# What every contrast in a run actually asked, and of which samples.
#
# The contrast grid used to be legible only by reading the labels: nothing
# said which design was fitted for a given row, how many samples were behind
# it, or which coefficients it combined, and the only prose explanation
# (interaction_question_table()) covered interaction rows alone. The manifest
# answers those questions in one table, written next to the results and shown
# at the top of the DESeq2 report page.

#' Samples a contrast's model is fitted on
#'
#' Not every contrast is fitted on every sample: a pooled one-vs-rest
#' comparison drops the levels it does not name, and a within-stratum
#' contrast is fitted on its stratum alone. A pairwise comparison of one
#' factor's levels, by contrast, is read off a model fitted on all samples.
#' @noRd
.contrast_n_samples <- function(spec, cd) {
    fm <- spec$fit_method[[1L]]
    # An interaction contrast names a stratum ("... within celltype = CB") but
    # is *fitted* on every sample -- the stratum is expressed in the
    # coefficients it combines, not by subsetting. That is the whole point of
    # preferring it to a within-stratum fit, so the count must say so.
    if (identical(fm, "coefs")) return(nrow(cd))
    keep <- rep(TRUE, nrow(cd))
    stratum <- contrast_stratum(spec)
    for (sv in names(stratum)) {
        if (!sv %in% names(cd)) next
        keep <- keep & as.character(cd[[sv]]) == stratum[[sv]]
    }
    tv <- contrast_test_vars(spec)
    comp <- spec$comparison[[1L]]
    pooled <- !is.na(comp) && grepl("+", comp, fixed = TRUE)
    if (pooled && length(tv) == 1L && tv %in% names(cd)) {
        sides <- strsplit(comp, "_vs_", fixed = TRUE)[[1L]]
        lv <- unlist(strsplit(sides, "+", fixed = TRUE), use.names = FALSE)
        keep <- keep & as.character(cd[[tv]]) %in% lv
    }
    sum(keep)
}

#' Design formula a contrast is fitted with
#' @noRd
.contrast_design <- function(spec, var, batch, dds = NULL) {
    fm <- spec$fit_method[[1L]]
    if (identical(fm, "lrt")) {
        full <- if ("full" %in% names(spec)) spec$full[[1L]] else NA_character_
        reduced <- if ("reduced" %in% names(spec)) spec$reduced[[1L]] else NA_character_
        if (is.na(full)) return(NA_character_)
        return(paste0(full, "  vs  ", reduced))
    }
    if (identical(fm, "precomputed")) return(NA_character_)
    tv <- contrast_test_vars(spec)
    f <- switch(
        fm,
        levels_additive = resolve_design(c(setdiff(var, tv), tv), batch = batch,
                                         include_interactions = FALSE),
        levels_subset = resolve_design(tv, batch = batch),
        coefs = resolve_design(var, batch = batch, include_interactions = TRUE),
        if (!is.null(dds)) DESeq2::design(dds) else resolve_design(var, batch = batch)
    )
    paste(deparse(f), collapse = "")
}

#' Tabulate what each contrast asks, how it is fitted, and of which samples
#'
#' One row per contrast, with the question in prose, the design formula
#' actually fitted, the number of samples behind it, the model coefficients it
#' combines (for interaction contrasts), and the directory its results are
#' written to. [run_diffexp_pipeline()] writes this as `contrasts.csv` beside
#' the per-contrast folders and the report shows it at the top of the DESeq2
#' page.
#'
#' @param specs A `contrast_specs` table.
#' @param dds Optional `DESeqDataSet`, for the sample counts.
#' @param var Grouping column name(s).
#' @param batch Optional covariate column name(s).
#' @return A data frame with columns `label`, `class`, `fit_method`,
#'   `question`, `design`, `n_samples`, `coefficients` and `path`.
#' @examples
#' \dontrun{
#' specs <- expand_contrast_specs(dds, var = c("tissue", "genotype"),
#'                                include_interactions = TRUE)
#' contrast_manifest(specs, dds, var = c("tissue", "genotype"))
#' }
#' @seealso [contrast_question()], [expand_contrast_specs()]
#' @export
contrast_manifest <- function(specs, dds = NULL, var = NULL, batch = NULL) {
    specs <- as_contrast_specs(specs)
    if (is.null(specs) || !nrow(specs)) {
        return(data.frame(
            label = character(0), class = character(0), fit_method = character(0),
            question = character(0), design = character(0),
            n_samples = integer(0), coefficients = character(0),
            path = character(0), stringsAsFactors = FALSE
        ))
    }
    cd <- if (is.null(dds)) {
        NULL
    } else {
        as.data.frame(SummarizedExperiment::colData(dds))
    }
    if (is.null(var)) {
        var <- unique(unlist(lapply(seq_len(nrow(specs)), function(i) {
            contrast_test_vars(specs[i, , drop = FALSE])
        })))
    }

    rows <- lapply(seq_len(nrow(specs)), function(i) {
        spec <- specs[i, , drop = FALSE]
        co <- contrast_coefs(spec)
        coefficients <- if (length(co$weights)) {
            paste0(names(co$weights), " \u00d7 ",
                   format(unname(co$weights), digits = 3, trim = TRUE),
                   collapse = ", ")
        } else if (length(co$pos)) {
            paste0(paste(co$pos, collapse = " + "),
                   if (length(co$neg)) paste0(" \u2212 ", paste(co$neg, collapse = " \u2212 ")) else "")
        } else {
            NA_character_
        }
        data.frame(
            label = spec$label[[1L]],
            class = spec$class[[1L]],
            fit_method = spec$fit_method[[1L]],
            question = spec$question[[1L]],
            design = .contrast_design(spec, var, batch, dds),
            n_samples = if (is.null(cd)) NA_integer_ else .contrast_n_samples(spec, cd),
            coefficients = coefficients,
            path = spec$path[[1L]],
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

#' Print a condensed contrast manifest to the console
#'
#' Shown before fitting starts, so a run says what it is about to compare
#' rather than only what it compared.
#' @noRd
.announce_contrasts <- function(manifest) {
    if (is.null(manifest) || !nrow(manifest)) return(invisible(NULL))
    message("Contrasts in this run (", nrow(manifest), "):")
    for (i in seq_len(nrow(manifest))) {
        n <- manifest$n_samples[[i]]
        message("  \u2022 ", manifest$label[[i]],
                if (!is.na(n)) paste0("  [n=", n, "]") else "",
                "\n      ", manifest$question[[i]])
    }
    invisible(NULL)
}
