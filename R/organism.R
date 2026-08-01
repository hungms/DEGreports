#' Detect organism (human or mouse) from a vector of gene symbols
#'
#' Guesses whether a set of gene symbols follows human (HGNC) or mouse (MGI)
#' naming conventions, so the package can pick the right `msigdbr` species
#' for GSEA without requiring the caller to specify it explicitly. HGNC
#' human gene symbols are conventionally all upper-case letters/digits
#' (e.g. `"TP53"`, `"HLA-A"`, `"IL6"`); MGI mouse gene symbols are
#' conventionally title-case -- a single leading capital followed by
#' lower-case letters/digits (e.g. `"Trp53"`, `"Il6"`). The organism whose
#' convention the larger share of `genes` matches is returned, provided
#' that share clears `min_fraction` and exceeds the other convention's
#' share; otherwise an error is raised (rather than guessing) since the
#' whole point of auto-detection is to avoid silently mislabelling gene
#' sets fetched for GSEA.
#'
#' @param genes Character vector of gene symbols (e.g. `rownames(dds)`, or
#'   the values of a diffexp results table's gene-identifier column).
#'   `NA`s and empty strings are dropped before classification.
#' @param min_fraction Minimum fraction of `genes` that must match the
#'   winning convention for detection to succeed. Default `0.5`.
#' @return `"human"` or `"mouse"`.
#' @examples
#' detect_organism(c("TP53", "BRCA1", "HLA-A", "IL6"))
#' detect_organism(c("Trp53", "Brca1", "Il6"))
#' @export
detect_organism <- function(genes, min_fraction = 0.5) {
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0) {
        stop("`genes` contains no usable (non-NA, non-empty) gene symbols to detect an organism from.",
             call. = FALSE)
    }

    # HGNC human convention: entirely upper-case letters/digits, optionally
    # with internal hyphens/dots/underscores (e.g. "HLA-A", "C1QA").
    human_re <- "^[A-Z][A-Z0-9]*([._-][A-Z0-9]+)*$"
    # MGI mouse convention: one leading capital, then lower-case letters/
    # digits, optionally with internal hyphens/dots/underscores (e.g.
    # "Trp53", "Gm1234").
    mouse_re <- "^[A-Z][a-z0-9]+([._-][A-Za-z0-9]+)*$"

    frac_human <- mean(grepl(human_re, genes))
    frac_mouse <- mean(grepl(mouse_re, genes))

    if (frac_human >= min_fraction && frac_human > frac_mouse) {
        return("human")
    }
    if (frac_mouse >= min_fraction && frac_mouse > frac_human) {
        return("mouse")
    }

    stop(
        "Could not confidently detect organism from gene symbols: ",
        sprintf("%.0f%%", frac_human * 100), " look human-style (all-caps, e.g. \"TP53\"), ",
        sprintf("%.0f%%", frac_mouse * 100), " look mouse-style (title-case, e.g. \"Trp53\"). ",
        "Detection needs standard HGNC/MGI gene *symbols* as gene identifiers (not Ensembl/Entrez ",
        "IDs and not a mix of organisms) -- check `rownames(dds)` / the gene-identifier column.",
        call. = FALSE
    )
}
