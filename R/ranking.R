#' Compute a gene ranking score from a DESeq2 results data frame
#'
#' Computes the ranking metric used to order genes for GSEA and for
#' "top DEG" tables, given a results data frame with (at minimum)
#' `log2FoldChange`, `pvalue`/`padj`, and `stat` columns (as produced by
#' [as.data.frame()] on a `DESeqResults` object).
#'
#' Supported metrics:
#' \describe{
#'   \item{`"pxfc"`}{`-log10(pvalue) * sign(log2FoldChange) * abs(log2FoldChange)`
#'     -- a combined significance-and-magnitude score. Ties are broken toward
#'     genes with larger fold-change when p-values are equal.}
#'   \item{`"padj"`}{`-log10(padj)`, signed by the direction of
#'     `log2FoldChange` (i.e. `-log10(padj) * sign(log2FoldChange)`).}
#'   \item{`"log2FoldChange"`}{The shrunken/raw log2 fold-change itself.}
#'   \item{`"stat"`}{The DESeq2 Wald statistic column, `stat`.}
#' }
#'
#' @param res_df A data frame with columns `log2FoldChange`, `pvalue`,
#'   `padj`, `stat` (as from `as.data.frame(DESeq2::results(...))`).
#' @param ranking One of `"pxfc"`, `"padj"`, `"log2FoldChange"`, `"stat"`.
#' @return A numeric vector, the same length as `nrow(res_df)`, giving the
#'   ranking score for each gene (higher = more "up"/significant in the
#'   test direction).
#' @examples
#' df <- data.frame(
#'   log2FoldChange = c(2, -1.5, 0.1),
#'   pvalue = c(0.001, 0.02, 0.9),
#'   padj = c(0.01, 0.08, 0.95),
#'   stat = c(5, -3, 0.2)
#' )
#' compute_ranking_score(df, "pxfc")
#' @export
compute_ranking_score <- function(res_df, ranking = "log2FoldChange") {
    ranking <- validate_ranking(ranking)
    required <- c("log2FoldChange", "pvalue", "padj", "stat")
    missing <- setdiff(required, colnames(res_df))
    if (length(missing) > 0) {
        stop("`res_df` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }

    score <- switch(ranking,
        pxfc = {
            p <- res_df$pvalue
            p[is.na(p)] <- 1
            fc <- res_df$log2FoldChange
            fc[is.na(fc)] <- 0
            -log10(p) * sign(fc) * abs(fc)
        },
        padj = {
            p <- res_df$padj
            p[is.na(p)] <- 1
            fc <- res_df$log2FoldChange
            fc[is.na(fc)] <- 0
            -log10(p) * sign(fc)
        },
        log2FoldChange = {
            fc <- res_df$log2FoldChange
            fc[is.na(fc)] <- 0
            fc
        },
        stat = {
            s <- res_df$stat
            s[is.na(s)] <- 0
            s
        }
    )

    # replace any remaining non-finite values (e.g. p == 0 -> Inf) with a
    # large but finite sentinel so downstream sorts/fgsea don't choke
    score[is.infinite(score) & score > 0] <- max(score[is.finite(score)], 0, na.rm = TRUE) + 1
    score[is.infinite(score) & score < 0] <- min(score[is.finite(score)], 0, na.rm = TRUE) - 1
    score[is.na(score)] <- 0
    score
}

#' Build a named, sorted rank vector for GSEA from a results data frame
#'
#' @param res_df A data frame as for [compute_ranking_score()], additionally
#'   containing a gene identifier column.
#' @param gene_col Name of the column holding gene identifiers to use as
#'   names on the returned vector (default `"gene"`).
#' @param ranking Ranking metric, passed to [compute_ranking_score()].
#' @return A named numeric vector sorted in decreasing order, with
#'   `NA`/duplicate gene identifiers removed (keeping the first, highest
#'   |score|, occurrence of any duplicate).
#' @export
build_gsea_ranks <- function(res_df, gene_col = "gene", ranking = "log2FoldChange") {
    if (!gene_col %in% colnames(res_df)) {
        stop("`gene_col` (\"", gene_col, "\") is not a column of `res_df`.", call. = FALSE)
    }
    score <- compute_ranking_score(res_df, ranking)
    genes <- res_df[[gene_col]]
    keep <- !is.na(genes) & !is.na(score)
    score <- score[keep]
    genes <- genes[keep]

    ord <- order(abs(score), decreasing = TRUE)
    score <- score[ord]
    genes <- genes[ord]
    dup <- duplicated(genes)
    score <- score[!dup]
    genes <- genes[!dup]

    ranks <- stats::setNames(score, genes)
    sort(ranks, decreasing = TRUE)
}
