#' Sum count rows that share a gene symbol
#'
#' Collapses a genes/transcripts × samples count matrix by summing every
#' set of rows whose symbol in `gene_sym_vec` is identical. Empty / `NA`
#' symbols are dropped.
#'
#' @param input Numeric/integer matrix or data frame (rows = features,
#'   columns = samples).
#' @param gene_sym_vec Character vector of gene symbols, length
#'   `nrow(input)`.
#' @param ... Ignored (accepted for compatibility with callers that pass
#'   extra args such as `normalize`).
#' @return A matrix with one row per unique non-empty symbol and the same
#'   columns as `input`. Row names are the gene symbols.
#' @export
summarize_genes <- function(input, gene_sym_vec, ...) {
    if (is.data.frame(input)) {
        mat <- as.matrix(input)
    } else {
        mat <- input
    }
    if (!is.matrix(mat)) {
        stop("`input` must be a matrix or data frame of counts.", call. = FALSE)
    }
    gene_sym_vec <- as.character(gene_sym_vec)
    if (length(gene_sym_vec) != nrow(mat)) {
        stop("`gene_sym_vec` must have length nrow(input) (", nrow(mat),
             "); got ", length(gene_sym_vec), ".", call. = FALSE)
    }
    ok <- !is.na(gene_sym_vec) & nzchar(gene_sym_vec)
    if (!any(ok)) {
        stop("No non-empty gene symbols available to summarize.", call. = FALSE)
    }
    n_drop <- sum(!ok)
    if (n_drop > 0L) {
        message("Dropping ", n_drop, " feature(s) with missing/empty gene symbols.")
    }
    mat <- mat[ok, , drop = FALSE]
    gene_sym_vec <- gene_sym_vec[ok]
    # Preserve integer counts when possible (DESeq2 prefers integers).
    if (is.integer(mat) || all(mat == floor(mat), na.rm = TRUE)) {
        storage.mode(mat) <- "integer"
    }
    out <- rowsum(mat, group = gene_sym_vec, reorder = FALSE)
    n_in <- length(gene_sym_vec)
    n_out <- nrow(out)
    if (n_out < n_in) {
        message("Summarized ", n_in, " feature(s) → ", n_out,
                " unique gene symbol(s).")
    }
    out
}

#' Collapse a `DESeqDataSet` to one row per gene symbol
#'
#' Sums count rows that share the same symbol in `rowData(dds)[[gene_sym_col]]`
#' and returns a new `DESeqDataSet` with the original design formula and
#' matching `colData`. Useful when the object is still transcript-level, or
#' when multiple feature IDs map to the same gene symbol.
#'
#' @param dds A `DESeqDataSet`.
#' @param gene_sym_col Column name in `rowData(dds)` holding gene symbols.
#'   Default `"gene"`.
#' @param ... Passed to [summarize_genes()] (currently ignored aside from
#'   accepting legacy args such as `normalize`).
#' @return A new `DESeqDataSet` with gene-symbol rownames and
#'   `rowData(dds.new)$gene` set to those rownames.
#' @examples
#' \dontrun{
#' rowData(dds)$gene <- gene_symbols
#' dds_gene <- summarize_genes_dds(dds)
#' }
#' @export
summarize_genes_dds <- function(dds, gene_sym_col = "gene", ...) {
    dds <- validate_dds(dds)
    gene_sym_col <- as.character(gene_sym_col)[1]
    rd <- SummarizedExperiment::rowData(dds)
    if (!gene_sym_col %in% colnames(rd)) {
        stop("`gene_sym_col` \"", gene_sym_col, "\" not found in rowData(dds). ",
             "Available columns: ", paste(colnames(rd), collapse = ", "),
             call. = FALSE)
    }
    gene_sym_vec <- as.character(rd[[gene_sym_col]])
    counts_mat <- SummarizedExperiment::assay(dds)
    # Prefer raw counts slot when present (DESeqDataSet assays$"counts").
    if ("counts" %in% SummarizedExperiment::assayNames(dds)) {
        counts_mat <- SummarizedExperiment::assay(dds, "counts")
    }
    df <- summarize_genes(
        input = counts_mat,
        gene_sym_vec = gene_sym_vec,
        ...
    )
    # Integer matrix required by DESeqDataSetFromMatrix.
    if (!is.integer(df)) {
        storage.mode(df) <- "integer"
    }

    cd <- as.data.frame(SummarizedExperiment::colData(dds),
                        stringsAsFactors = FALSE)
    # Align colData to the (possibly reordered) count columns.
    missing_cols <- setdiff(colnames(df), rownames(cd))
    if (length(missing_cols)) {
        stop("Summarized count columns not in colData(dds): ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    metadata <- cd[colnames(df), , drop = FALSE]

    dds_new <- DESeq2::DESeqDataSetFromMatrix(
        countData = df,
        colData = metadata,
        design = DESeq2::design(dds)
    )
    SummarizedExperiment::rowData(dds_new)$gene <- rownames(df)
    dds_new
}

#' Load an nf-core/rnaseq `DESeqDataSet` and collapse to gene symbols
#'
#' Loads a saved `dds` object from an `.RData` / `.rda` file produced by
#' [nf-core/rnaseq](https://nf-co.re/rnaseq), attaches gene symbols from a
#' three-column `tx2gene` TSV (column 2 must match `rownames(dds)`; column 3
#' is the gene symbol), and collapses features that share a symbol via
#' [summarize_genes_dds()].
#'
#' @param rdata Path to an `.RData` / `.rda` file containing an object named
#'   `dds` (`DESeqDataSet`).
#' @param tx2gene Path to a tab-separated mapping file **without a header**,
#'   with at least three columns: typically transcript ID, feature ID
#'   matching `rownames(dds)`, and gene symbol.
#' @return A gene-symbol-level `DESeqDataSet`.
#' @examples
#' \dontrun{
#' dds <- load_nfcore(
#'   rdata = "salmon.merged.gene_counts.rds.RData",
#'   tx2gene = "tx2gene.tsv"
#' )
#' }
#' @export
load_nfcore <- function(rdata, tx2gene) {
    rdata <- normalizePath(rdata, mustWork = TRUE)
    tx2gene <- normalizePath(tx2gene, mustWork = TRUE)

    env <- new.env(parent = emptyenv())
    loaded <- load(rdata, envir = env)
    if (!exists("dds", envir = env, inherits = FALSE)) {
        stop("No object named `dds` found in ", rdata,
             ". Objects present: ", paste(loaded, collapse = ", "),
             call. = FALSE)
    }
    dds <- env$dds
    if (!methods::is(dds, "DESeqDataSet")) {
        stop("`dds` loaded from ", rdata, " is not a DESeqDataSet ",
             "(got ", paste(class(dds), collapse = "/"), ").", call. = FALSE)
    }

    genes <- utils::read.table(
        tx2gene, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
        quote = "", comment.char = ""
    )
    if (ncol(genes) < 3L) {
        stop("`tx2gene` must have at least 3 columns ",
             "(feature id in column 2, gene symbol in column 3); got ",
             ncol(genes), " in ", tx2gene, call. = FALSE)
    }
    genes$V2 <- as.character(genes$V2)
    genes$V3 <- as.character(genes$V3)
    # One symbol per feature id (first occurrence wins).
    genes <- genes[!is.na(genes$V2) & nzchar(genes$V2), c("V2", "V3"), drop = FALSE]
    genes <- genes[!duplicated(genes$V2), , drop = FALSE]

    rownames(dds) <- as.character(rownames(dds))
    sym <- genes$V3[match(rownames(dds), genes$V2)]
    n_miss <- sum(is.na(sym))
    if (n_miss > 0L) {
        message(n_miss, " / ", nrow(dds),
                " feature(s) had no match in tx2gene column 2.")
    }
    if (all(is.na(sym) | !nzchar(sym))) {
        stop("None of the dds rownames matched tx2gene column 2.",
             call. = FALSE)
    }
    SummarizedExperiment::rowData(dds)$gene <- sym

    summarize_genes_dds(dds, gene_sym_col = "gene")
}
