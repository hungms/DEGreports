test_that("summarize_genes sums rows that share a symbol", {
    mat <- matrix(
        c(1L, 2L, 3L,
          4L, 5L, 6L,
          7L, 8L, 9L),
        nrow = 3L, byrow = TRUE,
        dimnames = list(c("t1", "t2", "t3"), c("s1", "s2", "s3"))
    )
    out <- summarize_genes(mat, c("A", "A", "B"))
    expect_equal(rownames(out), c("A", "B"))
    expect_equal(as.integer(out["A", ]), c(5L, 7L, 9L))
    expect_equal(as.integer(out["B", ]), c(7L, 8L, 9L))
})

test_that("summarize_genes drops empty symbols", {
    mat <- matrix(1:4, nrow = 2L, dimnames = list(c("t1", "t2"), c("s1", "s2")))
    out <- summarize_genes(mat, c("A", NA_character_))
    expect_equal(rownames(out), "A")
    expect_equal(nrow(out), 1L)
})

test_that("summarize_genes_dds rebuilds a DESeqDataSet by gene symbol", {
    counts <- matrix(
        c(10L, 20L,
          1L, 2L,
          3L, 4L),
        nrow = 3L, byrow = TRUE,
        dimnames = list(c("ENST1", "ENST2", "ENST3"), c("s1", "s2"))
    )
    coldata <- data.frame(
        condition = factor(c("A", "B")),
        row.names = c("s1", "s2")
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~condition)
    SummarizedExperiment::rowData(dds)$gene <- c("GeneX", "GeneX", "GeneY")

    dds2 <- summarize_genes_dds(dds)
    expect_s4_class(dds2, "DESeqDataSet")
    expect_equal(rownames(dds2), c("GeneX", "GeneY"))
    expect_equal(as.integer(SummarizedExperiment::assay(dds2)["GeneX", ]), c(11L, 22L))
    expect_equal(as.character(SummarizedExperiment::rowData(dds2)$gene),
                 c("GeneX", "GeneY"))
    expect_equal(DESeq2::design(dds2), DESeq2::design(dds))
})

test_that("load_nfcore loads RData, maps symbols, and summarizes", {
    counts <- matrix(
        c(5L, 6L,
          7L, 8L,
          1L, 1L),
        nrow = 3L, byrow = TRUE,
        dimnames = list(c("ENSG1", "ENSG2", "ENSG3"), c("s1", "s2"))
    )
    coldata <- data.frame(
        condition = factor(c("A", "B")),
        row.names = c("s1", "s2")
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~condition)

    tmp <- withr::local_tempdir()
    rdata <- file.path(tmp, "nfcore.RData")
    save(dds, file = rdata)

    tx2gene <- file.path(tmp, "tx2gene.tsv")
    utils::write.table(
        data.frame(
            V1 = c("tx1", "tx2", "tx3"),
            V2 = c("ENSG1", "ENSG2", "ENSG3"),
            V3 = c("Same", "Same", "Other"),
            stringsAsFactors = FALSE
        ),
        tx2gene, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )

    out <- load_nfcore(rdata, tx2gene)
    expect_s4_class(out, "DESeqDataSet")
    expect_equal(sort(rownames(out)), c("Other", "Same"))
    expect_equal(as.integer(SummarizedExperiment::assay(out)["Same", ]), c(12L, 14L))
})
