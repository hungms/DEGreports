# Tests for prepare_diffexp_from_deg() / run_deg() (DEG-table entry point).

.make_merged_deg <- function(n_genes = 40, seed = 1) {
    set.seed(seed)
    genes <- paste0("GENE", seq_len(n_genes))
    # Mix in a few real human symbols so detect_organism() succeeds.
    genes[seq_len(min(8L, n_genes))] <- c(
        "TP53", "BRCA1", "MYC", "EGFR", "STAT3", "JUN", "FOS", "IL6"
    )[seq_len(min(8L, n_genes))]
    comps <- c("groupA_vs_groupB", "groupA_vs_groupC", "groupB_vs_groupC")
    rows <- lapply(comps, function(comp) {
        data.frame(
            gene = genes,
            comparison = comp,
            log2FoldChange = rnorm(n_genes, sd = 1.2),
            padj = runif(n_genes, 0, 1),
            baseMean = rexp(n_genes, rate = 0.01) + 1,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

test_that("prepare_diffexp_from_deg() validates required columns and splits by comparison", {
    deg <- .make_merged_deg()
    out <- prepare_diffexp_from_deg(deg, alpha = 0.05, lfc_threshold = 0.5)

    expect_equal(out$source, "deg_table")
    expect_equal(out$gene_col, "gene")
    expect_setequal(names(out$by_comparison), unique(deg$comparison))
    expect_true(all(c("stat", "sig", "rank_score", "mean_cpm", "pvalue") %in%
                        colnames(out$results)))
    expect_equal(nrow(out$results), nrow(deg))
    expect_true(all(out$contrast_specs$comparison %in% names(out$by_comparison)))
})

test_that("prepare_diffexp_from_deg() remaps common aliases", {
    deg <- data.frame(
        Gene = c("TP53", "BRCA1"),
        contrast = c("A_vs_B", "A_vs_B"),
        log2FC = c(1.5, -2),
        FDR = c(0.01, 0.2),
        stringsAsFactors = FALSE
    )
    out <- prepare_diffexp_from_deg(deg)
    expect_true("gene" %in% colnames(out$results))
    expect_true("comparison" %in% colnames(out$results))
    expect_true("log2FoldChange" %in% colnames(out$results))
    expect_true("padj" %in% colnames(out$results))
    expect_equal(names(out$by_comparison), "A_vs_B")
})

test_that("prepare_diffexp_from_deg() errors clearly without comparison / padj", {
    expect_error(
        prepare_diffexp_from_deg(data.frame(gene = "TP53", log2FoldChange = 1)),
        "missing required column"
    )
    expect_error(
        prepare_diffexp_from_deg(data.frame(
            gene = "TP53", comparison = "A_vs_B", log2FoldChange = 1
        )),
        "padj"
    )
})

test_that("run_deg() writes deseq2 CSVs and multihit without counts", {
    skip_if_not_installed("msigdbr")
    skip_if_not_installed("fgsea")

    deg <- .make_merged_deg(n_genes = 30)
    # Force a few significant rows so multihit has DEG hits.
    deg$padj[1:5] <- 1e-4
    deg$log2FoldChange[1:5] <- c(2, -2, 1.5, -1.8, 2.2)

    out_dir <- withr::local_tempdir()
    # PPI always runs in the pipeline now; stub it out so the test does not
    # reach for the STRING download.
    local_mocked_bindings(run_ppi_pipeline = function(...) NULL)
    res <- suppressWarnings(run_deg(
        deg_table = deg,
        collections = "HALLMARK",
        run_gsea = TRUE,
        run_decoupler_tf = FALSE,
        run_viper = FALSE,
        run_cspa = FALSE,
        render_html = FALSE,
        save_dir = out_dir
    ))

    expect_true(file.exists(file.path(out_dir, "deseq2",
                                      "all_comparisons_diffexp.csv")))
    expect_true(dir.exists(file.path(out_dir, "deseq2", "groupA_vs_groupB")))
    expect_false(is.null(res$diffexp))
    expect_false(is.null(res$gsea))
    # Multihit may be empty of layers when TF/CSPA off, but should not error.
    expect_true(is.null(res$multihit) || is.list(res$multihit))
})

test_that("render_report() skips QC when dds is NULL", {
    skip_if_not_installed("rmarkdown")
    skip_if_not_installed("htmltools")
    skip_if_not_installed("plotly")
    skip_if_not_installed("DT")
    skip_if_not_installed("msigdbr")
    skip_if(!rmarkdown::pandoc_available(), "pandoc not available")

    deg <- .make_merged_deg(n_genes = 25)
    deg$padj[1:4] <- 1e-5
    deg$log2FoldChange[1:4] <- c(2, -2, 1.8, -1.5)
    diffexp_out <- prepare_diffexp_from_deg(deg)
    gsea_out <- run_gsea_pipeline(
        diffexp_out$results,
        collections = "HALLMARK",
        run_clusterprofiler = FALSE
    )

    out_dir <- withr::local_tempdir()
    paths <- render_report(
        dds = NULL,
        var = "comparison",
        distance_out = NULL,
        diffexp_out = diffexp_out,
        gsea_out = gsea_out,
        wgcna_out = NULL,
        run_args = list(analysis_mode = "deg_table"),
        save_dir = out_dir
    )

    expect_false("qc" %in% names(paths))
    expect_true(all(c("index", "deseq2", "clusterProfiler") %in% names(paths)))
    expect_true(all(file.exists(paths)))

    index_html <- paste(readLines(paths[["index"]]), collapse = "\n")
    expect_true(grepl("color:#999;\">QC</span>", index_html))
    expect_true(grepl("DEG table", index_html, ignore.case = TRUE))
})
