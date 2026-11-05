# Reproduce this DEGreports run (run_deseq2)
# Written 2026-11-05 10:27 by DEGreports 0.4.0.
# Run from the folder that holds index.html; results go to reproduced_run/.

# Load the package
library(DEGreports)

# Load the inputs saved with this run
inputs <- readRDS("reproduce/inputs.rds")

# The DESeqDataSet passed to run_deseq2(), before QC
dds <- inputs$dds

# Run QC, PCA, DESeq2, GSEA and the enabled modules, then render the HTML report
out <- run_deseq2(
    dds,
    var = "isotype",
    design = NULL,
    comparisons = NULL,
    batch = "batch",
    platform = "rnaseq",
    filter_method = "cpm",
    cpm_percentile = 0.05,
    remove_mt = TRUE,
    remove_xy = FALSE,
    remove_vdj = "BCR",
    ranking = "log2FoldChange",
    collections = c("HALLMARK", "GOBP", "GOCC", "GOMF", "KEGG", "REACTOME", "PID"),
    run_gsva = TRUE,
    run_wgcna = TRUE,
    run_decoupler_tf = TRUE,
    run_decoupler_ccc = TRUE,
    run_viper = TRUE,
    render_html = TRUE,
    parallel_render = TRUE,
    render_workers = 50L,
    save_dir = "reproduced_run"
)

# Open the report
# browseURL("reproduced_run/index.html")
