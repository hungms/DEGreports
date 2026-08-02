# Re-run deseq2reports on the deseq2pip BM IgA/IgG/IgM test dataset.
# Counts: deseq2pip/tests/rnaseq/pipeline/qc_results/dds_counts.txt
# Metadata taken from the companion dds_qc.rds (same samples/genes).
#
# Usage:
#   R_LIBS_USER="$HOME/Library/R/arm64/4.5/library" Rscript --vanilla run_bm_ig_pipeline.R

Sys.setenv(R_LIBS_USER = path.expand("~/Library/R/arm64/4.5/library"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

suppressPackageStartupMessages({
  library(DESeq2)
})
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

counts_path <- "/Users/hungm/The Francis Crick Dropbox/CaladoD/Matthew/git-repos/softwares/RNAseq/deseq2pip/tests/rnaseq/pipeline/qc_results/dds_counts.txt"
rds_path <- "/Users/hungm/The Francis Crick Dropbox/CaladoD/Matthew/git-repos/softwares/RNAseq/deseq2pip/tests/rnaseq/pipeline/qc_results/dds_qc.rds"
out_dir <- file.path(getwd(), "example_run")

message("Reading counts: ", counts_path)
counts <- as.matrix(utils::read.delim(counts_path, check.names = FALSE, row.names = 1))
storage.mode(counts) <- "integer"

message("Reading sample metadata from companion RDS: ", rds_path)
dds_ref <- readRDS(rds_path)
stopifnot(identical(colnames(counts), colnames(dds_ref)))
stopifnot(identical(rownames(counts), rownames(dds_ref)))

coldata <- as.data.frame(SummarizedExperiment::colData(dds_ref))
coldata <- coldata[, setdiff(colnames(coldata), "sizeFactor"), drop = FALSE]
coldata$Group2 <- factor(as.character(coldata$Group2), levels = c("IgM", "IgG", "IgA"))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~Group2
)

message("Samples: ", ncol(dds), " | Genes: ", nrow(dds))
message("Group2 levels: ", paste(levels(dds$Group2), collapse = ", "))
message("Detected organism preview: ", detect_organism(rownames(dds)))
message("Writing outputs to: ", out_dir)

out <- run_deseq2_pipeline(
  dds = dds,
  var = "Group2",
  design = NULL,
  filter_method = "cpm",
  cpm_percentile = 0.05,
  remove_mt = TRUE,
  remove_xy = FALSE,
  remove_vdj = "BCR",
  shrinkage = "ashr",
  lfc_threshold = 0.5,
  ranking = "pxfc",
  collections = c("HALLMARK", "GOBP", "GOCC", "GOMF", "KEGG_MEDICUS", "KEGG_LEGACY", "REACTOME", "PID"),
  run_wgcna = TRUE,
  wgcna_min_module_size = 100,
  run_decoupler = TRUE,
  run_viper = TRUE,
  render_html = TRUE,
  save_dir = out_dir
)

cat("\nPipeline finished.\n")
cat("Report index: ", out$report_path[["index"]], "\n", sep = "")
