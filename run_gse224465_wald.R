# Wald pipeline on example_data/GSE224465_counts.tsv (mouse GC B cells).
# Sample names encode genotype.celltype_replicate, e.g. WT.CB_R1, CK.CC_R3.
# Two-factor mode with include_interactions = TRUE: design ~ genotype * celltype
# and the full within- + interaction contrast grid.
#
# Usage:
#   mkdir -p example_run_GSE224465_wald/logs
#   Rscript --vanilla run_gse224465_wald.R \
#     > example_run_GSE224465_wald/logs/console.log 2>&1

.libPaths(c("/Users/hungm/Library/R/arm64/4.5/library", .libPaths()))

options(warn = 1)
options(error = function() {
  cat("\n*** ERROR — traceback follows ***\n", file = stderr())
  traceback(2)
  if (!interactive()) quit(save = "no", status = 1)
})

suppressPackageStartupMessages({
  library(DESeq2)
})
try(unloadNamespace("deseq2reports"), silent = TRUE)
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

root <- getwd()
counts_path <- file.path(root, "example_data", "GSE224465_counts.tsv")
out_dir <- file.path(root, "example_run_GSE224465_wald")
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

message("Reading counts: ", counts_path)
raw <- utils::read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene_name" %in% colnames(raw)) {
  stop("Expected a `gene_name` column in ", counts_path, call. = FALSE)
}

sample_cols <- setdiff(colnames(raw), c("", "Geneid", "gene_name"))
counts <- as.matrix(raw[, sample_cols, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- make.unique(as.character(raw$gene_name))

parse_sample <- function(nm) {
  parts <- strsplit(nm, "\\.", perl = TRUE)[[1]]
  if (length(parts) < 2) {
    stop("Cannot parse genotype.celltype from sample name: ", nm, call. = FALSE)
  }
  genotype <- parts[[1]]
  celltype <- sub("_R[0-9]+$", "", parts[[2]])
  if (!celltype %in% c("CB", "CC")) {
    stop("Unexpected celltype in sample name: ", nm, " (got ", celltype, ")",
         call. = FALSE)
  }
  if (!genotype %in% c("WT", "C", "K", "CK")) {
    stop("Unexpected genotype in sample name: ", nm, " (got ", genotype, ")",
         call. = FALSE)
  }
  list(genotype = genotype, celltype = celltype)
}

meta <- lapply(sample_cols, parse_sample)
coldata <- data.frame(
  genotype = factor(
    vapply(meta, `[[`, character(1), "genotype"),
    levels = c("WT", "C", "K", "CK")
  ),
  celltype = factor(
    vapply(meta, `[[`, character(1), "celltype"),
    levels = c("CB", "CC")
  ),
  row.names = sample_cols,
  check.names = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ genotype * celltype
)

message("Samples: ", ncol(dds), " | Genes: ", nrow(dds),
        " (", sum(duplicated(as.character(raw$gene_name))),
        " duplicate gene_name(s) made unique)")
message("genotype levels: ", paste(levels(dds$genotype), collapse = ", "))
message("celltype levels: ", paste(levels(dds$celltype), collapse = ", "))
message("Detected organism preview: ", detect_organism(rownames(dds)))
message("Writing outputs to: ", out_dir)

out <- run_deseq2_wald(
  dds = dds,
  var = c("genotype", "celltype"),
  include_interactions = TRUE,
  design = NULL,
  comparisons = NULL,
  filter_method = "cpm",
  cpm_percentile = 0.05,
  remove_mt = TRUE,
  remove_xy = FALSE,
  remove_vdj = "BCR",
  shrinkage = "ashr",
  alpha = 0.05,
  lfc_threshold = 0.5,
  ranking = "log2FoldChange",
  collections = c("HALLMARK", "GOBP", "GOCC", "GOMF", "KEGG_MEDICUS",
                  "KEGG_LEGACY", "REACTOME", "PID"),
  run_wgcna = TRUE,
  wgcna_min_module_size = 100,
  run_decoupler = TRUE,
  run_viper = TRUE,
  run_ppi = TRUE,
  run_cspa = TRUE,
  render_html = TRUE,
  save_dir = out_dir
)

cat("\nWald pipeline finished.\n")
cat("Report index: ", out$report_path[["index"]], "\n", sep = "")
