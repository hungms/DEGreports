# LRT pipeline on example_data/GSE189410_counts.txt (mouse BM plasma cells).
# Demo: treat isotype (IgM → IgG → IgA) as the discrete `time` axis
# (time-only LRT; no secondary `var`). For true genotype × time designs,
# pass var = "genotype", time = "day", include_interactions = TRUE.
#
# Usage:
#   mkdir -p example_run_GSE189410_lrt/logs
#   Rscript --vanilla run_gse189410_lrt.R \
#     > example_run_GSE189410_lrt/logs/console.log 2>&1

.libPaths(c("/Users/hungm/Library/R/arm64/4.5/library", .libPaths()))

# Surface failures immediately (no background; no silent continue).
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
counts_path <- file.path(root, "example_data", "GSE189410_counts.txt")
out_dir <- file.path(root, "example_run_GSE189410_lrt")
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

message("Reading counts: ", counts_path)
counts <- as.matrix(
  utils::read.delim(counts_path, check.names = FALSE, row.names = 1,
                    stringsAsFactors = FALSE)
)
storage.mode(counts) <- "integer"
rownames(counts) <- make.unique(rownames(counts))

parse_sample <- function(nm) {
  parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
  if (length(parts) < 3) {
    stop("Cannot parse tissue_isotype_rep from sample name: ", nm, call. = FALSE)
  }
  list(tissue = parts[[1]], isotype = parts[[2]])
}

sample_cols <- colnames(counts)
meta <- lapply(sample_cols, parse_sample)
coldata <- data.frame(
  tissue = factor(
    vapply(meta, `[[`, character(1), "tissue"),
    levels = c("BM", "Spleen")
  ),
  # Ordered levels for the LRT "time" axis (IgM -> IgG -> IgA).
  isotype = factor(
    vapply(meta, `[[`, character(1), "isotype"),
    levels = c("IgM", "IgG", "IgA")
  ),
  row.names = sample_cols,
  check.names = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ isotype
)

message("Samples: ", ncol(dds), " | Genes: ", nrow(dds))
message("isotype levels (time axis): ", paste(levels(dds$isotype), collapse = ", "))
message("Detected organism preview: ", detect_organism(rownames(dds)))
message("Writing outputs to: ", out_dir)

out <- run_deseq2_lrt(
  dds = dds,
  var = NULL,
  time = "isotype",
  include_interactions = FALSE,
  filter_method = "cpm",
  cpm_percentile = 0.05,
  remove_mt = TRUE,
  remove_xy = FALSE,
  remove_vdj = "BCR",
  alpha = 0.05,
  collections = c("HALLMARK", "GOBP", "GOCC", "GOMF", "KEGG_MEDICUS",
                  "KEGG_LEGACY", "REACTOME", "PID"),
  run_wgcna = TRUE,
  wgcna_min_module_size = 100,
  run_decoupler = TRUE,
  run_viper = TRUE,
  run_ppi = TRUE,
  run_cspa = TRUE,
  degpatterns_minc = 15,
  render_html = TRUE,
  save_dir = out_dir
)

cat("\nLRT pipeline finished.\n")
cat("Report index: ", out$report_path[["index"]], "\n", sep = "")
