# LRT pipeline on example_data/GSE189410_counts.txt (mouse BM plasma cells).
# Demo: treat isotype (IgM → IgG → IgA) as the discrete `time` axis
# (time-only LRT; no secondary `var`). Paired design: IgA/IgG/IgM from the
# same replicate are matched, so replicate (`batch` column) is passed as
# `pair` (full ~ batch + isotype, reduced ~ batch). For true genotype × time
# designs, pass var = "genotype", time = "day", include_interactions = TRUE.
#
# Usage (Slurm, all examples): sbatch run_examples.sbatch
# Single run:  <conda seurat5>/bin/Rscript run_gse189410_deseq2_lrt.R

# Project root on NEMO: outputs go to <root>/example_run_*/.
root <- "/nemo/lab/caladod/working/Matthew/git-repos/r-dev/DEGreports"
setwd(root)
# Use the project renv (a snapshot of the seurat5 conda library).
if (!nzchar(Sys.getenv("RENV_PROJECT"))) source(file.path(root, "renv", "activate.R"))
# Render workers = Slurm tasks (falls back to 4 outside Slurm).
n_workers <- as.integer(Sys.getenv("SLURM_NTASKS", "4"))

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
try(unloadNamespace("DEGreports"), silent = TRUE)
pkgload::load_all(root, export_all = FALSE, quiet = TRUE)

counts_path <- file.path(root, "example_data", "GSE189410_counts.txt")
out_dir <- file.path(root, "example_run_GSE189410_deseq2_lrt")
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

# Read counts; rownames are gene symbols.
counts <- as.matrix(
  utils::read.delim(counts_path, check.names = FALSE, row.names = 1,
                    stringsAsFactors = FALSE)
)
storage.mode(counts) <- "integer"
gene_names <- as.character(rownames(counts))
# Temporary unique rownames so DESeqDataSetFromMatrix accepts the matrix;
# collapse to gene symbols via summarize_genes_dds() below.
rownames(counts) <- make.unique(gene_names)

parse_sample <- function(nm) {
  parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
  if (length(parts) < 3) {
    stop("Cannot parse tissue_isotype_rep from sample name: ", nm, call. = FALSE)
  }
  rep <- parts[[3]]
  if (!grepl("^[0-9]+$", rep)) {
    stop("Unexpected replicate in sample name: ", nm, " (got ", rep, ")",
         call. = FALSE)
  }
  list(
    tissue = parts[[1]],
    isotype = parts[[2]],
    batch = paste0("rep", rep)
  )
}

sample_cols <- colnames(counts)
meta <- lapply(sample_cols, parse_sample)
batch_levels <- paste0(
  "rep",
  sort(unique(as.integer(sub("^rep", "", vapply(meta, `[[`, character(1), "batch")))))
)
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
  batch = factor(
    vapply(meta, `[[`, character(1), "batch"),
    levels = batch_levels
  ),
  row.names = sample_cols,
  check.names = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ isotype
)
# Sum duplicate gene symbols (no-op when already unique).
SummarizedExperiment::rowData(dds)$gene_name <- gene_names
dds <- summarize_genes_dds(dds, gene_sym_col = "gene_name")

out <- run_deseq2(
  mode = "lrt",
  dds = dds,
  var = NULL,
  time = "isotype",
  pair = "batch",
  include_interactions = FALSE,
  filter_method = "cpm",
  cpm_percentile = 0.05,
  remove_mt = TRUE,
  remove_xy = TRUE,
  remove_vdj = "BCR",
  collections = c("HALLMARK", "GOBP", "GOCC", "GOMF",
                  "KEGG", "REACTOME", "PID"),
  run_gsva = TRUE,
  run_wgcna = TRUE,
  run_decoupler_tf = TRUE,
  run_decoupler_ccc = TRUE,
  run_viper = TRUE,
  degpatterns_minc = 15,
  render_html = TRUE,
  parallel_render = TRUE,
  render_workers = n_workers,
  save_dir = out_dir
)

cat("\nLRT pipeline finished.\n")
cat("Report index: ", out$report_path[["index"]], "\n", sep = "")
