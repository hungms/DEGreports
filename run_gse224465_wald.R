# Wald pipeline on example_data/GSE224465_counts.tsv (mouse GC B cells).
# Sample names encode genotype.celltype_replicate, e.g. WT.CB_R1, CK.CC_R3.
#
# Genotype WT/C/K/CK is split into two binary mutation factors so the model is
# the three-way interaction:
#   ~ celltype * C_mutation * K_mutation
# where
#   celltype    = CB, CC
#   C_mutation  = noC (WT, K) / C (C, CK)
#   K_mutation  = noK (WT, C) / K (K, CK)
#
# The replicate number is used as batch: CB and CC of the same genotype and
# replicate (e.g. WT.CB_R1 / WT.CC_R1) come from the same animal/sort.
#
# Usage (Slurm, all three examples): sbatch run_examples.sbatch
# Single run:  <conda seurat5>/bin/Rscript run_gse224465_wald.R

# Project root on NEMO: outputs go to <root>/example_run_*/.
root <- "/nemo/lab/caladod/working/Matthew/git-repos/r-dev/DESeq2Reports"
setwd(root)
# Use the project renv (a snapshot of the seurat5 conda library).
if (!nzchar(Sys.getenv("RENV_PROJECT"))) source(file.path(root, "renv", "activate.R"))
# Render workers = Slurm tasks (falls back to 4 outside Slurm).
n_workers <- as.integer(Sys.getenv("SLURM_NTASKS", "4"))

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
pkgload::load_all(root, export_all = FALSE, quiet = TRUE)

counts_path <- file.path(root, "example_data", "GSE224465_counts.tsv")
out_dir <- file.path(root, "example_run_GSE224465_wald")
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

raw <- utils::read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene_name" %in% colnames(raw)) {
  stop("Expected a `gene_name` column in ", counts_path, call. = FALSE)
}

sample_cols <- setdiff(colnames(raw), c("", "Geneid", "gene_name"))
counts <- as.matrix(raw[, sample_cols, drop = FALSE])
storage.mode(counts) <- "integer"
gene_names <- as.character(raw$gene_name)
# Temporary unique Geneid rownames so DESeqDataSetFromMatrix accepts the matrix;
# collapse to gene_name symbols via summarize_genes_dds() below.
rownames(counts) <- make.unique(as.character(raw$Geneid))

parse_sample <- function(nm) {
  parts <- strsplit(nm, "\\.", perl = TRUE)[[1]]
  if (length(parts) < 2) {
    stop("Cannot parse genotype.celltype from sample name: ", nm, call. = FALSE)
  }
  genotype <- parts[[1]]
  replicate <- sub("^[^_]+_R", "", parts[[2]])
  if (!grepl("^[0-9]+$", replicate)) {
    stop("Cannot parse replicate from sample name: ", nm, call. = FALSE)
  }
  celltype <- sub("_R[0-9]+$", "", parts[[2]])
  if (!celltype %in% c("CB", "CC")) {
    stop("Unexpected celltype in sample name: ", nm, " (got ", celltype, ")",
         call. = FALSE)
  }
  if (!genotype %in% c("WT", "C", "K", "CK")) {
    stop("Unexpected genotype in sample name: ", nm, " (got ", genotype, ")",
         call. = FALSE)
  }
  list(
    genotype = genotype,
    celltype = celltype,
    replicate = replicate,
    # C mutation present in C and CK; K mutation present in K and CK.
    C_mutation = if (genotype %in% c("C", "CK")) "C" else "noC",
    K_mutation = if (genotype %in% c("K", "CK")) "K" else "noK",
    batch = paste0("rep", replicate)
  )
}

meta <- lapply(sample_cols, parse_sample)
batch_levels <- paste0(
  "rep",
  sort(unique(as.integer(vapply(meta, `[[`, character(1), "replicate"))))
)
coldata <- data.frame(
  celltype = factor(
    vapply(meta, `[[`, character(1), "celltype"),
    levels = c("CB", "CC")
  ),
  C_mutation = factor(
    vapply(meta, `[[`, character(1), "C_mutation"),
    levels = c("noC", "C")
  ),
  K_mutation = factor(
    vapply(meta, `[[`, character(1), "K_mutation"),
    levels = c("noK", "K")
  ),
  genotype = factor(
    vapply(meta, `[[`, character(1), "genotype"),
    levels = c("WT", "C", "K", "CK")
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
  design = ~ C_mutation * K_mutation
)
# Sum duplicate gene_name symbols (no-op when already unique).
SummarizedExperiment::rowData(dds)$gene_name <- gene_names
dds <- summarize_genes_dds(dds, gene_sym_col = "gene_name")


out <- run_deseq2_wald(
  dds = dds,
  var = c("celltype", "C_mutation", "K_mutation"),
  batch = "batch",
  include_interactions = TRUE,
  platform = "rnaseq",
  design = NULL,
  comparisons = NULL,
  filter_method = "cpm",
  cpm_percentile = 0.05,
  remove_mt = TRUE,
  remove_xy = FALSE,
  remove_vdj = "BCR",
  ranking = "log2FoldChange",
  collections = c("HALLMARK", "GOBP", "GOCC", "GOMF",
                  "KEGG", "REACTOME", "PID"),
  run_gsva = TRUE,
  run_wgcna = TRUE,
  run_decoupler_tf = TRUE,
  run_decoupler_ccc = TRUE,
  run_viper = TRUE,
  render_html = TRUE,
  parallel_render = TRUE,
  render_workers = n_workers,
  save_dir = out_dir
)

cat("\nWald pipeline finished.\n")
cat("Report index: ", out$report_path[["index"]], "\n", sep = "")
