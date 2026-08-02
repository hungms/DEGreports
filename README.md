# DESeq2Reports

A self-contained, modular R package for RNA-seq differential expression
analysis built on DESeq2 -- quality control, sample distance/PCA,
differential expression across pairwise and one-vs-rest comparisons, gene
set enrichment analysis (GSEA) against MSigDB, optional WGCNA, TF
activity, VIPER (bcellViper), and STRING PPI analysis, and a multi-page
HTML report -- all through a single one-command entry point or as
standalone modules.

## Install

```r
# install.packages("devtools")
devtools::install_local("path/to/deseq2reports")
```

## Quick start

```r
library(deseq2reports)
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)

out <- run_deseq2_pipeline(
  dds = dds,
  var = "condition",
  save_dir = "results"
)
```

Organism (human/mouse) is auto-detected from gene symbols -- there is no
`org` argument. By default this writes per-module folders under `results/`
(`qc/`, `distance/`, `counts/`, `deseq2/`, `clusterprofiler/`) plus `index.html`
and linked section pages (`qc/qc.html`, `deseq2/deseq2.html`, …). Optional
stages add `wgcna/`, `decoupler/`, and `viper/` when `run_wgcna = TRUE` /
`run_decoupler = TRUE` / `run_viper = TRUE`.

## Documentation

```r
vignette("deseq2reports", package = "deseq2reports")           # one-command + modular usage
vignette("deseq2reports-modules", package = "deseq2reports") # detailed module reference (HTML)
browseVignettes("deseq2reports")
```

Standalone modules: `run_qc()`, `run_distance_pca()`,
`run_diffexp_pipeline()`, `run_gsea_pipeline()`, `run_wgcna_pipeline()`,
`run_decoupler_pipeline()`, `run_viper_pipeline()`, `run_ppi_pipeline()`,
`render_report()`.

## Design notes

* **Self-contained**: only CRAN/Bioconductor dependencies -- no `Remotes:`.
* **Switchable arguments** for everything analysts commonly vary per run:
  count filtering (`filter_method = "cpm"` / `"count"` / `NULL`), batch-
  correction method for PCA, LFC shrinkage estimator
  (`ashr`/`apeglm`/`normal`/`none`), significance thresholds, DEG ranking
  metric (`pxfc`/`padj`/`log2FoldChange`/`stat`), and MSigDB collections.
* **Rank-deficiency guard**: `assert_full_rank_model_matrix()` checks the
  model matrix before fitting, so a comparison that would make the design
  singular fails immediately with an actionable message.
* **Exact-match level merging** for pooled one-vs-rest comparisons (never
  regex/substring matching), avoiding mis-mapping when one group name is a
  substring of another.
* CSV outputs throughout (no TSV), a multi-page self-contained HTML report
  per run, and a plain-text + JSON reproducibility log
  (`logs/run_<timestamp>.{log,json}`).

## License

MIT
