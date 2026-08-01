# deseq2reports

A self-contained, modular R package for RNA-seq differential expression
analysis built on DESeq2 -- quality control, sample distance/PCA,
differential expression across pairwise and one-vs-rest comparisons, gene
set enrichment analysis (GSEA) against MSigDB, optional WGCNA and TF
activity analysis, and a multi-page HTML report -- all through a single
one-command entry point or as standalone modules.

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
(`qc/`, `distance/`, `counts/`, `diffexp/`, `pathway/`) plus `index.html`
and linked section pages (`qc/qc.html`, `diffexp/diffexp.html`, …). Optional
stages add `wgcna/` and `tf_activity/` when `run_wgcna = TRUE` /
`run_decoupler = TRUE`.

See `vignette("deseq2reports")` for the full argument reference and the
standalone module functions (`run_qc()`, `run_distance_pca()`,
`run_diffexp_pipeline()`, `run_gsea_pipeline()`, `run_wgcna_pipeline()`,
`run_decoupler_pipeline()`, `render_report()`).

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
