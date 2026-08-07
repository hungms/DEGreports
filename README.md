# DESeq2Reports

A self-contained, modular R package for RNA-seq differential expression
analysis built on DESeq2 -- quality control, sample distance/PCA,
differential expression (Wald or LRT), GSEA / cluster ORA, optional WGCNA,
TF activity, VIPER (bcellViper), STRING PPI, CSPA surface markers, and a
multi-page HTML report -- through one-command entry points or standalone
modules.

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

# Group contrasts (Wald)
out <- run_deseq2_wald(
  dds = dds,
  var = "condition",
  save_dir = "results"
)

# Time-course / ordered trajectory (LRT + degPatterns)
# out <- run_deseq2_lrt(dds, var = NULL, time = "day", save_dir = "results_lrt")
```

(`run_deseq2_pipeline()` is a deprecated alias of `run_deseq2_wald()`.)

Organism (human/mouse) is auto-detected from gene symbols -- there is no
`org` argument. By default this writes per-module folders under `results/`
(`qc/`, `distance/`, `counts/`, `deseq2/`, `clusterprofiler/`) plus `index.html`
and linked section pages. Optional stages add `wgcna/`, `decoupler/`,
`viper/`, `ppi/`, and `multihit/` when enabled.

From nf-core/rnaseq outputs:

```r
dds <- load_nfcore("path/to/dds.RData", "path/to/tx2gene.tsv")
```

## Documentation

```r
vignette("deseq2reports", package = "deseq2reports")                 # one-command + modular usage
vignette("deseq2reports-modules", package = "deseq2reports")         # module reference
vignette("deseq2reports-analysis-flow", package = "deseq2reports")   # calculation flow
browseVignettes("deseq2reports")
```

Example runners in the source tree: `run_gse189410_wald.R`,
`run_gse224465_wald.R`, `run_gse189410_lrt.R`.

## Design notes

* **Self-contained**: only CRAN/Bioconductor dependencies -- no `Remotes:`.
* **Switchable arguments** for everything analysts commonly vary per run:
  count filtering (`filter_method = "cpm"` / `"count"` / `NULL`), batch-
  correction method for PCA, LFC shrinkage estimator
  (`ashr`/`apeglm`/`normal`/`none`), significance thresholds, DEG ranking
  metric (`log2FoldChange`/`pxfc`/`padj`/`stat`), and MSigDB collections.
* **Rank-deficiency guard**: `assert_full_rank_model_matrix()` checks the
  model matrix before fitting.
* **Exact-match level merging** for pooled one-vs-rest comparisons.
* CSV outputs throughout, multi-page HTML reports with relative `lib/`
  assets, and a plain-text + JSON reproducibility log.

## License

MIT
