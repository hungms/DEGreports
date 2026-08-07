# DESeq2Reports

A self-contained, modular R package for RNA-seq differential expression
analysis built on DESeq2 -- quality control, sample distance/PCA,
differential expression (Wald or LRT), GSEA / cluster ORA, optional WGCNA,
TF activity, VIPER (bcellViper), STRING PPI, CSPA surface markers, and a
multi-page HTML report -- through one-command entry points or standalone
modules.

## Install

```r
# install.packages(c("devtools", "BiocManager"))
devtools::install_github("hungms/deseq2reports")
```

Or from a local clone:

```r
devtools::install_local("path/to/deseq2reports")
```

**Reference data** (CSPA workbook + date-stamped CollecTRI networks) ships
under [`inst/extdata/`](https://github.com/hungms/deseq2reports/tree/main/inst/extdata)
and is loaded automatically after install. If a file is missing locally,
the pipeline downloads it from this GitHub repo into the user cache
(`tools::R_user_dir("deseq2reports", "data")`). See
[`inst/extdata/DATA_SOURCES.md`](https://github.com/hungms/deseq2reports/blob/main/inst/extdata/DATA_SOURCES.md).

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

## Citations and resources

Please cite the underlying methods when publishing analyses run with this
package. Key resources:

| Resource | Used for | Citation |
|---|---|---|
| **DESeq2** | Differential expression | Love MI, Huber W, Anders S. *Genome Biol* 2014. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8) |
| **fgsea** / **msigdbr** | GSEA + MSigDB gene sets | Korotkevich et al.; Liberzon et al. *Cell Syst* 2015 ([MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)) |
| **clusterProfiler** | Over-representation analysis (LRT clusters / WGCNA) | Wu T et al. *Innovation* 2021. [doi:10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141) |
| **WGCNA** | Co-expression modules | Langfelder P, Horvath S. *BMC Bioinformatics* 2008. [doi:10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559) |
| **decoupleR** + **CollecTRI** | TF activity (ULM) | Badia-i-Mompel et al. *Bioinformatics Advances* 2022 ([decoupleR](https://doi.org/10.1093/bioadv/vbac016)); Müller-Dott et al. *Nat Commun* 2023 ([CollecTRI](https://doi.org/10.1038/s41467-023-39123-y)). Networks retrieved via [OmniPath](https://omnipathdb.org/) |
| **VIPER** / **bcellViper** | Protein activity / B-cell regulon | Alvarez MJ et al. *Nat Genet* 2016. [doi:10.1038/ng.3593](https://doi.org/10.1038/ng.3593) |
| **STRINGdb** | PPI networks | Szklarczyk D et al. *Nucleic Acids Res* (STRING database) |
| **CSPA** | Cell-surface marker heatmaps | Bausch-Fluck D et al. *PLOS ONE* 2015. [doi:10.1371/journal.pone.0121314](https://doi.org/10.1371/journal.pone.0121314) (supplementary S3 file bundled as `inst/extdata/CSPA_pone.0121314.s003.xlsx`) |
| **DEGreport** | Temporal gene clustering (LRT) | Pantano L. DEGreport. Bioconductor |
| **ashr** / **apeglm** | LFC shrinkage (optional) | Stephens; Zhu et al. |

Bundled CollecTRI snapshots (`collectri_{human,mouse}_2026-08-06.csv`) are
OmniPath extracts dated **2026-08-06**, collapsed to `source` / `target` /
`mor` for ULM. Re-download and replace under `inst/extdata/` when updating
priors.

## License

MIT
