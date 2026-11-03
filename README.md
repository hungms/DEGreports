# DESeq2Reports

A self-contained, modular R package for RNA-seq differential expression
analysis built on DESeq2 -- quality control, sample distance/PCA,
differential expression (Wald or LRT; 1–3 factors), GSEA / cluster ORA,
STRING PPI, CSPA surface markers, optional GSVA / WGCNA / decoupler-tf /
decoupler-ccc / VIPER, cross-analysis gene hits, and a multi-page HTML
report -- through one-command entry points or standalone modules.

## Install

Prefer Bioconductor’s installer so CRAN **and** Bioconductor Imports
(WGCNA, VIPER, STRINGdb, DEGreport, clusterProfiler, …) are resolved
together:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hungms/deseq2reports", dependencies = TRUE)
```

Or with remotes / devtools (same idea — set `dependencies = TRUE`):

```r
# install.packages(c("devtools", "BiocManager"))
devtools::install_github("hungms/deseq2reports", dependencies = TRUE)
```

Or from a local clone:

```r
devtools::install_local("path/to/deseq2reports", dependencies = TRUE)
```

**System requirement:** the HTML report is knitted with pandoc
(**>= 1.12.3**), which is not an R package and so is not installed by any of
the commands above. RStudio bundles one, but an HPC batch job, an `Rscript`
invocation, or a container does not see it — check with
`rmarkdown::pandoc_available("1.12.3")`. Install it system-wide
(`module load pandoc`, `conda install -c conda-forge pandoc`, or
[pandoc.org](https://pandoc.org/installing.html)); with no admin rights,

```r
install.packages(c("pandoc", "gh"))
pandoc::pandoc_install()
```

puts a binary in your user data directory, which this package then activates
by itself in later sessions. Without pandoc every other stage still runs and
writes its tables and figures; only the report is skipped, with a warning in
the first seconds of the run.

Pipeline stages (including optional-toggle modules GSVA, WGCNA, decoupler-tf,
decoupler-ccc, and VIPER) declare their R packages as **Imports**, so they are
installed with the package. Suggests holds only what an analysis can run
without: vignette/test helpers (`knitr`, `testthat`, `pkgload`), the optional
rendering back-ends `future` / `future.apply` (`parallel_render = TRUE`) and
`targets` (`use_targets = TRUE`), and `liana`, whose Consensus resource
decoupler-ccc prefers for ligand–receptor pairs but falls back to OmniPath
without.
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

# Three-factor Wald (main effects + interactions, e.g.
# deseq2/C_vs_noC_x_K_vs_noK_in_CB/)
# out <- run_deseq2_wald(
#   dds, var = c("celltype", "C_mutation", "K_mutation"),
#   include_interactions = TRUE, save_dir = "results_3factor"
# )

# From a precomputed merged DEG table (no counts / no DESeqDataSet)
# out <- run_deseq2_deg(
#   deg_table = "all_contrasts_deg.csv",  # must have gene, comparison, log2FC, padj
#   save_dir = "results_from_deg"
# )
```

(`run_deseq2_pipeline()` is a deprecated alias of `run_deseq2_wald()`.)

**DEG table format** for `run_deseq2_deg()` / `prepare_diffexp_from_deg()`:
one stacked table with a `comparison` column distinguishing contrasts
(e.g. `groupA_vs_groupB`, `groupA_vs_groupC`). Required columns: `gene`,
`comparison`, `log2FoldChange`, and `padj` (or `pvalue`). Strongly
recommended: `stat` (for TF / optional GSEA ranking). Common aliases
(`log2FC`, `FDR`, `contrast`, …) are remapped automatically.

### Choosing the comparisons

Every run writes `deseq2/contrasts.csv` (and shows it at the top of the DESeq2
page): one row per contrast, with the question it asks, the design formula
fitted for it, the samples behind it, and the coefficients it combines.

```r
# what would run, before running it
specs <- expand_contrast_specs(dds, var = c("Tissue", "Isotype"))
contrast_manifest(specs, dds, var = c("Tissue", "Isotype"))

# a subset of the automatic grid
run_deseq2_wald(dds, var = c("Tissue", "Isotype"),
                comparisons = c("SP_vs_BM", "IgG_vs_IgM"), save_dir = "results")

# a contrast the grid does not generate
grid <- add_contrasts(
  specs,
  contrast_spec("earlyIso_vs_lateIso", test_var = "Isotype",
                comparison = "IgM+IgD_vs_IgG+IgA"),
  contrast_spec_from_coefs("SP_specific_IgG", weights = c(TissueSP.IsotypeIgG = 1))
)
run_deseq2_wald(dds, var = c("Tissue", "Isotype"), comparisons = grid,
                include_interactions = TRUE, save_dir = "results_custom")
```

`stratified = TRUE` adds within-stratum contrasts (each factor compared inside
one fixed level combination of the others, fitted on that subset alone), and
`opts = deseq2reports_opts(alpha = 0.01, ...)` changes the statistical and
per-module settings without dropping to module-by-module calls. See
`vignette("deseq2reports")`, section *Customising the comparisons*.

Organism (human/mouse) is auto-detected from gene symbols -- there is no
`org` argument. By default this writes per-module folders under `results/`
(`qc/` incl. `counts/` and `distance/`, `deseq2/`, `clusterprofiler/`) plus
`index.html` and linked section pages. Optional stages add `gsva/`, `wgcna/`,
`decoupler-tf/`, `decoupler-ccc/`, and `viper/` when enabled. PPI and CSPA
run by default (`deseq2/<contrast>/ppi/`, `deseq2/<contrast>/cspa/`;
module PPI under `wgcna/<Module_N>/ppi/` (LRT-cluster PPI under `ppi/`).

From nf-core/rnaseq outputs:

```r
dds <- load_nfcore("path/to/dds.RData", "path/to/tx2gene.tsv")
```

## HTML report

Nav: **Index → QC → DESeq2 → clusterProfiler → GSVA → WGCNA → decoupleR-TF → decoupleR-CCC → VIPER → Multihit**.

| Page | Highlights |
|---|---|
| **QC** | Sample QC (detected genes / total transcripts) → Library QC → PCA (incl. PC–covariate association) → similarity |
| **DESeq2** | Hub summary of significant up/down DEGs; per-contrast MA / volcano / tables |
| **clusterProfiler** | Hub summary of significant up/down pathways; per-contrast GSEA (or LRT ORA) |
| **GSVA** / **WGCNA** | Sample-level pathway scores; co-expression modules, eigengenes, module ORA / PPI |
| **decoupleR-TF** / **VIPER** | TF and regulator activity, per contrast or per LRT cluster |
| **decoupleR-CCC** | Ligand–receptor ULM; in LRT mode one tab per degPatterns cluster plus the pooled whole-test table |
| **Multihit** | Cross-analysis gene ranking (own page after VIPER) |

Multi-factor Wald runs always test each factor's main effect across all levels
of the others, in flat folders (`deseq2/C_vs_noC/`). With
`include_interactions = TRUE`, each pair of factors additionally gets a
difference-of-differences contrast — "does the K effect differ between C and
noC, within CB" — in `deseq2/C_vs_noC_x_K_vs_noK_in_CB/`, plus the three-way
`deseq2/CC_vs_CB_x_C_vs_noC_x_K_vs_noK/`. All of them come from one
`~ celltype * C_mutation * K_mutation` fit on every sample, and the hub page
tabulates the question each one answers.

Rendering a page per contrast is the slowest part of a large run, so
`render_report()` takes `parallel_render = TRUE` to knit contrast pages
concurrently across `render_workers` processes (**future.apply**). The Wald and
LRT entry points also accept `use_targets = TRUE`, which renders through
`render_report_targets()`: hub pages become **targets** so a re-render skips the
ones whose inputs have not changed, with contrast pages following afterwards.
Every page links its assets, figures and CSVs by relative path, so
the whole results tree can be moved or copied and still opened offline -- just
keep it together, since a page needs its sibling `lib/` and `<page>_figures/`
folders.

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
  Pipeline module packages are **Imports** (installed with the package).
* **Switchable arguments** for everything analysts commonly vary per run:
  count filtering (`filter_method = "cpm"` / `"count"` / `NULL`), assay
  platform (`platform = "rnaseq"` / `"microarray"`, which picks ComBat-seq
  or ComBat when `batch` is set), DEG ranking metric
  (`log2FoldChange`/`pxfc`/`padj`/`stat`), and MSigDB collections.
  Statistical and plotting thresholds -- LFC shrinkage, `alpha`, the LFC
  cutoff, WGCNA module size, and hub/TF/marker counts -- are fixed at the
  pipeline level so runs stay comparable; call the module functions
  directly to override them.
* **1–3 grouping factors** via `var`; each factor's main effect is always
  tested, and `include_interactions = TRUE` adds difference-of-differences
  interaction contrasts, each a combination of coefficients from one
  product-design fit. A factor with three or more levels is compared
  one-vs-rest (`CA_vs_CB+CC_x_K_vs_noK_in_C`), the pooled side weighting each
  level equally rather than by sample count.
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
| **GSVA** | Sample-level pathway scores | Hänzelmann S, Castelo R, Guinney J. *BMC Bioinformatics* 2013. [doi:10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7) |
| **LIANA** / **OmniPath** | Ligand–receptor pairs (decoupler-ccc) | Dimitrov D et al. *Nat Commun* 2022 ([LIANA](https://doi.org/10.1038/s41467-022-30755-0)); Türei D et al. *Nat Methods* 2021 ([OmniPath](https://doi.org/10.1038/s41592-021-01279-0)) |
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
