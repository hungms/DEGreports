# DEGreports (development version)

## Breaking changes

* `run_deseq2_wald()` and `run_deseq2_lrt()` are merged into
  `run_deseq2(mode = c("wald", "lrt"))`; `run_deseq2_deg()` is now `run_deg()`.
  The old names are removed.
* Arguments that belong to the other mode now raise an error (`design`,
  `comparisons`, `stratified` for Wald; `time`, `degpatterns_minc` for LRT).
* `ranking` defaults to `"log2FoldChange"` (Wald) or `"stat"` (LRT).
* Batch correction for plots and sample-level modules is now
  `limma::removeBatchEffect()` on the VST matrix (log2 expression for
  `run_limma()`), replacing ComBat-seq / ComBat. DE models are unchanged.
* `platform` is removed from every function (it only chose the ComBat
  variant). `sva` is no longer a dependency.
* The batch-corrected count export is now `counts_batch_corrected_vst.csv`
  (log scale), replacing `counts_combat_seq.csv` / `counts_combat_vst.csv`.

## New features

* `run_limma()` runs the pipeline on TPM/FPKM data (limma-trend on
  log2(x + 1)); `expr_to_dds()` wraps a `SummarizedExperiment` or `DGEList`.
* Every run writes `reproduce/reproduce_run.R` and `reproduce/inputs.rds`,
  linked from the index page.
* Multihit UpSet plots show only genes shared by two or more methods.
* **Paired designs.** New `pair` argument on `run_deseq2()` (Wald and LRT)
  and `run_limma()`: adds a blocking term (`~ [batch +] pair + var`) to
  every model and removes it from the sample-level matrix. Errors when a
  tested factor is constant within every unit, or when `batch` is nested in
  `pair`; warns about units seen in only one condition.
* `batch` accepts several columns in every module.

# DEGreports 0.4.0

Existing calls and outputs are unchanged, apart from the new
`deseq2/contrasts.csv`.

## New features

* **Contrast manifest.** Every run writes `deseq2/contrasts.csv` (question
  asked, design fitted, sample count, coefficients, output folder) and shows
  it on the DESeq2 and index pages. Built by `contrast_manifest()`.
* **Custom contrasts.** `contrast_spec()`, `contrast_spec_from_coefs()` and
  `add_contrasts()` add contrasts the automatic grid does not generate.
  `run_diffexp_coef()` is exported. A single factor accepts any level
  comparison directly, e.g. `comparisons = "d0+d1_vs_d7+d14"`.
* **More than three factors** are supported in the contrast grid.
* `stratified = TRUE` is available from the Wald entry point.
* **Adjustable settings.** `degreports_opts()` (passed as `opts`) sets
  shrinkage, `alpha`, the LFC cutoff and per-module top-N values. The report
  lists only the settings that differ from the defaults.

## Bug fixes

* Passing a contrast-spec table to `run_diffexp_pipeline()` was rejected.
* Factor levels containing `"_vs_"`, `"+"` or `/` now fail early with a clear
  message instead of producing wrong labels or folders.
* A failure in WGCNA, decoupleR-TF or cluster ORA no longer aborts the run;
  all optional modules now warn and continue.
* LRT progress counter stopped at `[7/9]`.

## Internal

* One stage driver runs all entry points (shared order, progress and
  failure policy).
* Contrast encodings consolidated in `R/contrast_spec.R`.

# DEGreports 0.3.0

## Breaking changes

* **Interactions are differences of differences.** `include_interactions =
  TRUE` now tests whether one factor's effect differs between levels of
  another (e.g. `C_vs_noC_x_K_vs_noK_in_CB`), fitted once on all samples.
  Labels are renamed (`interaction__TissueSpleen_x_IsotypeIgG` →
  `Spleen_vs_BM_x_IgG_vs_IgM`).
* **Main effects always tested.** With 2–3 factors, each factor is compared
  across all levels of the others (`~ [batch +] others + factor`).
  `include_interactions = FALSE` now means main effects only.
* Factors with 3+ levels are interacted one-vs-rest, with each pooled level
  weighted equally.
* `batch_method` is removed. Batch correction follows `platform`:
  `"rnaseq"` uses ComBat-seq, `"microarray"` uses ComBat.
* `run_distance_pca()` and `run_counts_export()` outputs renamed from
  `*_combat_seq` to `*_corrected` / `batch_corrected`.
* Statistical thresholds and module top-N values are fixed at pipeline level
  (`.pipeline_const`); `run_ppi` is removed (PPI always runs).

## New features

* decoupleR-CCC scores ligand–receptor pairs per degPatterns cluster in LRT
  mode, with one tab per cluster.
* DESeq2 hub page explains each interaction contrast as a question
  (`interaction_question_table()`).
* PCA, sample similarity and all sample-level modules share one
  batch-corrected matrix.
* PC–covariate association heatmap (top 15 PCs), before and after batch
  correction.
* `run_deseq2_deg()` runs the pipeline from a merged DEG table.
* Optional modules (default off): GSVA, WGCNA, decoupleR-TF, decoupleR-CCC,
  VIPER. PPI and CSPA default on.
* Parallel page rendering (`parallel_render`) and resumable rendering
  (`use_targets`).
* DESeq2 and clusterProfiler hub pages summarise significant up/down counts
  per contrast.
* Multihit has its own page; QC page gains a Sample QC table.
* `SystemRequirements: pandoc`; the package activates a
  `pandoc::pandoc_install()` copy automatically.

## Bug fixes

* Missing pandoc no longer discards the run; the report is skipped with an
  early warning.
* Older msigdbr versions (no `CP:KEGG_LEGACY`) no longer stop cluster ORA.
* An old plotly no longer crashes the run; static figures are used instead.
  Requires `plotly (>= 4.12.0)`.
* `filter_method = "cpm"` removed nothing on real data. Unexpressed genes are
  now dropped before the percentile cutoff.
* A failed contrast page is retried once and no longer aborts the report.
* Parallel renders no longer overwrite each other's intermediate files.
* Heatmaps were clipped; pages are now sized from content
  (`gtable_page_size()`).
* Report rendering could stall for 30+ minutes on inline images; figures are
  now written to `<page>_figures/`.
* Download PDFs are written to disk instead of embedded (pages ~⅓ smaller).
* Stratum links rendered as literal `**CB__noK**`.

## Report and layout

* `render_workers` defaults to `min(4, detectCores() - 1)`.
* Before/after tabs name the correction method (ComBat-seq / ComBat).
* QC figure sizes and PCA loadings layout improved.
* Count and PCA outputs moved to `qc/counts/` and `qc/distance/`.
* Collection alias `"KEGG_LEGACY"` renamed to `"KEGG"`; `"KEGG_MEDICUS"`
  dropped.

# DEGreports 0.2.0

* `run_deseq2_pipeline()` renamed to `run_deseq2_wald()`.
* New `run_deseq2_lrt()`: time-course LRT with degPatterns clusters, cluster
  ORA, and per-cluster CSPA, PPI, decoupleR, VIPER and multihit.
* `summarize_genes()`, `summarize_genes_dds()` and `load_nfcore()` collapse
  IDs to gene symbols (including nf-core/rnaseq output).
* Bundled reference data in `inst/extdata/` (CSPA, CollecTRI), with download
  fallback.
* All pipeline module packages moved to Imports.

# DEGreports (development notes prior to 0.2.0)

* New modules: WGCNA, VIPER (bcellViper), STRING PPI, CSPA surface markers,
  GSEA leading-edge overlap, multihit ranking.
* Multiple MSigDB collections per run (HALLMARK, GO, KEGG, REACTOME, PID).
* Organism auto-detected from gene symbols (`org` argument removed).
* `remove_vdj` drops BCR/TCR genes.
* Two-factor designs default to additive; `include_interactions = TRUE`
  adds interaction contrasts.
* Default `ranking` is `"log2FoldChange"`.
* HTML report split into hub and per-contrast pages, with shared `lib/`
  assets, download buttons on every figure and table, and no external
  MathJax.
* Plot fixes: volcano/MA colouring and axes, GSEA barcode panels, PCA sizing,
  interactive legends outside the panel.

# DEGreports 0.1.0

* Initial release: one-command DESeq2 pipeline (`run_deseq2_pipeline()`),
  QC, distance/PCA, differential expression and GSEA modules, CSV outputs and
  an HTML report.
