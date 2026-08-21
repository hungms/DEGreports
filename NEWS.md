# deseq2reports 0.3.0

* **Ligand–receptor pairs are now scored per degPatterns cluster in LRT mode.**
  decoupleR-CCC was the last LRT module still working off the whole-test
  signature, so its barplot ranked pairs by a statistic pooled over every
  temporal pattern and could not say which pattern a pair belonged to (while
  decoupleR-TF and VIPER had used per-cluster signatures for some time).
  `run_decoupler_ccc_pipeline()` takes a new `cluster_diffexp` argument — the
  same `lrt_clusters_diffexp_for_tf()` object those two modules consume — and
  runs one ULM per cluster, writing `<lrt>/cluster_<n>/ulm_lr.csv` and
  `ulm_lr_barplot.pdf` next to the existing whole-test outputs; the returned
  list gains `cluster_acts`, `cluster_barplots` and `clusters_by_test`.
  `run_deseq2_lrt()` passes it automatically. Because a cluster signature
  holds only that cluster's genes, `minsize = 2` means a pair is reported for
  a cluster only when **both** its ligand and its receptor are members —
  which is exactly the cluster attribution the whole-test fit could not make,
  and also why a cluster can legitimately score no pair at all (reported as a
  message, not a warning). On the LRT decoupleR-CCC pages the single
  whole-test barplot is replaced by one tab per cluster (barplot + table +
  CSV), with the pooled table kept below under "All pairs (whole test)" as
  the better-powered but pattern-blind view. Wald and DEG-table runs are
  unchanged.
* **Breaking (interactions are now differences of differences).**
  `include_interactions = TRUE` used to mean "compare inside each stratum":
  `K_vs_noK` refitted on the 6 samples of `CC__C`, and separately on `CC__noC`,
  leaving the reader to eyeball whether the two differed. It now generates the
  contrast that *asks* that question. For `var = c("celltype", "C_mutation",
  "K_mutation")` the 12 within-stratum contrasts are replaced by 7 interaction
  contrasts: one per factor pair per level of the third factor
  (`C_vs_noC_x_K_vs_noK_in_CB` — does the K effect differ between C and noC
  within CB, which is the same question as "does the C effect differ between K
  and noK within CB", so a pair gives one contrast and not two), plus the
  three-way `CC_vs_CB_x_C_vs_noC_x_K_vs_noK`. Each is fitted as a sum of
  coefficients from a **single** `DESeq()` of `~ celltype * C_mutation *
  K_mutation` on **all** samples, so the interactions share the dispersion
  estimates and sample size of the main effects instead of being estimated from
  a quarter of the data — and the fit happens once for the whole grid, where
  each two-factor interaction row used to refit the dataset for itself.
  Two-factor interaction labels are renamed accordingly, from
  `interaction__TissueSpleen_x_IsotypeIgG` to `Spleen_vs_BM_x_IgG_vs_IgM`
  (breaking for anyone reading those folders by name).
  `run_diffexp_pipeline()` now fails up front, naming the empty combinations,
  when the factorial has a cell with no samples, rather than letting the
  rank-deficiency check report it generically.
  The within-stratum grid remains implemented and tested behind
  `expand_contrast_specs(stratified = TRUE)` /
  `run_diffexp_pipeline(stratified = TRUE)`, which the pipeline entry points do
  not expose.
* **A factor with 3+ levels is interacted one-vs-rest, not pairwise.** With
  `celltype = CA/CB/CC`, the interaction grid now asks
  `CA_vs_CB+CC_x_K_vs_noK_in_C` ("does the K effect differ between CA and the
  average of CB and CC, within C") in place of the `CB_vs_CA`, `CC_vs_CA`,
  `CC_vs_CB` rows, matching how that factor's main effects already read and
  keeping the grid linear in the number of levels rather than quadratic (a
  4-level factor gives 4 rows per stratum, not 6). The pooled side is a *mean*
  of coefficients, each level weighted `1/(k - 1)` however many samples it has,
  which a summed `results(contrast = list(...))` cannot express: such rows
  instead carry a `coef_weights` map (`"name=weight;..."`) in the spec table,
  which `run_diffexp_coef()` turns into a numeric contrast vector over
  `resultsNames()`. Rows whose weights are all ±1 are unchanged and still go
  through the list contrast. Coefficient weights are accumulated rather than
  appended, so a term that cancels between the two sides (as the two-way
  interaction does inside a pooled three-way row) drops out instead of being
  double-counted.
* The DESeq2 hub page adds a **"What each interaction asks"** table, stating
  every interaction contrast as a question ("Does the K_mutation effect
  (`K_vs_noK`) differ between celltype CA and the average of CB and CC, within
  C_mutation = C?") together with what a positive `log2FoldChange` means for it,
  in both of its symmetric readings. Available programmatically as
  `interaction_question_table()`, built from the new `ix_terms` spec column.
* **Breaking (missing main effects):** with two or three grouping factors, the
  contrast grid only ever asked the *stratified* question. `var = c(v1, v2, v3)`
  gave `v1` compared inside each `v2 × v3` cell and nothing else, so the plain
  "does `v1` differ, across all levels of `v2` and `v3`" contrast -- usually the
  first thing anyone wants -- could not be requested at all, and
  `include_interactions` was ignored outright for three factors (it only changed
  the design formula). `expand_contrast_specs()` now always expands each
  factor's **main effect** (new spec `class = "main"`): its pairwise and, with
  ≥ 3 levels, one-vs-rest comparisons across every level of the other factors,
  e.g. `v1l1_vs_v1l2+v1l3` and `v2l1_vs_v2l2`. `include_interactions = TRUE`
  adds the interaction contrasts described below on top; `FALSE` now means main
  effects only, where it used to mean strata. A 3 × 2 × 2 design goes from 36
  contrasts, none of them marginal, to 8 main effects (+ 18 interactions with
  interactions on).
* `run_diffexp_pipeline()` fits a `"main"` row on **all** samples with the
  design rebuilt additively as `~ [batch +] other factors + tested factor`. The
  factors not under test are therefore adjusted for, rather than left in the
  residual (`~ tested factor` alone) or held at their reference levels (which
  is what an interaction design's own main-effect coefficient reports). Labels
  are flat (`deseq2/C_vs_noC/C_vs_noC_diffexp.csv`). Main effects lead the
  report's contrast list, and the DESeq2 / clusterProfiler / decoupleR-CCC hubs
  now head each group of links ("Main effects", "Interactions (difference of
  differences)") through the new `emit_contrast_nav()`, instead of filing
  everything but the strata under "Other".
* Fix: stratum links on the DESeq2 / clusterProfiler / decoupleR-CCC hubs
  rendered as literal `**CB__noK**`. Pandoc read the `__` inside the name as a
  strong emphasis of its own, which swallowed the closing `**` of the one
  around it; the underscores are now escaped.
* **Download PDFs are written to disk instead of embedded.** Every "Download
  PDF" button put its whole PDF into the page as base64 inside an `onclick`
  handler -- 2.6 MB of the ~7 MB of a DESeq2 contrast page went on two copies
  of the MA and volcano PDFs that `run_diffexp_pipeline()` had *already* written
  in the same folder. `download_plot_button()` /
  `download_grid_plot_button()` now write the PDF into the page's
  `<page>_figures/` directory and link it by relative path, the way
  `figure_img_tag()` already did for figures; called outside
  `render_report()` they still embed, so a standalone button stays
  self-contained. Contrast pages drop by roughly a third (measured: 0.74 MB
  where the same page was 1.1 MB), which is that much less for pandoc to
  convert, and the PDFs arrive with readable names
  (`A_vs_B_figures/A_vs_B_MA_plot.pdf`).
* **Fix (one bad page threw away the whole run):** a contrast page whose pandoc
  step failed in a `mclapply()` worker aborted `render_report()`, which aborted
  the pipeline -- past the point where every analysis had already been computed,
  so a ten-minute run ended with no index page, no `dds.rds`, and no `ppi/`
  cleanup because 2 of 15 sub-pages had not converted. Such a page is now
  **retried once, alone, in the main session**: the failures seen in practice
  (`pandoc document conversion failed with error 1`, i.e. an IO exception, on a
  machine with 90 GB free) come from four pandoc processes writing multi-MB
  documents into the same tree at once, not from the page's own content, and the
  retry both clears them and puts pandoc's own message on the console instead of
  losing it with the dead child. A page that fails the retry too now warns and
  leaves a gap in the report rather than taking the rest of it down with it.
* **Breaking:** `batch_method` is gone and limma batch correction is removed
  entirely (`limma` dropped from Imports). Batch correction is now chosen by
  the new `platform` argument on `run_deseq2_wald()` / `run_deseq2_lrt()` /
  `run_distance_pca()` / `run_counts_export()` and the sample-level modules:
  `"rnaseq"` (default) uses `sva::ComBat_seq()` on raw counts then VST,
  `"microarray"` uses `sva::ComBat()` on the VST matrix. Both protect the
  grouping variable. `platform` only affects batch correction -- the DESeq2
  model and gene filtering are unchanged.
* PCA, the QC sample-similarity heatmap, and the sample-level modules
  (WGCNA, GSVA, ULM, VIPER, CSPA) now share one corrected matrix, so the
  report can no longer show one correction while a module used another.
  Previously PCA used limma/ComBat on VST while the modules used ComBat-seq.
  This also halves the work: ComBat(-seq) runs once per pipeline instead of
  once for PCA and again for the similarity heatmap.
* **Breaking:** `run_distance_pca()` renames `similarity_plot_combat_seq` /
  `similarity_matrix_combat_seq` to `similarity_plot_corrected` /
  `similarity_matrix_corrected` and adds `batch_correction` (`"ComBat-seq"`,
  `"ComBat"`, or `"none"`). The saved files become
  `sample_similarity_{heatmap,matrix}_batch_corrected.{pdf,csv}`.
  `run_counts_export()` renames its `combat_seq` element to
  `batch_corrected` and writes `counts_combat_vst.csv` for microarray.
* **Breaking:** the following pipeline arguments were removed and are now
  fixed internally (see `.pipeline_const`): `shrinkage`, `alpha`,
  `lfc_threshold`, `gsva_n_heatmap`, `wgcna_power`,
  `wgcna_min_module_size`, `decoupler_n_tfs`, `decoupler_n_lr`,
  `viper_n_tfs`, `ppi_top_n_hubs`, `cspa_n_markers`. `run_ppi` is removed
  too -- STRING PPI always runs. The module functions still expose all of
  them, so call those directly to deviate. Note
  `wgcna_min_module_size` is fixed at `100` (the value the example runners
  used), not `run_wgcna_pipeline()`'s own default of `30`. The
  "Pipeline arguments" report table shrinks accordingly.
* New exported helpers `batch_correction_label()` and
  `expr_correction_label()` turn a module's `expr_correction` tag into
  prose for the report pages.
* When `batch` is set, WGCNA / GSVA / sample-level ULM and VIPER / CSPA
  heatmaps use the batch-corrected matrix (protecting the grouping
  variable). DESeq2, GSEA, and contrast-level TF/CCC/msVIPER still use raw
  counts with `batch` in the design.
* **Fix (apparent hang after WGCNA):** report figures are written to a
  `<page>_figures/` folder beside each HTML page instead of being inlined as
  base64 `data:` URIs. For a non-self-contained document,
  `rmarkdown:::html_document_base()`'s post-processor pushes every
  `<img src>` through `utils::URLdecode()`, which walks the string one byte
  at a time and regrows its result on each step -- quadratic, so free for a
  file path and ruinous for base64. The WGCNA page carried 17 MB of inline
  PNG across 83 figures (its largest alone costing ~14 min), so the pipeline
  sat silently in that post-processor for over half an hour *after* knitting
  and pandoc had both finished. `wgcna.html` drops from 19 MB to ~2 MB and the
  page renders in seconds. New exported `figure_img_tag()` replaces
  `htmltools::plotTag()` in the templates; outside a render it still falls
  back to inline base64.
* Render progress messages now report each page's size alongside its wall
  time, since size is what drives render cost.
* QC: the top PCA loadings figure is drawn wide enough to read its rotated
  gene labels, and is exempt from the package's square-panel convention -- a
  square panel is capped by the available panel *height*, so the extra width
  used to become whitespace rather than label spacing. Its download/saved PDF
  is portrait (9 × 12 in) rather than a squashed 20 × 5 in strip; new
  `pca_loadings_figure_size()` keeps the report and `run_distance_pca()` in
  agreement, and `emit_plot_with_download()` / `download_plot_button()` /
  `download_grid_plot_button()` gained `pdf_width` / `pdf_height` for figures
  whose on-page shape does not suit a printable page.
* QC: PC–covariate association heatmap uses short/wide cells (normal-aspect
  text, not CSS-stretched). Sample-similarity heatmap display width ×1.3
  (aspect preserved). Download-PDF pages match the on-page HTML box and scale
  the grob uniformly to fit (no cropping, no tiny inset). CSV and PDF
  download buttons share the same `btn-sm` size.
* **Fix (`filter_method = "cpm"` removed nothing):** the percentile was taken
  over the whole mean-CPM distribution and applied as `>= cutoff`. On real bulk
  RNA-seq a full annotation leaves roughly half the genes zero in every sample,
  so the quantile landed inside that tie at zero -- and since
  `log2(0 + 1) == 0`, every gene in the dataset satisfied `>= 0`. The default
  filter was therefore a silent no-op for any dataset with more than
  `cpm_percentile` unexpressed genes: GSE224465 kept all 54,860 genes and the
  report cheerfully reported "removed the bottom 5%". Unexpressed genes are now
  dropped first, then `cpm_percentile` applies to what remains (GSE224465:
  27,169 + 1,419 = 28,588 removed, 26,272 kept). `cpm_percentile = 0` drops only
  the unexpressed genes. The progress message and the QC page now spell out both
  parts, and the cutoff line on the count-distribution plot is a real boundary
  instead of sitting invisibly on the axis at zero.
* Docs: `run_qc()`'s Details said low-expression filtering was "off by default
  (`filter_method = NULL`)" while the signature has defaulted to `"cpm"` since
  0.2.0.
* **Fix (clipped heatmaps):** every `pheatmap` in the report was drawn on a page
  smaller than its own layout, and because `gtable` centres its layout, grid ate
  the figure from *both* ends at once -- the column dendrogram off the top and
  the sample labels off the bottom. Two causes: `plot_gsva_heatmap()` scaled the
  gtable's `heights` by 1.3, which for a `pheatmap` (whose matrix row and column
  are `npc`) makes the layout 1.3× *whatever page it is given*, so a 30% crop was
  guaranteed at every size; and the templates sized the page from
  `sum(gt$heights)`, which for the same reason just reports back whichever
  graphics device happened to be current when the chunk ran. New exported
  `gtable_page_size()` sizes the page from content instead -- fixed furniture
  plus one cell per row and column, with floors for wrapped row labels, the
  stacked legend ticks, and a title wider than the matrix column -- and
  `emit_heatmap_with_download()` draws through it, scaling to fit rather than
  clipping if a cap bites. All five heatmap sites (QC similarity and PC–covariate
  association, GSVA, VIPER, ULM, WGCNA eigengenes) and the saved PDFs
  (`save_pheatmap()`) go through it, so the GSVA heatmap now shows its dendrogram
  and sample names, and each is shown at print size rather than shrunk to 55-70%.
* QC: every before/after tabset (PCA scatter, PC–covariate association,
  sample similarity) and the PCA / scree / loadings plot titles now name the
  method -- **After ComBat-seq** or **After ComBat** -- instead of the generic
  "After batch correction", falling back to the generic wording only when no
  `batch` was supplied. The uncorrected side stays **Before batch correction**
  (it is the same VST matrix whichever method the corrected side used), which
  also retires the jargon tab label "Before (VST)".
* QC: the on-page PCA loadings figure shrinks from 16.4 × 5.5 to 7.2 × 3.6 in, so
  it sits inside the report's content column instead of overhanging it. The
  portrait 9 × 12 in download stays the readable copy of the gene labels.
* QC: dropped the "Adjust via the `cpm_percentile` argument" aside from the
  low-expression filter note.
* Fix: every report page now knits with its own private `intermediates_dir`.
  `rmarkdown::render()` defaulted this to the template's own folder, so
  per-contrast pages rendered with `parallel_render = TRUE` overwrote a single
  shared `report_<module>.knit.md` and every contrast HTML came out
  byte-identical. It also stopped multi-MB `.knit.md` files (19 MB for WGCNA)
  from being written into the installed/source package tree, which stalled
  rendering when that tree lives in a synced folder.
* `render_workers` now defaults to `min(4, detectCores() - 1)`. Forked
  contrast knits inherit the parent heap, and an unbounded fan-out pushed the
  machine into swap so the next hub page (WGCNA) appeared to hang.
* Report rendering drops heavy plot objects from memory after WGCNA / GSVA /
  QC pages so later pages (decoupleR-TF, VIPER, …) start promptly instead of
  stalling on serializing unused ggplot/ plotly handles. Progress messages
  now print `Finished HTML page: … (<elapsed>)` for every hub *and* contrast
  page, so a slow page is distinguishable from a stalled one.
* Top-level `ppi/` (LRT-cluster PPI + STRINGdb cache) is removed automatically
  at pipeline completion via `cleanup_ppi_dir()`; contrast PPI and WGCNA
  module PPI under `wgcna/<Module_N>/ppi/` are kept.
* Download PDF / CSV button rows use shared flex alignment (`.deseq2reports-dl-btn-wrap`)
  so mixed `downloadthis` and link buttons sit on the same baseline.
* Optional parallel HTML rendering: `parallel_render = TRUE` on
  `render_report()` / `run_deseq2_wald()` (via **future.apply**). Optional
  resumable rendering via `use_targets = TRUE` / `render_report_targets()`
  (requires Suggested **targets**).
* Output layout: count exports and sample distance/PCA artifacts live under
  `qc/counts/` and `qc/distance/` (with other QC outputs in `qc/`). WGCNA
  module STRING PPI is written to `wgcna/<Module_N>/ppi/`; LRT temporal-cluster
  PPI stays under `ppi/<cluster>/ppi/` until pipeline cleanup.
* Optional modules (default **off**): `run_gsva`, `run_wgcna`,
  `run_decoupler_tf`, `run_decoupler_ccc`, `run_viper`. STRING PPI and CSPA
  stay **on** (`run_ppi = TRUE`, `run_cspa = TRUE`). `run_deseq2_deg()`
  matches this (decoupler-tf / CCC / VIPER off unless requested; GSEA / PPI /
  CSPA on; GSVA needs counts so it is skipped).
* Optional **GSVA** HTML page (`run_gsva`): sample-level scores on VST
  (ComBat-seq then VST when `batch` is set), independent of clusterProfiler
  GSEA. HALLMARK is always scored and
  heatmapped in full. GO / KEGG / REACTOME / PID are scored in full (CSV)
  but heatmapped as the top 50 pathways associated with the grouping
  variable (Kruskal–Wallis on GSVA scores), not GSEA padj. Heatmaps are
  1.3× taller; pathway names wrap at 40 characters and row text shrinks
  when any plotted name wraps.
* Optional **decoupleR-CCC** module (`run_decoupler_ccc`): LIANA / OmniPath
  ligand–receptor ULM per comparison, hub + per-contrast HTML
  (`decoupler-ccc/decoupleR-CCC.html`), GSEA-style L–R barplots showing
  the **top 10 up and top 10 down** pairs. Distinct
  from CollecTRI TF activity (`run_decoupler_tf`, `run_decoupler_tf_pipeline()`,
  `render_decoupleR_tf_page()`, HTML `decoupler-tf/decoupleR-TF.html`).
  Multihit hit columns / plot factors are `decoupler-tf` and `decoupler-ccc`.
* Three-factor Wald contrasts use **nested folders**:
  `deseq2/C_vs_noC/CB__noK/` (short comparison + stratum levels only),
  with `contrast_file_stem()` for artifact prefixes. 1-/2-factor layouts
  stay flat. Pairwise direction follows factor level order (mutant vs WT).
* QC report: **Sample QC** table (detected genes + total transcripts per
  sample) before Library QC; former Sample QC block renamed to **PCA**.
* Renamed GSEA/ORA collection alias `"KEGG_LEGACY"` to `"KEGG"` (still
  fetched from msigdbr `CP:KEGG_LEGACY`). `"KEGG_LEGACY"` is remapped.
* Dropped `"KEGG_MEDICUS"` from supported GSEA/ORA collections (alias map,
  validation, pipeline defaults, examples, and docs). Use `"KEGG"`
  for KEGG pathways.
* clusterProfiler contrast pages: GSEA collections (**HALLMARK**, **GOBP**,
  …) are Bootstrap tabs again (same pattern as LRT cluster ORA).
* PCA: covariate association heatmap for the top 15 PCs vs `batch` /
  grouping `var` / `time` (LRT) (`pca_covariate_association()`,
  `plot_pca_covariate_heatmap()`), shown **before and after** batch
  correction (30% smaller) when `batch` is set.
* When `batch` is set, DEG / LRT models always include it
  (`~ batch + var` / `~ batch + time`); LRT also stamps that formula onto
  `design(dds)` so the report header matches the fitted models.
* Multihit moved to its own HTML page (`multihit/multihit.html`) after
  VIPER in the nav (no longer embedded on the index).
* DESeq2 hub: per-contrast significant up/down DEG counts and %
  (`summarize_deg_contrasts()`).
* clusterProfiler hub: per-contrast significant up/down pathway counts
  and % (`summarize_gsea_contrasts()`).
* New `run_deseq2_deg()` / `prepare_diffexp_from_deg()`: start from a
  **merged DEG table** with a `comparison` column (e.g.
  `groupA_vs_groupB`) and run GSEA, GSEA overlap, CSPA marker tables,
  contrast-level decoupleR ULM / msVIPER, DEG PPI, multihit ranking, and
  HTML (QC / WGCNA / sample heatmaps skipped; no counts required).
* `run_decoupler_tf_pipeline()` / `run_viper_pipeline()` accept `dds = NULL`
  for contrast-level scoring only.
* `render_report(dds = NULL)` skips the QC page and greys it out in the nav.
* Suggests gains the optional rendering back-ends `future` / `future.apply`
  (`parallel_render = TRUE`) and `targets` (`use_targets = TRUE` /
  `render_report_targets()`), plus `liana` for the ligand–receptor resource
  when no OmniPath download is available and `pkgload` for the test helpers.

# deseq2reports 0.2.0

* Rename `run_deseq2_pipeline()` → `run_deseq2_wald()` (deprecated alias kept).
* New `run_deseq2_lrt()` for time-course LRT analysis:
  - Nested full/reduced model grid via `expand_lrt_specs()`.
  - degPatterns temporal clustering (`run_degpatterns_pipeline()`).
  - Cluster ORA on the clusterProfiler page (ranked GSEA is not used).
  - CSPA top **10** genes per degPatterns cluster (`run_lrt_cspa_heatmap()`).
  - PPI plots per cluster in Bootstrap tabsets on the DESeq2 page.
  - decoupleR / VIPER score **per cluster** (`lrt_clusters_diffexp_for_tf()`);
    contrast-level TF / msVIPER sections are omitted in LRT HTML.
  - Multihit tables per cluster via `run_multihit_lrt_pipeline()`
    (layers: LRT DEG, WGCNA, decoupleR, VIPER) on the Multihit page.
* Loading helpers: `summarize_genes()`, `summarize_genes_dds()`, and
  `load_nfcore()` to collapse transcript/feature IDs to gene symbols
  (including nf-core/rnaseq `.RData` + `tx2gene` TSV).
* Packaged reference data under `inst/extdata/` (GitHub install default):
  - CSPA S3 workbook `CSPA_pone.0121314.s003.xlsx`
  - Date-stamped CollecTRI networks
    `collectri_{human,mouse}_2026-08-06.csv`
  - Auto-resolve via installed package path, else download from this
    GitHub repo into the user cache (`DATA_SOURCES.md`).
* Pipeline modules (WGCNA, decoupleR, VIPER, PPI, CSPA, LRT degPatterns,
  GSEA overlap / simplifyEnrichment, ashr/apeglm shrinkage, report
  widgets) are **Imports** so a normal install pulls them in; only
  `knitr` / `testthat` remain Suggests.
* Example runners: `run_gse189410_wald.R`, `run_gse224465_wald.R`,
  `run_gse189410_lrt.R`.
* Vignettes / README updated for nested 3-factor contrasts, Sample QC /
  PCA association, Multihit page, hub DEG/pathway summaries, and
  required pipeline dependencies.

# deseq2reports (development notes prior to 0.2.0)

* CSPA runs immediately after DESeq2 and writes under `deseq2/cspa/`.
* GSEA overlap writes under `clusterprofiler/<contrast>/gsea_overlap/`.
* DESeq2 and clusterProfiler HTML are split into a hub page plus one HTML
  page per contrast (lower peak memory while rendering).
* CSPA heatmaps on the DESeq2 page are ~65% smaller; per-contrast CSPA sits
  after PPI, and the combined CSPA heatmap is at the bottom of the page.
* GSEA leading-edge overlap no longer has a standalone HTML page. Jaccard
  similarity is clustered with simplifyEnrichment `binary_cut`, plotted
  without pathway name labels, with side word clouds of pathway-name tokens
  (Suggests `simplifyEnrichment` / ComplexHeatmap). The section lives at the
  bottom of the **clusterProfiler** page. The top-5-|NES| barplot was removed.
* Multi-hit scoring: GSEA contributes at most +1 per gene (any significant
  pathway leading edge); CSPA marker selection is an additional hit layer.
* GSEA overlap heatmaps remain capped at the top 40 pathways by `|NES|`
  (`max_heatmap_pathways`; full Jaccard CSV still written for all
  `padj < alpha` pathways).
* Multihit scoring via `run_multihit_pipeline()` / `run_multihit_lrt_pipeline()`
  (Wald: top-10 up/down barplots; LRT: per-cluster tables). Rendered on the
  Multihit HTML page as of 0.2.0.
* HTML report `lib/` assets are written next to each page (`qc/lib`,
  `deseq2/lib`, …) so rmarkdown 2.x can resolve them under `output_dir`.
* Removed unused `plot_sample_correlation()` and the opt-in
  `run_qc(add_correlation_heatmap / correlation_transform)` args. Sample
  correlation in the HTML report remains the Spearman heatmap from
  `plot_sample_similarity()` / `run_distance_pca()`.
* Two-factor mode defaults to an additive design (`~ A + B`) with
  within-factor contrasts only. Set `include_interactions = TRUE` on
  `run_deseq2_wald()` / `run_diffexp_pipeline()` / `resolve_design()` /
  `expand_contrast_specs()` to use `~ A * B` and expand interaction
  contrasts: reference interactions plus pairwise differences among
  non-reference interaction coefficients (e.g. CK vs C × cell type). HTML
  reports list these under **Interaction** as *Vs reference* and
  *Pairwise (non-reference)*.
* CSPA surface-marker heatmaps: `run_cspa_marker_pipeline()` intersects DEGs
  with the Cell Surface Protein Atlas (Bausch-Fluck et al. S3; bundled under
  `inst/extdata/`, human/mouse sheets, all categories), takes the top N
  markers per contrast, and plots row-scaled VST heatmaps across **all**
  samples (plus a combined union heatmap). Enabled by default via
  `run_cspa = TRUE` on `run_deseq2_wald()` (requires Suggests `readxl`).
* HTML reports use folder-bundled delivery (`self_contained = FALSE` +
  shared `lib/`) instead of single-file inlined HTML. Large tables show a
  truncated DT and link to the on-disk CSV via `link_output_file()` rather
  than embedding a second copy with `downloadthis`.
* Default `ranking` is now `"log2FoldChange"` (was `"pxfc"`; still supported).
* PPI network plots default to `max_plot_nodes = 50`.
* Per-contrast and per-module artifacts under `decoupler/`, `viper/`,
  and `wgcna/` (incl. per-module `ppi/` subfolders) are now written into
  subfolders matching the deseq2 contrast labels (raw) or sanitized module
  names (`Module_1/`). LRT-cluster PPI uses `ppi/<cluster>/ppi/`.
  Run-level summaries (sample activities, soft-threshold plots, STRING
  cache, etc.) stay at the module folder root.
* Module output folders are lowercase (`clusterprofiler/`, `decoupler/`)
  while comparison subfolders keep their original labels (e.g.
  `IgA_vs_IgM`). Report page titles/nav labels are unchanged.
* Report pages now show only the minimal `deseq2reports` module call needed
  to reproduce each section (not upstream package internals). New vignette
  `deseq2reports-modules` is a detailed HTML module reference; the intro
  vignette focuses on one-command vs modular usage.
* GSEA barcode panels use a single track per pathway coloured by NES sign
  (red = NES > 0, blue = NES < 0), instead of splitting member genes into
  up/down tracks. Pathway and VIPER barcode plot sizes are halved.
* WGCNA module-gene tables include the full module membership (sorted by
  |kME|), showing the top 10 rows by default instead of truncating HTML to
  100 genes. PPI degree tables likewise list every gene with ≥1 STRING
  edge among seeds/module members (plot still degree-capped).
* GSEA pathway page: under each NES barplot, a gene-enrichment barcode
  panel (`plot_gsea_barcode_panel()`) shows the same top/bottom pathways
  with names on the left, member-gene ticks along the ranked list, and
  padj + NES on the right.
* New VIPER / msVIPER module: `run_viper_pipeline()` scores the human
  B-cell regulon from **bcellViper** via `viper` / `msviper` (sample-level
  VST activities + contrast-level barcodes from DESeq2 Wald `stat`). Mouse
  genes are mapped to human orthologs with **babelgene**. Enabled with
  `run_viper = TRUE` on `run_deseq2_wald()`; writes a separate
  `viper/viper.html` page. Requires Suggests packages `viper`,
  `bcellViper`, and (for mouse) `babelgene`.
* New PPI connectivity analysis: `run_ppi_pipeline()` uses **STRINGdb**
  interactions and BioNERO-style interactive plotly networks (Louvain
  communities, degree-sized nodes, hover any gene, top hubs labelled),
  plus native STRINGdb network PNGs on disk. Enabled by default via
  `run_ppi = TRUE` on `run_deseq2_wald()`. Requires Suggests packages
  `STRINGdb` and `igraph`.
* Sample-similarity heatmap: size `plotTag` so device inches match the
  absolute pheatmap gtable (`width = inches*72`, `pixelratio = dpi/72`),
  fixing the empty band under the section title.
* Volcano y-axis is `-log10(padj)`; x-axis only clips *super*-extreme
  unshrunken |log2FC| outliers (99.5th pct, factor 3).
* GSEA barplot: `geom_col(width = 0.86)`, height −30%, NES-axis
  `expansion(mult = c(0.25, 0.25))`.
* WGCNA HTML page: show the 12 largest modules (full CSVs/PDFs still on
  disk) so self-contained pandoc does not OOM with dozens of modules.
* PCA scatter size increased 25% (interactive `base_size` 344); both axes
  use `expansion(mult = c(0.1, 0.1))`. PCA loadings plot width +20%.
* Volcano x-axis / colouring use unshrunken (MLE) `log2FoldChange`;
  results tables keep shrunk `log2FoldChange` plus
  `log2FoldChange_unshrunken`.
* QC HTML figures (except the sample-similarity heatmap) are half the
  previous display size and pinned to a fixed pixel CSS size so they no
  longer stretch with the browser window; library-size boxplot width still
  scales with sample count. Sample-similarity heatmap rasterization uses
  higher dpi + Retina `pixelratio` for sharper output.
* Volcano (and MA) colouring now defaults to `|log2FC| > 0.5` in addition
  to `padj < alpha`; grey `ns` points are drawn behind up/down points; the
  volcano x-axis uses a robust |log2FC| limit so a handful of extreme
  outliers no longer blow out the scale.
* GSEA: added `PID` collection; barplot height/width scale with the number
  of pathways shown; CSV / tables / barplot stars use fgsea's native BH
  `padj` (no separate Storey `qvalue` recalculation).
* New `remove_vdj` argument on `run_qc()` / `run_deseq2_wald()`
  (`NULL`/`FALSE`, `TRUE`/`"both"`, or `"BCR"` / `"TCR"`) to drop
  immunoglobulin / T-cell receptor gene symbols. Index page now reports
  VDJ removal and the full pipeline argument list.
* Interactive MA / volcano / PCA widgets (`make_interactive(square =
  TRUE)`) now place the colour legend in the right margin
  (`legend.x > 1` paper coordinates) outside a square data panel, rather
  than over the points. Static `plot_ma()` / `plot_volcano()` also set
  `legend.position = "right"` explicitly; `aspect.ratio = 1` continues
  to apply only to the data panel (legend and axis text are excluded).
* Fixed a rendering regression where every figure/table's corner download
  button(s), and sometimes the figure/widget itself, could render as
  literal escaped HTML text inside a `<pre><code>` block instead of live
  markup. `emit_widget()` now renders with
  `htmltools::doRenderTags(..., indent = FALSE)`.
* Removed the external `mathjax.rstudio.com` `<script>` dependency;
  each template's YAML header sets `mathjax: null`.
* Organism (human/mouse) is now auto-detected from gene symbols via
  `detect_organism()`, so `run_gsea_pipeline()` / `run_deseq2_wald()` no
  longer take an `org` argument.
* `run_gsea_pipeline()` / `collections` supports multiple MSigDB
  collections per run: `"HALLMARK"`, `"GOBP"`, `"GOCC"`, `"GOMF"`,
  `"KEGG"`, `"REACTOME"`, and `"PID"`.
* New WGCNA co-expression module: `run_wgcna_pipeline()` with plotting
  helpers for soft-threshold, module sizes, trait / eigengene heatmaps,
  and ORA dotplots.
* WGCNA modules are labelled `"Module 1"`, `"Module 2"`, … by descending
  gene count; leftover bucket is `"Unassigned"` (excluded from ORA).
* New sample-level WGCNA visualization: `plot_module_eigengene_heatmap()`.
* Every figure and table in the HTML report has corner download button(s)
  via `wrap_with_corner_button()` / `emit_plot_with_download()` /
  `emit_widget_with_download()`.
* `plot_ma()` uses log2 mean CPM on the x-axis (full-dataset); results
  include a `mean_cpm` column.
* `render_report()` writes cross-linked HTML pages (index, QC, DESeq2,
  clusterProfiler, WGCNA, decoupleR, VIPER) with `report_nav_links()`.

# deseq2reports 0.1.0

* Initial release: modular RNA-seq DESeq2 pipeline with one-command wrapper
  (`run_deseq2_pipeline()`), standalone QC / distance-PCA / differential
  expression / GSEA modules, CSV outputs, and an automated HTML report
  (`render_report()`).
