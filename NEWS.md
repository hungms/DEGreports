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
  - Index multihit tables per cluster via `run_multihit_lrt_pipeline()`
    (layers: LRT DEG, WGCNA, decoupleR, VIPER).
* Loading helpers: `summarize_genes()`, `summarize_genes_dds()`, and
  `load_nfcore()` to collapse transcript/feature IDs to gene symbols
  (including nf-core/rnaseq `.RData` + `tx2gene` TSV).
* Local CollecTRI priors: `example_data/collectri_{human,mouse}.csv` are
  preferred by `run_decoupler_pipeline()` before OmniPath download;
  optional `collectri_file=` override.
* Suggests `DEGreport` for temporal clustering.
* Example runners: `run_gse189410_wald.R`, `run_gse224465_wald.R`,
  `run_gse189410_lrt.R`.
* Vignettes updated for Wald / LRT entry points, loading helpers, and
  analysis-flow diagrams.

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
* Index page **Cross-analysis gene hits**: per-contrast stacked barplots of
  the top 10 genes across DEG, GSEA LE, VIPER, decoupleR, and CSPA
  (`run_multihit_pipeline()`).
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
  `ppi/`, and `wgcna/` are now written into subfolders matching the
  deseq2 contrast labels (raw) or sanitized module names (`Module_1/`).
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
  `"KEGG_MEDICUS"`, `"KEGG_LEGACY"`, `"REACTOME"`, and `"PID"`.
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
