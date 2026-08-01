# deseq2reports (development version)

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
  of pathways shown; CSV / tables include Storey `qvalue` (via the
  Suggests-only `qvalue` package, with BH fallback); barplot significance
  stars use `qvalue`.
* New `remove_vdj` argument on `run_qc()` / `run_deseq2_pipeline()`
  (`NULL`/`FALSE`, `TRUE`/`"both"`, or `"BCR"` / `"TCR"`) to drop
  immunoglobulin / T-cell receptor gene symbols. Index page now reports
  VDJ removal and the full pipeline argument list.
* Interactive MA / volcano / PCA widgets (`make_interactive(square =
  TRUE)`) now place the colour legend in the right margin
  (`legend.x > 1` paper coordinates) outside a square data panel, rather
  than over the points. The previous approach of shrinking
  `xaxis.domain` and pinning `legend.x = domain[2]` looked correct in
  `plotly_build()` but still rendered the legend inside the axes in the
  browser. Static `plot_ma()` / `plot_volcano()` also set
  `legend.position = "right"` explicitly; `aspect.ratio = 1` continues
  to apply only to the data panel (legend and axis text are excluded).
* Fixed a rendering regression where every figure/table's corner download
  button(s), and sometimes the figure/widget itself, could render as
  literal escaped HTML text inside a `<pre><code>` block instead of live
  markup. Root cause: `emit_widget()` serialized its htmltools tags with
  the default pretty-printer, which indents each nested `<div>` by two
  spaces per nesting level; `wrap_with_corner_button()`'s three levels of
  nesting pushed lines to four or more leading spaces, which Pandoc's
  Markdown parser interprets as an indented code block. `emit_widget()`
  now renders with `htmltools::doRenderTags(..., indent = FALSE)`, which
  never indents regardless of nesting depth.
* Removed the external `mathjax.rstudio.com` `<script>` dependency that
  `self_contained: true` alone does not strip from `rmarkdown::html_document`
  output. None of the 5 report pages use LaTeX math, so each template's
  YAML header now sets `mathjax: null`, making every rendered page truly
  offline-viewable with no remote script to be blocked by a
  network-restricted browser/proxy (which could otherwise present as
  broken-looking navigation between pages).
* Organism (human/mouse) is now auto-detected from gene symbols via the new
  `detect_organism()`, so `run_gsea_pipeline()`/`run_deseq2_pipeline()` no
  longer take an `org` argument.
* `run_gsea_pipeline()`/`collections` now supports multiple MSigDB
  collections per run: `"HALLMARK"` (renamed from `"H"`), `"GOBP"`,
  `"GOCC"`, `"GOMF"`, `"KEGG_MEDICUS"`, and `"REACTOME"`. GSEA runs
  independently for every requested collection against every comparison,
  and `run_deseq2_pipeline()`'s default now runs all six. GSEA barplot
  titles are now just the comparison name (the collection is already the
  surrounding tab/section header, so the old `"<comparison> -- <collection>"`
  title was redundant).
* New WGCNA co-expression module: `run_wgcna_pipeline()` (soft-threshold
  power selection via `pick_wgcna_power()`, module detection, module-trait
  correlation, hub-gene identification, and per-module over-representation
  analysis against `GOBP`/`GOCC`/`GOMF`/`KEGG_MEDICUS`/`REACTOME` by
  default), with plotting helpers `plot_soft_threshold()`,
  `plot_module_sizes()`, `plot_module_trait_heatmap()`,
  `plot_module_eigengene_heatmap()`, and `plot_module_ora_dotplot()`.
* `plot_pca()`, `plot_ma()`, and `plot_volcano()` now use tightened axis
  expansion (1.5% instead of 3-12%) so the dotted/dashed zero-reference
  lines touch the panel borders and dead space around the data is
  minimized. Panels remain square (`aspect.ratio = 1`). The `plot_ma()`
  up/down count annotations are now anchored at fixed relative (`Inf`)
  positions instead of the data range, so they no longer force extra
  headroom; `plot_volcano()` keeps a slightly larger top expansion (10%)
  so its count annotations clear repelled gene-symbol labels near the
  most significant points.
  Enabled in the pipeline via `run_wgcna = TRUE` and shown in the HTML
  report as a new "Co-expression modules (WGCNA)" section.
* WGCNA modules are now labelled `"Module 1"`, `"Module 2"`, ... by
  descending gene count (instead of WGCNA's internal colour names), with
  the leftover/unassigned bucket (WGCNA's `"grey"`) relabelled
  `"Unassigned"` and excluded from the numbering; `module_df` retains the
  original colour in a new `module_color` column (used by
  `plot_module_sizes()` for bar fill). `"Unassigned"` continues to be
  excluded from per-module over-representation analysis across every
  collection, since it is a heterogeneous leftover of genes WGCNA could
  not place in any real co-expression module, not a module with a
  coherent expression profile.
* New sample-level WGCNA visualization: `plot_module_eigengene_heatmap()`
  draws a `pheatmap` of each module's eigengene value per sample
  (modules x samples, both axes hierarchically clustered), styled like
  the existing sample-similarity heatmap. It complements the existing
  condition-level `plot_module_trait_heatmap()` and is included in
  `run_wgcna_pipeline()`'s output (`$module_eigengene_heatmap`) and the
  HTML report's WGCNA section.
* The report's per-module over-representation section is renamed
  "Per-module pathway enrichment" (from "Per-module HALLMARK
  over-representation") and now runs against every requested MSigDB
  collection (previously HALLMARK only), rendered as a dot plot
  (`plot_module_ora_dotplot()`: x = gene-set size, colour = -log10
  adjusted p-value) instead of a barplot, nested by collection then by
  module.
* Every figure and table in the HTML report now has its download
  button(s) positioned in the top-right corner of the figure/table
  itself (absolutely positioned, `top: 6px; right: 6px`), rather than
  as a separate element below it. New exported helpers:
  `capture_plot_tag()` (rasterizes a static `ggplot`/grid-drawable to a
  self-contained base64 `<img>` tag), `wrap_with_corner_button()` (the
  general-purpose positioning wrapper, with a `center` toggle for
  static figures vs. full-width responsive widgets), and the
  convenience wrappers `emit_plot_with_download()` (static figure + PDF
  button) and `emit_widget_with_download()` (interactive `plotly`/`DT`
  widget + one or more buttons) used throughout the report template.
* New download-button helpers -- `download_table_button()`,
  `download_plot_button()`, `download_grid_plot_button()` (built on
  `downloadthis::download_this()`/`download_file()`) -- now accompany every
  table and static figure in the HTML report with a self-contained
  "Download CSV"/"Download PDF" button, so the fully offline report still
  offers per-table/per-figure downloads.
* `plot_ma()` now plots log2 mean CPM (computed across the full,
  unfiltered dataset) on the x-axis instead of DESeq2's `baseMean`, which
  is comparison-specific under one-vs-rest contrasts and so not
  comparable across comparisons; `run_diffexp()`'s output gained a
  `mean_cpm` column accordingly.
* Polished `plot_pca()`/`plot_ma()`/`plot_volcano()`: tightened axis
  padding below ggplot2's default expansion, added on-plot annotation of
  up-/down-regulated gene counts to both `plot_ma()` and `plot_volcano()`,
  and restyled `plot_volcano()` to match the visual conventions of the
  `EnhancedVolcano` package (point borders, dashed threshold lines,
  symmetric x-axis).
* Removed prose describing visualization/rendering mechanics (gridlines,
  aspect ratios, etc.) from the report template(s); interpretive text
  (hover instructions, significance-star legends) is retained.
* `render_report()` now writes 5 separate, cross-linked HTML pages instead
  of one monolithic file: `index.html` (run summary and links to the
  other 4), `qc.html`, `diffexp.html`, `pathway_enrichment.html`, and
  `wgcna.html` (omitted -- with its nav-bar link greyed out on the other
  4 pages rather than removed -- when `wgcna_out = NULL`). Every page is
  still rendered with `self_contained = TRUE`, so each opens standalone
  offline; the cross-links are plain relative `<a href>`s, so the 5 (or
  4) files just need to stay together in one folder. `render_report()`
  now returns a named character vector of the rendered paths (`"index"`,
  `"qc"`, `"diffexp"`, `"pathway_enrichment"`, `"wgcna"`) instead of a
  single path, and no longer takes a `report_name` argument (each page
  has a fixed output filename). The bundled
  `inst/rmarkdown/report_template.Rmd` is replaced by 5 templates:
  `report_index.Rmd`, `report_qc.Rmd`, `report_diffexp.Rmd`,
  `report_pathway_enrichment.Rmd`, `report_wgcna.Rmd`. New exported
  helper `report_nav_links()` renders the shared navigation bar reused
  identically at the top of all 5 pages.

# deseq2reports 0.1.0

* Initial release: modular RNA-seq DESeq2 pipeline with one-command wrapper
  (`run_deseq2_pipeline()`), standalone QC / distance-PCA / differential
  expression / GSEA modules, CSV outputs, and an automated HTML report
  (`render_report()`).
