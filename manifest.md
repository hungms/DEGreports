# Claude progress log — DEGreports

Running notes for picking work back up in a later session. Newest entries first.
Release-level detail lives in `NEWS.md`; this file records *state, decisions and loose ends*.

---

## Repo / workflow facts

- GitHub: https://github.com/hungms/DEGreports (branch `main`).
- Working copy on NEMO: `/nemo/lab/caladod/working/Matthew/git-repos/r-dev/DEGreports`
  is a git checkout of `main` (since 2026-11-03) — commit and push from here.
  Files on NEMO are `rwxrwx` while the repo stores `644`; if `git status` lists
  every file as modified, run `git config core.fileMode false`.
- `inst/lib/` (rendered htmlwidget deps) and `.DS_Store` are generated junk and
  gitignored; `inst/doc/` (built vignettes) is gitignored too — vignettes live in
  `vignettes/*.Rmd` only.
- `CODE_WALKTHROUGH.md` lists every function in pipeline order (generated from the
  call graph 2026-11-03; regenerate or hand-update if functions move).
  `manifest.md` and `CODE_WALKTHROUGH.md` are in `.Rbuildignore`.
- `gh` 2.103.3 installed in `~/.local/bin/gh` (2026-11-03); `gh auth setup-git` is
  configured, so plain `git push` over HTTPS works. A fine-grained PAT without
  Contents: write on this repo gives a 403 — re-run `gh auth login` (web browser) if so.
- **Running tests:** no R on the default PATH. Use the conda R
  `/nemo/lab/caladod/working/Matthew/.conda/envs/seurat5/bin/Rscript` (R 4.3.2, all
  deps present; DESCRIPTION now requires plotly >= 4.10.3, which it has). Run tests
  on a copy of the tree — some tests write output into the working dir.
- **Example report hosting:** the GSE189410 Wald report is on the orphan branch
  `gh-pages` (built in a scratch dir from `example_run_GSE189410_deseq2_wald/`,
  minus `*.rds` and `logs/`, plus `.nojekyll`). README button links to
  https://hungms.github.io/DEGreports/. Re-run the example → rebuild and push the
  branch by hand.
- **Branching (from 2026-11-05):** `main` = verified releases only; `devel` =
  integration; new work on `feature/<topic>` / `fix/<topic>` off `devel`, merged
  back with `--no-ff`. Merge `devel` → `main` (and tag `vX.Y.Z`) only after tests
  pass, the example report renders, and NEWS + this file are updated. No force-push.
- **Conventions skill:** `~/.claude/skills/r-pipeline-package/` (in
  `/nemo/lab/caladod/working/Matthew/.claude/`) records this package's layout,
  style, README/NEWS templates and release steps for future packages.

## Current version: 0.4.0

### 2026-11-05 — imputed sex, per-variable PCA tabs, DESCRIPTION text (branch `feature/impute-sex`)

- **Sex imputation** (R/sex.R): `impute_sample_sex()` wraps
  `biowrenchR::impute_sex()` (k-means, 2 clusters, on log2 CPM of X-escape + Y
  genes). Always run in `run_qc()`, which is **reordered**: expression filter →
  sex → MT → X/Y → VDJ (user asked: after CPM filter, before remove_xy). Stored as
  `colData(dds)$imputed_sex`; `qc_imputed_sex.csv/.pdf`; QC metadata `$sex`.
  Split trusted only if *both* marker sets differ >= 1 log2 in the expected
  direction (a one-set gap was fooled by noise in a single-sex simulation).
  Fails softly (status `failed`) for non-symbol IDs; biowrenchR "markers disagree"
  is reported as `undetermined`.
- **GSE189410 is single-sex (female):** Xist ~5 log2CPM in every sample, only 2 Y
  genes in the matrix -> `undetermined`, sex not modelled. Expected.
- **Model:** silent stage `model` after QC calls `.model_covariates()`: adds
  `imputed_sex` to `c(batch, pair)` only if the additive (and product, with
  interactions) design stays full rank. Later stages use `st$model$covars`, so sex
  is also removed by removeBatchEffect. Stratified subset fits drop covariates not
  estimable in the subset (`.estimable_covariates()`). Decision shown in the report
  args (`imputed_sex`) and the QC sex section.
- **PCA:** `run_distance_pca(colour_vars=)` -> `pca_by_uncorrected` /
  `pca_by_corrected`; QC page: Before tab then After tab, each a Bootstrap tabset
  Group / Pair / Sex / Batch (`emit_pca_tabs()`).
- **DESCRIPTION:** Title + Description rewritten to match README; `biowrenchR` in
  Imports with `Remotes: hungms/biowrenchR` (breaks the "no Remotes" convention by
  user request). Installed into renv with `remotes::install_github(dependencies =
  FALSE)` (renv::install failed upgrading dbplyr); recorded in renv.lock.
- `remove_xy = TRUE` without `gene_chr` now uses `biowrenchR::xy_genes(org)`
  (before, it warned and skipped). User set `remove_xy = TRUE` and reduced Slurm
  resources (30 tasks, 50G, no --nodes) in the runners; committed as-is.
- Tests: new test-sex.R; qc-steps order test updated. Full suite 294 tests, same 7
  pre-existing failures. E2E render on simulated paired two-sex data OK (sex
  correctly not modelled: nested in donor).

### 2026-11-05 — `pair` argument + removeBatchEffect replaces ComBat (branch `feature/pair-design`)

- **`pair`** on `run_deseq2()` (Wald + LRT) and `run_limma()`. Internally
  `covars <- c(batch, pair)` is passed as `batch` to every module, so pair enters
  every design (`~ batch + pair + var`, LRT full and reduced) and the sample-level
  correction. `validate_pair()` (R/validate_args.R): errors if a tested factor is
  constant within every unit (between-subject factor; nested designs not supported),
  or if a `batch` column is constant within every unit; warns for units in one
  condition only.
- **Correction:** user chose to switch *all* ComBat-seq/ComBat to
  `limma::removeBatchEffect()` on VST (sum-to-zero covariates, `design = ~protect`),
  in `rbe_transformed_matrix()` / `batch_corrected_matrix()` (R/distance.R).
  Tags are now `rbe_vst` / `rbe_log2expr`. `platform` removed everywhere; `sva`
  dropped from Imports. Count export: `counts_batch_corrected_vst.csv`.
- `validate_batch()` accepts several columns.
- Tests: new `test-pair.R`; ComBat/platform tests rewritten. Full suite: 283 tests,
  7 failures, all **pre-existing on devel** (test-comparisons x2, test-diffexp-plots
  x2, test-square-panels x2, test-wgcna x1 — plot/ORA-dotplot expectations).
- End-to-end check on simulated paired data: design `~donor + condition`, PCA
  corrected, export written.
- Not done: `CODE_WALKTHROUGH.md` still lists the removed combat_* /
  validate_platform functions (file is untracked).
- PCA shapes points by `pair` when set (unless `shape_var` given); `plot_pca()`
  uses up to 25 distinct shapes, colour only beyond that (merged f931feb).
- Example runners now use `pair = "batch"` (replicate is the matched unit; user
  confirmed). Job 63189435 cancelled; resubmitted on `devel` @ e491a78 as Slurm
  array **63195356** (0 = Wald, 1 = LRT). Logs in `example_run_logs/` (folder had been deleted; recreated).
  Outputs overwrite `example_run_GSE189410_deseq2_{wald,lrt}/`. Runners load the
  working tree, so keep it on `devel` until both finish; then rebuild `gh-pages`.

### 2026-11-05 — examples trimmed, README rewrite, licence, hosted report (pushed)

- Commits `b57521a` and `017ef79` on `main`. `b57521a` also committed the earlier
  "uncommitted" entries below (rename, `run_deseq2(mode=)`, `run_limma()`,
  `reproduce/`, multihit UpSet, parallel-render tmpfile fix, renv/Slurm setup).
- **Examples:** removed GSE224465 and GSE116971 (data + runners). Only
  `run_gse189410_deseq2_wald.R` and `run_gse189410_deseq2_lrt.R` remain;
  `run_examples.sbatch` is now `--array=0-1`. The GSE116971 entry below is
  historical. Local output dirs `example_run_GSE224465_deseq2_wald/` and
  `example_run_GSE116971_deseq2_wald/` were left on disk (gitignored).
- **README:** restructured to Why DEGreports (current limitations → approach →
  additional features → example-report button) / Analyses table / Install (with
  `renv.lock` exact versions) / Quick start / Citations / License. Removed the
  design-notes, documentation and contrast-customisation sections (the vignette
  covers them). User prefers professional, succinct wording.
- **Licence:** added `LICENSE.md` (full MIT text, in `.Rbuildignore`); CRAN
  `LICENSE` stub kept. `inst/extdata/DATA_SOURCES.md` now says CSPA (CC BY 4.0)
  and CollecTRI (OmniPath terms) keep their own licences.
- `renv.lock` committed; `.Rprofile` and `renv/` deliberately not.
- **NEWS.md** condensed from 688 to ~170 lines: one-line bullets under Breaking
  changes / New features / Bug fixes per release. User finds long NEWS hard to read;
  root-cause detail belongs in commits or here.


### 2026-11-04 — fix: random page failures under `parallel_render` (uncommitted)

- Symptom: in every example run ~1 contrast page per module failed with
  `pandoc document conversion failed with error 1`; the sequential retry always
  rescued it. pandoc's stderr: `/tmp/.../rmarkdown-str<hex>.html: withBinaryFile:
  does not exist`.
- Cause: rmarkdown (2.30) writes includes as `tempdir()/rmarkdown-str<hex>.html`
  and `clean_tmpfiles()` deletes **all** files matching that pattern when a render
  ends. `mclapply()` forks share one `tempdir()`, so one worker finishing deleted
  another worker's includes mid-render.
- Fix: `.isolate_rmarkdown_tmpfiles()` (`R/report.R`), called at the start of
  each forked worker, sets rmarkdown's `tmpfile_pattern` to
  `rmarkdown-str-p<pid>-` (via `assignInNamespace`, affects only that child).
  The sequential retry stays as a safety net.
- Verified with a standalone repro on ncpu (64 renders, 32 forks):
  35/64 failed without the fix, 0/64 with it. `test-render-report.R` passes.
- Note: compute nodes cannot see the login node's `/tmp`; put Slurm test files
  on `/nemo`.

### 2026-11-04 — new example: GSE116971 (human ASC pop2/3/5), Wald (uncommitted)

- Data: `example_data/GSE116971_Normalized_ASC_Expression_Data.csv.gz` from GEO
  (17 samples: 6 pop2, 6 pop3, 5 pop5 across 6 patients PB48/51/103/1131/1201/1438;
  PB1131 has no pop5). GEO has **only log10-normalized expression** (confirmed:
  Ig/MT genes at 10^4–10^5, 10^x sums ≈ 2M/sample); raw reads are SRA SRP153114.
- User chose **pseudo-counts**: `round(10^x)` → DESeq2. Approximate p-values
  (library sizes / count noise lost) — stated in the script header. Re-quantify
  from SRA (nf-core/rnaseq) if a definitive analysis is needed.
- User chose design `~ patient + population` (`var = "population"`,
  `batch = "patient"`) — populations compared within individuals.
- `run_gse116971_wald.R` (same pattern as the other examples; human, remove_vdj
  = "BCR", all optional modules on). Added as array task 3 in
  `run_examples.sbatch` (`--array=0-3`). First run: job 62795534.

### 2026-11-04 — multi-hit: explicit hit thresholds + UpSet plots (uncommitted)

Why: users read the CSPA column as "is a surface protein", but it meant "in that
contrast's top-25 CSPA heatmap markers" (e.g. L1cam: rank 21 in IgA_vs_IgM → hit,
rank 37 in IgA_vs_IgG → no hit). Each layer also used a different, undocumented rule.

New hit definitions (Wald + DEG-table multi-hit), one place: exported
`multihit_hit_sets()` in `R/multihit.R`, used by both the ranking and the UpSet:
- DEG: `padj < alpha (0.05)` & `|log2FC| > lfc_threshold (0.5)` (from `opts`).
- CSPA: every DEG above that is a CSPA surface protein (no longer the top-N
  heatmap list). `run_cspa_marker_pipeline()` / DEG-mode CSPA now return
  `cspa_symbols` for this.
- decoupleR-TF: `padj < 0.05` & `|score| > 1` (`tf_min_score`).
- decoupleR-CCC: ligand/receptor of pairs with `padj < 0.05` & `|score| > 1`
  (`ccc_min_score`; `.ccc_sig_genes(min_score =)`).
- VIPER: `padj < 0.05` & `|NES| > 2` (`viper_min_nes`).
- GSEA leading edge: unchanged (not in the user's spec) — any significant pathway.
- The TF/VIPER/CCC cutoffs also apply in LRT multi-hit (shared `.sig_tf_symbols()`
  / `.ccc_sig_genes()`).

Pipeline / report changes:
- `R/decoupler_tf.R`: contrast-level ULM now has `padj` (BH within each
  comparison); written to `ulm_contrast.csv`.
- `run_multihit_pipeline()` gains `lfc_threshold`, `tf_min_score`,
  `ccc_min_score`, `viper_min_nes`; returns `hit_sets`, `upsets`, `thresholds`;
  writes `<contrast>_hit_sets.csv` and `<contrast>_upset.pdf` per contrast.
- New exported `plot_multihit_upset()` (ComplexHeatmap `UpSet`, wrapped as a
  ggplot via `annotation_custom` so existing save/emit helpers work).
- `report_multihit.Rmd`: intro lists the thresholds; new **Method overlap
  (UpSet)** section with the same contrast tabs, placed before **Top genes per
  contrast** (bar plots). Tab layout factored into `emit_mh_sections(emit)`.
- Tests updated for the new rules (`test-multihit.R`, `test-decoupler-ccc.R`) +
  new threshold/UpSet test; both files pass.
- Docs: regenerated with renv's roxygen2 7.3.3, then reverted pure tool noise
  (older knitr/htmltools link targets, NAMESPACE importFrom layout from 8.1.0);
  kept real changes only.
- Verification rerun of all three examples: Slurm job 62793475.

### 2026-11-04 — pinned to the seurat5 conda env + renv + Slurm examples (uncommitted)

- **Target environment = seurat5** (`/nemo/lab/caladod/working/Matthew/.conda/envs/seurat5`,
  R 4.3.2, plotly 4.10.3, ggplot2 3.5.2, DESeq2 1.40.2, msigdbr 7.5.1). Do **not**
  upgrade packages there — adapt the package instead.
- `DESCRIPTION`: `plotly (>= 4.12.0)` → `(>= 4.10.3)`. plotly < 4.10.4 errors in
  `ggplotly()` on any plot with a legend under ggplot2 ≥ 3.5; new internal
  `.ggplotly_compat()` (`R/distance.R`) retries with the ggplot legend hidden and
  plotly's own legend on, so interactive PCA/MA/volcano still work.
- `tests/testthat/helper-simulate.R`: msigdbr call works with both `category =`
  (≤ 7.x) and `collection =` (≥ 10).
- **renv**: `renv.lock` (276 pkgs) hydrated from seurat5 — exact versions, no
  downloads. Cache `/nemo/lab/caladod/working/Matthew/.renv-cache` (set in the
  gitignored `.Renviron`; home quota is tiny). Project library lives in
  `~/.cache/R/renv/library/DESeq2Reports-*/` as symlinks into that cache. Always run
  with the seurat5 `Rscript` (R 4.3.2) from the project dir so `.Rprofile` activates renv.
- **Examples**: `run_gse*.R` now write to `<repo>/example_run_*/`, activate renv,
  use `render_workers = SLURM_NTASKS`. `run_examples.sbatch` = 3-task array on
  `ncpu`, 50 tasks / 100 GB each, `module load Pandoc/2.19` (no pandoc in seurat5).
  First submission: job 62773855 (2026-11-04).
- Tests on seurat5 versions: 941 pass / 7 fail — all 7 stale tests, not versions:
  `test-comparisons.R` ×2 (pre-0.4.0 return type), MA/volcano arrow position ×2
  (arrows moved to the right side), ORA dotplot ×3 (expansion now `c(0.02, 0.12)`;
  fixture lacks `overlap`).

### 2026-11-03 — dead-code / junk cleanup (uncommitted at time of writing)

Removed (all verified unused by parsing R/, report Rmds, vignettes, tests):
- `R/report.R`: the old `future.apply` parallel-render path —
  `.as_isolated_future_fun`, `.parallel_hub_worker`, `.parallel_contrast_worker`,
  `.parallel_render_packages`, `.with_parallel_plan`, `.parallel_render_apply` —
  plus `.degreports_src_dir`. Contrast pages already render via
  `parallel::mclapply()`. Deleted `tests/testthat/test-parallel-render.R`.
- `R/contrast_spec.R`: `.contrast_spec_fields()`; `R/contrasts.R`: `.contrast_up_to_module()`.
- DESCRIPTION: `wordcloud` (Imports; never used), `future.apply` and `pkgload`
  (Suggests; only the dead path used them). `future` kept — `render_report_targets()` uses it.
- Docs that said contrast pages render via `future.apply` now say `mclapply`
  (README, modules vignette, `targets.R`).
- Junk files: committed test-output spill `tests/testthat/group*/`, built
  `vignettes/*.html`, stale `inst/doc/` (main vignette 84 lines behind source).

Tests (conda R 4.3.2 + plotly 4.12.1): before 893 pass / 4 fail / 7 error / 9 skip;
after 887 / 4 / 7 / 9 — identical failure set (−6 passes = deleted parallel test).

Deliberately kept: exported functions not called internally (`render_*_page`,
`summarize_deg_contrasts`, `summarize_gsea_contrasts`, `lrt_sig_diffexp_for_tf`,
`emit_contrast_nav`, `interaction_question_table`) — public API; deprecated alias
`run_deseq2_pipeline()`; `run_gse*.R` example scripts (already `.Rbuildignore`d).

### 2026-11-03 — pushed 0.4.0 working-copy edits to GitHub

Previous remote HEAD: `d560316 fix pandoc` (2026-08-25, v0.3.0).
Local edits dated 2026-11-02/03 were pushed as one commit. Summary (full text in `NEWS.md`):

- **Contrasts are now a first-class object** — new `R/contrast_spec.R`:
  `contrast_spec()`, `contrast_spec_from_coefs()`, `add_contrasts()`,
  `as_contrast_specs()`, `validate_contrast_specs()`, accessors
  `contrast_coefs()`, `contrast_stratum()`, `contrast_test_vars()`,
  `contrast_path()`, `contrast_question()`, `validate_level_names()`.
  The four ad-hoc string encodings are now in one place; the fitting dispatch uses a
  `fit_method` field.
- **Contrast manifest** — new `R/contrast_manifest.R`: `contrast_manifest()`; every
  run writes `deseq2/contrasts.csv` (question, design, n samples, coefficients, folder),
  shown on the DESeq2 page + index, printed to console; returned as `manifest` from
  `run_diffexp_pipeline()`.
- **Custom contrasts** beyond the automatic grid via `comparisons =` (spec table);
  `run_diffexp_coef()` is now exported. Single-factor `var` accepts labels like
  `"d0+d1_vs_d7+d14"` directly.
- **Fix:** `run_diffexp_pipeline()` can now actually take a contrast-spec table
  (`validate_comparisons()` used to reject it).
- **Grid no longer capped at 3 factors** (removed `.expand_contrast_specs_two/_three`).
- `stratified = TRUE` reachable from `run_deseq2_wald()`.
- **Fix:** optional modules (CSPA, etc.) can no longer destroy a finished analysis.
- **Shared orchestration** — new `R/pipeline_driver.R`; the three entry points no
  longer duplicate it.
- **Settable analysis settings** — new `R/options.R`: `degreports_opts()`.
- Unsafe factor levels (`_vs_`, `+`, path separators) fail at grid build.
- Docs regenerated (`man/`, `NAMESPACE`), README + vignette updated.

### Open / to check next session

- Fix the 7 pre-existing test failures (comparisons, diffexp-plots arrows,
  square-panels / wgcna ORA dotplot).
- Re-run both GSE189410 examples after the removeBatchEffect switch and rebuild
  `gh-pages` (the hosted report still shows ComBat-seq).
- Enable GitHub Pages (Settings → Pages → `gh-pages` / root) if not done; check
  the README button renders.
- `DESCRIPTION` `Description:` still names the removed `run_deseq2_wald()` /
  `run_deseq2_lrt()` / `run_deseq2_deg()`; update to `run_deseq2()` / `run_limma()`
  / `run_deg()`.
- Test spill is untracked in `inst/rmarkdown/` and `tests/testthat/` (PDF/CSV and
  `group*_vs_*` dirs): delete and/or gitignore; fix the tests that write there.
- Confirm with the Crick whether the institute should be the copyright holder.
- `CODE_WALKTHROUGH.md` is untracked: decide whether to commit it.

- [ ] 11 pre-existing test failures (present before the cleanup):
  - `test-comparisons.R` ×2 — stale: expect `validate_comparisons()` to return a
    character vector; since 0.4.0 it returns a `contrast_specs` table. Update tests.
  - (fixed 2026-11-04) msigdbr ×4 — test helper used the msigdbr ≥ 10 API.
  - plot assertions ×5 — stale tests (see 2026-11-04 entry), update them.
- [ ] Decide whether to drop unused exported functions (breaking API change).
- [ ] No new tests yet for `contrast_spec()`, `contrast_manifest()`,
      `degreports_opts()`, `pipeline_driver.R`, >3-factor grids.
- [ ] Run `devtools::check()` before tagging a 0.4.0 release.
