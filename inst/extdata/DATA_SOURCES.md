# Packaged reference data (`inst/extdata/`)

These files ship with the GitHub package and are the **default** inputs for
CSPA and decoupleR CollecTRI when you install via:

```r
devtools::install_github("hungms/DEGreports")
```

If a file is missing locally (e.g. a partial install), helpers download it
from this repository’s `main` branch into the user cache
(`tools::R_user_dir("DEGreports", "data")`).

| File | Retrieved | Source |
|---|---|---|
| `CSPA_pone.0121314.s003.xlsx` | paper supplement (2015) | Cell Surface Protein Atlas S3 ([Bausch-Fluck et al., *PLOS ONE* 2015](https://doi.org/10.1371/journal.pone.0121314)) |
| `collectri_human_2026-08-06.csv` | 2026-08-06 | [CollecTRI](https://doi.org/10.1038/s41467-023-39123-y) via OmniPath (`resources=CollecTRI`, taxid 9606), collapsed to `source` / `target` / `mor` |
| `collectri_mouse_2026-08-06.csv` | 2026-08-06 | Same, taxid 10090 |

Override at call time with `cspa_file=` or `collectri_file=`, or point
`DEGREPORTS_GITHUB_RAW` at another raw GitHub base URL.
