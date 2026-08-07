#' Bundled extdata filenames (shipped under `inst/extdata/`)
#'
#' CollecTRI snapshots are date-stamped (`YYYY-MM-DD` = OmniPath download day).
#' Update these when replacing the CSVs and keep
#' [inst/extdata/DATA_SOURCES.md] in sync.
#' @noRd
.extdata_catalog <- function() {
    list(
        cspa = "CSPA_pone.0121314.s003.xlsx",
        collectri_human = "collectri_human_2026-08-06.csv",
        collectri_mouse = "collectri_mouse_2026-08-06.csv"
    )
}

#' GitHub raw base URL for packaged extdata
#' @noRd
.extdata_github_raw_base <- function() {
    Sys.getenv(
        "DESEQ2REPORTS_GITHUB_RAW",
        unset = "https://raw.githubusercontent.com/hungms/deseq2reports/main"
    )
}

#' User cache directory for downloaded extdata
#' @noRd
.extdata_cache_dir <- function() {
    path <- tools::R_user_dir("deseq2reports", which = "data")
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

#' Candidate directories that may contain `inst/extdata` files
#' @noRd
.extdata_search_dirs <- function() {
    pkg_path <- tryCatch(find.package("deseq2reports"), error = function(e) "")
    unique(c(
        system.file("extdata", package = "deseq2reports"),
        if (nzchar(pkg_path)) {
            c(
                file.path(pkg_path, "extdata"),
                file.path(pkg_path, "inst", "extdata"),
                file.path(pkg_path, "example_data")
            )
        },
        file.path("inst", "extdata"),
        file.path("example_data"),
        .extdata_cache_dir()
    ))
}

#' Locate an extdata file, optionally downloading from the GitHub repo
#'
#' Search order: installed package `extdata/` → source-tree
#' `inst/extdata/` → user cache → download from
#' `https://raw.githubusercontent.com/hungms/deseq2reports/main/inst/extdata/`
#' into the user cache.
#'
#' @param filename Basename under `inst/extdata/` (e.g. dated CollecTRI CSV).
#' @param download If `TRUE` (default), fetch from GitHub when missing locally.
#' @return Absolute path to the file.
#' @noRd
.resolve_extdata_file <- function(filename, download = TRUE) {
    filename <- as.character(filename)[1]
    if (is.na(filename) || !nzchar(filename)) {
        stop("`filename` must be a non-empty character string.", call. = FALSE)
    }

    for (dir in .extdata_search_dirs()) {
        if (!nzchar(dir) || !dir.exists(dir)) next
        cand <- file.path(dir, filename)
        if (file.exists(cand)) {
            return(normalizePath(cand, mustWork = TRUE))
        }
    }

    if (!isTRUE(download)) {
        stop("Bundled file \"", filename, "\" not found locally. ",
             "Install deseq2reports from GitHub or place the file under ",
             "inst/extdata/.", call. = FALSE)
    }

    url <- paste0(
        .extdata_github_raw_base(),
        "/inst/extdata/", filename
    )
    dest <- file.path(.extdata_cache_dir(), filename)
    message("Downloading ", filename, " from GitHub (",
            .extdata_github_raw_base(), ") ...")
    status <- tryCatch(
        utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE),
        error = function(e) e
    )
    if (inherits(status, "error") || !identical(status, 0L) || !file.exists(dest) ||
            file.info(dest)$size < 10) {
        if (file.exists(dest)) unlink(dest)
        stop("Failed to download \"", filename, "\" from ", url, ". ",
             "Install deseq2reports from GitHub ",
             "(devtools::install_github(\"hungms/deseq2reports\")) ",
             "or pass an explicit local file path.",
             if (inherits(status, "error")) {
                 paste0(" Error: ", conditionMessage(status))
             } else {
                 ""
             },
             call. = FALSE)
    }
    normalizePath(dest, mustWork = TRUE)
}

#' Resolve the current CollecTRI CSV for an organism
#'
#' Prefers the catalogued date-stamped file shipped in the GitHub package,
#' then any newer `collectri_<org>_YYYY-MM-DD.csv` found on disk, then
#' undated `collectri_<org>.csv`.
#' @noRd
.collectri_default_path <- function(organism = c("human", "mouse"),
                                    download = TRUE) {
    organism <- match.arg(organism)
    catalog <- .extdata_catalog()
    preferred <- catalog[[paste0("collectri_", organism)]]

    # Prefer an already-present file (catalog name or any dated / undated).
    found <- character(0)
    pat <- paste0("^collectri_", organism, "(_[0-9]{4}-[0-9]{2}-[0-9]{2})?\\.csv$")
    for (dir in .extdata_search_dirs()) {
        if (!nzchar(dir) || !dir.exists(dir)) next
        hits <- list.files(dir, pattern = pat, full.names = TRUE)
        found <- c(found, hits)
    }
    found <- unique(found[file.exists(found)])
    if (length(found)) {
        # Newest date stamp wins; undated sorts first and loses to dated.
        bn <- basename(found)
        dates <- sub(
            paste0("^collectri_", organism, "_([0-9]{4}-[0-9]{2}-[0-9]{2})\\.csv$"),
            "\\1", bn
        )
        dates[!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", dates)] <- "0000-00-00"
        found <- found[order(dates, decreasing = TRUE)]
        return(normalizePath(found[[1]], mustWork = TRUE))
    }

    # Nothing local: fetch the catalogued snapshot from GitHub.
    .resolve_extdata_file(preferred, download = download)
}

#' Resolve the bundled CSPA Excel workbook
#' @noRd
.cspa_default_path <- function(download = TRUE) {
    .resolve_extdata_file(.extdata_catalog()$cspa, download = download)
}
