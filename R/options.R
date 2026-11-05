# The analysis settings the one-command entry points used to hard-code.
#
# They were fixed (in `.pipeline_const`) so that every run is comparable and
# the report's "Pipeline arguments" table stays about study design rather than
# tuning knobs. That reasoning still holds -- which is why the defaults here
# are exactly the old constants, and why the report lists only the ones a run
# actually changed -- but fixing them meant an analyst who needed a different
# `alpha` had to abandon the entry point and wire ten module calls by hand.

#' Analysis settings for the pipeline entry points
#'
#' Collects the statistical and per-module settings that
#' [run_deseq2()], [`run_deseq2(mode = "lrt")`][run_deseq2] and [run_deg()] apply. The
#' defaults are the values those functions used to hard-code, so
#' `opts = degreports_opts()` reproduces previous behaviour exactly; pass a
#' modified object to deviate deliberately.
#'
#' Only the settings that differ from the defaults are listed in the report's
#' "Pipeline arguments" table, so a standard run still reads as a statement of
#' study design.
#'
#' @param shrinkage LFC shrinkage estimator: `"ashr"` (default), `"apeglm"`,
#'   `"normal"` or `"none"`.
#' @param alpha Adjusted-p significance threshold. Default `0.05`.
#' @param lfc_threshold Absolute log2 fold-change a gene must exceed to count
#'   as significant. Default `0.5`.
#' @param gsva_n_heatmap Pathways shown per GSVA heatmap. Default `50`.
#' @param wgcna_min_module_size Minimum WGCNA module size. Default `100`.
#' @param decoupler_n_tfs Transcription factors shown per decoupleR plot.
#'   Default `25`.
#' @param decoupler_n_lr Ligand-receptor pairs shown per plot. Default `10`.
#' @param viper_n_tfs Regulators shown per VIPER plot. Default `10`.
#' @param ppi_top_n_hubs Hub genes labelled per PPI network. Default `10`.
#' @param cspa_n_markers Surface markers per contrast (Wald / DEG-table
#'   modes). Default `25`.
#' @param cspa_n_markers_lrt Surface markers per degPatterns cluster (LRT
#'   mode, where there are many more groups). Default `10`.
#' @return A validated list of class `degreports_opts`.
#' @examples
#' opts <- degreports_opts(alpha = 0.01, shrinkage = "apeglm")
#' opts$alpha
#' @seealso [run_deseq2()]
#' @export
degreports_opts <- function(shrinkage = "ashr",
                               alpha = 0.05,
                               lfc_threshold = 0.5,
                               gsva_n_heatmap = 50L,
                               wgcna_min_module_size = 100L,
                               decoupler_n_tfs = 25L,
                               decoupler_n_lr = 10L,
                               viper_n_tfs = 10L,
                               ppi_top_n_hubs = 10L,
                               cspa_n_markers = 25L,
                               cspa_n_markers_lrt = 10L) {

    count <- function(x, nm) {
        if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1) {
            stop("`", nm, "` must be a single positive number.", call. = FALSE)
        }
        as.integer(x)
    }
    out <- list(
        shrinkage = validate_shrinkage(shrinkage),
        alpha = validate_alpha(alpha),
        lfc_threshold = validate_lfc_threshold(lfc_threshold),
        gsva_n_heatmap = count(gsva_n_heatmap, "gsva_n_heatmap"),
        wgcna_min_module_size = count(wgcna_min_module_size, "wgcna_min_module_size"),
        decoupler_n_tfs = count(decoupler_n_tfs, "decoupler_n_tfs"),
        decoupler_n_lr = count(decoupler_n_lr, "decoupler_n_lr"),
        viper_n_tfs = count(viper_n_tfs, "viper_n_tfs"),
        ppi_top_n_hubs = count(ppi_top_n_hubs, "ppi_top_n_hubs"),
        cspa_n_markers = count(cspa_n_markers, "cspa_n_markers"),
        cspa_n_markers_lrt = count(cspa_n_markers_lrt, "cspa_n_markers_lrt")
    )
    class(out) <- c("degreports_opts", "list")
    out
}

#' Coerce / validate an `opts` argument
#'
#' Accepts the object, a plain list of overrides, or `NULL`.
#' @noRd
.as_opts <- function(opts) {
    if (is.null(opts)) return(degreports_opts())
    if (inherits(opts, "degreports_opts")) return(opts)
    if (is.list(opts)) {
        unknown <- setdiff(names(opts), names(degreports_opts()))
        if (length(unknown)) {
            stop("Unknown setting(s) in `opts`: ", paste(unknown, collapse = ", "),
                 ". See degreports_opts().", call. = FALSE)
        }
        return(do.call(degreports_opts, opts))
    }
    stop("`opts` must be a degreports_opts() object, a named list, or `NULL`.",
         call. = FALSE)
}

#' Settings a run changed from the defaults
#'
#' Keeps the report's argument table about study design: a run that touched
#' nothing contributes nothing.
#' @return Named list of the differing settings (empty when all are default).
#' @noRd
.opts_diff <- function(opts) {
    def <- degreports_opts()
    changed <- vapply(names(def), function(nm) !identical(def[[nm]], opts[[nm]]),
                      logical(1))
    if (!any(changed)) return(list())
    out <- opts[names(def)[changed]]
    stats::setNames(out, paste0("opts$", names(out)))
}

#' @export
print.degreports_opts <- function(x, ...) {
    def <- degreports_opts()
    cat("<degreports_opts>\n")
    for (nm in names(def)) {
        marker <- if (identical(def[[nm]], x[[nm]])) "" else "   (changed)"
        cat(sprintf("  %-22s %s%s\n", nm, format(x[[nm]]), marker))
    }
    invisible(x)
}
