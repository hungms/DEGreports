# Shared geometry helpers used across the package's plotting functions to
# keep every figure's data panel (the plotting area itself -- excluding
# axis tick text, axis titles, and legends) an exact square, with zero
# extra padding on continuous axes.

#' Zero data-padding for a continuous ggplot2 scale
#'
#' Returns a `ggplot2::expansion()` with both the multiplicative and
#' additive terms at exactly `0`, so the plotted data range extends
#' edge-to-edge in the rendered panel (no blank margin between the
#' outermost point/bar/line and the panel border). Used as
#' `scale_x_continuous(expand = zero_expansion())` /
#' `scale_y_continuous(expand = zero_expansion())` on every continuous
#' axis across the package's static ggplot figures.
#'
#' Discrete/categorical axes are deliberately left alone -- ggplot2's
#' small built-in discrete expansion (`add = 0.6` on each side) stays in
#' place wherever a plot has a discrete x/y aesthetic, since without it
#' the outermost bar/box/tile would be sliced exactly in half by the
#' panel border.
#'
#' @return A `ggplot2::expansion()` object usable as the `expand =`
#'   argument of `scale_*_continuous()`.
#' @noRd
zero_expansion <- function() {
    ggplot2::expansion(mult = 0, add = 0)
}

#' `coord_cartesian()` with clipping disabled, for use alongside
#' [zero_expansion()]
#'
#' Zero data-padding means the outermost point/line/box sits with its
#' *data coordinate* exactly on the panel border -- but most geoms have
#' rendered extent beyond that coordinate (a point marker's radius, a
#' line's stroke width, a boxplot whisker/outlier). `ggplot2`'s default
#' `coord_cartesian(clip = "on")` hard-clips any drawing at the panel
#' border, so with zero expansion that extent gets sliced off (e.g. a
#' boundary point rendered as a half-circle) -- this is a *distinct*
#' defect from panel squareness/shape and shows up as data seemingly
#' "cut off" at the plot edges. Passing `clip = "off"` (what this helper
#' sets) turns off that hard clip so the full marker/line/box renders,
#' still anchored at the same zero-padded data coordinate.
#'
#' Use this in place of a bare `ggplot2::coord_cartesian()` call on every
#' static ggplot in the package that (a) uses [zero_expansion()] on a
#' continuous axis and (b) draws a geom with rendered extent beyond its
#' data coordinate (points, lines, boxplots, or an outlined/stroked
#' `geom_col()`). Extra arguments (e.g. `ylim =`) are passed straight
#' through to `coord_cartesian()`, so an existing
#' `coord_cartesian(ylim = ...)` call can be swapped in place for
#' `unclipped_coord(ylim = ...)`.
#'
#' For a `coord_flip()`-based plot, use [unclipped_coord_flip()] instead
#' (a plot has exactly one coordinate system, so the two are mutually
#' exclusive).
#'
#' @param ... Passed to `ggplot2::coord_cartesian()` (e.g. `xlim =`,
#'   `ylim =`).
#' @return A `ggplot2::coord_cartesian()` object with `clip = "off"`.
#' @noRd
unclipped_coord <- function(...) {
    ggplot2::coord_cartesian(..., clip = "off")
}

#' `coord_flip()` with clipping disabled -- the flipped-axes counterpart
#' of [unclipped_coord()]
#'
#' @param ... Passed to `ggplot2::coord_flip()`.
#' @return A `ggplot2::coord_flip()` object with `clip = "off"`.
#' @noRd
unclipped_coord_flip <- function(...) {
    ggplot2::coord_flip(..., clip = "off")
}

#' Force a ggplot2 panel to a square physical aspect ratio
#'
#' `ggplot2::theme(aspect.ratio = 1)` fixes the *drawn panel rectangle*
#' (the data plotting area, excluding axis text/titles/legend) to be
#' square regardless of the data's own x/y range -- the mechanism
#' already used by [plot_pca()], [plot_ma()], and [plot_volcano()],
#' extended here to every other static ggplot in the package. It works
#' identically for faceted plots (each facet panel becomes square) and
#' after `coord_flip()`.
#'
#' @return A `ggplot2::theme()` object usable as `+ square_panel_theme()`.
#' @noRd
square_panel_theme <- function() {
    ggplot2::theme(aspect.ratio = 1)
}

