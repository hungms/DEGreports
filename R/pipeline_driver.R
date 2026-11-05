# One driver for the three one-command entry points.
#
# `run_deseq2()`, `run_deseq2(mode = "lrt")` and `run_deg()` differ in how
# they test for differential expression and in a handful of mode-specific
# stages; everything else -- the order of the modules, the `[3/11] ...`
# progress counter, what happens when a module fails -- was written out three
# times. Describing a run as a list of stages puts that shared behaviour in
# one place, so adding or reordering a module no longer means editing a step
# counter in three files, and no module can quietly end up with a different
# failure policy from its neighbours.

#' Describe one pipeline stage
#'
#' @param name Element of the run state this stage's value is stored under
#'   (and therefore the name later stages read it by). `NULL` for a stage run
#'   purely for its side effects.
#' @param label Human-readable stage name, used for the progress message and
#'   for any failure warning.
#' @param run Function of the run state (a named list of earlier stage
#'   values), returning this stage's value.
#' @param enabled Logical; `FALSE` skips the stage entirely, and skipped
#'   stages do not count towards the progress total.
#' @param optional Logical. `TRUE` (the default for every analysis module)
#'   catches an error, warns, records `NULL`, and carries on -- a failing
#'   co-expression network should not destroy the differential expression that
#'   already ran. `FALSE` is for the stages nothing downstream can do without
#'   (QC, the DE test itself), where failing fast is the honest outcome.
#' @return A `pipeline_stage` list.
#' @noRd
#' @param silent Logical; a preparatory step that derives an input for later
#'   stages rather than producing an output of its own. It runs in order and
#'   obeys the same failure policy, but is neither announced nor counted, so
#'   the progress total still reflects the analyses the user asked for.
#' @noRd
.stage <- function(name, label, run, enabled = TRUE, optional = TRUE,
                   silent = FALSE) {
    stopifnot(is.function(run))
    structure(
        list(name = name, label = label, run = run,
             enabled = isTRUE(enabled), optional = isTRUE(optional),
             silent = isTRUE(silent)),
        class = "pipeline_stage"
    )
}

#' Run a list of stages, threading the run state through them
#'
#' @param stages List of [.stage()] descriptors, in execution order.
#' @param state Named list of starting values (typically `list(dds = dds)`).
#' @return The run state after every enabled stage has run. A stage that was
#'   skipped or that failed leaves `NULL` under its name, so callers can treat
#'   "not requested" and "did not work" the same way.
#' @noRd
.run_stages <- function(stages, state = list()) {
    active <- Filter(function(s) isTRUE(s$enabled), stages)
    n <- sum(!vapply(active, function(s) isTRUE(s$silent), logical(1)))
    shown <- 0L
    # Names of skipped stages still exist in the state, as NULL, so the report
    # and the return value do not have to test for missing elements.
    for (s in stages) {
        if (!is.null(s$name) && !isTRUE(s$enabled) && is.null(state[[s$name]])) {
            state[s$name] <- list(NULL)
        }
    }
    for (i in seq_along(active)) {
        s <- active[[i]]
        if (!isTRUE(s$silent)) {
            shown <- shown + 1L
            message("[", shown, "/", n, "] ", s$label, "...")
        }
        value <- if (isTRUE(s$optional)) {
            tryCatch(s$run(state), error = function(e) {
                warning(s$label, " failed (", conditionMessage(e),
                        "); continuing without it.", call. = FALSE)
                NULL
            })
        } else {
            s$run(state)
        }
        if (!is.null(s$name)) state[s$name] <- list(value)
    }
    state
}
