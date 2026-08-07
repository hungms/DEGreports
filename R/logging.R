#' Run an expression while logging its call, console output, and session info
#'
#' Captures the calling expression (`call`), all messages/warnings printed
#' during evaluation, and a `sessionInfo()` snapshot, and writes them as a
#' timestamped plain-text log plus a machine-readable JSON manifest to
#' `save_dir/logs/`. This gives a per-run audit trail without depending on
#' `renv` (so it works even for users who do not use `renv` in their
#' project).
#'
#' @param call The call to record, typically `match.call()` from the caller.
#' @param save_dir Directory under which a `logs/` subdirectory will be created.
#' @param expr An expression (wrapped in `quote()`) to evaluate in the
#'   calling function's environment (`parent.frame()`).
#' @param params Optional named list of resolved parameter values to record
#'   in the JSON manifest (e.g. `var`, `comparisons`, `min_count`).
#' @param envir Environment in which to evaluate `expr`. Default
#'   `parent.frame()` -- i.e. the environment of the function that called
#'   `log_output()`, so `expr` can see that function's local variables and
#'   arguments.
#' @return The value of evaluating `expr`.
#' @export
log_output <- function(call, save_dir, expr, params = list(), envir = parent.frame()) {
    save_dir <- validate_save_dir(save_dir)
    log_dir <- file.path(save_dir, "logs")
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_path <- file.path(log_dir, paste0("run_", timestamp, ".log"))
    manifest_path <- file.path(log_dir, paste0("run_", timestamp, ".json"))

    con <- file(log_path, open = "wt")
    on.exit(close(con), add = TRUE)

    writeLines(paste0("deseq2reports run log -- ", Sys.time()), con)
    writeLines(paste0("Call: ", paste(deparse(call), collapse = " ")), con)
    writeLines("", con)

    # Stream messages/warnings into the log as they arrive (so a crash still
    # leaves a useful trail), and also echo them to the console.
    write_log_line <- function(line) {
        writeLines(line, con)
        flush(con)
    }
    result <- withCallingHandlers(
        eval(expr, envir = envir),
        message = function(m) {
            write_log_line(paste0("[message] ", conditionMessage(m)))
            # Echo to stderr without re-entering this handler via message().
            cat(conditionMessage(m), "\n", file = stderr(), sep = "")
            invokeRestart("muffleMessage")
        },
        warning = function(w) {
            write_log_line(paste0("[warning] ", conditionMessage(w)))
            cat("[warning] ", conditionMessage(w), "\n", file = stderr(), sep = "")
            invokeRestart("muffleWarning")
        }
    )

    writeLines("", con)
    writeLines("--- sessionInfo() ---", con)
    writeLines(utils::capture.output(utils::sessionInfo()), con)

    manifest <- list(
        timestamp = timestamp,
        call = paste(deparse(call), collapse = " "),
        package_version = as.character(utils::packageVersion("deseq2reports")),
        r_version = as.character(getRversion()),
        params = params
    )
    writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"), manifest_path)

    result
}
