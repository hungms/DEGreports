# A runnable record of each pipeline run.
#
# Every entry point writes `reproduce/reproduce_run.R` under `save_dir`: the
# same entry-point call with only the arguments the user set, plus the input
# object (and any argument too big to read inline) in `reproduce/inputs.rds`,
# so the script runs as is. The index page links and shows it.

#' Folder, script and data names for the reproducible script (relative to
#' `save_dir`)
#' @noRd
.repro_dir <- "reproduce"
.repro_script <- "reproduce_run.R"
.repro_inputs <- "inputs.rds"

#' Deparse a value as R code, or `NULL` if it is too big to read inline
#' @noRd
.repro_deparse <- function(x, max_lines = 8L) {
    if (methods::is(x, "S4") || is.data.frame(x) || is.environment(x) ||
            is.function(x)) {
        return(NULL)
    }
    if (inherits(x, "formula")) {
        return(paste(deparse(x, width.cutoff = 500L), collapse = " "))
    }
    code <- deparse(x, width.cutoff = 70L)
    if (length(code) > max_lines) return(NULL)
    paste(code, collapse = "\n")
}

#' Write the reproducible script for a pipeline run
#'
#' @param fn Entry-point name, e.g. `"run_deseq2"`.
#' @param data_arg Name of the entry point's data argument (`"dds"`,
#'   `"object"`, `"deg_table"`).
#' @param data_var Variable name the script gives the input.
#' @param data_comment Comment line describing the input.
#' @param args Named list of every formal's value at call time.
#' @param call The matched call, used to tell which arguments the user set.
#' @param save_dir Run output folder.
#' @param run_comment One-line comment printed above the entry-point call.
#' @return Invisibly, the script path; `NULL` (with a warning) if writing
#'   failed -- a missing script never stops the analysis.
#' @noRd
.write_repro_script <- function(fn, data_arg, data_var, data_comment, args, call,
                                save_dir, run_comment) {
    tryCatch(.write_repro_script_impl(fn, data_arg, data_var, data_comment, args,
                                      call, save_dir, run_comment),
             error = function(e) {
                 warning("Could not write the reproducible R script: ",
                         conditionMessage(e), call. = FALSE)
                 NULL
             })
}

#' @noRd
.write_repro_script_impl <- function(fn, data_arg, data_var, data_comment, args,
                                     call, save_dir, run_comment) {
    out_dir <- file.path(validate_save_dir(save_dir), .repro_dir)
    dir.create(out_dir, showWarnings = FALSE)
    unlink(file.path(out_dir, .repro_inputs))

    # Only what the user set: defaults are left to the function.
    given <- setdiff(intersect(names(as.list(call)), names(args)),
                     c(data_arg, "save_dir"))

    inputs <- list()
    data_value <- args[[data_arg]]
    data_line <- if (is.character(data_value) && length(data_value) == 1L &&
                         file.exists(data_value)) {
        # A file path: point at the file rather than copying it.
        sprintf("%s <- %s", data_var,
                deparse(normalizePath(data_value, mustWork = TRUE)))
    } else {
        inputs[[data_var]] <- data_value
        sprintf("%s <- inputs$%s", data_var, data_var)
    }

    arg_lines <- character(0)
    opts_lines <- character(0)
    for (nm in given) {
        val <- args[[nm]]
        if (identical(nm, "opts")) {
            diff <- .opts_diff(.as_opts(val))
            if (!length(diff)) next
            names(diff) <- sub("^opts\\$", "", names(diff))
            inner <- vapply(names(diff), function(k) {
                paste0(k, " = ", .repro_deparse(diff[[k]]))
            }, character(1))
            opts_lines <- c(
                "# Analysis settings changed from the defaults",
                paste0("opts <- degreports_opts(",
                       paste(inner, collapse = ", "), ")"),
                ""
            )
            arg_lines <- c(arg_lines, "opts = opts")
            next
        }
        code <- .repro_deparse(val)
        if (is.null(code)) {
            inputs[[nm]] <- val
            code <- paste0("inputs$", nm)
        }
        arg_lines <- c(arg_lines, paste0(nm, " = ", code))
    }
    arg_lines <- c(arg_lines, "save_dir = \"reproduced_run\"")
    arg_lines <- gsub("\n", "\n    ", arg_lines, fixed = TRUE)

    if (length(inputs)) saveRDS(inputs, file.path(out_dir, .repro_inputs))

    lines <- c(
        paste0("# Reproduce this DEGreports run (", fn, ")"),
        paste0("# Written ", format(Sys.time(), "%Y-%m-%d %H:%M"),
               " by DEGreports ", utils::packageVersion("DEGreports"), "."),
        "# Run from the folder that holds index.html; results go to reproduced_run/.",
        "",
        "# Load the package",
        "library(DEGreports)",
        "",
        if (length(inputs)) c(
            "# Load the inputs saved with this run",
            sprintf("inputs <- readRDS(\"%s/%s\")", .repro_dir, .repro_inputs),
            ""
        ),
        data_comment,
        data_line,
        "",
        opts_lines,
        paste0("# ", run_comment),
        paste0("out <- ", fn, "("),
        paste0("    ", c(data_var, arg_lines),
               c(rep(",", length(arg_lines)), "")),
        ")",
        "",
        "# Open the report",
        "# browseURL(\"reproduced_run/index.html\")"
    )
    path <- file.path(out_dir, .repro_script)
    writeLines(lines, path)
    message("Reproducible R script: ", path)
    invisible(path)
}

#' Script text and link for the index page, or `NULL` when there is none
#' @noRd
.repro_for_report <- function(save_dir) {
    path <- file.path(save_dir, .repro_dir, .repro_script)
    if (!file.exists(path)) return(NULL)
    list(href = paste0(.repro_dir, "/", .repro_script),
         code = readLines(path, warn = FALSE))
}
