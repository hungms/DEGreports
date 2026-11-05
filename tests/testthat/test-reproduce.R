# Reproducible R script written with every run.

test_that(".write_repro_script keeps only user-set arguments and saves inputs", {
    save_dir <- withr::local_tempdir()
    f <- function(dds, var, alpha_dummy = 1, run_gsva = FALSE, gene_chr = NULL,
                  opts = degreports_opts(), save_dir = getwd()) {
        call <- match.call()
        .write_repro_script(
            "run_deseq2", data_arg = "dds", data_var = "dds",
            data_comment = "# input", args = mget(names(formals()), environment()),
            call = call, save_dir = save_dir, run_comment = "Run it")
    }
    gc <- data.frame(gene = paste0("g", 1:50), chr = "1")
    path <- suppressMessages(f(
        dds = list(a = 1), var = c("A", "B"), run_gsva = TRUE, gene_chr = gc,
        opts = degreports_opts(alpha = 0.1), save_dir = save_dir))
    expect_true(file.exists(path))
    code <- readLines(path)
    expect_silent(parse(text = code))
    expect_true(any(grepl('var = c("A", "B"),', code, fixed = TRUE)))
    expect_true(any(grepl("run_gsva = TRUE,", code, fixed = TRUE)))
    expect_true(any(grepl("gene_chr = inputs$gene_chr,", code, fixed = TRUE)))
    expect_true(any(grepl("opts <- degreports_opts(alpha = 0.1)", code, fixed = TRUE)))
    expect_false(any(grepl("alpha_dummy", code, fixed = TRUE)))
    expect_true(any(grepl('save_dir = "reproduced_run"', code, fixed = TRUE)))
    # Every code block has a comment above it.
    expect_true(all(c("# Load the package", "# Run it") %in% code))

    inputs <- readRDS(file.path(save_dir, "reproduce", "inputs.rds"))
    expect_equal(inputs$dds, list(a = 1))
    expect_equal(inputs$gene_chr, gc)

    rep <- .repro_for_report(save_dir)
    expect_identical(rep$href, "reproduce/reproduce_run.R")
    expect_identical(rep$code, code)
})

test_that("a file-path input is referenced, not copied", {
    save_dir <- withr::local_tempdir()
    csv <- file.path(save_dir, "deg.csv")
    writeLines("gene,log2FoldChange", csv)
    f <- function(deg_table, save_dir = getwd()) {
        .write_repro_script(
            "run_deg", data_arg = "deg_table", data_var = "deg_table",
            data_comment = "# input", args = mget(names(formals()), environment()),
            call = match.call(), save_dir = save_dir, run_comment = "Run it")
    }
    code <- readLines(suppressMessages(f(csv, save_dir = save_dir)))
    expect_true(any(grepl(normalizePath(csv), code, fixed = TRUE)))
    expect_false(file.exists(file.path(save_dir, "reproduce", "inputs.rds")))
    expect_false(any(grepl("readRDS", code, fixed = TRUE)))
})

test_that("the script written by run_limma reproduces the run", {
    dds <- simulate_dds(n_per_group = 3)
    x <- DESeq2::counts(dds) * stats::runif(nrow(dds), 0.5, 1.5)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(tpm = sweep(x, 2, colSums(x), "/") * 1e6),
        colData = SummarizedExperiment::colData(dds))
    save_dir <- withr::local_tempdir()
    suppressWarnings(suppressMessages(run_limma(
        se, var = "condition", render_html = FALSE, run_cspa = FALSE,
        collections = "HALLMARK", save_dir = save_dir)))
    code <- readLines(file.path(save_dir, "reproduce", "reproduce_run.R"))
    expect_true(any(grepl("run_limma(", code, fixed = TRUE)))
    expect_true(any(grepl("render_html = FALSE,", code, fixed = TRUE)))

    # Rerun with the package already loaded (library() is a no-op here).
    withr::with_dir(save_dir, suppressWarnings(suppressMessages(
        eval(parse(text = code[!grepl("^library\\(", code)]), envir = new.env())
    )))
    rel <- file.path("deseq2", "Treated_vs_Control", "Treated_vs_Control_diffexp.csv")
    expect_equal(utils::read.csv(file.path(save_dir, "reproduced_run", rel)),
                 utils::read.csv(file.path(save_dir, rel)))
})
