# run_deseq2(): mode dispatch and mode-specific arguments.

test_that("run_deseq2 rejects arguments that belong to the other mode", {
    dds <- simulate_dds()
    dir <- withr::local_tempdir()
    expect_error(run_deseq2(dds, var = "condition", time = "day", save_dir = dir),
                 "`time` only applies to mode = \"lrt\"")
    expect_error(run_deseq2(dds, mode = "lrt", time = "condition",
                            comparisons = "Treated_vs_Control", stratified = TRUE,
                            save_dir = dir),
                 "`comparisons`, `stratified` only apply to mode = \"wald\"")
    expect_error(run_deseq2(dds, mode = "lrt", save_dir = dir),
                 "`time` is required")
    expect_error(run_deseq2(dds, var = "condition", mode = "bogus", save_dir = dir),
                 "mode")
})

test_that("the old entry-point names are gone", {
    ns <- asNamespace("DEGreports")
    for (old in c("run_deseq2_wald", "run_deseq2_lrt", "run_deseq2_pipeline",
                  "run_deseq2_deg", "run_deseq2_limma")) {
        expect_false(exists(old, envir = ns, inherits = FALSE), info = old)
    }
})
