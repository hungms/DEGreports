test_that("isolating a worker function does not use baseenv() or mutate the original", {
    skip_if_not_installed("future")
    skip_if_not_installed("future.apply")

    orig_env <- environment(deseq2reports:::.parallel_hub_worker)
    isolated <- deseq2reports:::.as_isolated_future_fun(deseq2reports:::.parallel_hub_worker)
    expect_false(identical(environment(isolated), baseenv()))
    expect_identical(parent.env(environment(isolated)), emptyenv())
    expect_identical(environment(deseq2reports:::.parallel_hub_worker), orig_env)
})

test_that("isolated hub worker can be sent to a future without parent-env cycles", {
    skip_if_not_installed("future")
    skip_if_not_installed("future.apply")
    skip_if_not_installed("pkgload")

    pkg_dir <- deseq2reports:::.deseq2reports_src_dir()
    skip_if(is.null(pkg_dir), "not a pkgload/source checkout")

    tmp <- withr::local_tempfile(fileext = ".rds")
    saveRDS(list(ok = TRUE), tmp)

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = 2)

    out <- future.apply::future_lapply(
        "probe",
        deseq2reports:::.as_isolated_future_fun(function(key, shared_rds, pkg_dir) {
            pkgload::load_all(pkg_dir, export_all = FALSE, quiet = TRUE)
            list(
                key = key,
                ok = isTRUE(base::readRDS(shared_rds)$ok),
                has_helper = base::exists(
                    ".render_one_hub_page",
                    envir = base::asNamespace("deseq2reports"),
                    inherits = FALSE
                )
            )
        }),
        shared_rds = tmp,
        pkg_dir = pkg_dir,
        future.packages = "pkgload",
        future.globals = FALSE,
        future.seed = TRUE
    )

    expect_equal(out[[1]]$key, "probe")
    expect_true(out[[1]]$ok)
    expect_true(out[[1]]$has_helper)
})
