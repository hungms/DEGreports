test_that("extdata catalog resolves dated CollecTRI and CSPA from inst/extdata", {
    human <- .collectri_default_path("human", download = FALSE)
    mouse <- .collectri_default_path("mouse", download = FALSE)
    cspa <- .cspa_default_path(download = FALSE)
    expect_true(file.exists(human))
    expect_true(file.exists(mouse))
    expect_true(file.exists(cspa))
    expect_match(basename(human), "^collectri_human_20[0-9]{2}-[0-9]{2}-[0-9]{2}\\.csv$")
    expect_match(basename(mouse), "^collectri_mouse_20[0-9]{2}-[0-9]{2}-[0-9]{2}\\.csv$")
    expect_equal(basename(cspa), "CSPA_pone.0121314.s003.xlsx")
})

test_that("summarize_genes still works after extdata helpers load", {
    # Smoke: helpers did not break package load / other exports
    expect_true(is.function(summarize_genes))
    expect_true(is.function(load_cspa_genes))
})
