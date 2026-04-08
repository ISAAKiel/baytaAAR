test_that("diagnostics.max.min() throws an error", {
  expect_error(diagnostics.max.min(10))
})

test_that("diagnostics.max.min() throws an error when x not of class
          diagnostic_summary", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
  class(bay.ta_compare_diag) <- class(bay.ta_compare_diag)[-1]
  expect_error(diagnostics.max.min(bay.ta_compare_diag))
})

test_that("diagnostics.max.min() produces correct output", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
  expect_snapshot(diagnostics.max.min(bay.ta_compare_diag))
})
