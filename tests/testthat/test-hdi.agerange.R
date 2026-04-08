test_that("hdi.agerange() throws an error", {
  expect_error(hdi.agerange(10))
})

test_that("hdi.agerange() throws an error when x not of class
          diagnostic_summary", {
            bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
            bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
            class(bay.ta_compare_diag) <- class(bay.ta_compare_diag)[-1]
            expect_error(hdi.agerange(bay.ta_compare_diag))
          })

test_that("hdi.agerange() produces correct output", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
  expect_snapshot(hdi.agerange(bay.ta_compare_diag))
})
