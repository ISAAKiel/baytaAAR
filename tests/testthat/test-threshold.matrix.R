test_that("threshold.matrix() throws an error", {
  expect_error(threshold.matrix(10))
})

test_that("threshold.matrix() throws an error when x not of class
          diagnostic_summary", {
            bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
            bay.ta_compare_mcmc_list <- threshold.chains(bay.ta_compare)
            bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare_mcmc_list)
            class(bay.ta_compare_diag) <- class(bay.ta_compare_diag)[-1]
            expect_error(threshold.matrix(bay.ta_compare_diag))
          })

test_that("threshold.chains() produces correct output", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  bay.ta_compare_mcmc_list <- threshold.chains(bay.ta_compare)
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare_mcmc_list)
  expect_equal(threshold.matrix(bay.ta_compare_diag),
               matrix(c(14.80124860, 25.57174185, 33.37472069, 46.46155074,
                        63.97373410, 59.65893684, 89.83167260), nrow = 1))
})
