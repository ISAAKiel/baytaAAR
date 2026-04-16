test_that("threshold.chains() throws an error", {
  expect_error(threshold.chains(10), "Assertion on 'mcmc_list' failed: Must inherit from class 'mcmc.list', but has class 'numeric'.")
})

test_that("threshold.chains() throws an error when x not of class mcmc.list", {
            bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
            class(bay.ta_compare) <- class(bay.ta_compare)[-1]
            expect_error(threshold.chains(bay.ta_compare), "Assertion on 'mcmc_list' failed: Must inherit from class 'mcmc.list', but has class 'list'.")
          })

test_that("threshold.chains() produces correct output", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  bay.ta_thresholds <- threshold.chains(bay.ta_compare)
  expect_snapshot(threshold.chains(bay.ta_compare))
})
