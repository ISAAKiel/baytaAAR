test_that("corr.mat.mean() throws an error", {
  expect_error(corr.mat.mean(10))
})

test_that("corr.mat.mean() throws an error when x not of class
          mcmc.list", {
            bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
            class(bay.ta_compare) <- class(bay.ta_compare)[-1]
            expect_error(corr.mat.mean(bay.ta_compare))
          })

test_that("corr.mat.mean() produces correct output", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_snapshot(corr.mat.mean(bay.ta_compare))
})
