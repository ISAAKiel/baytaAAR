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
  test_matrix <- matrix(c(1.000000000,  0.498422732,  0.144545407,  0.424000111,
                          -0.051273130, 0.498422732,  1.000000000,  0.360682651,
                          0.608044545,  0.197598667,  0.144545407,  0.360682651,
                          1.000000000,  0.446595062, -0.008130608,  0.424000111,
                          0.608044545,  0.446595062,  1.000000000,  0.202390885,
                          -0.051273130,  0.197598667, -0.008130608, 0.202390885,
                          1.000000000), ncol = 5)
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  bay.ta_corr.mat.mean <- corr.mat.mean(bay.ta_compare)
  expect_equal(corr.mat.mean(bay.ta_compare), test_matrix)
})
