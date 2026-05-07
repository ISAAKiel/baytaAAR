test_that("bay.ta() throws an error", {
  expect_error(bay.ta(framework = "Test", method = as.matrix(sorsum_as[,2])))
})

test_that("bay.ta() with NIMBLE produced sensible output", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  skip_if_testcoverage()
  skip_if_on_ci_and_mac()
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  bay.ta_output <- bay.ta(framework = "NIMBLE",
                          method = sorsum_as[,2],
                          numSavedSteps = 1000,
                          seed = 1234)
  expect_equal(bay.ta_compare, bay.ta_output)
})

test_that("bay.ta() with JAGS produced sensible output", {
  skip_on_cran()
  skip_if_no_jags()
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_jags.Rds"))
  bay.ta_output <- bay.ta(framework = "JAGS",
                          method = sorsum_as[,2],
                          numSavedSteps = 1000,
                          seed = 1234)
  expect_equal(bay.ta_compare, bay.ta_output)
})

test_that("bay.ta() with NIMBLE and mnorm produced sensible output", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  skip_if_testcoverage()
  skip_if_on_ci_and_mac()
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  bay.ta_output <- bay.ta(framework = "NIMBLE",
                          algorithm = "mnorm",
                          method = spitalfields[,c(2:6)],
                          parameters = c("b", "a", "beta0", "beta",
                                         "thresh", "age.s", "Ustar"),
                          numSavedSteps = 1000,
                          seed = 1234)
  expect_equal(bay.ta_compare, bay.ta_output)
})
