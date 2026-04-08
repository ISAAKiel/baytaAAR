test_that("sequential.binom.test() throws an error", {
  expect_error(sequential.binom.test(10))
})

test_that("sequential.binom.test() throws an error when x not of class mcmc.list", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  class(bay.ta_compare) <- class(bay.ta_compare)[-1]
  expect_error(sequential.binom.test(bay.ta_compare))
})

test_that("sequential.binom.test() throws an error when no known age data is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(sequential.binom.test(bay.ta_compare))
})

test_that("sequential.binom.test() throws an error when only NAs for known age data is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(sequential.binom.test(bay.ta_compare, known_age = c(NA, NA)))
})

test_that("sequential.binom.test() throws an error when only NA for HDImass is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(sequential.binom.test(bay.ta_compare, known_age = spitalfields$Age,
                                  HDImass = NA))
})

test_that("sequential.binom.test() throws an error when only non-unique values for HDImass are supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(sequential.binom.test(bay.ta_compare, known_age = spitalfields$Age,
                                     HDImass = c(0.5, 0.6, 0.5)))
})

test_that("sequential.binom.test() produces correct output with default HDImass value", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_snapshot(sequential.binom.test(bay.ta_compare,
                                        known_age = spitalfields$Age))
})

test_that("sequential.binom.test() produces correct output with several HDImass values", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_snapshot(sequential.binom.test(bay.ta_compare,
                                        known_age = spitalfields$Age,
                                        HDImass = c(0.75, 0.95)))
})
