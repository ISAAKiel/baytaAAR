test_that("diagnostic.summary() throws an error", {
  expect_error(diagnostic.summary(10))
})

test_that("diagnostic.summary() throws an error when HDI > 1", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  expect_error(diagnostic.summary(bay.ta_compare, HDImass = 1.1))
})

test_that("diagnostic.summary() throws an error when gelman_diag not T/F", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  expect_error(diagnostic.summary(bay.ta_compare, gelman_diag = 1))
})

test_that("diagnostic.summary() produces correct output", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  expect_snapshot(diagnostic.summary(bay.ta_compare))
})
