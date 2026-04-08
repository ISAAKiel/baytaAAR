test_that("prob.cat() throws an error", {
  expect_error(prob.cat(10))
})

test_that("prob.cat() throws an error when x not of class mcmc.list", {
            bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
            class(bay.ta_compare) <- class(bay.ta_compare)[-1]
            expect_error(threshold.chains(bay.ta_compare))
          })

test_that("prob.cat() throws an error when no original data is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare))
})

test_that("prob.cat() throws an error when no data.frame is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare, df_orig = as.matrix(spitalfields)))
})

test_that("prob.cat() throws an error when no grouping name is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare, df_orig = spitalfields))
})

test_that("prob.cat() throws an error when no mode is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare, df_orig = spitalfields,
                        group_col = "Sex"))
})

test_that("prob.cat() throws an error when wrong mode is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare, df_orig = spitalfields,
                        group_col = "Sex", mode = "median"))
})

test_that("prob.cat() throws an error when wrong age-modifier is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare, df_orig = spitalfields,
                        group_col = "Sex", mode = "mean", age_identifier = "age"))
})

test_that("prob.cat() produces correct output for mean", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_snapshot(prob.cat(bay.ta_compare, df_orig = spitalfields,
                           group_col = 7, mode = "mean"))
})

test_that("prob.cat() produces correct output for mode = summed", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_snapshot(prob.cat(bay.ta_compare, df_orig = spitalfields,
                           group_col = 7, mode = "summed"))
})
