test_that("prob.cat() throws an error", {
  expect_error(prob.cat(10))
})

test_that("prob.cat() throws an error when x not of class mcmc.list", {
            bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
            class(bay.ta_compare) <- class(bay.ta_compare)[-1]
            expect_error(prob.cat(bay.ta_compare))
          })

test_that("prob.cat() throws an error when no grouping vector is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare))
})

test_that("prob.cat() throws an error when no mode is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare, group_vec = spitalfields$Sex))
})

test_that("prob.cat() throws an error when wrong mode is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare, group_vec = spitalfields$Sex,
                        mode = "median"))
})

test_that("prob.cat() throws an error when wrong age-modifier is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(prob.cat(bay.ta_compare, group_vec = spitalfields$Sex,
                        mode = "mean", age_identifier = "age"))
})

test_that("prob.cat() produces correct output for mean", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_snapshot(prob.cat(bay.ta_compare, group_vec = spitalfields$Sex,
                           mode = "mean"))
})

test_that("prob.cat() produces correct output for mode = summed", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_snapshot(prob.cat(bay.ta_compare, group_vec = spitalfields$Sex,
                           mode = "summed"))
})
