test_that("age.comp.summarize() throws an error", {
  expect_error(age.comp.summarize(10))
})

test_that("age.comp.summarize() throws an error when x not of class mcmc.list", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  class(bay.ta_compare) <- class(bay.ta_compare)[-1]
  expect_error(age.comp.summarize(bay.ta_compare))
})

test_that("age.comp.summarize() throws an error when no known age data is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(age.comp.summarize(bay.ta_compare))
})

test_that("age.comp.summarize() throws an error when only NAs for known age data is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(age.comp.summarize(bay.ta_compare, known_age = c(NA, NA)))
})

test_that("age.comp.summarize() throws an error when wrong mode is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_error(age.comp.summarize(bay.ta_compare, known_age = spitalfields$Age,
                        mean_choice = "average"))
})

test_that("age.comp.summarize() produces correct output", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  expect_snapshot(age.comp.summarize(bay.ta_compare, known_age = spitalfields$Age))
})

test_that("age.comp.summarize() produces correct output wit reduced dataset", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  spitalfields_Age <- spitalfields$Age
  spitalfields_Age[c(10,32,56,110,138)] <- NA
  expect_snapshot(age.comp.summarize(bay.ta_compare, known_age = spitalfields_Age))
})
