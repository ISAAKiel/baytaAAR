test_that("age.comp.plot() throws an error when x not of class diagnostic_summary", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
  class(bay.ta_compare_diag) <- class(bay.ta_compare_diag)[-1]
  expect_error(age.comp.plot(bay.ta_compare_diag))
})

test_that("age.comp.plot() throws an error when no known age data is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
  expect_error(age.comp.plot(bay.ta_compare_diag))
})

test_that("age.comp.plot() throws an error when only NAs for known age data is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
  expect_error(age.comp.plot(bay.ta_compare_diag, known_age = c(NA, NA)))
})

test_that("age.comp.plot() throws an error when wrong mode is supplied", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
  expect_error(age.comp.plot(bay.ta_compare_diag, known_age = spitalfields$Age,
                                  mean_choice = "average"))
})

test_that("age.comp.plot() produces correct output", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
  vdiffr::expect_doppelganger("spitalfields_age_comparison",
                              suppressWarnings(age.comp.plot(bay.ta_compare_diag,
                                    known_age = spitalfields$Age)))
})

test_that("age.comp.plot() produces correct output wit reduced dataset", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  bay.ta_compare_diag <- diagnostic.summary(bay.ta_compare)
  spitalfields_Age <- spitalfields$Age
  spitalfields_Age[c(10,32,56,110,138)] <- NA
  vdiffr::expect_doppelganger("spitalfields_age_comparison_missing_data",
                              suppressWarnings(age.comp.plot(bay.ta_compare_diag,
                                    known_age = spitalfields_Age)))
})
