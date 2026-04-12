test_that("HDIofMCMC() throws an error", {
  expect_error(HDIofMCMC("ten"))
})

test_that("HDIofMCMC() produces correct output with default HDImass value", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  expect_equal(HDIofMCMC(sampleVec = round(as.matrix(bay.ta_compare)[,"b"],6)),
               round(c(0.020043, 0.080820),6))
})
