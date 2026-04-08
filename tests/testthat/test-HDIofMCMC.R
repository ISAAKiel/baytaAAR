test_that("HDIofMCMC() throws an error", {
  expect_error(HDIofMCMC("ten"))
})

test_that("HDIofMCMC() produces correct output with default HDImass value", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  expect_equal(HDIofMCMC(sampleVec = round(as.matrix(bay.ta_compare)[,"b"],8)),
               round(c(0.02009131, 0.08394821),8))
})
