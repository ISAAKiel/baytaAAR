test_that("summarizePost() produces correct output with default HDImass value", {
  bay.ta_compare <- readRDS(test_path("fixtures", "sorsum_res_nimble.Rds"))
  expect_snapshot(summarizePost(as.matrix(bay.ta_compare)[,"b"]))
})
