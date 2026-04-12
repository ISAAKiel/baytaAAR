test_that("tmnlp() throws an error", {
  expect_error(tmnlp(10))
})

test_that("tmnlp() produces correct output", {
  bay.ta_compare <- readRDS(test_path("fixtures", "spitalfields_res.Rds"))
  x_mcmcMat = as.matrix(bay.ta_compare, chains=TRUE)
  ages <- x_mcmcMat[,grep("^age.s\\[",colnames(x_mcmcMat))]
  expect_equal(tmnlp(spitalfields$Age, ages), 4.4610)
})
