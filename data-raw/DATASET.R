## pre-computed MCMC output for vignette on Spitalfields

spitalfields_res <- bay.ta(
  framework = "NIMBLE",
  algorithm = "mnorm",
  multicore = F,
  method = spitalfields[,c(2:6)],
  minimum_age = 16,
  parameters = c("b", "a", "beta0", "beta", "thresh", "age.s", "Ustar"),
  thinSteps = 200,
  numSavedSteps = 500,
  seed = 331
)

usethis::use_data(spitalfields_res, internal = TRUE, overwrite = TRUE,
                  compress = "xz")
