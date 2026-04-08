# the following script creates the mcmc objects used for testing
sorsum_res_nimble <- bay.ta(framework = "NIMBLE",
                            method = as.matrix(sorsum_as[,2]),
                            numSavedSteps = 1000,
                            seed = 1234)
saveRDS(sorsum_res_nimble,
     file = file.path("./tests/testthat/fixtures", "sorsum_res_nimble.Rds"))

sorsum_res_jags <- bay.ta(framework = "JAGS", method = as.matrix(sorsum_as[,2]),
                          numSavedSteps = 1000,
                          seed = 1234)
saveRDS(sorsum_res_jags,
        file = file.path("./tests/testthat/fixtures", "sorsum_res_jags.Rds"))

spitalfields_res <- bay.ta(framework = "NIMBLE",
                           algorithm = "mnorm",
                            method = as.matrix(spitalfields[,c(2:6)]),
                           parameters = c("b", "a", "beta0", "beta",
                                          "thresh", "age.s", "Ustar"),
                            numSavedSteps = 1000,
                            seed = 1234)
saveRDS(spitalfields_res,
        file = file.path("./tests/testthat/fixtures", "spitalfields_res.Rds"))
