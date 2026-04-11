#' Bayesian Transition Analysis with JAGS or NIMBLE

#' @param numSteps number of steps
#'
#' @rdname bay.ta
#'
#' @export
#'
bay.ta.nimble <- function(
    algorithm,
    method,
    parameters,
    eta = 1,
    gomp_b = NA,
    error_sd = NA,
    minimum_age = 15,
    maximum_age = 100,
    burnInSteps = 2000,
    nChains = 3,
    thinSteps = 1,
    numSteps = 10000,
    seed = FALSE
){
  n_methods = ncol(method)
  Ntotal = nrow(method)
  nYlevels <- NULL
  thresh <- NULL
  for (i in 1:n_methods) nYlevels[i] <-
    as.numeric(max(stats::na.omit(method[,i])))
  nthresh <- nYlevels-1
  thresh_init <- matrix(NA, n_methods, max(2,nthresh))
  for (i in 1:n_methods) {
    for (j in 2: max(2,nthresh)) {
      thresh_init[i,j] <- j - 0.5
    } }
  thresh <- matrix(NA, n_methods, max(2,nthresh))
  for (i in 1:n_methods) thresh[i,1] <-  0.5

  thresh_k <- thresh_init
  for (i in 1:n_methods) thresh_k[i,1] <-  0.5

  y_init <- matrix(NA, nrow = Ntotal, ncol = n_methods)

  for (j in 1:n_methods) {
    for (i in 1:Ntotal) {
      if (is.na(method[i,j])) {
        y_init[i,j] <- sample(1:nYlevels[j], 1)
      }
    }
  }

  ystar_init <- matrix(NA, nrow = Ntotal, ncol = n_methods)
  for (j in 1:n_methods) {
    for (i in 1:Ntotal) {
      k <- if (!is.na(method[i,j])) method[i,j] else y_init[i,j]
      ystar_init[i,j] <- k - stats::runif(1, -0.2, 0.2)
    }
  }

  if(!is.na(gomp_b)) {
    gomp_b_beg <- gomp_b - 0.001
    gomp_b_end <- gomp_b + 0.001
  } else {
    gomp_b_beg <- 0.02
    gomp_b_end <- 0.1
  }

  # Generate values for Gompertz alpha if minimum age is not 15
  gomp_a0 <- gomp.a0(minimum_age = minimum_age)

  initsList <- function(){
    init_list <- list(
      y = y_init - 1,
      thresh = thresh_init,
      ystar = ystar_init - 1,
      beta = stats::runif(n_methods, 0.5, 1),
      beta0 = stats::runif(n_methods, -10, -3),
      age = stats::runif(Ntotal, 20, 40),
      b = stats::runif(1, 0.02, 0.1)
    )
    if(algorithm == "mnorm") init_list[["Ustar"]] <- diag(rep(1, n_methods))
    return(init_list)
  }

  constantList = list(
    Ntotal = Ntotal,
    nthresh = nthresh,
    n_methods = n_methods,
    minimum_age = minimum_age,
    maximum_age = maximum_age,
    gomp_b_beg = gomp_b_beg,
    gomp_b_end = gomp_b_end,
    gomp_a0_m = gomp_a0[1],
    gomp_a0_ic = gomp_a0[2]
  )

  dataList = list(y = method - 1,
                  thresh = thresh,
                  thresh_k = thresh_k
  )

  bay_ta_common <- nimble::nimbleCode({
    for ( i in 1:Ntotal ) {
      for (j in 1:n_methods) {
        mu[i,j] <- beta0[j] + beta[j] * log_age[i]
      }
      log_age[i] <- log(age.s[i])
      age[i] ~ T(dgomp(b, a), 0, maximum_age - minimum_age)
      age.s[i] <- age[i] + minimum_age
    }
    for (m in 1 : n_methods) {
      beta[m] ~ T(dnorm(0, 1/10^2), 0, )
      beta0[m] ~ T(dnorm( 0 , 1/10^2 ), , 0)
    }
    b  ~ dunif(gomp_b_beg, gomp_b_end)
    a <- exp(gomp_a0_m * b + gomp_a0_ic)
  })

  bay_ta_norm <- nimble::nimbleCode({
    for ( i in 1:Ntotal ) {
      for (j in 1:n_methods) {
      ystar[i,j] ~  dnorm( mu[i,j], 1)
    } } } )

  bay_ta_mnorm <- nimble::nimbleCode({
    for ( i in 1:Ntotal ) {
      ystar[i,1:n_methods] ~
        dmnorm( mu[i,1:n_methods], cholesky = Ustar[1:n_methods,
                                                    1:n_methods],prec_param = 0)
    }
    Ustar[1:n_methods,1:n_methods] ~ dlkj_corr_cholesky(eta, n_methods)
    } )

  bay_ta_error <- nimble::nimbleCode({
    for ( i in 1:Ntotal ) {
      age.s_c[i] ~ T(dnorm(age.s[i], 1/(error_sd)^2), minimum_age, maximum_age)
  }})

    bay_ta_nthresh_single <- nimble::nimbleCode({
    for (m in 1 : n_methods) {
      for ( k in 2: max(2, nthresh) ) {
        thresh[m,k] ~  T(dnorm(thresh_k[m,k], 1/10^2), thresh[m,k-1], )
      }
    }
      for ( i in 1:Ntotal ) {
        for (j in 1:n_methods) {
          y[i,j] ~ dinterval(ystar[i,j],thresh[j,1:nthresh])
        }
      }
    })

      bay_ta_nthresh <- nimble::nimbleCode({
        for (m in 1 : n_methods) {
          for ( k in 2: max(2, nthresh[m]) ) {
            thresh[m,k] ~  T(dnorm(thresh_k[m,k], 1/10^2), thresh[m,k-1], )
          }
        }
        for ( i in 1:Ntotal ) {
          for (j in 1:n_methods) {
            y[i,j] ~ dinterval(ystar[i,j],thresh[j,1:nthresh[j]])
          }
      } })

  if(algorithm == "mnorm") {
    constantList[["eta"]] <- eta
    bay_ta <- nimble::nimbleCode(bay_ta_common, bay_ta_mnorm)
  } else { # simple ordinal regression
    bay_ta <- nimble::nimbleCode(bay_ta_common, bay_ta_norm)
  }

  if(length(nthresh) == 1) {
    bay_ta <- nimble::nimbleCode(bay_ta, bay_ta_nthresh_single)
  } else {
    bay_ta <- nimble::nimbleCode(bay_ta, bay_ta_nthresh)
  }

  if(!is.na(error_sd)) {
    constantList[["error_sd"]] <- error_sd
    bay_ta <- nimble::nimbleCode(bay_ta, bay_ta_error)
  }

      model <- nimble::nimbleModel(
        code = bay_ta,
        constants = constantList,
        data = dataList,
        inits = initsList(),
        check = TRUE
      )

      cmodel <- nimble::compileNimble(model, showCompilerOutput = TRUE)
      conf <- nimble::configureMCMC(model, monitors = parameters)
      mcmc <- nimble::buildMCMC(conf)
      cmcmc <- nimble::compileNimble(mcmc, project = model)

      samples <- nimble::runMCMC(
        cmcmc,
        niter = numSteps,
        nburnin = burnInSteps,
        thin = thinSteps,
        nchains = nChains,
        samplesAsCodaMCMC = TRUE,
        setSeed = seed
      )
  return(samples)
}
