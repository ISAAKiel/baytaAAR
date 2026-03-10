#' Nimble bay.ta
#'
#' @param method matrix
#' @param parameters parameters
#' @param gomp_b Gompertz
#' @param minimum_age min age
#' @param maximum_age maxage
#' @param burnInSteps burn-in
#' @param nChains number of chains
#' @param thinSteps thinning
#' @param numSteps number of steps
#' @param seed seed
#'
#' @return
#' coda object
#'
#' @export
#' @noRd
#'
#' @examples
#' NULL
#'
bay.ta.nimble.norm <- function(
    method,
    parameters,
    gomp_b = NA,
    minimum_age = 15,
    maximum_age = 100,
    burnInSteps = 2000,
    nChains = 3,
    thinSteps = 1,
    numSteps = 10000,
    seed = FALSE
){
  library(nimble)

  n_methods = ncol(method)
  Ntotal = nrow(method)
  nYlevels <- NULL
  thresh <- NULL
  for (i in 1:n_methods) nYlevels[i] <- as.numeric(max(na.omit(method[,i])))
  nthresh <- nYlevels-1
  thresh_init <- matrix(NA, n_methods, max(2,nthresh))
  for (i in 1:n_methods) {
    #if ( nthresh[i] > 1) {
      for (j in 2: max(2,nthresh)) {
        thresh_init[i,j] <- j - 0.5
      } } #}#
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
      ystar_init[i,j] <- k - runif(1, -0.2, 0.2)
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
      beta = runif(n_methods, 0.5, 1),
      beta0 = runif(n_methods, -10, -3),
      age = runif(Ntotal, 20, 40),
      b = runif(1, 0.02, 0.1)
    )
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

  bay_ta <- nimbleCode({
    for ( i in 1:Ntotal ) {
      for (j in 1:n_methods) {
        y[i,j] ~ dinterval(ystar[i,j],thresh[j,1:nthresh[j]])
        mu[i,j] <- beta0[j] + beta[j] * log_age[i]
        ystar[i,j] ~  dnorm( mu[i,j], 1)
      }
      log_age[i] <- log(age.s[i])
      age[i] ~ T(dgomp(b, a), 0, maximum_age - minimum_age)
      age.s[i] <- age[i] + minimum_age
    }

    for (m in 1 : n_methods) {
      beta[m] ~ T(dnorm(0, 1/10^2), 0, )
      beta0[m] ~ T(dnorm( 0 , 1/10^2 ), , 0)
      for ( k in 2: max(2, nthresh[m]) ) {
        thresh[m,k] ~  T(dnorm(thresh_k[m,k], 1/10^2), thresh[m,k-1], )
      }
    }
    b  ~ dunif(gomp_b_beg, gomp_b_end)
    a <- exp(gomp_a0_m * b + gomp_a0_ic)
  })

  bayta_model <- nimbleModel(
    code = bay_ta,
    name = "bayta",
    constants = constantList,
    data = dataList,
    inits = initsList()
  )
  bayta_conf <- configureMCMC(bayta_model, onlySlice = FALSE)
  bayta_conf$addMonitors(parameters)
  comp_model <- compileNimble(bayta_model)
  bayta_MCMC <- buildMCMC(bayta_conf)
  comp_bayta_MCMC <- compileNimble(bayta_MCMC)
  samples <- runMCMC(comp_bayta_MCMC,
          niter = numSteps,
          nburnin = burnInSteps,
          thin = thinSteps,
          nchains = nChains,
          setSeed = seed,
          samplesAsCodaMCMC = TRUE
          )
  return(samples)
}
