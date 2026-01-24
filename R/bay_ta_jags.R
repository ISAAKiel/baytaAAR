#' Bayesian Transition Analysis with JAGS
#'
#' This function implements a version of the Bayesian Transition Analysis with
#' MCMCM with JAGS. It presupposes an installed version of JAGS which it
#' interfaces with runjags (##).
#'
#' @param method a matrix with traits, starting from 1.
#'
#' @param parameters a vector with parameters to monitor.
#'
#' @param gomp_b numeric. Optional Gompertz parameter. Default: NA.
#'
#' @param minimum_age numeric. Minimum age. Default: 15.
#'
#' @param maximum_age numeric. Maximum age. Default: 100.
#'
#' @param error_sd numeric. Optional calibration term for estimation of age.
#' Default: 7.5.
#'
#' @param adaptSteps integer. Number of adaptation steps. Default: 2000.
#'
#' @param burnInSteps integer. Number of adaptation steps. Default: 3000.
#'
#' @param thinSteps integer. Number of thinning steps. Default: 1.
#'
#' @param numSavedSteps integer. Number of saved steps. Default: 10000.
#'
#' @param nChains integer. Number of chains. Default: 3.
#'
#' @param runjagsMethod string. Mode to run runjags, options: "rjags",
#' "rjparallel", "parallel". Default: "rjags".
#'
#' @param silent.jags TRUE/FALSE Silent mode to run jags. Default: FALSE.
#'
#' @param silent.runjags TRUE/FALSE Silent mode to run runjags. Default: FALSE.
#'
#'
#' @return
#' a coda object.
#'
#' #@examples
#'
#' result <- bay.ta.jags(method_matrix, c("beta", "age"))
#'
#' @rdname bay.ta.jags
#' @export

bay.ta.jags <- function(
    method,
    parameters,
    gomp_b = NA,
    minimum_age = 15,
    maximum_age = 100,
    error_sd = 7.5,
    adaptSteps = 2000,
    burnInSteps = 3000,
    runjagsMethod ="rjags",
    nChains = 3,
    thinSteps = 1,
    numSavedSteps = 10000,
    silent.jags = F,
    silent.runjags = F) {

  checkmate::assertMatrix(method)
  checkmate::assertChoice(runjagsMethod, c("rjags", "rjparallel", "parallel"))
  checkmate::assertCount(nChains, positive = TRUE)

  n_methods <- ncol(method)
  Ntotal <- nrow(method)
  nYlevels <- c()
  for (i in 1:n_methods) nYlevels[i] <- as.numeric(max(na.omit(method[,i])))
  nthresh <- nYlevels-1
  thresh_init <- matrix(NA,n_methods, max(nthresh),)
  for (i in 1:n_methods) {
    if ( nthresh[i] > 1) {
      for (j in 2: (nthresh[i] )) {
        thresh_init[i,j] <- j - 0.5
      } } }
  thresh <- matrix(NA,n_methods, max(nthresh),)
  for (i in 1:n_methods) thresh[i,1] <-  0.5

  ones <- rep(1,Ntotal)
  C <- 100000

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
    RNG_list <- c("base::Wichmann-Hill",
                  "base::Marsaglia-Multicarry",
                  "base::Super-Duper",
                  "base::Mersenne-Twister")
    init_list <- list(
      .RNG.name = sample(RNG_list, 1),
      .RNG.seed = sample(1:1e+06, 1),
      thresh = thresh_init,
      ystar = method - runif(1, 0.8, 1.2),
      beta = runif(n_methods, 0.5, 1),
      beta0 = runif(n_methods, -10, -3),
      age = runif(Ntotal, 20, 40)
    )
    return(init_list)
  }
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    Ntotal = Ntotal,
    y = method - 1,
    n_methods = n_methods,
    nthresh = nthresh,
    thresh = thresh,
    minimum_age = minimum_age,
    maximum_age = maximum_age,
    ones = ones,
    C = C,
    error_sd = error_sd,
    gomp_b_beg = gomp_b_beg,
    gomp_b_end = gomp_b_end,
    gomp_a0_m = gomp_a0[1],
    gomp_a0_ic = gomp_a0[2]
  )

  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      for (j in 1:n_methods) {
        y[i,j] ~ dinterval(ystar[i,j],thresh_sort[j,1:nthresh[j]])
        mu[i,j] <- beta0[j] + beta[j] * log_age[i]
        ystar[i,j] ~  dnorm( mu[i,j], 1)
      }
      log_age[i] <- log(age.s[i])
      age[i] ~ dunif(0, maximum_age - minimum_age)
      age.s[i] <- age[i] + minimum_age
      spy[i] <- a * exp(b * age[i]) * exp(-a/b * (exp(b * age[i]) - 1)) / C # implementing Gompertz probability density
      ones[i] ~ dbern( spy[i]  )
      age.s_c[i] ~ dnorm(age.s[i], 1/(error_sd)^2) T(minimum_age, maximum_age)
      }
    for (m in 1 : n_methods) {
        beta[m] ~ dnorm(0, 1/10^2) T(0,)
        beta0[m] ~ dnorm( 0 , 1/10^2 )T(,0)
        thresh_sort[m,1:nthresh[m]] <- sort(thresh[m,1:nthresh[m]])
        for ( k in 1: (nthresh[m] ) ) {
          thresh_age_log[m,k] <- ( thresh[m,k] - beta0[m] ) / (beta[m] + 0.001) # adding 0.001 to prevent division by 0
          thresh_age[m,k] <- exp(thresh_age_log[m,k])
        }
      for ( k in 2: (nthresh[m]) ) {
            thresh[m,k] ~ dnorm(k - 0.5, 1/10^2) T(0.5,)
         }
    }
    b  ~ dunif(gomp_b_beg, gomp_b_end)
    a <- exp(gomp_a0_m * b + gomp_a0_ic)
    M <- 1 / b * log (b/a) + minimum_age
  }
  " # close quote for modelString

  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  runjags.options(
    silent.jags = silent.jags,
    silent.runjags = silent.runjags
    )

  runJagsOut <- run.jags( method = runjagsMethod,
                          model=modelString ,
                          monitor=parameters ,
                          data=dataList ,
                          inits=initsList ,
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps ,
                          sample=ceiling(numSavedSteps/nChains) ,
                          modules = c("dic","glm"),
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  return(codaSamples)
}
