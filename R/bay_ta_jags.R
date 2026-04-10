#' Bayesian Transition Analysis with JAGS or NIMBLE
#'
#' This function implements a version of the Bayesian Transition Analysis with
#' MCMCM with JAGS. It presupposes an installed version of JAGS which it
#' interfaces with `runjags` (##).
#'
#' @param runjagsMethod string. Mode to run `runjags`, options: "rjags",
#' "rjparallel", "parallel". Default: "rjags".
#'
#' @rdname bay.ta
#'
#' @export

bay.ta.jags <- function(
    method,
    parameters,
    gomp_b = NA,
    minimum_age = 15,
    maximum_age = 100,
    error_sd = NA,
    adaptSteps = 2000,
    burnInSteps = 3000,
    runjagsMethod ="rjags",
    nChains = 3,
    thinSteps = 1,
    numSavedSteps = 10000,
    silent.jags = F,
    silent.runjags = F,
    seed = seed) {

  checkmate::assertMatrix(method)
  checkmate::assertChoice(runjagsMethod, c("rjags", "rjparallel", "parallel"))
  checkmate::assertCount(nChains, positive = TRUE)

  if(is.na(error_sd)) error_sd =  1

  set.seed(seed)

  n_methods <- ncol(method)
  Ntotal <- nrow(method)
  nYlevels <- c()
  nthresh <- c()
  for (i in 1:n_methods) {
    nYlevels[i] <- as.numeric(max(stats::na.omit(method[,i])))
    nthresh[i] <- max(nYlevels[i]-1, 1)
  }
  thresh_init <- matrix(NA, n_methods, max(nthresh))
  for (i in 1:n_methods) {
    if ( nthresh[i] > 1) {
      for (j in 2: (nthresh[i] )) {
        thresh_init[i,j] <- j - 0.5
      } } }
  thresh <- matrix(NA, n_methods, max(nthresh))
  for (i in 1:n_methods) thresh[i,1] <-  0.5

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
      .RNG.seed = as.integer(sample(1:1e+06, 1)),
      thresh = thresh_init,
      ystar = ystar_init - 1,
      beta = stats::runif(n_methods, 0.5, 1),
      beta0 = stats::runif(n_methods, -10, -3),
      age = stats::runif(Ntotal, 20, 40)
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
      for ( k in 2: (nthresh[m]) ) {
            thresh[m,k] ~ dnorm(k - 0.5, 1/10^2) T(0.5,)
         }
    }
    b  ~ dunif(gomp_b_beg, gomp_b_end)
    a <- exp(gomp_a0_m * b + gomp_a0_ic)
  }
  " # close quote for modelString

  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  runjags::runjags.options(
    silent.jags = silent.jags,
    silent.runjags = silent.runjags
    )

  runJagsOut <- runjags::run.jags( method = runjagsMethod,
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
  codaSamples = coda::as.mcmc.list( runJagsOut )
  return(codaSamples)
}
