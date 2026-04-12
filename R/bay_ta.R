#' Bayesian Transition Analysis with JAGS or NIMBLE
#'
#' \code{bay.ta()} is a wrapper for the functions \code{bay.ta.jags()} and
#' \code{bay.ta.nimble()}. Both implement the Bayesian Transition Analysis with
#' MCMC. NIMBLE allows the user to run models with multinormal ordered regression,
#' also with parallel clusters. In this respect, however, JAGS tends to be more
#' stable. The latter presupposes, however, that you have installed JAGS outside
#' of R.
#'
#' @details blabla
#'
#' @param framework character string. Either \code{JAGS} \code{NIMBLE}. Default:
#'   \code{NIMBLE}.
#' @param algorithm character string. Either \code{norm} for 'simple' ordered
#'   regression or \code{mnorm} for multinormal ordered regression. Default:
#'   \code{norm}.
#' @param multicore TRUE/FALSE. If TRUE each chain is assigned to a dedicated
#'   core. Default: FALSE.
#' @param seed integer. Random number for reproducibility. In parallel
#'   processing, each cluster automatically gets different seeds. If no seed is
#'   specified, the value is set to today's date as integer.
#' @param method matrix of integers. Ordinal trait(s) for age estimation.
#' @param eta numeric. Prior for the Cholesky factor of the LKJ distribution,
#'   must be > 0. Only used for multinormal ordered regression for the correlation
#'   matrix. 1 implies equal correlations, lower values assume stronger
#'   correlations. Default: 1.
#' @param gomp_b numeric. Optional prior for parameter Gompertz beta. Default:
#'   NA.
#' @param error_sd numeric.
#' @param minimum_age numeric. Minimum age for Gompertz distribution. Default:
#'   15.
#' @param maximum_age numeric. Maximum age for Gompertz distribution. Default:
#'   100.
#' @param parameters vector of character strings. Parameters to monitor.
#' @param nChains integer. Number of chains. Default: 3.
#' @param adaptSteps integer. Number of adaptation steps, ignored when
#'   \code{framework} is set to "NIMBLE". Default: 2000.
#' @param burnInSteps integer. Number of steps for burn-in. Default: 3000.
#' @param thinSteps integer. Thinning, i. e. which \emph{i}th step should be
#'   saved. Default: 1 (no thinning).
#' @param numSavedSteps integer. Number of saved steps. Default: 10000. The
#'   total number of steps equals \code{thinSteps × numSavedSteps}.
#' @param silent.jags TRUE/FALSE Silent mode to run JAGS. Default: FALSE.
#'   Ignored when \code{framework} is set to "NIMBLE".
#' @param silent.runjags TRUE/FALSE Silent mode to run runjags. Default: FALSE.
#'   Ignored when \code{framework} is set to "NIMBLE".
#'
#'
#' @return A list of MCMC chains of class \code{coda::mcmc.list}.
#'
#' @export
#'
#' @examplesIf interactive()
#'
#'   # select Sorsum data with auricular surface after Lovejoy et al. 1985 and
#'   # convert to matrix
#'   sorsum <- as.matrix(sorsum_as[,2])
#'
#'   # example with default settings sorsum_res <- bay.ta(method = sorsum)
#'
#'   # example with framework JAGS sorsum_res <- bay.ta(framework = "JAGS",
#'   method = sorsum)
#'
#'   # example with framework JAGS and multiple cores (parallel computing)
#'   sorsum_res <- bay.ta(framework = "JAGS", multicore = TRUE, method = sorsum)
#'
#'   # example with 10,000 saved iterations and a thinning of 10 (= 100,000
#'   # iterations)
#'   sorsum_res <- bay.ta(method = sorsum, numSavedSteps = 10000, thin = 10)
#'
#'   # select Spitalfields data with multiple traits and convert to matrix
#'   spitalfields_traits <- as.matrix(spitalfields[,c(2:6)])
#'
#'   # example with multinormal likelihood
#'   spitalfields_res <- bay.ta(framework = "NIMBLE", algorithm = "mnorm",
#'   method = spitalfields_traits)
#'
bay.ta  <- function(
    framework = "NIMBLE",
    algorithm = "norm",
    multicore = FALSE,
    seed = as.integer(format(Sys.Date(), "%Y%m%d")),
    method,
    eta = 1,
    gomp_b = NA,
    error_sd = NA,
    minimum_age = 15,
    maximum_age = 100,
    parameters = c("b", "a", "beta0", "beta", "thresh", "age.s"),
    nChains = 3L,
    adaptSteps = 2000,
    burnInSteps = 3000,
    thinSteps = 1,
    numSavedSteps = 10000,
    silent.jags = F,
    silent.runjags = F) {

  checkmate::assertChoice(framework, c("JAGS", "NIMBLE"))
  checkmate::assertChoice(algorithm, c("norm", "mnorm"))
  checkmate::assertMatrix(method)
  checkmate::assertLogical(multicore)
  checkmate::assertCount(seed)
  checkmate::assertNumeric(eta, lower = 0)
  checkmate::assertNumeric(gomp_b, lower = 0.01, upper = 0.15)
  checkmate::assertNumeric(error_sd, lower = 0)
  checkmate::assertNumeric(minimum_age, lower = 0)
  checkmate::assertNumeric(maximum_age, lower = minimum_age)
  checkmate::assertSubset(parameters, empty.ok = FALSE,
                          choices = c("b", "a", "beta0", "beta",
                                      "thresh", "age.s", "age.s_c", "Ustar"))
  checkmate::assertCount(nChains, positive = TRUE)
  checkmate::assertNumeric(adaptSteps, lower = 0 )
  checkmate::assertNumeric(burnInSteps, lower = 0 )
  checkmate::assertNumeric(thinSteps, lower = 1)
  checkmate::assertNumeric(numSavedSteps)
  checkmate::assertLogical(silent.jags)
  checkmate::assertLogical(silent.runjags)
  available_cores <- parallel::detectCores(logical = FALSE)
  if(framework == "JAGS" & algorithm == "mnorm") {
    stop(message("Framework is set to JAGS. However, JAGS currently does not
                 support multinormal ordinal regression.\n"))
  }
  if(multicore == TRUE & nChains > available_cores - 1) {
    stop(message(paste0("
    Your machine has only ", available_cores, " physical cores. To maintain usability,
    it is advisable to reduce the number of concurring chains or to set multicore to FALSE.\n")))
  }
  if(algorithm == "mnorm" & ncol(method) < 2) {
    stop(message("With multinormal ordinal regression, there need to be two traits or more.\n"))
  }
  if(multicore == TRUE) {
    runjagsMethod <- "parallel"
  } else {
    runjagsMethod <- "rjags"
  }

  start_time <- Sys.time()
  cat("Starting Time:", format(start_time, "%d %b %Y %X"), "\n")

  # single core or multi core
  if(framework == "JAGS") { # simple ordinal probit regression with JAGS,
    # either single or multicore processing
    results <- bay.ta.jags(
      method = method,
      gomp_b = gomp_b,
      error_sd = error_sd,
      minimum_age = minimum_age,
      maximum_age = maximum_age,
      parameters = parameters,
      nChains = nChains,
      runjagsMethod = runjagsMethod,
      adaptSteps = adaptSteps,
      burnInSteps = burnInSteps,
      thinSteps = thinSteps,
      numSavedSteps = numSavedSteps,
      seed = seed,
      silent.jags = silent.jags,
      silent.runjags = silent.runjags
    )
  } else {
    if(multicore == FALSE) { # single core
      if(length(seed) > 0) {
        nimble.set.seed <- TRUE
        set.seed(seed)
      } else {
        nimble.set.seed <- FALSE
      }
      results <- bay.ta.nimble ( # multinormal ordinal probit regression
        algorithm = algorithm,
        method = method,
        gomp_b = gomp_b,
        error_sd = error_sd,
        eta = eta,
        minimum_age = minimum_age,
        maximum_age = maximum_age,
        parameters = parameters,
        nChains = nChains,
        burnInSteps = burnInSteps,
        thinSteps = thinSteps,
        numSteps = burnInSteps + ceiling(numSavedSteps * thinSteps / nChains),
        seed = nimble.set.seed
      )
    } else { # multicore

    this_cluster <- parallel::makeCluster(nChains)
    on.exit(parallel::stopCluster(this_cluster))
    nChains_ <- nChains

    # Create argument sets
    shared_args <- list(
      algorithm = algorithm,
      method = method,
      gomp_b = gomp_b,
      error_sd = error_sd,
      eta = eta,
      minimum_age = minimum_age,
      maximum_age = maximum_age,
      parameters = parameters,
      burnInSteps = burnInSteps,
      thinSteps = thinSteps,
      numSteps = burnInSteps + ceiling(numSavedSteps * thinSteps / nChains_),
      nChains = 1
    )

    # Worker function
    worker_fun <- function(i, args_nimble, seed) {
      args_nimble$seed <- seed + i
      do.call(bay.ta.nimble, args_nimble)
    }

    # Export needed functions and objects
    parallel::clusterExport(
      this_cluster,
      varlist = c("bay.ta.nimble","dgomp", "pgomp",
                  "qgomp", "rgomp","gomp.a0", "shared_args", "seed",
                  "worker_fun"),
      envir = environment()
    )
    parallel::clusterEvalQ(this_cluster, library("nimble"))

    #parallel::clusterExport(this_cluster, varlist = "worker_fun", envir = environment())

    results <- parallel::parLapply(
      cl = this_cluster,
      X = 1:nChains,
      fun = function(i) worker_fun(i, shared_args, seed)
    )
    }
  }
  # Calculate time difference
  time_diff <- difftime(Sys.time(), start_time, units = "secs")
  time_secs <- as.numeric(time_diff)

  # Convert to readable format
  if (time_secs < 60) {
    cat("Execution Time:", round(time_secs, 2), "seconds\n")
  } else if (time_secs < 3600) {
    cat("Execution Time:", round(time_secs / 60, 2), "minutes\n")
  } else if (time_secs < 86400) {
    cat("Execution Time:", round(time_secs / 3600, 2), "hours\n")
  } else {
    cat("Execution Time:", round(time_secs / 86400, 2), "days\n")
  }

  results <- coda::as.mcmc.list(lapply(results, function(chain) {
    s <- as.matrix(chain)
    a <- s[, "a"]
    b <- s[, "b"]
    M <- (1 / b) * log(b / a) + minimum_age
    coda::as.mcmc(cbind(M = M, s))
  }))
  return(results)
}
