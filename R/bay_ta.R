#' Bayesian Transition Analysis with JAGS or NIMBLE
#'
#' \code{bay.ta()} implements latent trait analysis within a Bayesian Markov
#' Chain Monte Carlo (MCMC) framework. It is intended to estimate the
#' age-of-death of adult individuals for whom one or several ordinal traits
#' have been assessed. It produces probability densities for the individual
#' ages but also for the respective population as a whole. \code{bay.ta()} has
#' been introduced and tested in the paper ###.
#' \code{bay.ta()} is a wrapper for the functions \code{bay.ta.jags()} and
#' \code{bay.ta.nimble()}. NIMBLE allows the user to run models with multinormal
#' ordered regression, also with parallel clusters. In this respect, however,
#' JAGS tends to be more stable. The latter presupposes, however, that you have
#' installed JAGS outside of R.
#'
#' @section Data requirements:
#' As input, \code{bay.ta()} assumes a \code{matrix} of trait expressions. In
#' its simplest form, this may contain only one column with a single trait. NAs
#' are allowed but neither must all entries in any of the rows be \code{NA} nor
#' can this be the case for one or several of the columns. \code{bay.ta} will
#' reject to run in such cases, and the offending rows or columns need to be
#' removed from analysis. Please see the article on Chelsea 'Old church' for an
#' example how this can be accomplished.
#' The levels of all traits must start at \code{1}. Binary traits are possible.
#' Mixing of levels like \code{1.5} as short-cut for a trait-expression between
#' \code{1} and \code{2}, however, should be an absolute no-go as this would
#' violate basic principles of ordinal scaling. Thus, for such cases a decision
#' for one of the neighboring levels has to be made or they need to be set to
#' \code{NA}.
#' The nodes (= rows of the matrix) do not have to be fully observed for the
#' multinormal model to run because with
#' \href{https://r-nimble.org/release-notes.html#february-14-2026-weve-released-version-1.4.1}{NIMBLE vers. 1.4.1.},
#' the NIMBLE team introduced a sampler for only partly observed multivariate
#' normal random variables.
#'
#' @param framework character string. Either \code{JAGS} or \code{NIMBLE}.
#'  Default: \code{NIMBLE}.
#' @param algorithm character string. Either \code{norm} for 'simple' ordered
#'   regression or \code{mnorm} for multinormal ordered regression. Default:
#'   \code{norm}.
#' @param multicore \code{TRUE/FALSE}. If \code{TRUE} each chain is assigned to
#'  a dedicated core. Default: \code{FALSE}.
#' @param seed integer. Random number for reproducibility. In parallel
#'   processing, each cluster automatically gets different seeds. If no seed is
#'   specified, the value is set to today's date as integer.
#' @param method matrix of integers. Ordinal trait(s) for age estimation.
#' @param eta numeric. Parameter for the LKJ distribution, must be > 0. Only
#'  used for multinormal ordered regression for the correlation matrix.
#'  \code{1} implies equal correlations, lower values assume stronger
#'   correlations. Default: \code{1}.
#' @param gomp_b numeric. Optional prior for parameter Gompertz beta. Default:
#'   \code{NA}.
#' @param error_sd numeric. Optional error parameter for age estimates. Default:
#'   \code{NA}.
#' @param minimum_age numeric. Minimum age for Gompertz distribution. Default:
#'   \code{15}.
#' @param maximum_age numeric. Maximum age for Gompertz distribution. Default:
#'   \code{100}.
#' @param parameters vector of character strings. Parameters to monitor.
#' @param nChains integer. Number of chains. Default: \code{3}.
#' @param adaptSteps integer. Number of adaptation steps, ignored when
#'   \code{framework} is set to \code{NIMBLE}. Default: \code{2000}.
#' @param burnInSteps integer. Number of steps for burn-in. Default: \code{3000}.
#' @param thinSteps integer. Thinning, i. e. which \emph{i}th step should be
#'   saved. Default: \code{1} (no thinning).
#' @param numSavedSteps integer. Number of saved steps. Default: \code{10000}.
#'  The total number of steps equals \code{thinSteps × numSavedSteps}.
#' @param silent.jags TRUE/FALSE Silent mode to run JAGS. Default: \code{FALSE}.
#'   Ignored when \code{framework} is set to \code{NIMBLE}.
#' @param silent.runjags TRUE/FALSE Silent mode to run runjags. Default:
#'  \code{FALSE}. Ignored when \code{framework} is set to \code{NIMBLE}.
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
#'   # example with default settings
#'   sorsum_res <- bay.ta(method = sorsum)
#'
#'   # example with framework JAGS
#'   sorsum_res <- bay.ta(framework = "JAGS", method = sorsum)
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
#'   # example with multinormal likelihood, please be patient
#'   spitalfields_res <- bay.ta(falgorithm = "mnorm",
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
    nChains = 3,
    adaptSteps = 2000,
    burnInSteps = 3000,
    thinSteps = 1,
    numSavedSteps = 10000,
    silent.jags = F,
    silent.runjags = F) {

  checkmate::assertChoice(framework, c("JAGS", "NIMBLE"))
  checkmate::assertChoice(algorithm, c("norm", "mnorm"))
  checkmate::assertMatrix(method)
  for (i in 1:ncol(method)) {
    checkmate::assertIntegerish(method[,i], all.missing = FALSE, lower = 1) }
  for (i in 1:nrow(method)) {
    checkmate::assertIntegerish(method[i,], all.missing = FALSE, lower = 1) }
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
  checkmate::assertCount(adaptSteps )
  checkmate::assertCount(burnInSteps)
  checkmate::assertCount(thinSteps, positive = TRUE)
  checkmate::assertCount(numSavedSteps)
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
      varlist = c("bay.ta.nimble", "gomp.a0", "shared_args", "seed",
                  "worker_fun"), #"dgomp", "pgomp", "qgomp", "rgomp"
      envir = environment()
    )
    parallel::clusterEvalQ(this_cluster, library("nimble"))

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
