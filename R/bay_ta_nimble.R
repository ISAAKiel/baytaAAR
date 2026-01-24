#' Bayesian Transition Analysis with NIMBLE
#'
#' This is a wrapper for the functions \code{bay.ta.nimble.norm} and
#' \code{bay.ta.nimble.mnorm}. Both implement the Bayesian Transition Analysis
#' with MCMCM with the framework NIMBLE and allow the user to run models with
#' 'simple' ordered regression or multinormal ordered regression, also with
#' parallel clusters.
#'
#' @param algorithm character string. Either \code{norm} for 'simple' ordered
#' regression or \code{mnorm} for multinormal ordered regression. Default:
#' \code{norm}.
#' @param cluster_n integer. Number of clusters for parallel processing chains,
#' maximally the number of cores -1. Default: 1.
#' @param seed integer. Random number for reproducibility. In parallel
#' processing, each cluster automatically gets different seeds.
#' @param method matrix of integers. Ordinal trait(s) for age estimation.
#' @param eta double. Prior for the Cholesky factor, must be > 0. Only used for
#' multinormal ordered regression for the correlation matrix. 1 implies equal
#' correlations, lower values assume higher correlations. Default: 1.
#' The correlation matrix stems from a LKJ distribution. We implemented it
#' according to the nimble manual
#' (https://r-nimble.org/manual/cha-writing-models.html#lkj-distribution-for-correlation-matrices).
#' @param gomp_b double. Optional prior for parameter Gompertz beta. Default:
#' NA.
#' @param minimum_age double. Minimum age for Gompertz distribution. Default:
#' 15.
#' @param maximum_age double. Maximum age for Gompertz distribution. Default:
#' 100.
#' @param parameters vector of character strings. Parameters to monitor.
#' @param nChains integer. Number of chains. Default: 4.
#' @param burnInSteps integer. Number of steps for burn-in. Default: 2000.
#' @param thinSteps integer. Thinning, i. e. which ith step should be saved.
#' Default: 1 (no thinning).
#' @param numSteps integer. Number of steps to run the model. Default: 10000.
#'
#' @return
#' A coda object.
#'
#' @export
#'
#' @examples
#' NULL
#'
#'
bay.ta.nimble  <- function(
    algorithm = "norm",
    cluster_n = 1,
    seed = FALSE,
    method,
    eta = 1,
    gomp_b = NA,
    minimum_age = 15,
    maximum_age = 100,
    parameters,
    nChains = 3,
    burnInSteps = 2000,
    thinSteps = 1,
    numSteps = 10000) {

  checkmate::assertChoice(algorithm, c("norm", "mnorm"))
  checkmate::assertCount(cluster_n, positive = TRUE)
  checkmate::assertMatrix(method)
  available_cores <- parallel::detectCores(logical = FALSE)
  if(cluster_n > available_cores - 1) {
    stop(message(paste0("
    Your machine has only ", available_cores, " physical cores. To maintain usability,
    it is advisable to choose one cluster less than this number.\n")))
  }
  if(algorithm == "mnorm" & ncol(method) < 2) {
    stop(message("With multinormal ordinal regression, there need to be two traits or more.\n"))
  }
  if(algorithm == "mnorm" & anyNA(method)) {
    stop(message("
    With multinormal ordinal regression, all nodes must be observed or unobserved,
    traits with partially NAs (missing data) are not allowed. You have 4 options:\n
    1) Delete all columns which contain NAs.\r
    2) Delete all rows which contain NAs.\r
    3) Impute the missing data manually.\r
    4) Run a multiple ordinal regression (algorithm = \"norm\") instead.\n"))
  }

  start_time <- Sys.time()
  cat("Starting Time:", format(start_time, "%d %b %Y %X"), "\n")

  # single core or multi core
  if(cluster_n == 1) { # single core
    if(algorithm == "norm") { # simple ordinal probit regression
      results <- bay.ta.nimble.norm (
        method = method,
        gomp_b = gomp_b,
        minimum_age = minimum_age,
        maximum_age = maximum_age,
        parameters = parameters,
        nChains = nChains,
        burnInSteps = burnInSteps,
        thinSteps = thinSteps,
        numSteps = ceiling(numSteps / nChains)
      )
    } else {
      results <- bay.ta.nimble.mnorm ( # multinormal ordinal probit regression
        method = method,
        gomp_b = gomp_b,
        eta = eta,
        minimum_age = minimum_age,
        maximum_age = maximum_age,
        parameters = parameters,
        nChains = nChains,
        burnInSteps = burnInSteps,
        thinSteps = thinSteps,
        numSteps = ceiling(numSteps / nChains)
      )
    }

  } else { # multicore

    this_cluster <- parallel::makeCluster(cluster_n)
    on.exit(parallel::stopCluster(this_cluster))

    # Create argument sets
    shared_args_norm <- list(
      method = method,
      gomp_b = gomp_b,
      minimum_age = minimum_age,
      maximum_age = maximum_age,
      parameters = parameters,
      burnInSteps = burnInSteps,
      thinSteps = thinSteps,
      numSteps = ceiling(numSteps / cluster_n),
      nChains = 1
    )

    shared_args_mnorm <- c(shared_args_norm, eta = eta)

    # Export needed functions and objects
    parallel::clusterExport(
      this_cluster,
      varlist = c("bay.ta.nimble.norm", "bay.ta.nimble.mnorm","dgomp", "pgomp",
                  "qgomp", "rgomp","gomp.a0", "shared_args_norm",
                  "shared_args_mnorm", "algorithm", "seed"),
      envir = environment()
    )

    # Worker function
    worker_fun <- function(i, algorithm, args_norm, args_mnorm, seed) {
      if (algorithm == "norm") {
        args_norm$seed <- seed + i
        do.call(bay.ta.nimble.norm, args_norm)
      } else {
        args_mnorm$seed <- seed + i
        do.call(bay.ta.nimble.mnorm, args_mnorm)
      }
    }

    parallel::clusterExport(this_cluster, varlist = "worker_fun", envir = environment())

    results <- parallel::parLapply(
      cl = this_cluster,
      X = 1:cluster_n,
      fun = function(i) worker_fun(i, algorithm, shared_args_norm, shared_args_mnorm, seed)
    )
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
  return(results)
}
