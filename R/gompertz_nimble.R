#' @title Gompertz probability distribution
#'
#' @description
#' The Gompertz distribution and Gompertz function for internal use in
#' \code{NIMBLE} models
#'
#' @param x x
#' @param n n
#' @param b b
#' @param a a
#' @param log log
#' @param q q
#' @param p p
#' @param lower.tail lower.tail
#' @param log.p log.p
#'
#' @return depending on the kind of function, either
#'
#' @noRd
NULL

#' rgomp
#' @noRd
rgomp <- nimble::nimbleFunction(
  run = function(n = integer(0), b = double(0), a = double(0)) {
    returnType(double(0))

    # Check for valid inputs
    if (n <= 0) nimStop("Number of samples must be positive")
    if (b <= 0 | a <= 0) nimStop("Invalid parameters for Gompertz distribution")

    # Generate a single Gompertz-distributed random value
    u <- runif(1)  # Uniform random number
    x <- (1 / b) * log(1 - (b / a) * log(1 - u))  # Inverse CDF sampling

    return(x)
  }
)

#' dgomp
#' @noRd
dgomp <- nimble::nimbleFunction(
  run = function(x = double(0), b = double(0), a = double(0),
                 log = integer(0, default = 0)) {  # log.p is now an argument
    returnType(double(0))

    # Compute the density of the Gompertz distribution
    log_density <- log(a) + b * x - (a / b) * (exp(b * x) - 1)

    # Return the log density if log.p = TRUE
    if(log) return(log_density)
    else return(exp(log_density))
  }
)

#' pgomp
#' @noRd
pgomp <- nimble::nimbleFunction(
  run = function(q = double(0), b = double(0), a = double(0),
                 lower.tail = logical(0), log.p = integer(0, default = 0)) {
    returnType(double(0))

    # Compute CDF of Gompertz distribution
    logp <- -a / b * (exp(b * q) - 1)

    if (!lower.tail) {
      p <- exp(logp)
      if (log.p) return(log(p))
      return(p)
    } else {
      p <- 1 - exp(logp)
      if (!log.p) return(p)
      else return(log(p))
    }
  }
)


#' qgomp
#' @noRd
qgomp <- nimble::nimbleFunction(
  run = function(p = double(0), b = double(0), a = double(0),
                 lower.tail = logical(0),
                 log.p = integer(0, default = 0)) {
    returnType(double(0))

    # Adjust for log.p argument
    if (log.p) {
      p <- exp(p)  # If log.p = TRUE, convert back from log scale
    }

    # Adjust for lower.tail argument
    if (!lower.tail) {
      p <- 1 - p  # If lower.tail = FALSE, reverse the probability
    }

    # Solve for the quantile (inverse CDF)
    return( (1/b) * log(1 - (b/a) * log(1 - p)) )
  }
)



# Register the functions in Nimble
nimble::registerDistributions(list(
  dgomp = list(
    BUGSdist = "dgomp(b, a)",  # Include log.p in BUGSdist
    Rdist = "dgomp(b, a)",    # Include log.p in Rdist
    pqAvail = TRUE,
    range = c(0, Inf)
  )
))
