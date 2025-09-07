# gomp
#' The Gompertz distribution and Gompertz function for internal use in
#' \code{nimble} models

#' @rdname rgomp
#' @export
#' @noRd
rgomp <- nimbleFunction(
  run = function(n = integer(0), b = double(0), a = double(0)) {
    returnType(double(0))

    # Check for valid inputs
    if (n <= 0) stop("Number of samples must be positive")
    if (b <= 0 | a <= 0) stop("Invalid parameters for Gompertz distribution")

    # Generate a single Gompertz-distributed random value
    u <- runif(1)  # Uniform random number
    x <- (1 / b) * log(1 - (b / a) * log(1 - u))  # Inverse CDF sampling

    return(x)
  }
)

#' @rdname dgomp
#' @export
#' @noRd
dgomp <- nimbleFunction(
  run = function(x = double(0), b = double(0), a = double(0),
                 log = integer(0, default = 0)) {  # log.p is now an argument
    returnType(double(0))

    # Compute the density of the Gompertz distribution
    density_val <- a * exp(b * x) * exp(-a / b * (exp(b * x) - 1))

    # Return the log density if log.p = TRUE
    if(log) return(log(density_val))
    else return(density_val)
  }
)

#' @rdname pgomp
#' @export
#' @noRd
pgomp <- nimbleFunction(
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


#' @rdname qgomp
#' @export
#' @noRd
qgomp <- nimbleFunction(
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
    return(-log(1 - p) / b + log(a) / b)
  }
)



# Register the functions in Nimble
registerDistributions(list(
  dgomp = list(
    BUGSdist = "dgomp(b, a)",  # Include log.p in BUGSdist
    Rdist = "dgomp(b, a)",    # Include log.p in Rdist
    pqAvail = TRUE,
    range = c(0, Inf)
  )
))
