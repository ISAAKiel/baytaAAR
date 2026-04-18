utils::globalVariables(c("beta0", "age", "gomp_a0_m", "gomp_a0_ic", "b", "a",
                          "id", "y", "chosen_mean", "HDIhigh", "HDIlow", "runif") )

#' @import nimble
NULL

ignore_unused_imports <- function() {
  Rdpack::get_usage()
}
