utils::globalVariables(c("beta0", "age", "gomp_a0_m", "gomp_a0_ic", "b", "a",
                         "HDIlow", "HDIhigh", "id", "y", "chosen_mean",
                         "runif","nimStop") )

#' @import nimble
NULL

ignore_unused_imports <- function() {
  Rdpack::get_usage()
}
