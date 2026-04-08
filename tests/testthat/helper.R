skip_if_no_jags <- function() {
  jags_installed <- tryCatch({
    rjags::jags.version()
    TRUE
  }, error = function(e) FALSE)
  if (!jags_installed) {
    skip("JAGS not available")
  }
}
