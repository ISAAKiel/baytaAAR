skip_if_no_jags <- function() {
  jags_installed <- tryCatch({
    rjags::jags.version()
    TRUE
  }, error = function(e) FALSE)
  if (!jags_installed) {
    skip("JAGS not available")
  }
}

skip_if_testcoverage <- function() {
  is_cov <- Sys.getenv("TEST_COVERAGE") == "true"
  if (is_cov) {
    testthat::skip("Skipping during coverage run")
  }
}

skip_if_on_ci_and_mac <- function() {
  is_gha <- Sys.getenv("GITHUB_ACTIONS") == "true"
  is_mac <- Sys.info()[["sysname"]] == "Darwin"

  if (is_gha && is_mac) {
    testthat::skip("Skipping on GitHub Actions macOS (NIMBLE segfault)")
  }
}
