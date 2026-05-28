# Extract thresholds

A convenience function to extract mean thresholds values from the output
of
[`diagnostic.summary()`](https://isaakiel.github.io/baytaAAR/reference/diagnostic.summary.md)
which in turn was derived from a
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html)
computed with
[`threshold.chains()`](https://isaakiel.github.io/baytaAAR/reference/threshold.chains.md)

## Usage

``` r
threshold.matrix(x, mean_choice = "Mode")
```

## Arguments

- x:

  output from function
  [`diagnostic.summary()`](https://isaakiel.github.io/baytaAAR/reference/diagnostic.summary.md)

- mean_choice:

  a character string of either "Mean", "Median" or "Mode". Default:
  "Mode".

## Value

A matrix with threshold values of traits. The number of rows corresponds
to the number of traits, and the number of columns to the maximum number
of levels of one of the traits.

## Examples

``` r
if (FALSE) { # interactive()
# select Sorsum data with auricular surface after Lovejoy et al. 1985
sorsum <- sorsum_as[,2]

# example with default settings, please be patient
sorsum_res <- bay.ta(method = sorsum)

# compute threshold chains
threshold_chains <- threshold.chains(sorsum_res)

# compute summary diagnostics
threshold_diag <- diagnostic.summary(threshold_chains)

# extract threshold matrix (for sorsum only 1 row)
threshold.matrix(threshold_diag)
}
```
