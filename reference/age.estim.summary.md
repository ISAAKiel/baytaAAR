# Summary of age estimates

Convenience function to quickly extract the age-related estimates from
the result of the function `diagnostics.summary()`.

## Usage

``` r
age.estim.summary(x, age_identifier = "age.s")
```

## Arguments

- x:

  output from function `diagnostics.summary()`

- age_identifier:

  a character string of either "age.s" or "age.s_c" to select the
  uncalibrated or calibrated age estimates. Default: "age.s".

## Value

A data.frame with mean, median and mode as well as the HDI ranges as
specified in the output of `diagnostics.summary()` as columns and the
following rows:

- `b` Gompertz parameter \\\beta\\.

- `a` Gompertz parameter \\\alpha\\.

- `M` Modal age, derived from the Gompertz parameters \\\alpha\\ and
  \\\beta\\ according to the equation (1 / \\\beta\\) \* log(\\\beta\\ /
  \\\alpha\\) + minimum_age.

- `age_mean` Mean ages.

- `hdi_diff` *Highest density intervals*.

## Examples

``` r
if (FALSE) { # interactive()
# select Sorsum data with auricular surface after Lovejoy et al. 1985
sorsum <- sorsum_as[,2]

# example with default settings, please be a little bit patient
sorsum_res <- bay.ta(method = sorsum)

sorsum_diag <- diagnostic.summary(sorsum_res)

# show summary of age-related estimates
age.estim.summary(sorsum_diag)
}
```
