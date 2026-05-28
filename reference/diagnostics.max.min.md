# Maximum and minimum diagnostic values

Convenience function to quickly extract maximum and mininum diagnostics
values of the function `diagnostics.summary()` over all parameters. The
maximum values of the PSRF should be below 1.1 while the minumum ESS
should be above 10,000. If either of this is not the case, consider to
increase the length of the chains, i. e. the number of iterations.

## Usage

``` r
diagnostics.max.min(x)
```

## Arguments

- x:

  output from function `diagnostics.summary()`

## Value

A data.frame with one row and the following numeric columns:

- `PSRF_max` Maximum value of the *potential scale reduction factor*.

- `PSRF_upper_max` Mximum value of the upper limit of the
  0.95-confidence interval of the PSRF.

- `ESS_min` Minimum of the *effective sample size*.

## Examples

``` r
if (FALSE) { # interactive()
# select Sorsum data with auricular surface after Lovejoy et al. 1985
sorsum <- sorsum_as[,2]

# example with default settings, please be atient
sorsum_res <- bay.ta(method = sorsum)

# compute diagnostics of the MCMC samples
sorsum_diag <- diagnostic.summary(sorsum_res)

# show maximum and minimum values
diagnostics.max.min(sorsum_diag)
}
```
