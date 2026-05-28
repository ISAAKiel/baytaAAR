# Quality measures of age estimation

Comparison of estimated age and known age-at-death with the help of
several Goodness-of-fit measures. For most of the measures smaller is
better. The only exception is *corrPearson* where larger is better.

## Usage

``` r
age.comp.summary(
  mcmc_list,
  known_age,
  mean_choice = "Mode",
  age_identifier = "age.s",
  ...
)
```

## Arguments

- mcmc_list:

  MCMC output from coda chains.

- known_age:

  a vector of known age-at-death. NAs are allowed and those entries will
  subsequently be ignored.

- mean_choice:

  a character string of either "Mean", "Median" or "Mode". Default:
  "Mode".

- age_identifier:

  a character string of either "age.s" or "age.s_c" to select the
  uncalibrated or calibrated age estimates. Default: "age.s".

- ...:

  Arguments passed on to
  [`diagnostic.summary`](https://isaakiel.github.io/baytaAAR/reference/diagnostic.summary.md)

  `HDImass`

  :   numeric. Value within 0 and 1. Default = 0.95.

  `gelman_diag`

  :   logical. If TRUE, the Gelman-Rubin diagnostics for computing the
      PSRF is invoked. Default: TRUE.

## Value

A data.frame with one row and eight columns with age estimation quality
parameters as follows:

- `Bias` Arithmetic mean of the difference between known and estimated
  age.

- `corrPearson` Correlation of known and estimated age.

- `corr_p` p-value of the correlation of known and estimated age.

- `Residual_slope` Slope of the regression line of the difference
  between known and estimated age.

- `Inaccuracy` Arithmetic mean of the absolute difference between known
  and estimated age.

- `RMSE` *Root mean square error* of known and estimated age.

- `TMNLP` *Test mean log posterior*, a local evaluation of the
  probability density at the point of known age.

- `CRPS` *Continuous ranked probability score*, a global evaluation of
  the probability density at the point of known age.

## Examples

``` r
if (FALSE) { # interactive()

  # select Spitalfields data with multiple traits
  spitalfields_traits <- spitalfields[,c(2:6)]

  # example with multinormal likelihood, please be patient
  spitalfields_res <- bay.ta(framework = "NIMBLE", algorithm = "mnorm",
  method = spitalfields_traits)

  # compute age summary statistics
  age.comp.summary(spitalfields_res, known_age = spitalfields$Age)
}
```
