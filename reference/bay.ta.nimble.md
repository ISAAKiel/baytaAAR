# Bayesian Transition Analysis with NIMBLE

Bayesian Transition Analysis with NIMBLE

## Usage

``` r
bay.ta.nimble(
  algorithm,
  method,
  parameters,
  eta = 1,
  gomp_b = NA,
  error_sd = NA,
  minimum_age = 15,
  maximum_age = 100,
  burnInSteps = 2000,
  nChains = 3,
  thinSteps = 1,
  numSteps = 10000,
  seed = FALSE
)
```

## Arguments

- algorithm:

  character string. Either `norm` for 'simple' ordered regression or
  `mnorm` for multinormal ordered regression. Default: `norm`.

- method:

  matrix of integers, converted to matrix if not already matrix. Ordinal
  trait(s) for age estimation.

- parameters:

  vector of character strings. Parameters to monitor.

- eta:

  numeric. Parameter for the LKJ distribution, must be \> 0. Only used
  for multinormal ordered regression for the correlation matrix. `1`
  implies equal correlations, lower values assume stronger correlations.
  Default: `1`.

- gomp_b:

  numeric. Optional prior for parameter Gompertz beta. Default: `NA`.

- error_sd:

  numeric. Optional error parameter for age estimates. Default: `NA`.

- minimum_age:

  numeric. Minimum age for Gompertz distribution. Default: `15`.

- maximum_age:

  numeric. Maximum age for Gompertz distribution. Default: `100`.

- burnInSteps:

  integer. Number of steps for burn-in. Default: `3000`.

- nChains:

  integer. Number of chains. Default: `3`.

- thinSteps:

  integer. Thinning, i. e. which *i*th step should be saved. Default:
  `1` (no thinning).

- numSteps:

  number of steps

- seed:

  integer. Random number for reproducibility. In parallel processing,
  each cluster automatically gets different seeds. If no seed is
  specified, the value is set to today's date as integer.

## Value

A list of MCMC chains of class
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html).

## Examples

``` r
if (FALSE) { # interactive()

  # select Sorsum data with auricular surface after Lovejoy et al. 1985 and
  # convert to matrix
  sorsum <- as.matrix(sorsum_as[,2])

  # example with default settings
  sorsum_res <- bay.ta(method = sorsum)
}
```
