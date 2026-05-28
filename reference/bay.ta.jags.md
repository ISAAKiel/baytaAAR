# Bayesian Transition Analysis with JAGS

Bayesian Transition Analysis with JAGS

## Usage

``` r
bay.ta.jags(
  method,
  parameters,
  gomp_b = NA,
  minimum_age = 15,
  maximum_age = 100,
  error_sd = NA,
  adaptSteps = 2000,
  burnInSteps = 3000,
  runjagsMethod = "rjags",
  nChains = 3,
  thinSteps = 1,
  numSavedSteps = 10000,
  silent.jags = F,
  silent.runjags = F,
  seed = seed
)
```

## Arguments

- method:

  matrix of integers, converted to matrix if not already matrix. Ordinal
  trait(s) for age estimation.

- parameters:

  vector of character strings. Parameters to monitor.

- gomp_b:

  numeric. Optional prior for parameter Gompertz beta. Default: `NA`.

- minimum_age:

  numeric. Minimum age for Gompertz distribution. Default: `15`.

- maximum_age:

  numeric. Maximum age for Gompertz distribution. Default: `100`.

- error_sd:

  numeric. Optional error parameter for age estimates. Default: `NA`.

- adaptSteps:

  integer. Number of adaptation steps, ignored when `framework` is set
  to `NIMBLE`. Default: `2000`.

- burnInSteps:

  integer. Number of steps for burn-in. Default: `3000`.

- runjagsMethod:

  string. Mode to run \`runjags\`, options: "rjags", "rjparallel",
  "parallel". Default: "rjags".

- nChains:

  integer. Number of chains. Default: `3`.

- thinSteps:

  integer. Thinning, i. e. which *i*th step should be saved. Default:
  `1` (no thinning).

- numSavedSteps:

  integer. Number of saved steps. Default: `10000`. The total number of
  steps equals `thinSteps × numSavedSteps`.

- silent.jags:

  TRUE/FALSE Silent mode to run JAGS. Default: `FALSE`. Ignored when
  `framework` is set to `NIMBLE`.

- silent.runjags:

  TRUE/FALSE Silent mode to run runjags. Default: `FALSE`. Ignored when
  `framework` is set to `NIMBLE`.

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

  # example with framework JAGS
  sorsum_res <- bay.ta(framework = "JAGS", method = sorsum)
}
```
