# Bayesian Transition Analysis with JAGS or NIMBLE

`bay.ta()` implements latent trait analysis within a Bayesian Markov
Chain Monte Carlo (MCMC) framework. It is intended to estimate the
age-of-death of adult individuals for whom one or several ordinal traits
have been assessed. It produces probability densities for the individual
ages but also for the respective population as a whole. `bay.ta()` has
been introduced and tested by Müller-Scheeßel et al. (2026).

## Usage

``` r
bay.ta(
  framework = "NIMBLE",
  algorithm = "norm",
  multicore = FALSE,
  seed = as.integer(format(Sys.Date(), "%Y%m%d")),
  method,
  eta = 1,
  gomp_b = NA,
  error_sd = NA,
  minimum_age = 15,
  maximum_age = 100,
  parameters = c("b", "a", "beta0", "beta", "thresh", "age.s"),
  nChains = 3,
  adaptSteps = 2000,
  burnInSteps = 3000,
  thinSteps = 1,
  numSavedSteps = 10000,
  silent.jags = F,
  silent.runjags = F
)
```

## Arguments

- framework:

  character string. Either `JAGS` or `NIMBLE`. Default: `NIMBLE`.

- algorithm:

  character string. Either `norm` for 'simple' ordered regression or
  `mnorm` for multinormal ordered regression. Default: `norm`.

- multicore:

  `TRUE/FALSE`. If `TRUE` each chain is assigned to a dedicated core.
  Default: `FALSE`.

- seed:

  integer. Random number for reproducibility. In parallel processing,
  each cluster automatically gets different seeds. If no seed is
  specified, the value is set to today's date as integer.

- method:

  matrix of integers, converted to matrix if not already matrix. Ordinal
  trait(s) for age estimation.

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

- parameters:

  vector of character strings. Parameters to monitor.

- nChains:

  integer. Number of chains. Default: `3`.

- adaptSteps:

  integer. Number of adaptation steps, ignored when `framework` is set
  to `NIMBLE`. Default: `2000`.

- burnInSteps:

  integer. Number of steps for burn-in. Default: `3000`.

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

## Value

A list of MCMC chains of class
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html).

## Details

`bay.ta()` is a wrapper for the functions
[`bay.ta.jags()`](https://isaakiel.github.io/baytaAAR/reference/bay.ta.jags.md)
and
[`bay.ta.nimble()`](https://isaakiel.github.io/baytaAAR/reference/bay.ta.nimble.md).
NIMBLE allows the user to run models with multinormal ordered
regression, also with parallel clusters. In this respect, however, JAGS
tends to be more stable. The latter presupposes, however, that you have
installed JAGS outside of R.

## Data requirements

As input, `bay.ta()` assumes a `matrix` of trait expressions. In its
simplest form, this may contain only one column with a single trait. NAs
are allowed but neither must all entries in any of the rows be `NA` nor
can this be the case for one or several of the columns. `bay.ta` will
reject to run in such cases, and the offending rows or columns need to
be removed from analysis. Please see the article on Chelsea 'Old church'
for an example how this can be accomplished. The levels of all traits
must start at `1`. Binary traits are possible. Mixing of levels like
`1.5` as short-cut for a trait-expression between `1` and `2`, however,
should be an absolute no-go as this would violate basic principles of
ordinal scaling. Thus, for such cases a decision for one of the
neighboring levels has to be made or they need to be set to `NA`. The
nodes (= rows of the matrix) do not have to be fully observed for the
multinormal model to run because with [NIMBLE vers.
1.4.1.](https://r-nimble.org/release-notes.html#february-14-2026-weve-released-version-1.4.1),
the NIMBLE team introduced a sampler for only partly observed
multivariate normal random variables.

## References

Müller-Scheeßel N, Rinne C, Fuchs K (2026). “A Fully Bayesian Approach
to Adult Skeletal Age Estimation: Multivariate Latent Trait Modeling
with Markov Chain Monte Carlo Sampling.” *American Journal of Biological
Anthropology*.

## Examples

``` r
if (FALSE) { # interactive()

  # select Sorsum data with auricular surface after Lovejoy et al. 1985
  sorsum <- sorsum_as[,2]

  # example with default settings
  sorsum_res <- bay.ta(method = sorsum)

  # example with framework JAGS
  sorsum_res <- bay.ta(framework = "JAGS", method = sorsum)

  # example with framework JAGS and multiple cores (parallel computing)
  sorsum_res <- bay.ta(framework = "JAGS", multicore = TRUE, method = sorsum)

  # example with 10,000 saved iterations and a thinning of 10 (= 100,000
  # iterations)
  sorsum_res <- bay.ta(method = sorsum, numSavedSteps = 10000, thin = 10)

  # select Spitalfields data with multiple traits
  spitalfields_traits <- spitalfields[,c(2:6)]

  # example with multinormal likelihood, please be patient
  spitalfields_res <- bay.ta(falgorithm = "mnorm",
  method = spitalfields_traits)
}
```
