# Diagnostic summary of MCMC samples

Summarising diagnostics from a
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html), partly
derived from *Kruschke 2015*.

## Usage

``` r
diagnostic.summary(mcmc_list, HDImass = 0.95, gelman_diag = TRUE)
```

## Arguments

- mcmc_list:

  MCMC output from coda chains.

- HDImass:

  numeric. Value within 0 and 1. Default = 0.95.

- gelman_diag:

  logical. If TRUE, the Gelman-Rubin diagnostics for computing the PSRF
  is invoked. Default: TRUE.

## Value

A data.frame of class `diagnostic_summary` with the row names according
to the parameters to be monitored and the following numeric columns:

- `PSRF Point est.` *Potential scale reduction factor* (= Gelman-Rubin
  statistic), a measure of the mixing of chains.

- `PSRF Upper C.I.` The upper limit of the 0.95-confidence interval of
  the PSRF.

- `Mean` Arithmetic mean of the estimates.

- `Median` Median of the estimates.

- `Mode` Mode of the estimates.

- `ESS` *Effective sample size*, a control of autocorrelation.

- `MCSE` *Monte Carlo standard error*.

- `HDImass` Credibility level of the *highest density interval*.

- `HDIlow` Start of the *highest density interval*.

- `HDIhigh` End of the *highest density interval*.

## Details

Because the first threshold is fixed, the Gelman-Rubin multivariate PSRF
will always throw an error, so this is automatically set to `FALSE`. If
the gelman diagnostics still produce an error, deactivate `gelman_diag`
altogether by setting it to `FALSE`, too.

## References

Kruschke JK (2015). *Doing Bayesian data analysis: a tutorial with R,
JAGS, and Stan*. Academic Press, Amsterdam.

## Examples

``` r
if (FALSE) { # interactive()
# select Sorsum data with auricular surface after Lovejoy et al. 1985
sorsum <- sorsum_as[,2]

# example with default settings, please be patient
sorsum_res <- bay.ta(method = sorsum)

# compute diagnostics of the MCMC samples
sorsum_diag <- diagnostic.summary(sorsum_res)

# show first rows
head(sorsum_diag)
}
```
