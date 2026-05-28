# Compute thresholds for chains

The computation of the thresholds on the log- and the age-scale is done
outside of the MCMC simulation to spare the computation cost and the
memory. The function returns a
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) which
can be further processed.

## Usage

``` r
threshold.chains(mcmc_list)
```

## Arguments

- mcmc_list:

  MCMC output from coda chains.

## Value

A [`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) for
threshold values of traits on the age scale.

## Examples

``` r
if (FALSE) { # interactive()
# select Sorsum data with auricular surface after Lovejoy et al. 1985
sorsum <- sorsum_as[,2]

# example with default settings, please be patient
sorsum_res <- bay.ta(method = sorsum)

# compute threshold chains
threshold_chains <- threshold.chains(sorsum_res)
}
```
