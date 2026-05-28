# Extract correlation matrix from Cholesky factor

As the LKJ prior for the correlation matrix uses the Cholesky
decomposition of the correlation matrix, getting the correlation indices
from the coda chains is less straightforward than it seems. It involves
taking the cross product from the resulting coda estimates.

## Usage

``` r
corr.mat.mean(mcmc_list)
```

## Arguments

- mcmc_list:

  MCMC output from coda chains.

## Value

A symmetric matrix with correlations between traits. The number of rows
and columns corresponds to the number of traits.

## Examples

``` r
if (FALSE) { # interactive()

  # select Spitalfields data with multiple traits
  spitalfields_traits <- spitalfields[,c(2:6)]

  # example with multinormal likelihood, please be patient
  spitalfields_res <- bay.ta(algorithm = "mnorm",
  method = spitalfields_traits)

  # compute correlation matrix
  corr.mat.mean(spitalfields_res)
}
```
