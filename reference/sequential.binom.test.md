# Sequential cumulative binomial test

The *cumulative binomial* test asserts if the expected coverage, i.e.
the percentage of known ages within the *highest density intervals*, is
within the confidence interval of the realized coverage. This wrapper
function allows to run this test sequentially, i.e. with a sequence of
expected coverage levels, at a confidence level of 0.95.

## Usage

``` r
sequential.binom.test(
  mcmc_list,
  known_age,
  HDImass = 0.95,
  age_identifier = "age.s"
)
```

## Arguments

- mcmc_list:

  MCMC output from coda chains.

- known_age:

  a vector of known age-at-death. NAs are allowed and those entries will
  subsequently be ignored.

- HDImass:

  a numeric or a vector with the probability range.

- age_identifier:

  a character string of either "age.s" or "age.s_c" to select the
  uncalibrated or calibrated age estimates. Default: "age.s".

## Value

A dataframe with the number of rows equaling the length of the parameter
`HDImass` and six columns as follows:

- `coverage` Expected coverage.

- `n_in` Absolute number of known ages within the *highest density
  intervals*.

- `perc` Realized coverage.

- `CI_low` Lower limit of the confidence interval for the realized
  coverage.

- `CI_up` Upper limit of the confidence interval for the realized
  coverage.

- `p_value` p-value of the binomial test. If significant, the expected
  coverage is outside of the confidence intervals of the realized
  coverage.

## Examples

``` r
if (FALSE) { # interactive()

  # select Spitalfields data with multiple traits and convert to matrix
  spitalfields_traits <- spitalfields[,c(2:6)]

  # example with multinormal likelihood, please be patient
  spitalfields_res <- bay.ta(algorithm = "mnorm",
  method = spitalfields_traits)

  # compute sequential binomial tests at expected probability levels 0.75
  # and 0.95
  sequential.binom.test(spitalfields_res, known_age = spitalfields$Age,
  HDImass = c(0.75, 0.95))
}
```
