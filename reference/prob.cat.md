# Summed or mean probability densities per category

Summing or averaging probability densities per category. The resulting
data.frames can be used, for example, to produce illustrative diagrams.
See the vignettes for some examples.

## Usage

``` r
prob.cat(
  mcmc_list,
  age_identifier = "age.s",
  group_vec,
  mode = c("mean", "summed")
)
```

## Arguments

- mcmc_list:

  MCMC output from coda chains.

- age_identifier:

  a character string of either "age.s" or "age.s_c" to select the
  uncalibrated or calibrated age estimates. Default: "age.s".

- group_vec:

  a vector specifying the grouping category.

- mode:

  a string specifying the resulting data.frame of summed probabilities
  or mean probabilities per category. Either `mean` or `summed`.

## Value

A data.frame with either probability summed by category or mean per
category.

## Examples

``` r
if (FALSE) { # interactive()

  # select Spitalfields data with multiple traits
  spitalfields_traits <- spitalfields[,c(2:6)]

  # example with multinormal likelihood, please be patient
  spitalfields_res <- bay.ta(framework = "NIMBLE", algorithm = "mnorm",
  method = spitalfields_traits)

  # compute averaging probabilities per category Sex
  prob_cat_mean <- prob.cat(spitalfields_res, group_vec = spitalfields$Sex,
  mode = "mean")

  # compute summed probabilities per category Sex
  prob_cat_summed <- prob.cat(spitalfields_res, group_vec = spitalfields$Sex,
  mode = "summed")
}
```
