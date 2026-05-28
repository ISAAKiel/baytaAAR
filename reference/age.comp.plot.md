# Plots of quality measures of age estimation

Visualisation of the difference between estimated and known age with the
help of a combination of plots.

## Usage

``` r
age.comp.plot(
  x,
  age_identifier = "age.s",
  known_age,
  mean_choice = "Mode",
  hdi_color = c("chartreuse4", "coral2")
)
```

## Arguments

- x:

  output from the function
  [`diagnostic.summary()`](https://isaakiel.github.io/baytaAAR/reference/diagnostic.summary.md).

- age_identifier:

  a character string of either "age.s" or "age.s_c" to select the
  uncalibrated or calibrated age estimates. Default: "age.s".

- known_age:

  a vector of known age-at-death. NAs are allowed and those entries will
  subsequently be ignored.

- mean_choice:

  a character string of either "Mean", "Median" or "Mode". Default:
  "Mode".

- hdi_color:

  a character vector of exactly two entries with color values to
  differentiate estimated ages within the HDI from those outside the
  HDI. Default: c("chartreuse4", "coral2")

## Value

A ggplot object with 2 x 2 single plots, showing:

- Top left:

  Comparison of estimated *highest density intervals* with known ages,
  color1 = age within HDI, color2 = age outside HDI, individuals ordered
  according to known age-at-death.

- Top right:

  Comparison of the density of known ages with a Gompertz function
  derived from the arithmetic mean of the estimated population
  parameters \\\alpha\\ and \\\beta\\.

- Bottom left:

  Scatter plot of known and estimated ages with regression line in blue.
  The dotted line marks perfect equivalence.

- Bottom right:

  Slope of the regression line from the left bottom image (cf.
  goodness-of-fit measure `Residual_slope` from the function
  [`age.comp.summary()`](https://isaakiel.github.io/baytaAAR/reference/age.comp.summary.md)).

## Examples

``` r
if (FALSE) { # interactive()

  # select Spitalfields data with multiple traits
  spitalfields_traits <- spitalfields[,c(2:6)]

  # example with multinormal likelihood, please be patient
  spitalfields_res <- bay.ta(algorithm = "mnorm",
  method = spitalfields_traits)

  # compute age summary statistics
  age.comp.plot(spitalfields_res, known_age = spitalfields$Age)
}
```
