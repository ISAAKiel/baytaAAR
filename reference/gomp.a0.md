# gomp.a0

Internal function for generating starting values for the Gompertz model
if the starting age is not 15 years. Not run if the minimum age is
actually 15. The original formula derives from *Sasaki and Kondo 2016*.

## Usage

``` r
gomp.a0(sampling = 1e+05, b_min = 0.02, b_max = 0.1, minimum_age = 15)
```

## Arguments

- sampling:

  integer. Number of sampling steps. Default: 100000.

- b_min:

  numeric. Minimum of Gompertz \\\beta\\ parameter. Default: 0.02.

- b_max:

  numeric. Maximum of Gompertz \\\beta\\ parameter. Default: 0.1.

- minimum_age:

  numeric. Minimum age in years. Default: 15.

## Value

A vector with coefficients for generating \\\alpha\\ from \\\beta\\, the
parameters of the Gompertz function.

## References

Sasaki T, Kondo O (2016). “An informative prior probability distribution
of the gompertz parameters for bayesian approaches in paleodemography.”
*American Journal of Physical Anthropology*, **159**(3), 523–533.
