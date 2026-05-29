# Data input

As input,
[`bay.ta()`](https://isaakiel.github.io/baytaAAR/reference/bay.ta.md)
expects a `matrix` of trait expressions. In its simplest form, this may
contain only one column with a single trait as in the following example
from the Neolithic gallery grave at Sorsum/Germany (“auricular
surface”):

``` r

data(sorsum_as, package = "baytaAAR")
head(sorsum_as)
#>     id auricular_surface
#> 1  700                 4
#> 4  704                 5
#> 9  710                 8
#> 10 711                 4
#> 16 717                 6
#> 17 719                 3
```

In this case, the data are stored as a `data.frame`. If this has not
been done beforehand, the input is internally converted by
[`bay.ta()`](https://isaakiel.github.io/baytaAAR/reference/bay.ta.md) to
a `matrix`. This is necessary even if – as here – there is only one
trait.

Another example is provided with the data from Spitalfields which, among
others variables, contains five columns with traits.

``` r

data(spitalfields, package = "baytaAAR")
head(spitalfields[,2:6])
#>   TO ST MI MA AP
#> 1  1  1  1  1  2
#> 2  2  1  1  1  1
#> 3  1  2  1  1  1
#> 4  3  2  2  1  2
#> 5  1  1  3  1  1
#> 6  2  2  2  1  1
```

NAs are allowed but neither rows nor columns may consist entirely of
`NA` values. `bay.ta` will refuse to run in such cases, and the
offending rows or columns need to be removed before analysis. Please see
the
[`vignette("articles/worked_example")`](https://isaakiel.github.io/baytaAAR/articles/worked_example.md)
on Chelsea ‘Old church’ for an example of how this can be done.

The levels of all traits must start at `1`. Binary traits are possible.
Mixed levels like `1.5` as a short-cut for a trait expression between
`1` and `2`, however, should not be used as this would violate basic
principles of ordinal scaling. Thus, in such cases one of the
neighboring levels must be chosen or they need to be set to `NA`.

In contrast to what has been published before \[@ref_299013\], the nodes
(i.e., rows of the matrix) do not need to be fully observed for the
multinormal model to run, partly because, with [NIMBLE vers.
1.4.1](https://r-nimble.org/release-notes.html#february-14-2026-weve-released-version-1.4.1),
the NIMBLE team introduced a sampler for partially observed multivariate
normal variables.
