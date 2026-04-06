[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Coverage Status](https://img.shields.io/codecov/c/github/ISAAKiel/baytaAAR/master.svg)](https://app.codecov.io/github/ISAAKiel/baytaAAR?branch=master)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/baytaAAR)](https://CRAN.R-project.org/package=baytaAAR)
[![](https://cranlogs.r-pkg.org/badges/baytaAAR)](https://CRAN.R-project.org/package=baytaAAR)
[![](https://cranlogs.r-pkg.org/badges/grand-total/baytaAAR)](https://CRAN.R-project.org/package=baytaAAR)
[![license](https://img.shields.io/badge/license-GPL%203-B50B82.svg)](https://www.r-project.org/Licenses/GPL-3)

<img src="man/figures/logo.pdf" align="right" height="139" alt="" />

baytaAAR
-------

`baytaAAR` provides Bayesian age estimation for bioarchaeological skeletal  data using ordinal probit regression models implemented in JAGS and NIMBLE. The package is designed to handle multiple ordinal traits of adult individuals and incorporates a Gompertz prior on age to reflect population-level mortality trends. It accounts for observational error and supports full customization of model parameters and MCMC settings. Ideal for analyzing age-at-death in archaeological samples through a probabilistic and reproducible framework.

For further information, please have a look at the Documentation with several vignettes.


How to cite this package
------------

You can cite this package like this "we calculated archaeological life tables with the mortAAR R package
(Müller-Scheeßel et al. 2026)". Here is the full bibliographic reference to include in your reference list (don't forget to update the 'last accessed' date):

> N. Müller-Scheeßel, K. Fuchs, C. Rinne (2026). baytaAAR: Bayesian age estimations of adults (vers. 0.1.0). <https://doi.org/10.32614/CRAN.package.baytaAAR>.


Installation
------------

`baytaAAR` is available on [CRAN](https://CRAN.R-project.org/package=baytaAAR) and can be installed through `install.packages("baytaAAR")`. You can also install the development version with:

```
if(!require('remotes')) install.packages('remotes')
remotes::install_github('ISAAKiel/baytaAAR', build_vignettes = TRUE)
```

Licence
-------

`baytaAAR` is released under the [GNU General Public Licence, version 3](https://www.r-project.org/Licenses/GPL-3). Comments and feedback are welcome, as are code contributions.
