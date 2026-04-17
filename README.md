[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test coverage](https://codecov.io/gh/nmueller18/baytaAAR/graph/badge.svg)](https://app.codecov.io/gh/nmueller18/baytaAAR)[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/baytaAAR)](https://CRAN.R-project.org/package=baytaAAR)
[![CRAN status](https://img.shields.io/cran/v/baytaAAR)](https://CRAN.R-project.org/package=baytaAAR)
[![](https://cranlogs.r-pkg.org/badges/grand-total/baytaAAR)](https://CRAN.R-project.org/package=baytaAAR)
[![license](https://img.shields.io/badge/license-GPL%203-B50B82.svg)](https://www.r-project.org/Licenses/GPL-3)

baytaAAR
-------

<img src="man/figures/logo.png" align="right" height="139" alt="" />

`baytaAAR` provides Bayesian age estimation for bioarchaeological skeletal  data of human adults using ordinal probit regression models implemented in JAGS and NIMBLE. The package is designed to handle multiple ordinal traits of adult individuals and incorporates a Gompertz prior on age to reflect population-level mortality trends. It accounts for estimation uncertainties and supports full customization of model parameters and MCMC settings.

For further information, please have a look at the Documentation with several vignettes.


How to cite this package
------------

You can cite this package like this "we estimated age-at-death with the bataAAR R package (Müller-Scheeßel et al. 2026)". Here is the full bibliographic reference to include in your reference list (don't forget to update the 'last accessed' date):

> N. Müller-Scheeßel, K. Fuchs, C. Rinne (2026). baytaAAR: Bayesian age estimations of adults (vers. 0.1.0). <https://doi.org/10.32614/CRAN.package.baytaAAR>.


Installation
------------

`baytaAAR` is not yet available on CRAN but you can also the development version with:

```
if(!require('remotes')) install.packages('remotes')
remotes::install_github('ISAAKiel/baytaAAR', build_vignettes = TRUE)
```

Licence
-------

`baytaAAR` is released under the [GNU General Public Licence, version 3](https://www.r-project.org/Licenses/GPL-3). Comments and feedback are welcome, as are code contributions.
