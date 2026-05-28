## baytaAAR

\![baytaAAR logo is the ISAAK logo with Gompertz
curve\]![](reference/figures/logo.png)

`baytaAAR` provides Bayesian age estimation for bioarchaeological
skeletal data of human adults using ordinal probit regression models
implemented in JAGS and NIMBLE. The package is designed to handle
multiple ordinal traits and incorporates a Gompertz prior on age to
reflect population-level mortality trends. It accounts for estimation
uncertainties and supports full customization of model parameters and
MCMC settings.

For further information, please have a look at the Documentation with
several vignettes.

## How to cite this package

You can cite this package like this “we estimated age-at-death of human
adults with the baytaAAR R-package (Müller-Scheeßel et al. 2026)”. Here
is the full bibliographic reference to include in your reference list
(don’t forget to update the ‘last accessed’ date):

> N. Müller-Scheeßel, K. Fuchs, C. Rinne (2026). baytaAAR: Bayesian age
> estimations of human adults (vers. 1.0.0).
> <https://github.com/ISAAKiel/baytaAAR>.

## Installation

`baytaAAR` is not yet available on CRAN but you can install the
development version with:

    if(!require('remotes')) install.packages('remotes')
    remotes::install_github('ISAAKiel/baytaAAR', build_vignettes = TRUE)

## Licence

`baytaAAR` is released under the [GNU General Public Licence, version
3](https://www.r-project.org/Licenses/GPL-3). Comments and feedback are
welcome, as are code contributions.
