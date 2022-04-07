# bmstdr: Bayesian Modeling of Spatio-Temporal Data with R

## Author: <a href="https://www.sujitsahu.com/">Sujit K. Sahu </a>

### Date: April 07, 2022

# Introduction:

This is the github page for the `R` package
*[bmstdr](https://CRAN.R-project.org/package=bmstdr)*. This package
facilitates Bayesian modeling of both point referenced and areal unit
data with or without temporal replications. Three main functions in the
package: `Bspatial` for spatial only point referenced data, `Bsptime`
for spatio-temporal point reference data and `Bcartime` for areal unit
data, which may also vary in time, perform the main modeling and
validation tasks. Computations and inference in a Bayesian modeling
framework are done using popular `R` software packages such as
*[spBayes](https://CRAN.R-project.org/package=spBayes)*,
*[spTimer](https://CRAN.R-project.org/package=spTimer)*,
*[spTDyn](https://CRAN.R-project.org/package=spTDyn)*,
*[CARBayes](https://CRAN.R-project.org/package=CARBayes)*,
*[CARBayesST](https://CRAN.R-project.org/package=CARBayesST)* and also
code written using computing platforms *INLA* and
*[rstan](https://CRAN.R-project.org/package=rstan)*.

<p>

Point referenced data are modeled using the Gaussian error distribution
only but a top level generalized linear model is used for areal data
modeling. The user of
*[bmstdr](https://CRAN.R-project.org/package=bmstdr)* is afforded the
flexibility to choose an appropriate package and is also free to name
the rows of their input data frame for validation purposes. The package
incorporates a range of prior distributions allowable in the nominated
packages with default hyper-parameter values. The package allows quick
comparison of models using both model choice criteria, such as DIC and
WAIC, and facilitates K-fold cross-validation without much programming
effort. Familiar diagnostic plots and model fit exploration using the S3
methods such as `summary`, `residuals` and `plot` are included so that a
beginner user confident in model fitting using the base `R` function
`lm` can quickly learn to analyzing data by fitting a range of
appropriate spatial and spatio-temporal models. This vignette
illustrates the package using five built-in data sets. Three of these
are on point referenced data on air pollution and temperature at the
deep ocean and the other two are areal unit data sets on Covid-19
mortality in England.
