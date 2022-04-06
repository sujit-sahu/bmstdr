# bmstdr
title: "bmstdr: Bayesian Modeling of Spatio-Temporal Data with R"
author: <a href="https://www.sujitsahu.com/">Sujit K. Sahu </a>
Abstract: This is the R package bmstdr. The package facilitates Bayesian modeling of both point referenced and areal unit spatial 
and spatio-temporal data. Three main functions in the package: Bspatial and Bsptime for spatial and spatio-temporal point 
referenced data respectively and Bcartime for areal unit data, which may also vary in time, perform the main modeling 
and validation tasks. Computations and inference in a Bayesian modeling framework   
are done using packages such as spBayes, spTimer, spTDyn, CARBayes, CARBayesST and also bespoke 
code written in INLA and rstan. The user of bmstdr is afforded the flexibility to choose 
an appropriate package and is also free to name the rows of their input data frame for validation purposes. 
The package is illustrated with five data sets: three on point referenced data on air pollution, 
temperature at the deep ocean and two are areal unit data sets on Covid-19 mortality in England.


  Point reference data are modeled using the Gaussian error distribution only but a top 
  level generalized linear model is used for areal data modeling. The package allows 
  quick comparison of models using model choice criteria, such as DIC and WAIC, 
  and facilitates K-fold cross-validation without much programming effort. 
  Familiar diagnostic plots and model fit exploration using the S3 methods such as summary, 
  residuals and plot are included so that a beginner user confident in model fitting using the 
  base R function lm can quickly learn to analyzing data by fitting a range of 
  appropriate spatial and spatio-temporal models.
