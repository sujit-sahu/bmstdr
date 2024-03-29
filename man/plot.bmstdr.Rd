% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{plot.bmstdr}
\alias{plot.bmstdr}
\title{Plot method for  bmstdr objects.}
\usage{
\method{plot}{bmstdr}(x, segments = TRUE, ...)
}
\arguments{
\item{x}{A bmstdr model fit object.}

\item{segments}{TRUE or FALSE. It decides whether to draw the prediction intervals
as line segments.}

\item{...}{Any other additional arguments.}
}
\value{
It plots the observed values on the original scale 
against the predictions and the 95\% prediction intervals if validation has been 
performed. It then plots the residuals against fitted values. It then applies 
plotting method to the model fitted object as returned by the chosen named package. 
There is no return value.
}
\description{
Plot method for  bmstdr objects.
}
