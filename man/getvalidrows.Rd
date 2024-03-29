% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/newfunctions.R
\name{getvalidrows}
\alias{getvalidrows}
\title{Returns a vector of row numbers for validation.}
\usage{
getvalidrows(sn, tn, valids, validt = NULL, allt = FALSE)
}
\arguments{
\item{sn}{The total number of spatial locations.}

\item{tn}{The total number of time points in each location.}

\item{valids}{A vector of site numbers in (1:sn) to be used for validation.}

\item{validt}{A vector of time points in (1:tn) to be used for validation.}

\item{allt}{Whether all the time points should be used for validation.}
}
\value{
Integer vector providing the row numbers of the data frame for validation.
Output of this function is suitable as the argument \code{validrows} for the
\code{bmstdr} model fitting functions \code{Bsptime, Bcartime}.
}
\description{
Returns a vector of row numbers for validation.
}
\examples{
{
# To validate at site numbers 1, 5, and 10 at 31 randomly selected
# time points for the nysptime data set we issue the following commands
set.seed(44)
vt <- sample(62, 31)
vrows <- getvalidrows(sn=28, tn=62, valids=c(1, 5, 10), validt=vt)
# To validate at sites 1 and 2 at all time points
vrows <- getvalidrows(sn=28, tn=62, valids=c(1, 2), allt=TRUE)
}
}
