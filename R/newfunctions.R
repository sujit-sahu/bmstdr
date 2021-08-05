#' Banerjee, Carlin and Gelfand (2015) Mat'ern covariance function
#' @param d is the distance
#' @param phi is the rate of decay
#' @param nu is the smoothness parameter
#' @return Returns the Mat'ern covariance for distance object d
#' @export
materncov <- function(d, phi, nu) {
  u <- sqrt(2*nu) * d * phi
  covmat <- nu * log(u) + log(besselK(u, nu)) - lgamma(nu) - (nu-1) * log(2)
  diag(covmat) <- 0
  exp(covmat)
}
#' Banerjee et al Mat'ern covariance function
#' @param d is the distance
#' @param phi is the rate of decay
#' @param nu is the smoothness parameter
#' @return Returns the Mat'ern correlation function for distance object d
#' @export
maternfun <- function(d, phi, nu) {
  u <- sqrt(2*nu) * d * phi
  covvec <- nu * log(u) + log(besselK(u, nu)) - lgamma(nu) - (nu-1) * log(2)
  covvec[d==0] <- 1
  exp(covvec)
}
#' This function is used to delete values outside the state of New York
#' @param xyz A list containing the x values, y values and
#' interpolated z values at each x and y pair.
#' @return Returns the input but with NA placed in z values corresponding to
#' the locations whose x-y coordinates are outside the land boundary of the
#' USA.
#' @export
fnc.delete.map.XYZ <- function(xyz){
  x <- xyz$x
  y <- xyz$y
  z <- xyz$z
  xy <- expand.grid(x, y)
  eus<-(maps::map.where(database="state", x=xy[,1], y=xy[,2]))
  dummy <- rep(0, length(xy[,1]))
  eastUS <- NULL
  eastUS <- data.frame(lon=xy[,1],lat=xy[,2],state=eus,dummy=dummy)
  eastUS[!is.na(eastUS[,3]),4] <- 1
  eastUS[eastUS[,3]=="pennsylvania" & !is.na(eastUS[,3]),4]<-0
  eastUS[eastUS[,3]=="new jersey" & !is.na(eastUS[,3]),4]<-0
  eastUS[eastUS[,3]=="connecticut" & !is.na(eastUS[,3]),4]<-0
  eastUS[eastUS[,3]=="massachusetts:main" & !is.na(eastUS[,3]),4]<-0
  eastUS[eastUS[,3]=="new hampshire" & !is.na(eastUS[,3]),4]<-0
  eastUS[eastUS[,3]=="vermont" & !is.na(eastUS[,3]),4]<-0
  a <- eastUS[, 4]
  z <- as.vector(xyz$z)
  z[!a] <- NA
  z <- matrix(z, nrow = length(xyz$x))
  xyz$z <- z
  xyz
}
##

#' Draws a time series (ribbon) plot by combining fitted and predicted values
#' @param yobs A vector of the observed values
#' @param ylow A vector of the lower limits of the fitted or predicted values
#' @param ymed A vector of fitted or predicted values
#' @param yup A vector of the upper limits of the fitted or predicted values
#' @param misst An integer vector which contains the indices of the predicted
#' values
#' @return A ribbon plot,  ggplot2 object,  which shows observed values
#' in red color and open circle, predicted values in blue color and
#' filled circle. 
#' @export
#'
fig11.13.plot <- function(yobs, ylow, ymed, yup, misst) {
  tn <- length(yobs)
  yr <- range(c(yobs, ylow, yup), na.rm=T)
  ymiss <- rep(1, tn)
  ymiss[misst] <- 19 # which values were set aside
  adt <- data.frame(x=1:tn, yobs=yobs, ymed= ymed, ylow=ylow, yup=yup, ymiss=ymiss)
  adt$ymiss <- factor(adt$ymiss, levels = c("1", "19"),
                      labels = c("Observation", "Prediction"))
  obs <- adt[-misst, ]
  preds <- adt[misst, ]
  # library(ggplot2)
  p <- ggplot() +
    geom_point(data = obs, aes(x =x, y = yobs, shape=ymiss), col="red", size = 3) +
    geom_point(data = preds, aes(x =x, y = ymed, shape=ymiss), col ="blue", size = 3) +
    # geom_line(data = obs, aes(x =x, y = yobs), col ="red", size = 1) +
    scale_shape_manual(values = c(1, 19), guides("")) + # Suppress the guides
    geom_line(data = adt, aes(x =x, y = ymed), col ="blue", size = 1) +
    geom_ribbon(data = adt, aes(x =x, ymin =ylow, ymax = yup), alpha = 0.2, color = "grey50") +
    # labs(y = "O3 8hr max", x = "Days") +
    theme(legend.position = c(0.9, 0.9))
  p
}
#' Returns a vector of row numbers for validation.
#' @param sn The total number of spatial locations.
#' @param tn The total number of time points in each location.
#' @param valids A vector of site numbers in (1:sn) to be used for validation.
#' @param validt A vector of time points in (1:tn) to be used for validation.
#' @param allt Whether all the time points should be used for validation.
#' @return Integer vector providing the row numbers of the data frame for validation.
#' Output of this function is suitable as the argument \code{validrows} for the
#' \code{bmstdr} model fitting functions \code{Bsptime, Bcartime}.
#' @examples{
#' # To validate at site numbers 1, 5, and 10 at 31 randomly selected
#' # time points for the nysptime data set we issue the following commands
#' set.seed(44)
#' vt <- sample(62, 31)
#' vrows <- getvalidrows(sn=28, tn=62, valids=c(1, 5, 10), validt=vt)
#' # To validate at sites 1 and 2 at all time points
#' vrows <- getvalidrows(sn=28, tn=62, valids=c(1, 2), allt=TRUE)
#' }
#' @export
getvalidrows <- function(sn, tn, valids, validt=NULL, allt=FALSE) {
  # Assumes data are sorted first by site and then by time
  if (allt) {
    validt <- 1:tn
  } else {
    k <- length(validt)
    if (k==0) stop("Need to provide validation times \n or set allt=T")
  }
  allrows <- matrix(1:(sn*tn), nrow=sn, byrow=TRUE)
  vrows <- sort(c(allrows[valids, validt]))
  # checking
  # dats <- data.frame(sites=rep(1:sn, each=tn), times=rep(1:tn, each=sn), valflag=0)
  # dats$valflag[vrows] <- 1
  # print(dats)
  vrows
}

#' Calculate the hit and false alarm rates
#' @param thresh Threshold value
#' @param yobs A vector of observations, may include missing values
#' @param ypred Predictions
#' @return  A list containing the calculated hit and false alarm rates
#' @export
hitandfalsealarm <- function(thresh, yobs, ypred) {
  # thresh is threshold value
  # yobs : observaytions may include missing
  # ypred: predictions
  dat <- na.omit(data.frame(yobs=yobs, ypred=ypred))
  m <- nrow(dat)
  dat <- dat - thresh
  k <- sign(dat[, 1 ] * dat[, 2])
  htrate <- 100* length(k[k>= 0 ])/m
  ksub <- dat[which(dat[, 1] < 0), ]
  frate <- 100* length(ksub[ksub[,2]>0, 1])/m
  list(hitrate=htrate, falsealarm=frate)
}
