#' @import spTimer
#' @import spBayes
#' @import rstan
#' @import CARBayes
#' @import CARBayesST
#' @import ggplot2
#' @import MCMCpack
#' @import Rcpp
#' @import methods 
#' @import graphics
#' @import stats 
#' @importFrom mnormt dmnorm
#' @importFrom utils combn
#' @importFrom utils head
#' @importFrom rstan sampling 
#' @importFrom Rdpack reprompt
#' @importFrom inlabru bru_safe_inla
#' @importFrom ggpubr ggarrange
#' @useDynLib bmstdr
NULL
#if(getRversion() >= "2.15.1")  utils::globalVariables(c("."), add=FALSE)
utils::globalVariables(c("nyspatial", "nysptime",  "ydata", "fitvals", "residvals", "up", "low", "Ntrials"))
utils::globalVariables(c("distance", "variogram", "preds", "inornot", "Time", "s.index", "x", "y", "f2"))
NULL

.onLoad <-
  function(libname, pkgname)
  {
    library.dynam(pkgname, pkgname, lib.loc=libname)
  }


.onAttach <-
  function(libname, pkgname)
  {
    ## figureout the version automatically
    library(help=bmstdr)$info[[1]] -> version
    version <- version[pmatch("Version",version)]
    um <- strsplit(version," ")[[1]]
    version <- um[nchar(um)>0][2]
    packageStartupMessage("\n## bmstdr version: ", version," \n")
  }

#' Observed against predicted plot 
#' @param yobs A vector containing the actual observations 
#' @param predsums A data frame containing predictive summary 
#' statistics with the same number of rows as the length of the vector yobs. 
#' The data frame must have columns named as meanpred, medianpred, sd, low and up. 
#' Ideally this argument should be the output of the command 
#' \code{\link{get_validation_summaries}}.  
#' @param segments Logical: whether to draw line segments for the prediction intervals. 
#' @param summarystat Can take one of two values "median" (default) or "mean" 
#' indicating which one to use for the plot.   
#' @param plotit  Logical scalar value: whether to plot the predictions against the observed values.
#' @return Draws a plot only after removing the missing observations.  It also returns a list of two ggplot2 
#' objects: (i) a plot with intervals drawn \code{pwithseg} and (ii) a plot without the segments drawn: 
#' \code{pwithoutseg} and (iii) a simple plot not showing the range of the prediction intervals.    
#' @examples 
#' set.seed(4)
#' vrows <- sample(nrow(nysptime), 100)
#' M1 <- Bsptime(model="lm", formula=y8hrmax~xmaxtemp+xwdsp+xrh, data=nysptime, 
#' validrows=vrows, scale.transform = "SQRT")
#' psums <-  get_validation_summaries(M1$valpreds)
#' oplots <- obs_v_pred_plot(yobs=M1$yobs_preds$y8hrmax, predsum=psums)
#' names(oplots)
#' plot(oplots$pwithoutseg)
#' plot(oplots$pwithseg)
#' @export
obs_v_pred_plot <- function(yobs, predsums, segments=TRUE, summarystat="median", plotit=TRUE) {
## yobs is r by 1
## predsums is r by 4 data frame where the four columns are mean, sd, up and low
#
  if (!is.vector(yobs)) {
    stop("The yobs argument must be a vector\n")
  }
  if (length(yobs) != nrow(predsums)) {
    stop("Number of observed data is not the same as the number of predictions\n")
  }
  needs <- c("meanpred", "medianpred", "low", "up")
  pnames <- colnames(predsums)
  # pnames <- c("meanpred", "medianpred", "low", "up")
  a <- match(x=needs, table=pnames)
  k <- sum(is.na(a))
  if (k>0) stop("Some required prediction summaries are missing from obs_v_pred_plot")
  
  adat <- data.frame(yobs=as.numeric(yobs), predsums)
  adat <- adat[!is.na(adat$yobs), ]
  
  if (summarystat=="median") {
    adat$preds <- adat$medianpred
  } else { 
    adat$preds <- adat$meanpred
  }
  adat$inornot <- 1
  adat$inornot[adat$yobs >= adat$low & adat$yobs <= adat$up] <- 8
  coverage <- round(100*length(adat$inornot[adat$inornot=="8"])/length(adat$inornot), 1)
  
  # adat$inornot <- factor(adat$inornot)
  adat$cols <-"red4"
  adat$cols[adat$inornot>1] <-"grey1"
  adat$inornot <- factor(adat$inornot, levels=c("1", "8"), labels=c("out", "in"))
  
  # First draw a simple plot 
  yr <- range(c(adat$yobs, adat$preds))
  x1 <- seq(from=yr[1], to=yr[2], length=100)
  ddat <- data.frame(x=x1, y=x1)
  p1 <- ggplot() + 
    xlim(yr) + 
    ylim(yr) + 
    geom_point(data=adat, aes(x=yobs, y=preds, shape=inornot), col=adat$cols,  size=1) + 
    #geom_abline(intercept=0, slope=1, col="blue") +
    geom_line(data=ddat, aes(x=x, y=y), col="blue") + 
    labs(x="Observation", y="Prediction", title=paste("Coverage percentage=", coverage)) + 
    theme(legend.position=c(0.05, 0.9)) 
  # plot(p1)
  
 yfullr <- range(c(adat$yobs, adat$preds, adat$low, adat$up))
 arrow <- arrow(length = unit(0.03, "npc"))
 x1 <- seq(from=yfullr[1], to=yfullr[2], length=100)
 ddat <- data.frame(x=x1, y=x1)
 
 p2 <- ggplot() + 
  xlim(yfullr) + 
  ylim(yfullr) + 
  geom_point(data=adat, aes(x=yobs, y=preds, shape=inornot),  col=adat$cols, size=3) + 
  geom_line(data=ddat, aes(x=x, y=y), col="blue") + 
  #geom_abline(intercept=0, slope=1, col="blue") +
  scale_shape_manual(values=c(1,8), guide = guide_legend(reverse=TRUE)) +
  scale_color_manual(values=c("red4", "grey1")) +
  labs(x="Observation", y="Prediction", title=paste("Coverage percentage=", coverage)) + 
  theme(legend.position=c(0.05, 0.9)) 
 
  p3 <- p2 + geom_segment(data=adat, aes(x=yobs, y=preds, xend=yobs, yend=up), col=adat$cols,  linetype=1, arrow=arrow) +
    geom_segment(data=adat, aes(x=yobs, y=preds, xend=yobs, yend=low), col=adat$cols, linetype=1, arrow=arrow) 
  if (segments) { 
    if (plotit) plot(p3)    
  } else { 
    if (plotit) plot(p2)
  }
 return(list(pwithseg=p3, pwithoutseg=p2, pordinary=p1))
}
#'  Obtains parameter estimates from MCMC samples 
#' @param samps A matrix of N by p samples for the p parameters 
#' @param level Desired confidence level - defaults to 95\%. 
#' @return A data frame containing four columns: mean, sd, low (er limit), 
#' and up (per limit) for the p parameters.   
#' @examples 
#' samps <- matrix(rnorm(10000), ncol= 10 )
#' dim(samps)
#' a <- get_parameter_estimates(samps)
#' a
#' b <- get_parameter_estimates(samps, level=98)
#' b
#' @export
get_parameter_estimates <- function(samps, level=95) {
## samps must be N (mcmc) by k parameters
## Returns the parameter estimates from the samples
  lowcut <- (1 - level/100)/2
  upcut <- 1 - lowcut
  means <- apply(samps, 2, mean)
  sds <- apply(samps, 2, sd)
  low <- apply(samps, 2, quantile, probs=lowcut)
  up <- apply(samps, 2, quantile, probs=upcut)
  paramstable <- data.frame(mean=means, sd=sds, low=low, up=up)
 paramstable
}

#'  Obtains suitable validation summary statistics from MCMC samples 
#'  obtained for  validation.   
#' @param samps A matrix of N by p samples for the p parameters 
#' @param level Desired confidence level - defaults to 95\%. 
#' @return A data frame containing five columns: meanpred, 
#' medianpred, sdpred, low (er limit), 
#' and up (per limit) for the p parameters.   
#' @examples 
#' set.seed(4)
#' vrows <- sample(nrow(nysptime), 100)
#' M1 <- Bsptime(model="lm", formula=y8hrmax~xmaxtemp+xwdsp+xrh, data=nysptime, 
#' validrows=vrows, scale.transform = "SQRT")
#' samps<- M1$valpreds
#' valsums <- get_validation_summaries(samps)
#' head(valsums)  
#' @export
get_validation_summaries <- function(samps, level=95) {
  ## samps must be N (mcmc) by k parameters
  ## Returns the parameter estimates from the samples
  lowcut <- (1 - level/100)/2
  upcut <- 1 - lowcut
  means <- apply(samps, 2, mean)
  sds <- apply(samps, 2, sd)
  medians <- apply(samps, 2, median)
  low <- apply(samps, 2, quantile, probs=lowcut)
  up <- apply(samps, 2, quantile, probs=upcut)
  psums <- data.frame(meanpred=means,  sdpred=sds, medianpred=medians, low=low, up=up)
  psums
}
#' Calculates the four validation statistics: RMSE, MAE, CRPS and coverage
#' given the observed values and MCMC iterates. 
#' @param yval A vector containing n observed values of the response 
#' variable. 
#' @param yits A n by N matrix of predictive samples from the 
#' n observations contained in yval.   
#' @param level The nominal coverage level, defaults to 95\%. 
#' @param summarystat Summary statistics to use to calculate the validation 
#' predictions from the samples. It should be a function like mean or 
#' median which can be calculated by R. The default is mean. 
#' @return A list giving the rmse, mae, crps and coverage. 
#' @examples
#' set.seed(4)
#' vrows <- sample(nrow(nysptime), 100)
#' M1 <- Bsptime(model="lm", formula=y8hrmax~xmaxtemp+xwdsp+xrh, data=nysptime, 
#' validrows=vrows, scale.transform = "SQRT")
#' valstats <- calculate_validation_statistics(M1$yobs_preds$y8hrmax, 
#' yits=t(M1$valpreds))
#' unlist(valstats)
#' @export
calculate_validation_statistics <- function(yval, yits, level=95, summarystat="mean"){
  ## yval is the actual n observations
  ## yits is the mcmc samples with dim n by N iteration

  if (!is.vector(yval)) {
    stop("The yobs argument must be a vector\n")
  }
   if (length(yval) != nrow(yits)) {
     cat(length(yval), "length and dim", dim(yits))
    stop("Number of observed data is not the same as the number of prediction variables\n")
   }
  low <- (1 - level/100)/2
  up <- 1 - low
  yval <- as.numeric(yval)
  # if (summarystat="mean") meanpred <- apply(X=yits, 1, mean)
  # else meanpred <- apply(X=yits, 1, median)
  meanpred <- apply(X=yits, 1, summarystat)
  rmse  <- sqrt(mean((yval - meanpred)^2, na.rm=TRUE))
  mae <- mean(abs(yval - meanpred), na.rm=TRUE)
  zup <- apply(yits, 1, quantile, probs=up)
  zlow <- apply(yits, 1, quantile, probs=low)
  cvg <- cal_cvg(vdaty=yval, ylow=zlow, yup=zup)
  # crpsres <- crpsR(yval, yits) ## crps can accept NA's
   tmp <- cbind(yval,yits)
   tmp <- na.omit(tmp)
   crpsres <- crpscpp(tmp) # Use C++ to calculate CRPS
   a <- list(rmse=rmse, mae=mae, crps =crpsres, cvg=cvg)
   results <- list(stats=a)
   results
}


#' Calculates and plots the variogram cloud and an estimated variogram.
#' @param formula Its a formula argument for the response and the coordinates. 
#' @param coordtype Type of coordinates: utm, lonlat or plain with utm 
#' (supplied in meters) as the default. Distance will be calculated in units of kilometer
#' if this argument is either utm or lonlat. Euclidean distance will be calculated 
#' if this is given as the third type plain.  If  distance in meter is to be calculated 
#' then coordtype should be passed on as plain although the coords are supplied in UTM. 
#' @param data A data frame containing the response and the co-ordinates
#' @param nbins Number of bins for the variogram. Default is 30. 
#' @return A list containing:
#'   \itemize{
#'    \item cloud - A data frame containing the variogram cloud. 
#'    This contains pairs of all the data locations, distance 
#'    between the locations and the variogram value for the pair. 
#'    \item variogram  A data frame containing the variogram values in 
#'    each bin.   
#'    \item cloudplot A ggplot2 object of the plot of the  variogram cloud. 
#'    \item variogramplot A ggplot2 object of the plot of the  binned variogram
#'     values. 
#'     }
#' @examples 
#' a <- bmstdr_variogram(data=nyspatial, formula = yo3~utmx + utmy, 
#' coordtype="utm", nb=50)
#' names(a)
#' if (require(ggpubr)) ggarrange(a$cloudplot, a$variogramplot, nrow=1, ncol=2)
#' @export
bmstdr_variogram <- function(formula=yo3 ~ utmx + utmy, coordtype="utm", data=nyspatial, nbins=30)
{
  
  if (!is.data.frame(data)) stop("Need a data frame in the data argument")
  if (!inherits(formula, "formula")) stop("Need a valid formula")
  if (!is.null(coordtype)) coordtype <- match.arg(coordtype, 
                                                  choices=c("utm", "lonlat", "plain"))
  
  
  X <- model.matrix(formula, data=data)
  a <- model.frame(formula=formula, data=data)
  y <- model.response(a)
  z <- data.frame(y=y, X)
  z <- na.omit(z)
  y <- z$y
  
  xnames <- all.vars(formula)[-1]
  kutmx <- z[,xnames[1]]
  kutmy <- z[,xnames[2]]
  
  n <- length(y)
  nc2 <- t(combn(n, 2))
  k <- nrow(nc2)
  bigmat <- matrix(0, nrow=k, ncol=6)
  bigmat[, 1] <- kutmx[nc2[,1]]
  bigmat[, 2] <- kutmy[nc2[,1]]
  bigmat[, 3] <- kutmx[nc2[,2]]
  bigmat[, 4] <- kutmy[nc2[,2]]
  bigmat[, 5] <- as.vector(apply(bigmat[, 1:4],  1,  vdist, coordtype)) #distances in kilometers
  bigmat[, 6] <- ((y[nc2[,1]] - y[nc2[,2]])^2)/2.0 # variogram cloud
  
  colnames(bigmat) <- c("p1.x", "p1.y", "p2.x", "p2.y", "distance", "variogram" )
  
  
  a <- cut(bigmat[,5], nbins)
  mvar <- as.vector(tapply(bigmat[,6], a, mean))
  mdis <- as.vector(tapply(bigmat[,5], a, mean))
  z <- data.frame(distance = mdis, variogram= mvar)
  z <- na.omit(z)
  bigmat <- as.data.frame(bigmat)
  bigmat <- na.omit(bigmat)
  p1 <- ggplot(data=bigmat) + 
    geom_point(aes(x=distance, y=variogram), shape=8)
  plot(p1)
  
  p2 <- ggplot(data=z) + 
    geom_point(aes(x=distance, y=variogram), shape=8) + 
    geom_smooth(method = "loess", se=FALSE, aes(x=distance, y=variogram))
  plot(p2)
  list(cloud=bigmat, variogram=z, cloudplot=p1, variogramplot=p2)
}


#' Grid search method for choosing phi
#' Calculates the validation statistics using the spatial model with a given range of values of
#' the decay parameter phi.
#' @inheritParams Bspatial
#' @param formula An object of class "formula" (or one that can be coerced to that class): 
#' a symbolic description of the model to be fitted.
#' @param data The data frame for which the model formula is to be fitted. 
#' If a spatial model is to be fitted then the data 
#' frame should contain two columns containing the locations 
#' of the coordinates. See the coords argument below.
#' @param  coordtype Type of coordinates: utm, lonlat or plain with utm 
#' (supplied in meters) as the default. 
#' Distance will be calculated in units of kilometer 
#' if this argument is either utm or lonlat. 
#' Euclidean distance will be calculated if this is 
#' given as the third type plain. If distance in meter is to be 
#' calculated then coordtype should be passed on as plain although 
#' the coords are supplied in UTM.
#' @param coords 	A vector of size two identifying the two column 
#' numbers of the data frame to take as coordinates. Or this can 
#' be given as a matrix of number of sites by 2 
#' providing the coordinates of all the data locations.
#' @param scale.transform	Transformation of the response variable. 
#' It can take three values: SQRT, LOG or NONE.
#' @param phis A vector values of phi
#' @param s A vector giving the validation sites
#' @param verbose Logical. Should it print progress? 
#' @param ... Any additional parameter that may be passed to \code{Bspatial}
#' @return  A data frame giving the phi values and the corresponding 
#' validation statistics
#' @examples
#' \donttest{
#' a <- phichoice_sp(formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial, 
#' coordtype="utm", coords=4:5, 
#' phis=seq(from=0.1, to=1, by=0.4), scale.transform="NONE", 
#' s=c(8,11,12,14,18,21,24,28), N=20, burn.in=10, verbose=TRUE)
#' }
#' @export
phichoice_sp <- function(formula, data, coordtype, coords, phis, scale.transform, 
                         s, N, burn.in, verbose=TRUE, ...) {
  n <- length(phis)
  res <- matrix(NA, nrow=n+2, ncol=4)
  d <- Bspatial(model="lm", formula=formula, data=data, 
                coordtype=coordtype, coords=coords, validrows =s, 
                mchoice=FALSE, verbose=verbose, N=N, burn.in=burn.in, ...)
  res[n+2, ] <- unlist(d$stats)
  b <- Bspatial(model="spat",   formula=formula, data=data,
                coordtype=coordtype, coords=coords, validrows =s,  
                mchoice=FALSE, verbose=verbose, N=N, burn.in=burn.in, ...)
  res[n+1, ] <- unlist(b$stats)
  a <- c(phis, b$phi, 0)
  
  # pb <- txtProgressBar(min = 0, max = n, style = 3)   # set progress bar
  for (i in 1:n) {
    if (verbose) cat("Now doing ", i, "to go to ", n, "\n")
    phi <- phis[i]
    b <- Bspatial(model="spat", formula=formula, data=data, 
                  coordtype=coordtype, coords=coords, validrows=s, 
                  verbose=verbose, mchoice=FALSE, phi=phi, N=N, burn.in=burn.in, ...)
    res[i, ] <- unlist(b$stats)
    # setTxtProgressBar(pb, i)
  }
  # close(pb)
  results <- data.frame(phis=a,  res)
  dimnames(results)[[2]] <- c("phis", "rmse", "mae", "crps", "cvg")
  a <- results[order(results$rmse), ]
  a
}

# Grid search method for choosing phi(s) and phi(t)
#' Calculates the validation statistics using the spatial model with a given range of values of
#' the decay parameter phi.
#' @inheritParams Bsptime
#' @param formula An object of class "formula" (or one that can be coerced to that class): 
#' a symbolic description of the model to be fitted.
#' @param data The data frame for which the model formula is to be fitted. 
#' If a spatial model is to be fitted then the data 
#' frame should contain two columns containing the locations 
#' of the coordinates. See the coords argument below.
#' @param  coordtype Type of coordinates: utm, lonlat or plain with utm 
#' (supplied in meters) as the default. 
#' Distance will be calculated in units of kilometer 
#' if this argument is either utm or lonlat. 
#' Euclidean distance will be calculated if this is 
#' given as the third type plain. If distance in meter is to be 
#' calculated then coordtype should be passed on as plain although 
#' the coords are supplied in UTM.
#' @param coords 	A vector of size two identifying the two column 
#' numbers of the data frame to take as coordinates. Or this can 
#' be given as a matrix of number of sites by 2 
#' providing the coordinates of all the data locations.
#' @param scale.transform	Transformation of the response variable. 
#' It can take three values: SQRT, LOG or NONE.
#' @param phis A vector values of phi for spatial decay 
#' @param phit A vector values of phi for temporal decay
#' @param valids A vector giving the validation sites
#' @param verbose Logical. Should progress be printed? 
#' @return  A data frame giving the phi values and the corresponding validation statistics
#' @examples
#' \donttest{
#' a <-  phichoicep(formula=y8hrmax ~ xmaxtemp+xwdsp+xrh, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT",  
#' phis=c(0.001,  0.005, 0.025, 0.125, 0.625), phit=c(0.05, 0.25, 1.25, 6.25), 
#' valids=c(8,11,12,14,18,21,24,28), N=20, burn.in=10, verbose=TRUE)
#' }
#' @export 
phichoicep <- function(formula, data, coordtype, coords, scale.transform, 
                       phis, phit, valids, N, burn.in, verbose=TRUE) {
  a <- expand.grid(phis, phit)
  n <- length(a[,1])
  res <- matrix(NA, nrow=n+2, ncol=4)
  vrows <-  which(nysptime$s.index%in% valids)
  f2 <- y8hrmax ~ xmaxtemp+xwdsp+xrh
  b <- Bsptime(model="separable",formula=formula, data=data,
               coordtype=coordtype, coords=coords, validrows=vrows,
               scale.transform = scale.transform, N=N, burn.in=burn.in, 
               mchoice=FALSE, verbose=verbose)
  a[n+1, ] <- c(b$phi.s, b$phi.t)
  res[n+1, ] <- unlist(b$stats)
  
  b <- Bsptime(model="lm", formula=formula, data=data, coordtype=coordtype, coords=coords, validrows=vrows,
               scale.transform = scale.transform, N=N, burn.in=burn.in, 
               mchoice=FALSE, verbose=verbose)
  a[n+2, ] <- c(0, 0)
  res[n+2, ] <- unlist(b$stats)
  
  for (i in 1:n) {
    if (verbose) cat("Now doing ", i, "to go to ", n, "\n")
    phi.s <- a[i,1]
    phi.t <- a[i,2]
    b <- Bsptime(model="separable", formula=formula, data=data,
                 validrows=vrows, coordtype=coordtype, coords=coords,
                 N=N, burn.in=burn.in, 
                 verbose=verbose, mchoice=FALSE, phi.s=phi.s, phi.t=phi.t,
                 scale.transform =scale.transform)
    res[i, ] <- unlist(b$stats)
    
  }
  results <- data.frame(phis=a[,1], phit=a[,2], res)
  dimnames(results)[[2]] <- c("phis", "phit", "rmse", "mae", "crps", "cvg")
  a <- results[order(results$rmse), ]
  a
}

## Not exporteds 


## #' Find 2.5\%, 50\% and 97.5\% quantiles
## #' @param x a vector
## #' @return median and the 95\% credible interval
## #' @examples
## #' quant(rnorm(1000))
## #' quant(runif(10000))
quant <- function(x){
  quantile(x, prob=c(0.025, 0.5, 0.975))
}

## 
## #' Get X and y from formula and data 
## #' @param formula A model formula
## #' @param data A dataframe
## #' @return Returns a list containing the design matrix X and the 
## #' response data vector y.   
## #' @export
getXy <- function(formula, data) { 
  old <- options()  # Save old options
  on.exit(options(old)) # On exit restore old options
  options(na.action='na.pass')
  X <- model.matrix(formula, data=data)
  a <- model.frame(formula=formula, data=data)
  y <- as.vector(model.response(a))
  list(X=X, y=y)
}



## Evaluate quadratic form a^T V a
## @param a Vector of the quadratic form
## @param V Matrix of the quadratic form 
## @return Returns the value of the  qudratic form as a scalar
quadform <- function(a, V) {
  u <- t(a) %*% V %*% a
  u[1,1]
}

cal_cvg <- function(vdaty, ylow, yup) {
  ## Remove the NA's
  u <- data.frame(vdaty=vdaty, low=ylow, up=yup)
  u <- na.omit(u)
  u$inornot <- 0
  u$inornot[u$vdaty >= u$low & u$vdaty <= u$up] <- 1
  100*sum(u$inornot)/length(u$inornot)
}

##
cal_valstats_from_summary <- function(yval, pred_summary, nsample, level=95) {
  
  ## yval is the actual n observations
  ## pred_summary is summary  (n x 4) matrix of means, sds and intervals
  
  rmse  <- sqrt(mean((yval - pred_summary$mean)^2, na.rm=TRUE))
  mae <- mean(abs(yval - pred_summary$mean), na.rm=TRUE)
  cvg <- cal_cvg(vdaty=yval, ylow=pred_summary$low, yup=pred_summary$up)
  
  ## now sampling to get the crps
  k <- nsample*length(yval)
 
  # yits <- matrix(rnorm(k, mean=0, sd=1), ncol=nsample)
  # use the t_8 distribution as an approximation for the predictive distribution 
  yits <- matrix(rt(n=k, df=8), ncol=nsample)
  sdmat <- matrix(rep(pred_summary$sd, nsample),  ncol=nsample)
  means <- matrix(rep(pred_summary$mean, nsample), ncol=nsample)
  yits <- yits*sdmat + means
  # yits r by nsample
   #crpsvalue  <- crps(yval, yits)
  tmp <- cbind(yval,yits)
  tmp <- na.omit(tmp)
  crpsvalue <- crpscpp(tmp) # Use C++ to calculate CRPS
  a <- list(rmse=rmse, mae=mae, crps=crpsvalue, cvg=cvg)
  b <- list(stats=a, samples=yits)
  b
}

#' Calculate DIC function.
#' Has two arguments: (1) log full likelihood at thetahat and (2) vector of log-likelihood at the theta samples
#' Calculate the DIC criteria values
#' @param loglikatthetahat Log of the likelihood function at theta hat (Bayes). It is a scalar value.
#' @param logliks A vector of log likelihood values at the theta samples
#' @return a list containing four values pdic, pdicalt, dic and dicalt
## #' @export
calculate_dic <- function(loglikatthetahat, logliks) {
  if (!is.vector(logliks)) {
    stop("The logliks argument must be a vector\n")
  }
  expected_log_lik <- mean(logliks)
 # expected_log_lik
  pdic <- 2 * (loglikatthetahat -  expected_log_lik)
  pdicalt <- 2 * var(logliks)
  dicorig <- -2 * loglikatthetahat + 2*pdic
  dicalt <- -2 * loglikatthetahat + 2*pdicalt
  dicresults <- list(pdic=pdic, pdicalt=pdicalt, dicorig=dicorig, dicalt=dicalt)
  dicresults
}

#' Calculate WAIC  function.
#' Has the sole argument component wise likelihood at thetahat samples.
#' v must be a matrix
#' Calculate the WAIC criteria values
#' @param v must be a N (MCMC) by n (data) sample size matrix
#' of the log-likelihood evaluations
#'
#' @return a list containing four values p_waic, p_waic alt, waic and waic_alt
## #' @export
calculate_waic <- function(v) {
  if (!is.matrix(v)) {
    stop("The argument to calculate WAIC must be a 
         N (MCMC) by n (data) sample size matrix\n")
  }
  
  var_log_liks <-  apply(v, 2, var) ## Estimate of Var log( f(y_i| theta) from MCMC
  pwaic2 <- sum(var_log_liks)
  logmin <- min(v)  ## Remember the minimum
  liks <- exp(v - logmin)  # Evaluate likelihood values upto the constant to avoid numerical problems
  av_liks <- apply(liks, 2, mean) ## Estimate of f(y_i| y) from MCMC
  log_av_liks <- log(av_liks) + logmin
  sum_log_av_liks <-  sum(log_av_liks)
  expected_log_liks <- apply(v, 2, mean) ## Estimate of Expected log f(y_i| theta) from MCMC
  pwaic1 <- 2 * sum_log_av_liks - 2 * sum(expected_log_liks)
  waic1 <- -2 * sum_log_av_liks + 2 * pwaic1
  waic2 <- -2 * sum_log_av_liks + 2 * pwaic2
  waic_results <- list(pwaic1=pwaic1, pwaic2=pwaic2, waic1=waic1, waic2=waic2)
  waic_results
}

## Replaced by Rcpp function 
crpsR <- function(xval, xits){
  ## xval is the actual n observations
  ## xits is the mcmc samples with dim n x N iteration
  tmp <- cbind(xval,xits)
  tmp <- na.omit(tmp)
  its <- length(tmp[1,])-1
  fnc <- function(p,n){
    out<-NULL
    out1<-NULL
    for(i in 1:n){
      out[i]<-abs(p[i+1]-p[1])
      #out1[i]<-abs(p[i+1]-sample(p[-c(1,i)],1))*abs(p[i+1]-sample(p[-c(1,i)],1))
      out1[i]<-sum(abs(p[i+1]-p[-1]))
    }
    out<-sum(out)/n
    out1<-sum(out1)/(2*n*n)
    out<-out-out1
    out
  }
  tmp <- apply(tmp,1,fnc,its)
  tmp <- mean(tmp)
  tmp
}

pred_samples_lm <- function(pars, Xmat) { ## likelihood value for a given theta sigma2 as first and second components of pars
  p <- ncol(Xmat)
  n <- nrow(Xmat)
  fitmeans <- Xmat %*% as.vector(pars[1:p])
  unlist(fitmeans+ sqrt(pars[p+1]) * rnorm(n, mean=0, sd=1))
}
pred_samples_sp <- function(pars, Xmat, H) { ##
  p <- ncol(Xmat)
  fitmeans <- Xmat %*% as.vector(pars[1:p])
  u <- mnormt::rmnorm(n=1, mean=fitmeans, varcov=pars[p+1]*H)
  u
}

pred_samples_sptime <- function(pars, y, Xmat, Sinv, Tinv) {
  ## log likelihood value for a given beta sigma2 for WAIC
  p <- ncol(Xmat)
  sn <- nrow(Sinv)
  tn <- nrow(Tinv)
  fitmeans <- as.vector(Xmat %*% as.vector(pars[1:p]))

  qii_inv <- 1/diag(Sinv)
  Qsdiag_inv <- diag(qii_inv, sn, sn)
  Qsdiag_Hinv <-  Qsdiag_inv %*% Sinv

  qtt_inv <- 1/diag(Tinv)
  Qtdiag_inv <- diag(qtt_inv, tn, tn)
  Qtdiag_Tinv <-  Qtdiag_inv %*% Tinv

  err <- y - fitmeans
  materrbs <- matrix(err, byrow=TRUE, ncol=tn) ## This is sn by tn
  temp1 <- Qsdiag_Hinv %*% materrbs
  temp <- temp1 %*% Qtdiag_Tinv ## This is sn by
  cmu <- y - as.vector(t(temp))  ## paralls with y and is the conditional mean
  csd <- sqrt(pars[p+1]*kronecker(qii_inv, qtt_inv))
  cmu + rnorm(sn*tn) * csd
}

pred_samples_sptime2 <- function(pars, Xmat, Ls, Lt) { ## Ls and Lt are Cholesky factorisation of Sigma_s and Sigma_t
  p <- ncol(Xmat)
  fitmeans <- Xmat %*% as.vector(pars[1:p]) ## nT by 1
  sn <- nrow(Ls)
  tn <- nrow(Lt)
  u <- matrix(rnorm(n=sn*tn, sd=sqrt(pars[p+1])), nrow=sn)
  v <- Ls %*% u %*% t(Lt)
  u <- as.vector(t(v)) + fitmeans
  u
}


log_full_likelihood <- function(pars, y, Xmat) { ## log likelihood value for a given parameter value pars
  ## pars is (p+1) by 1 with first p regression components and the last component is sigma^2
  ## Assumes we have y and the X matrix in the current environment
  p <- ncol(Xmat)
  fitmeans <- Xmat %*% as.vector(pars[1:p])
  logdens <- unlist(lapply(y - fitmeans, dnorm, mean=0, sd=sqrt(pars[p+1]), log=TRUE))
  sum(logdens)
}

log_likelihoods_lm <- function(pars, y, Xmat) { ## likelihood values for a given theta sigma2 as first and second components of pars
  ## output is n by 1 vector of log likelihood functions
  p <- ncol(Xmat)
  fitmeans <- Xmat %*% as.vector(pars[1:p])
  logdens <- unlist(lapply(y-fitmeans, dnorm, mean=0, sd=sqrt(pars[p+1]), log=TRUE))
  logdens
}

log_full_likelihood_sp <- function(pars, y, Xmat, H) { ## log likelihood value for a given beta sigma2
  p <- ncol(Xmat)
  fitmeans <- as.vector(Xmat %*% as.vector(pars[1:p]))
  logden <- mnormt::dmnorm(y, mean=fitmeans, varcov =pars[p+1]*H, log=TRUE)
  logden
}


log_likelihoods_sp <- function(pars, y, Xmat, Hinv) {
  ## log likelihood value for just one data point
  ## for WAIC
  p <- ncol(Xmat)
  n <- nrow(Xmat)
  fitmeans <- as.vector(Xmat %*% as.vector(pars[1:p]))
  qii_inv <- 1/diag(Hinv)
  Qdiag_inv <- diag(qii_inv, n, n)
  Qdiag_Hinv <-  Qdiag_inv %*% Hinv
  aivec <-  Qdiag_Hinv %*% fitmeans
  Amat <- diag(1, n, n) - Qdiag_Hinv
  Ay  <- Amat %*% y
  cmu <- Ay + aivec
  csd <- sqrt(pars[p+1]*qii_inv)
  dnorm(y, mean=cmu, sd=csd, log=TRUE)
}


log_likelihoods_sptime <- function(pars, y, Xmat, Sinv, Tinv) {
  ## log likelihood value for a given beta sigma2 for WAIC
  p <- ncol(Xmat)
  sn <- nrow(Sinv)
  tn <- nrow(Tinv)
  fitmeans <- as.vector(Xmat %*% as.vector(pars[1:p]))

  qii_inv <- 1/diag(Sinv)
  Qsdiag_inv <- diag(qii_inv, sn, sn)
  Qsdiag_Hinv <-  Qsdiag_inv %*% Sinv

  qtt_inv <- 1/diag(Tinv)
  Qtdiag_inv <- diag(qtt_inv, tn, tn)
  Qtdiag_Tinv <-  Qtdiag_inv %*% Tinv

  err <- y - fitmeans
  materrbs <- matrix(err, byrow=TRUE, ncol=tn) ## This is sn by tn
  temp1 <- Qsdiag_Hinv %*% materrbs
  temp <- temp1 %*% Qtdiag_Tinv ## This is sn by
  cmu <- y - as.vector(t(temp))  ## paralls with y and is the conditional mean
  csd <- sqrt(pars[p+1]*kronecker(qii_inv, qtt_inv))
  dnorm(y, mean=cmu, sd=csd, log=TRUE)
}



log_full_likelihood_sptime <- function(pars, y, Xmat, Sinv, Tinv, log_detSinv, log_detTinv) {
  ## log full likelihood value for a given beta sigma2 for DIC
  p <- ncol(Xmat)
  sn <- nrow(Sinv)
  tn <- nrow(Tinv)
  fitmeans <- as.vector(Xmat %*% as.vector(pars[1:p]))
  err <- y - fitmeans
  Dmat <- matrix(err, ncol=tn, byrow=TRUE) # sn by tn matrix
  u1 <- Sinv %*% Dmat
  u2 <- Dmat %*% Tinv
  qform <- sum(u1*u2)
  logden <- -0.5 * (sn*tn*log(2 * pi * pars[p+1]) - sn * log_detTinv - tn * log_detSinv + qform/pars[p+1])
  logden
}


## convert seconds into min. hour. and day from spTimer
##
fancy.time<-function(t)
{
  ## convert seconds into min. hour. and day from spTimer
  if(t < 60){
    t <- round(t ,2)
    tt <- paste(t ," - Sec.")
    # message(paste("##\n# Total time taken::",t,"Sec.\n##\n"))
  }
  #
  if(t < (60*60) && t >= 60){
    t1 <- as.integer(t/60)
    t <- round(t-t1*60,2)
    tt <- paste(t1," - Mins.",t," - Sec.")
    # message(paste("##\n# Total time taken::",t1,"Min.",t,"Sec.\n##\n"))
  }
  #
  if(t < (60*60*24) && t >= (60*60)){
    t2 <- as.integer(t/(60*60))
    t <- t-t2*60*60
    t1 <- as.integer(t/60)
    t <- round(t-t1*60,2)
    tt <- paste(t2," - Hour/s.",t1," - Mins.",t," - Sec.")
    # message(paste("##\n# Total time taken::",t2,"Hour/s.",t1,"Min.",t,"Sec.\n##\n"))
  }
  #
  if(t >= (60*60*24)){
    t3 <- as.integer(t/(60*60*24))
    t <- t-t3*60*60*24
    t2 <- as.integer(t/(60*60))
    t <- t-t2*60*60
    t1 <- as.integer(t/60)
    t <- round(t-t1*60,2)
    tt <- paste(t3," - Day/s.",t2," - Hour/s.",t1," - Mins.",t," - Sec.")
    # message(paste("##\n# Total time taken::",t3,"Day/s.",t2,"Hour/s.",t1,"Mins.",t,"Sec.\n##\n"))
  }
  tt <- paste("##\n# Total time taken::", tt, sep=" ")
  tt
}


## Provides the likelihood evaluations for calculating DIC and WAIC
## Provides the conditional log-likelihoods from the marginal
## model. These are used to calculate WAIC.
## Inputs are: y (on the modelling scale), sn, tn, and the stanfitted object
## The output is a list of
## (i)   log full likelihood at theta hat
## (ii)  N (MCMC) dimensional vector of log full likelihood at N theta samples
## (iii) matrix of N (MCMC) by n (observed data points)
# #' @export
logliks_from_gp_marginal_stanfit <- function(y, X, sn, tn,  distmat, stanfit) {

  listofdraws <- rstan::extract(stanfit)
  phi <- listofdraws$phi
  itmax <- length(phi)
  # xbeta <- listofdraws$xbmodel
  
  beta <- listofdraws$beta # N by p 
  xbeta <- beta %*% t(X) # N by n 
  
  tau_sq <- listofdraws$tau_sq
  sigma_sq <- listofdraws$sigma_sq
  zmiss <- listofdraws$z_miss
  ntobs <- length(y[!is.na(y)])
  loglik <- matrix(NA, nrow=itmax, ncol=ntobs)
  yrep <- matrix(NA, nrow=itmax, ncol=ntobs)
  log_full_like_vec <- numeric()

  for (it in 1:itmax) {
    # it <- 1
    sigma2 <- sigma_sq[it]
    tau2   <- tau_sq[it]
    phi_it <- phi[it]
    yimputed  <- y
    yimputed[is.na(y)] <- zmiss[it, ]
    ymat <- matrix(yimputed, byrow=TRUE, ncol=tn)
    Sigma <- diag(tau2, nrow=sn, ncol=sn) + sigma2 * exp(-phi_it * distmat)
    Qmat <- solve(Sigma)
    meanvec <- xbeta[it, ]
    meanmat <- matrix(meanvec, byrow=TRUE, ncol=tn)
    meanmult <- diag(1/diag(Qmat), nrow=sn, ncol=sn) %*% Qmat
    
    condmean <- matrix(NA, nrow=sn, ncol=tn)
    condvar <-  matrix(NA, nrow=sn, ncol=tn)
    errs  <- ymat - meanmat
    cvar <-  1/diag(Qmat)
    tmp <- meanmult %*% errs
    cmean <- ymat - tmp 
    
    udens <- apply(errs, 2,  mnormt::dmnorm, varcov=Sigma, log=TRUE)
    log_full_like_vec[it] <- sum(udens)
    
    condmean_vec <- as.vector(t(cmean))
    condvar_vec <-  rep(cvar, each=tn)
    yobs <- y[!is.na(y)]
    ymean <-   condmean_vec[!is.na(y)]
    yvar <- condvar_vec[!is.na(y)]
    loglik[it, ] <- dnorm(yobs, mean=ymean, sd=sqrt(yvar), log=TRUE)
    yrep[it, ] <- ymean + rnorm(ntobs) * sqrt(yvar)
    
  }
  # print(calculate_waic(loglik))
  
  ## to calculate log full likelihood at theta hat
  ## calculate theta hat first
  
  sigma2 <- mean(sigma_sq)
  tau2   <- mean(tau_sq)
  phi_mean <- mean(phi)
  zmissmean <- apply(zmiss, 2, mean)
  yimputed  <- y
  yimputed[is.na(y)] <- zmissmean
  ymat <- matrix(yimputed, byrow=TRUE, ncol=tn)
  Sigma <- diag(tau2, nrow=sn, ncol=sn) + sigma2 * exp(-phi_mean * distmat)
  meanxbeta <-  apply(xbeta, 2, mean)
  meanmat <- matrix(meanxbeta, byrow=TRUE, ncol=tn)
  errs  <- ymat - meanmat
  tmp <- meanmult %*% errs
  udens <- apply(errs, 2,  mnormt::dmnorm, varcov=Sigma, log=TRUE)
  log_full_like_at_thetahat <- sum(udens)
  

  yrepmeans <- as.vector(apply(yrep, 2, mean))
  yrepvars <- as.vector(apply(yrep, 2, var))
  yobs <- y[!is.na(y)]
  gof <-   sum((yobs-yrepmeans)^2)
  penalty <- sum(yrepvars)
  pmcc <- list(gof=gof, penalty=penalty, pmcc=gof+penalty)

  list(log_full_like_at_thetahat=log_full_like_at_thetahat,  log_full_like_vec=log_full_like_vec,
       loglik=loglik, pmcc=pmcc)
}



## Needs more work
# #' @export
logliks_from_full_gp_spTimer <- function(gpfit) {
  
  y <- as.vector(gpfit$Y)
  
  #u <- gpfit$Y
  #u[1:5, 1] <- 1:5
  #b <- as.vector(u)
  
  if (gpfit$scale.transform == "SQRT") y <- sqrt(y)
  if (gpfit$scale.transform == "LOG") y <- log(y)
  
  X <- gpfit$X
  sn <- gpfit$n
  tn <- sum(gpfit$T)
  distmat <- gpfit$Distance.matrix
  
  phi <- gpfit$phip
  itmax <- length(phi)
  betasamp <- gpfit$betap  # p by itmax
  
  xbeta <- t(X %*% betasamp)  # itmax by nT
  
  tau_sq <- gpfit$sig2ep
  sigma_sq <- gpfit$sig2etap
  #
  nT <- sn*tn
  missing_flag <- rep(0, nT)
  missing_flag[is.na(y)] <- 1
  ntmiss <- sum(missing_flag)
  ntobs <- nT - ntmiss
  
  ##
  ovalues <- gpfit$op
  sig2eps <-  gpfit$sig2ep
  omissing <- ovalues[missing_flag>0, ]
  sige <- sqrt(sig2eps)
  
  sigemat <- matrix(rep(sige, each=ntmiss), byrow=FALSE, ncol=itmax)
  a <- matrix(rnorm(ntmiss*itmax), nrow=ntmiss, ncol=itmax)
  yits <- omissing + a * sigemat
  zmiss <- t(yits)
  ###
  
  loglik <- matrix(NA, nrow=itmax, ncol=ntobs)
  log_full_like_vec <- numeric()
  
  for (it in 1:itmax) {
     # it <- 1
    sigma2 <- sigma_sq[it]
    tau2   <- tau_sq[it]
    phi_it <- phi[it]
    yimputed  <- y
    yimputed[is.na(y)] <- zmiss[it, ]
    ymat <- matrix(yimputed, byrow=TRUE, ncol=tn)
    Sigma <- diag(tau2, nrow=sn, ncol=sn) + sigma2 * exp(-phi_it * distmat)
    Qmat <- solve(Sigma)
    meanvec <- xbeta[it, ]
    meanmat <- matrix(meanvec, byrow=TRUE, ncol=tn)
    meanmult <- diag(1/diag(Qmat), nrow=sn, ncol=sn) %*% Qmat
    
    condmean <- matrix(NA, nrow=sn, ncol=tn)
    condvar <-  matrix(NA, nrow=sn, ncol=tn)
    errs  <- ymat - meanmat
    cvar <-  1/diag(Qmat)
    tmp <- meanmult %*% errs
    cmean <- ymat - tmp 
    
    udens <- apply(errs, 2,  mnormt::dmnorm, varcov=Sigma, log=TRUE)
    log_full_like_vec[it] <- sum(udens)

    condmean_vec <- as.vector(t(cmean))
    condvar_vec <-  rep(cvar, each=tn)
    yobs <- y[!is.na(y)]
    ymean <-   condmean_vec[!is.na(y)]
    yvar <- condvar_vec[!is.na(y)]
    loglik[it, ] <- dnorm(yobs, mean=ymean, sd=sqrt(yvar), log=TRUE)
    
    
  }
  # print(calculate_waic(loglik))
  
  ## to calculate log full likelihood at theta hat
  ## calculate theta hat first
  
  sigma2 <- mean(sigma_sq)
  tau2   <- mean(tau_sq)
  phi_mean <- mean(phi)
  zmissmean <- apply(zmiss, 2, mean)
  yimputed  <- y
  yimputed[is.na(y)] <- zmissmean
  ymat <- matrix(yimputed, byrow=TRUE, ncol=tn)
  Sigma <- diag(tau2, nrow=sn, ncol=sn) + sigma2 * exp(-phi_mean * distmat)
  meanxbeta <-  apply(xbeta, 2, mean)
  meanmat <- matrix(meanxbeta, byrow=TRUE, ncol=tn)
  errs  <- ymat - meanmat
  tmp <- meanmult %*% errs
  udens <- apply(errs, 2,  mnormt::dmnorm, varcov=Sigma, log=TRUE)
  log_full_like_at_thetahat <- sum(udens)
  

  v <- list(log_full_like_at_thetahat=log_full_like_at_thetahat,  log_full_like_vec=log_full_like_vec, loglik=loglik)
}



# #' @export
logliks_from_full_AR_spTimer <- function(gpfit) {

  y <- as.vector(gpfit$Y)

  if (gpfit$scale.transform == "SQRT") y <- sqrt(y)
  if (gpfit$scale.transform == "LOG") y <- log(y)

  X <- gpfit$X
  sn <- gpfit$n
  tn <- sum(gpfit$T)
  distmat <- gpfit$Distance.matrix

  phi <- gpfit$phip
  rho <- gpfit$rhop
  itmax <- length(phi)
  betasamp <- gpfit$betap  # p by itmax

  xbeta <- t(X %*% betasamp)  # itmax by nT

  tau_sq <- gpfit$sig2ep
  sigma_sq <- gpfit$sig2etap
  #
  nT <- sn*tn
  missing_flag <- rep(0, nT)
  missing_flag[is.na(y)] <- 1
  ntmiss <- sum(missing_flag)
  ntobs <- nT - ntmiss

  ##
  ovalues <- gpfit$op
  sig2eps <-  gpfit$sig2ep
  omissing <- ovalues[missing_flag>0, ]
  sige <- sqrt(sig2eps)

  sigemat <- matrix(rep(sige, each=ntmiss), byrow=FALSE, ncol=itmax)
  a <- matrix(rnorm(ntmiss*itmax), nrow=ntmiss, ncol=itmax)
  yits <- omissing + a * sigemat
  zmiss <- t(yits)
  ###

  loglik <- matrix(NA, nrow=itmax, ncol=ntobs)
  log_full_like_vec <- numeric()

  for (it in 1:itmax) {
    #   it <- 1
    sigma2 <- sigma_sq[it]
    tau2   <- tau_sq[it]
    phi_it <- phi[it]
    yimputed  <- y
    yimputed[is.na(y)] <- zmiss[it, ]
    ymat <- matrix(yimputed, byrow=TRUE, ncol=tn)
    Sigma <- diag(tau2, nrow=sn, ncol=sn) # + sigma2 * exp(-phi_it * distmat)
    Qmat <- solve(Sigma)
    ##
    omat <-  matrix(ovalues[, it], byrow=TRUE, ncol=tn)
    meanvec <- xbeta[it, ]
    meanmat <- matrix(meanvec, byrow=TRUE, ncol=tn)
    meanmult <- diag(1/diag(Qmat), nrow=sn, ncol=sn) %*% Qmat

    condmean <- matrix(NA, nrow=sn, ncol=tn)
    condvar <-  matrix(NA, nrow=sn, ncol=tn)
    u <- 0.0
    for (k in 1:tn) {
      if (k==1) oprev <- rep(0, sn)
      else  oprev <- omat[, k-1]

      meanvec <- meanmat[,k] + rho[it] * oprev

      yvec <- ymat[, k]  # sn by 1
      condmean[, k] <- yvec - meanmult %*% (yvec - meanvec)
      condvar[, k] <- 1/diag(Qmat)
      logden_contr <- mnormt::dmnorm(yvec, mean=meanvec, varcov =Sigma, log=TRUE)
      u <- u +  logden_contr
    }
    log_full_like_vec[it] <- u

    condmean_vec <- as.vector(t(condmean))
    condvar_vec <-  as.vector(t(condvar))
    yobs <- y[!is.na(y)]
    ymean <-   condmean_vec[!is.na(y)]
    yvar <- condvar_vec[!is.na(y)]
    loglik[it, ] <- dnorm(yobs, mean=ymean, sd=sqrt(yvar), log=TRUE)
  }
  # print(calculate_waic(loglik))

  ## to calculate log full likelihood at theta hat
  ## calculate theta hat first

  sigma2 <- mean(sigma_sq)
  tau2   <- mean(tau_sq)
  phi_mean <- mean(phi)
  zmissmean <- apply(zmiss, 2, mean)
  yimputed  <- y
  yimputed[is.na(y)] <- zmissmean
  ymat <- matrix(yimputed, byrow=TRUE, ncol=tn)
  Sigma <- diag(tau2, nrow=sn, ncol=sn) # + sigma2 * exp(-phi_mean * distmat)
  meanxbeta <-  apply(xbeta, 2, mean)
  meanmat <- matrix(meanxbeta, byrow=TRUE, ncol=tn)

  u <- 0.0
  for (k in 1:tn) {
    meanvec <- meanmat[,k]
    yvec <- ymat[, k]  # sn by 1
    logden_contr <- mnormt::dmnorm(yvec, mean=meanvec, varcov =Sigma, log=TRUE)
    u <- u +  logden_contr
  }
  log_full_like_at_thetahat <- u

  v <- list(log_full_like_at_thetahat=log_full_like_at_thetahat,  log_full_like_vec=log_full_like_vec, loglik=loglik)

  print(calculate_dic(v$log_full_like_at_thetahat, v$log_full_like_vec))
  print(calculate_waic(v$loglik))
  v
}

# 
vdist <- function(points, coordtype)
{
  point1 <- c(points[1], points[2])
  point2 <- c(points[3], points[4])
  if (coordtype=="lonlat") d <- as.numeric(geodetic.distance(point1, point2))
  else d <- dist(rbind(point1, point2))
  if (coordtype=="utm") d <- d/1000
  d
}
# 
row_dists <- function(u, coordtype) { 
  ## u has the rows of the form lon/lat, lon/lat or utmx/utmy utmx  
  k <- length(u)/2 - 1
  p0 <- u[1:2]
  v <- matrix(u[-(1:2)], byrow=TRUE, ncol=2)
  d <- rep(NA, k)
  for (i in 1:k) { 
    pi <- v[i, ]
    if (coordtype=="lonlat") { 
      d[i] <- geodetic.distance(p0, pi)
    } else { 
      d[i]  <- dist(rbind(p0, pi))
      if (coordtype=="utm")  d[i] <- d[i]/1000
    }
  }
  d
}



#
geodetic.distance <- function(point1, point2)
{
  #The following program computes the distance on the surface
  #of the earth between two points point1 and point2.
  #Both the points are of the form (Longitude, Latitude)
  R <- 6371
  p1rad <- point1 * pi/180
  p2rad <- point2 * pi/180
  d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))
  d <- acos(d)
  as.numeric(R*d)
}
## #' Calculates distance (in kilometer) matrix from either a matrix of long-lat 
## #' pairs or UTM X-Y pairs  
## #' @param coords A two column matrix of coordinates 
## #' @param coordtype Type of coordinates: utm, lonlat or plain with utm 
## #' (supplied in meters) as the default. Distance will be calculated in kilometer
## #' if this argument is either utm or lonlat. Euclidean distance will be calculated 
## #' if this is given as the third type plain.  If  distance in meter is to be calculated 
## #' then coordtype should be passed on as plain although the coords are supplied in UTM. 
## #' @export
dist_mat <- function(coords, coordtype) { 
  n <- nrow(coords)
  if (coordtype=="lonlat") { 
    a <- 1:n 
    b <- 1:n 
    u <- expand.grid(a=a, b=b)
    v <- u[u$a<u$b, ] 
    w <- cbind(coords[v$a, ], coords[v$b, ])
    dvec <- apply(w, 1, vdist, coordtype="lonlat")
    v$dist <- as.vector(dvec)
    vlow <- data.frame(a=v$b, b=v$a, dist=v$dist)
    vdiag <- data.frame(a=1:n, b=1:n, dist=0)
    dall <- rbind(v, vlow, vdiag)
    u$index <- 1:nrow(u)
    dord <- merge(u, dall)
    dall <- dord[order(dord$index), ]
    d <- matrix(dall$dist, n, n)
  } else { 
   d  <- as.matrix(dist(as.matrix(coords)))
   if (coordtype=="utm")  d <- d/1000
  }
  as.matrix(d)
}

# 
g.dist <- function(points)
{
  point1 <- c(points[1], points[2])
  point2 <- c(points[3], points[4])
  #The following program computes the distance on the surface
  #The argument points is of the form (Long, Lat, Long, Lat)
  R <- 6371
  p1rad <- point1 * pi/180
  p2rad <- point2 * pi/180
  d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))
  d <- acos(d)
  as.numeric(R*d)
}
dist_mat_loop <- function(coords, coordtype) { 
  m <- nrow(coords)
  d <- matrix(0, nrow=m, ncol=m)
  if (coordtype=="lonlat") { 
    for (i in 1:(m-1)) { 
      for (j in (i+1):m) {
        # message("i=", i, "j=", j, "\n")
        d[i, j] <-   as.numeric(geodetic.distance(coords[i,], coords[j,]))
        d[j, i] <- d[i, j]
      }
    } 
  } else { 
    d  <- as.matrix(dist(as.matrix(coords)))
    if (coordtype=="utm")  d <- d/1000
  }
  as.matrix(d)
}
## Apply the truncation function. From the spTimer package.  
## @param Y Vector of values on which the truncation will be applied. 
## @param at The value at which truncation occurred. 
## @param lambda The lambda value for the transformation. 
## It is 2 for square-root etc. 
## @param both  Whether truncation at both end points. 
## @return A vector of same length as the input vector Y 
## after applying truncation. 
## @export
truncated.fnc<-function(Y, at=0, lambda=NULL, both=FALSE){
  #
  #
  # at is left-tailed
  #
  if(is.null(lambda)){
    stop("Error: define truncation parameter lambda properly using list ")
  }
  if(is.null(at)){
    stop("Error: define truncation point properly using list ")
  }
  if(at < 0){
    stop("Error: currently truncation point only can take value >= zero ")
  }
  zm <- cbind(Y,Y-at,1)
  zm[zm[, 2] <= 0, 3] <- 0
  zm[zm[, 3] == 0, 2] <- - extraDistr::rhnorm(nrow(zm[zm[,3]==0,]),sd(zm[zm[, 3] != 0, 1],na.rm=TRUE))
  ck <- min(zm[,2],na.rm=TRUE)
  zm[,2] <- zm[,2]-ck
  zm[,2] <- zm[,2]^(1/lambda)
  #zm[, 2] <- zm[, 2]^(1/lambda)
  #zm[zm[, 3] == 0, 2] <- - rhnorm(nrow(zm[zm[,3]==0,]),sd(zm[zm[, 3] != 0, 1]^(1/lambda)))
  #
  if (both == TRUE) {
    list(zm=zm,at2=ck)
  }
  else {
    list(zm=c(zm[,2]),at2=ck)
  }
  #
}
##from spTimer 
## for truncated model
##
## Reverse the truncated function for a vector of values
## @param Y Vector of values on which the truncation will be reversed. 
## @param at The value at which truncation occurred. 
## @param lambda The lambda value for the transformation.  It is 2 for square-root etc. 
## @param at2 is the second at parameter. 
## @return A vector of same length as the input vector Y. 
## @export
reverse.truncated.fnc<-function(Y, at=0, lambda=NULL, at2=0){
  #
  # at is left-tailed
  #
  if(is.null(lambda)){
    stop("Error: define truncation parameter lambda properly using list ")
  }
  #zm <- Y
  #zm[zm <= 0]<- 0
  #zm <- zm^(lambda)
  #zm <- zm + at
  zm <- Y
  zm <- zm^(lambda)
  zm <- zm + at2
  zm[zm <= 0]<- 0
  zm <- zm + at
  zm
  #
}
##
## probability below threshold (for truncated)
##from spTimer 
prob.below.threshold <- function(out, at){
  # out is the N x nItr MCMC samples
  fnc<-function(x, at){
    length(x[x <= at])/length(x)
  }   
  apply(out,1,fnc,at=at)
}
##
##
##
