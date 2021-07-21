#' Bayesian regression model fitting for point referenced spatio-temporal data. 
#' Calculates parameter estimates, validation statistics, and 
#' estimated values of several Bayesian model choice criteria. 
#' @param formula An object of class "formula" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param data The data frame for which the model formula is to be fitted.
#' The data frame should be in long format having one row for each location and  time
#' combination. The data frame must be ordered by time within each site, and should
#' optionally have a column, named s.index,  providing the site indices.
#' Thus the data,  with n sites and T times within each site, should be
#' organised in the order: (s1, t1), (s1, t2), ... (s1, T), ... (sn, t1), ... (sn, T). 
#' The data frame should also contain two columns giving the coordinates of the
#' locations for spatio temporal model fitting.  
#' @param package Which package is to be used in model fitting? Currently available 
#' packages are:
#' \itemize{  
#' \item{"spBayes" }{The model implemented is the dynamic spatio-temporal 
#' model fitted using the \code{\link{spDynLM}} function in the \code{spBayes} package.}  
#' \item{"stan" }{The model implemented is the marginal independent GP model.}
#' \item{"inla" }{The only model implemented is the AR model.}
#' \item{"spTimer" }{All possible models in this package can be fitted.}
#' \item{"sptDyn" }{All possible models in this package can be fitted.}
#' \item{"none". } {In this case case, the argument  \code{model} must be 
#' specified either as "lm" or "seperable". See below}
#' }
#' @param model The model to be fitted. This argument is passed to the fitting package. 
#' In case the package is none, then it can be either "lm" or "seperable". 
#' The "lm" option is for an independent error regression model 
#' while the other option fits a seperable model without any nugget effect.
#' The separable model fitting method cannot handle missing data. All missing data points
#' in the response variable will be replaced by the grand mean of the available observations.  
#' @param coordtype Type of coordinates: utm, lonlat or plain with utm 
#' (supplied in meters) as the default. Distance will be calculated in units of kilometer
#' if this argument is either utm or lonlat. Euclidean distance will be calculated 
#' if this is given as the third type plain.  If  distance in meter is to be calculated 
#' then coordtype should be passed on as plain although the coords are supplied in UTM. 
#' @param coords A vector of size two identifying the two column numbers 
#' of the data frame to take as coordinates. 
#' Or this can be given as a  matrix of number of sites by 2 providing the coordinates of all the
#' data locations. 
#' @param validrows A vector of row numbers of the supplied data frame 
#' which should be used for validation. When the model is "seperable" this argument 
#' must include all the time points for the sites to be validated.  Otherwise, the 
#' user is allowed to select the row numbers of the data frame validation as they wish. 
#' The default NULL value instructs that validation will not be performed. 
#' \itemize{ 
#' \item{lm model:} {If package is "none" and the model is "lm", this argument is 
#' a vector of row indices of the data frame which should be used for validation. }
#' \item{separable model:} {If package is "none" and the model is "seprable" this argument is 
#'  a vector of site indices which should be used for validation. The "seperable" 
#'  model  does not allow some sites to be used for both fitting and validation.
#'  Thus it is not possible to validate at selected time points using the separable model.}
#'  }
#' @param scale.transform Transformation of the response variable. It can take three values: SQRT, LOG or NONE.
#' Default value is "NONE". 
#' @param prior.beta0 A scalar value or a vector providing the prior mean for beta parameters.
#' @param prior.M Prior precision value (or matrix) for beta.  Defaults to a diagonal 
#' matrix with diagonal values 10^(-4).
#' @param prior.sigma2 Shape and scale parameter value for the gamma prior on 1/sigma^2, the precision. 
#' @param prior.tau2 Shape and scale parameter value for the gamma prior on tau^2, the nugget effect. 
#' @param prior.sigma.eta Shape and scale parameter value for the inverse gamma prior 
#' distribution for sigma^2 eta; only used in the spBayes package. 
#' @param phi.s Only used if the model is "seperable". The value of the 
#'  fixed spatial decay parameter for the exponential covariance function.
#' If this is not provided then a value is chosen which corresponds to an effective range
#' which is the maximum distance between the data locations.
#' @param phi.t  Only used if the model is "seperable". 
#' The fixed decay parameter for the exponential covariance function in the temporal domain.
#' If this is not provided then a value is chosen which corresponds to an effective temporal
#' range which is the maximum time of the data set.
#' @param prior.phi Specifies the prior distribution for \eqn{\phi} only when 
#' package is one of Stan, Sptimer or spTDyn.  Distribution options 
#' unform specified by "Unif" and gamma specified by "Gamm" have been implemented in 
#' both Stan and spTimer. Additionally a half-Cauchy prior distribution specified as "Cauchy"
#' has been implemented in Stan. In the case of spTimer the uniform distribution is discrete 
#' while in the case of Stan the uniform distribution is continuous. In the case of 
#' spTimer the option "FIXED" can be used to keep the value fixed. In that case the fixed value 
#' can be given by by a scalar value as the argument \code{prior.phi.param} below or it can be left 
#' unspecified in which case the  fixed value of \eqn{\phi} is chosen as 3/maximum distance between the 
#' data locations.  The  "FIXED" option is not available for the Stan package. 
#' @param prior.phi.param Lower and upper limits of the uniform prior distribution for
#' phi the spatial decay parameter. For the default uniform distribution the values correspond
#' to an effective range that is between 1\% and 100\% of the maximum distance
#' between the data locations. For the Gamma distribution the default values are 2 and 1
#' and for the Cauchy distribution the default values are 0, 1 which specifies
#' a half-Cauchy distribution in \eqn{(0, \infty)}.  
#' @param phi.tuning Only relevant for spTimer and spTDyn models. 
#' Tuning prameter fo sampling phi. See the help file for spT.Gibbs
#' @param phi.npoints  Only relevant for spTimer and spTDyn models. 
#' Number of points for the discrete uniform prior distribution on phi. See the help file for spT.Gibbs
#' @param g_size   Only relevant for GPP models fitted by either spTimer or spTDyn. 
#' The grid size c(m, n) for the knots for the GPP model. A square grid is assumed 
#' if this is passed on as a scalar. This does not need to be given if knots.coords is given instead.
#' @param knots.coords Only relevant for GPP models fitted by either spTimer or spTDyn. 
#' Optional two column matrix of UTM-X and UTM-Y coordinates of the knots in kilometers.
#' It is preferable to specify the g_size parameter instead.
#' @param time.data	Defining the segments of the time-series set up using the function \code{\link{spT.time}}.
#' Only used with the spTimer package. 
#' @param truncation.para	Provides truncation parameter lambda and truncation point "at" using list.
#' Only used with the spTimer package for a truncated model. 
#' @param newcoords	The locations of the prediction sites in similar format to the \code{coords} argument, 
#' only required if fit and predictions are to be performed simultaneously. 
#' If omitted, no predictions will be performed.
#' @param newdata	 The covariate values at the prediction sites specified by \code{newcoords}. 
#' This should have same space-time structure as the original data frame.
#' @param annual.aggrn This provides the options for calculating annual summary 
#' statistics by aggregating different time segments (e.g., annual mean). 
#' Currently implemented options are: "NONE", "ave" and "an4th", 
#' where "ave" = annual average, "an4th"= annual 4th highest. 
#' Only applicable if spT.time inputs more than one segment and 
#' when fit and predict are done simultaneously. Only used in the spTimer package. 
## @param fitted.values Only relevant for the spTimer package, 
##  this option provides calculating fitted values and corresponding 
##  sd in the original scale. Currently implemented options are: 
## "ORIGINAL" and "TRANSFORMED". Only applicable if scale.transform inputs "SQRT" or "LOG". 
## Note that the PMCC (model validation criteria) values are changed accordingly.
#' @param rhotp  Only relevant for models fitted by  spTDyn. 
#' Initial value for the rho parameters in the temporal dynamic model.
#' The default is rhotp=0 for which  parameters are sampled from the full conditional distribution
#' via MCMC with initial value 0.
#' If rhotp=1,parameters are not sampled and fixed at value 1.
#' @param prior.range A length 2 vector, with (range0, Prange) specifying 
#' that \eqn{P(\rho < \rho_0)=p_{\rho}}, where \eqn{\rho} is the spatial range of 
#' the random field. If Prange is NA, then range0 is used as a fixed range value. 
#' If this parameter is unspecified then range0=0.90 * maximum distance 
#' and Prange =0.95. If instead a single value is specified then the range is set at the single value.
#' @param prior.sigma A length 2 vector, with (sigma0, Psigma) specifying 
#' that \eqn{P(\sigma > \sigma_0)=p_{\sigma}}, where \eqn{\sigma} is the marginal 
#' standard deviation of the field. If Psigma is NA, then sigma0 is taken as the fixed value of 
#' this parameter.
#' @param offset Only used in INLA based modelling.  Offset parameter. See documentation for \code{inla.mesh.2d}.
#' @param max.edge Only used in INLA based modelling. See documentation for \code{inla.mesh.2d}.
#' @param cov.model Model for the covariance function. Only relevant for the spBayes, spTimer and the spTDyn packages.  Default is the exponential model. 
#' See the documentation for \code{\link{spLM}} in the package spBayes. 
#' @param tol.dist	Minimum separation distance between any two locations out of those specified by 
#' coords, knots.coords and pred.coords. The default is 0.005. The program 
#' will exit if the minimum distance is less than the non-zero specified value. 
#' This will ensure non-singularity of the covariance matrices.
#' @param N MCMC sample size. 
#' @param burn.in How many initial iterations to discard. 
#' Only relevant for MCMC based model fitting, i.e., when package is spBayes or Stan.  
#' @param rseed Random number seed that controls the starting point for the random numer stream.
#' A set value is required to help reproduce the results.
#' @param n.report  How many times to report in MCMC progress. This is only used when the package is spBayes or spTimer. 
#' @param no.chains Number of parallel chains to run in Stan. 
#' @param ad.delta Adaptive delta controlling the behavior of Stan during fitting. 
#' @param t.depth Maximum allowed tree depth in the fitting process of Stan. 
#' @param s.size step size in the fitting process of Stan. 
#' @param plotit  Logical scalar value: whether to plot the predictions against the observed values.
#' @param verbose Logical scalar value: whether to print various estimates and statistics.
#' @param mchoice Logical scalar value: whether model choice statistics should be calculated.
#' @param ... Any additional arguments that may be passed on to the fitting package. 
#' @return A list containing:
#'   \itemize{
#'    \item params - A table of parameter estimates 
#'    \item fit  -   The fitted model object. This is present only if a 
#'    named package, e.g.  \code{spTimer}  has been used.
#'    \item max.d  -   Maximum distance between data locations. 
#'    This is in unit of kilometers unless the \code{coordtype} argument 
#'    is set as \code{plain}.     
#'    \item fitteds  -   A vector of fitted values.   
#'     \item mchoice -   Calculated model choice statistics if 
#'     those have been requested by the input argument \code{mchoice=T}.
#'      Not all model fits will contain  all the model choice statistics. 
#'     \item stats   -   The four validation statistics: 
#'     rmse, mae, crps and coverage.  This is present only if 
#'     model validation has been performed. 
#'    \item yobs_preds  -    A data frame containing the 
#'    validation rows of the model fitting data  frame. 
#'    The last five columns of this data frame contains 
#'    the validation prediction summaries: mean, sd, 
#'    median, and 95\% prediction interval. 
#'    This is present only if model validation has been performed. 
#'    \item valpreds  -   A matrix containing the MCMC samples of 
#'    the validation predictions. 
#'    The dimension of this matrix is the number of validations times 
#'    the number of retained 
#'    MCMC samples. This is present only if model validation has 
#'    been performed.  
#'    \item residuals   -   A vector of residual values.  
#'    \item sn   -   The number of data locations used in fitting.  
#'    \item tn   -   The number of time points used in fitting. 
#'    \item phi.s, phi.t   -   Adopted value of the spatial and temporal 
#'    decay parameters if those were fixed during model fitting.  
#'    \item{prior.phi}  -   If present this contains the name of 
#'    the prior distribution for 
#'    the spatial decay parameter \eqn{phi} used to fit the 
#'    model. 
#'    \item prior.phi.param  -    If present this contains the values of 
#'    the hyperparameters of the prior distribution for the spatial 
#'    decay parameter \eqn{phi}.  
#'    \item prior.range   -    Present only if the \code{INLA} package 
#'    has been used in model fitting.  This contains the values of 
#'    the hyperparameters of the prior distribution for the range.  
#'    \item logliks   -    A list containing the log-likelihood 
#'    values used in calculation of the model choice statistics 
#'    if those have been requested in the first place. 
#'    \item knots.coords  -   The locations of the knots if the 
#'    model has been fitted using the GPP method.  
#'    \item formula   -   The input formula for the regression part of 
#'    the model.  
#'     \item scale.transform  -   The transformation adopted by the 
#'     input argument with the same name.  
#'    \item package   -   The name of the package used for model fitting.  
#'    \item model   -   The name of the fitted model.   
#'    \item call  -   The command used to call the model fitting function.  
#'    \item computation.time  -   Computation time required 
#'    to run the model fitting.  
#' }
#' @examples
#' library(bmstdr)
#' library(databmstdr)
#' library(ggplot2)
#' library(ggsn)
#' head(nysptime)
#' M1 <- Bsptime(model="lm", data=nysptime,  formula=y8hrmax~xmaxtemp+xwdsp+xrh)
#' names(M1)
#' M2 <- Bsptime(model="separable", data=nysptime, 
#' formula=y8hrmax~xmaxtemp+xwdsp+xrh, coordtype="utm", coords=4:5) 
#' names(M2)
#' \dontrun{
#' f2 <- y8hrmax~xmaxtemp+xwdsp+xrh
#' M1 <- Bsptime(model="lm", formula=f2, data=nysptime, 
#' scale.transform = "SQRT", N=5000)
#' summary(M1)
#' M1$tn
#' M1$sn
#' plot(M1)
#' a <- residuals(M1, numbers=list(sn=28, tn=62))
#' ## Spatio-temporal model fitting
#' M2 <- Bsptime(model="separable", formula=f2, data=nysptime, 
#' scale.transform = "SQRT", coordtype="utm", coords=4:5,  N=5000)
#' names(M2)
#' summary(M2)
#' M2$phi.s
#' M2$phi.t
#' valids <-  c(8,11,12,14,18,21,24,28)
#'vrows <-  which(nysptime$s.index%in% valids)
#'# a <- phichoicep(valids=c(8, 11)) 
#'# optimal phi.s=0.005, phi.t=0.05
#'## Re-run to see
#'M2.1 <- Bsptime(model="separable",  formula=f2, data=nysptime, 
#'validrows=vrows, coordtype="utm", coords=4:5,
#'phi.s=0.005, phi.t=0.05, scale.transform = "SQRT")
#'summary(M2.1)
#'plot(M2.1)
#'plot(M2.1, segments = F)
#'M1 <- Bsptime(model="lm", formula=f2, data=nysptime, 
#'validrows=vrows, scale.transform = "SQRT")
#'summary(M1)
#'M3 <- Bsptime(package="spTimer", formula=f2, data=nysptime, 
#'coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#'n.report = 5, N=5000)
#'summary(M3)
#'plot(M3)
#'nymap <- map_data(database="state",regions="new york")
#'head(nymap)
#'s <- c(1, 5, 10)
#'fcoords <- nyspatial[-s, c("Longitude", "Latitude")]
#'vcoords <- nyspatial[s,  c("Longitude", "Latitude")]
#'library(tidyr)
#'label <- tibble(long = -76.1, lat = 41.5,
#'label = "25 fitted (circles) & 3  \n  validation (numbered) sites")
#'vsites3 <- ggplot() +
#'geom_polygon(data=nymap, aes(x=long, y=lat, group=group),
#'color="black", size = 0.6, fill=NA) +
#'geom_point(data =fcoords, aes(x=Longitude,y=Latitude)) +
#'geom_text(data=vcoords, aes(x=Longitude,y=Latitude, label=s), col=4) +
#'labs( title= "28 air pollution monitoring sites in New York", x="Longitude", y = "Latitude") +
#'geom_text(aes(label=label, x=long, y=lat), data = label, vjust = "top", hjust = "right")  +
#' geom_rect(mapping=aes(xmin=-80.2, xmax=-77.3, ymin=41, ymax=41.6), color="black", fill=NA) + 
#' geom_rect(mapping=aes(xmin=-78.7, xmax=-75.8, ymin=41, ymax=41.6), color="black", fill=NA) + 
#' scalebar(data =nymap, dist = 100, location = "bottomleft", transform=T, dist_unit = "km",
#' st.dist = .05, st.size = 5, height = .06, st.bottom=T, model="WGS84") +
#' north(data=nymap, location="topleft", symbol=12) 
#' vsites3
#' set.seed(44)
#' tn <- 62
#' valids <- c(1, 5, 10)
#' validt <- sort(sample(1:62, size=31))
#' vrows <- getvalidrows(sn=28, tn=62, valids=valids, validt=validt)
#' ymat <- matrix(nysptime$y8hrmax, byrow=T, ncol=tn)
#' yholdout <- ymat[valids, validt]
#' M31 <- Bsptime(package="spTimer",formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, validrows=vrows, model="GP",  
#' scale.transform = "NONE",  n.report = 5, N=5000)
#' summary(M31)
#' modfit <- M31$fit
#' ## Extract the fits for the validation sites
#' fitall <- data.frame(modfit$fitted)
#' head(fitall)
#' fitall$s.index <- rep(1:28, each=tn)
#' library(spTimer)
#' vdat <- spT.subset(data=nysptime, var.name=c("s.index"), s=valids)
#' fitvalid <- spT.subset(data=fitall, var.name=c("s.index"), s=valids)
#' head(fitvalid)
#' fitvalid$low <- fitvalid$Mean - 1.96 * fitvalid$SD
#' fitvalid$up <- fitvalid$Mean + 1.96 * fitvalid$SD
#' fitvalid$yobs <- sqrt(vdat$y8hrmax)
#' fitvalid$yobs <- vdat$y8hrmax
#' yobs <- matrix(fitvalid$yobs, byrow=T, ncol=tn)
#' y.valids.low <- matrix(fitvalid$low, byrow=T, ncol=tn)
#' y.valids.med <- matrix(fitvalid$Mean, byrow=T, ncol=tn)
#' y.valids.up <- matrix(fitvalid$up, byrow=T, ncol=tn)
#' p1 <- fig11.13.plot(yobs[1, ], y.valids.low[1, ], y.valids.med[1, ], y.valids.up[1, ], misst=validt)
#' p1 <- p1 + ggtitle("Validation for Site 1")
#' p1
#' p2 <- fig11.13.plot(yobs[2, ], y.valids.low[2, ], y.valids.med[2, ], y.valids.up[2, ], misst=validt)
#' p2 <- p2 + ggtitle("Validation for Site 5")
#' p3 <- fig11.13.plot(yobs[3, ], y.valids.low[3, ], y.valids.med[3, ], y.valids.up[3, ], misst=validt)
#' p3 <- p3 + ggtitle("Validation for Site 10")
#' library(ggpubr)
#' ggarrange(p1, p2, p3, common.legend = TRUE, legend = "top", nrow = 3, ncol = 1)
#' ## Using stan 
#' a <- Bsptime(package="stan",formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, N=60, burn.in=10, verbose = F)
#' summary(a)
#' # Prediction map and a map of the standard deviation  
#' post <- M3$fit
#' library(spTimer)
#' gpred <- predict(post, newdata=gridnysptime, newcoords=~Longitude+Latitude)
#' meanmat <- post$op
#' sig2eps <-  post$sig2ep
#' dim(meanmat)
#' sige <- sqrt(sig2eps)
#' itmax <- ncol(meanmat)
#' nT <- nrow(nysptime)
#' sigemat <- matrix(rep(sige, each=nT), byrow=F, ncol=itmax)
#' ##
#' a <- matrix(rnorm(nT*itmax), nrow=nT, ncol=itmax)
#' ypreds <- meanmat + a * sigemat
#' if (scale.transform == "SQRT")  
#' ypreds <-  (ypreds)^2
#' dim(ypreds)
#' sitemeans <- function(a, sn=100, tn=62) { 
#' u <- matrix(a, nrow=sn, ncol=tn, byrow=T)
#' b <- apply(u, 1, mean)
#' as.vector(b)
#' }
#' v <- apply(ypreds, 2, sitemeans, sn=28)
#' a <- get_parameter_estimates(t(v)) 
#' fits <- data.frame(nyspatial[, 1:5], a) ## all the data to plot
#' fits
#' names(gpred)
#' gpred$predN
#' length(gpred$Mean)
#' dim(gpred$pred.samples)
#' 
#' u <- gpred$pred.samples
#' v <- apply(u, 2, sitemeans)
#' a <- get_parameter_estimates(t(v)) 
#' b <- data.frame(gridnyspatial[, 1:5], a) ## all the data to plot
#' b <- rbind(b, fits)
#' coord <- nyspatial[, c("Longitude","Latitude")]
#' library(akima)
#' xo <- seq(from=min(coord$Longitude)-0.5, to = max(coord$Longitude)+0.8, length=200)
#' yo <- seq(from=min(coord$Latitude)-0.25, to = max(coord$Latitude)+0.8, length=200)
#' surf <- interp(b$Longitude, b$Latitude, b$mean,  xo=xo, yo=yo)
#' names(surf) 
#' v <- fnc.delete.map.XYZ(xyz=surf)
#' interp1 <- data.frame(long = v$x, v$z )
#' names(interp1)[1:length(v$y)+1] <- v$y
#' library(tidyr)
#' interp1 <- gather(interp1,key = lat,value =Predicted,-long,convert = TRUE)
#' library(ggplot2)
#' nymap <- map_data(database="state",regions="new york")
#' mappath <- cbind(nymap$long, nymap$lat)
#' zr <- range(interp1$Predicted, na.rm=T)
#' com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#colour palette
#' P <- ggplot() +  
#' geom_raster(data=interp1, aes(x = long, y = lat,fill = Predicted)) +
#' geom_polygon(data=nymap, aes(x=long, y=lat, group=group), color="black", size = 0.6, fill=NA) + 
#' geom_point(data=coord, aes(x=Longitude,y=Latitude))  +
#' stat_contour(data=na.omit(interp1), aes(x = long, y = lat,z = Predicted), 
#' colour = "black", binwidth =2) +
#' scale_fill_gradientn(colours=com, na.value="gray95", limits=zr) +
#' theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#' ggsn::scalebar(data =interp1, dist = 100, location = "bottomleft", transform=T, dist_unit = "km",
#' st.dist = .05, st.size = 5, height = .06, st.bottom=T, model="WGS84") +
#' ggsn::north(data=interp1, location="topleft", symbol=12) +
#' labs(title= "Predicted map of average ozone air pollution in New York",
#'  x="Longitude", y = "Latitude", size=2.5) 
#' P
#' # Repeat the above for standard deviation 
#' surf <- interp(b$Longitude, b$Latitude, b$sd,  xo=xo, yo=yo)
#' names(surf)
#' v <- fnc.delete.map.XYZ(xyz=surf)
#' interp1 <- data.frame(long = v$x, v$z )
#' names(interp1)[1:length(v$y)+1] <- v$y
#' interp1 <- gather(interp1,key = lat,value =sd,-long,convert = TRUE)
#' nymap <- map_data(database="state",regions="new york")
#' mappath <- cbind(nymap$long, nymap$lat)
#' zr <- range(interp1$sd, na.rm=T)
#' Psd <- ggplot() +  
#' geom_raster(data=interp1, aes(x = long, y = lat,fill = sd)) +
#' geom_polygon(data=nymap, aes(x=long, y=lat, group=group), color="black", size = 0.6, fill=NA) + 
#' geom_point(data=coord, aes(x=Longitude,y=Latitude))  +
#' stat_contour(data=na.omit(interp1), aes(x = long, y = lat,z = sd), 
#' colour = "black", binwidth =0.1) +  
#' scale_fill_gradientn(colours=com, na.value="gray95", limits=zr) +
#' theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#' ggsn::scalebar(data =interp1, dist = 100, location = "bottomleft", transform=T, dist_unit = "km",
#' st.dist = .05, st.size = 5, height = .06, st.bottom=T, model="WGS84") +
#' ggsn::north(data=interp1, location="topleft", symbol=12) +
#' labs(title= "Sd map of predicted air pollution in New York", 
#' x="Longitude", y = "Latitude", size=2.5) 
#' Psd
#' 
#' 
#' ## Refit M1, M2, M3 with mchoice=T
#' M1.c <- Bsptime(model="lm", formula=f2, data=nysptime, 
#' scale.transform = "SQRT", mchoice=T)
#' M2.c <- Bsptime(model="separable",  formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, 
#' phi.s=0.005, phi.t=0.05, scale.transform = "SQRT", mchoice=T)
#' M3.c <- Bsptime(package="spTimer", model="GP", formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", n.report=5,
#'  mchoice=T, N=5000)
#' M4.c <- Bsptime(package="stan",formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", N=1500, burn.in=500,
#'  mchoice=T, verbose = F)
#'  aresid <- residuals(M4.c)
#'  valids <- c(8,11,12,14,18,21,24,28)
#'  vrows <- getvalidrows(sn=28, tn=62, valids=valids, allt=T)
#'  M1.v <- Bsptime(model="lm", formula=f2, data=nysptime, 
#'  scale.transform = "SQRT", validrows=vrows)
#'  M2.v <- Bsptime(model="separable",  formula=f2, data=nysptime, 
#'  coordtype="utm", coords=4:5, phi.s=0.005, phi.t=0.05, 
#'  scale.transform = "SQRT", validrows=vrows)
#'  M3.v <- Bsptime(package="spTimer", model="GP", formula=f2, data=nysptime, 
#'  coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#'  mchoice=T, validrows=vrows, n.report =5, N=5000)
#'  M4.v <- Bsptime(package="stan",formula=f2, data=nysptime, 
#'  coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#'  N=1500, burn.in=500, mchoice=T, validrows=vrows, verbose = F)
#'  plot(M4.v, segments = F)
#'  tablemchoice <- cbind(M1.c$mchoice, M2.c$mchoice, M3.c$mchoice, M4.c$mchoice)
#'  round(tablemchoice, 2)
#'  tablemvalid <- cbind(M1.v$stats, M2.v$stats, M3.v$stats, M4.v$stats)
#'  tablemvalid
#'  plot(M3.v, segments=F)
#'  names(M3.v)
#'  names(M3.v$yobs_preds)
#'  obs <- M3.v$yobs_preds$y8hrmax
#'  dim(M3.v$yobs_preds)
#'  predsums <- get_validation_summaries(M3.v$valpreds)
#'  head(predsums)
#'  pall <- obs_v_pred_plot(obs, predsums, segments=F)
#'  names(pall)
#'  M5 <- Bsptime(package="spTimer", model="AR", formula=f2, data=nysptime,
#'  coordtype="utm", coords=4:5, scale.transform = "SQRT", mchoice=T,  
#'  n.report=5, N=5000)
#'  summary(M5)
#'  a <- residuals(M5)
#'  M6 <- Bsptime(package="inla", model="AR", formula=f2, data=nysptime,
#'  coordtype="utm", coords=4:5, scale.transform = "SQRT", mchoice=T)
#'  # 5 minutes 
#'  summary(M6)
#'  valids <- c(8,11,12,14,18,21,24,28)
#'  vrows <- getvalidrows(sn=28, tn=62, valids=valids, allt=T)
#'  M5.v <- Bsptime(package="spTimer", model="AR", formula=f2, data=nysptime, 
#'  coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#'  validrows=vrows, mchoice=T, n.report=5, N=5000)
#'  M6.v <- Bsptime(package="inla", model="AR", formula=f2, data=nysptime, 
#'  coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#'  validrows=vrows, mchoice=T)
#'  a1 <- c(M5$mchoice, unlist(M5.v$stats))
#'  b <- as.list(M6$mchoice)
#'  b
#'  b1 <- c(b$gof, b$penalty, b$pmcc,  unlist(M6.v$stats))
#'  b1 <- c(M6$mchoice$gof, M6$mchoice$penalty, M6$mchoice$pmcc,  unlist(M6.v$stats))
#'  table78 <- rbind(a1, b1)
#'  round(table78, 2)
#'  ## sptDyn model fitting 
#'  library(spTDyn)
#'  ## spatio-temporally varying
#'  f3 <- y8hrmax~ xmaxtemp + sp(xmaxtemp)+ tp(xwdsp) + xrh
#'  M7 <- Bsptime(package="sptDyn", model="GP", formula=f3, data=nysptime,
#'  coordtype="utm", coords=4:5, scale.transform = "SQRT", mchoice=T,  
#'  N=5000, n.report=2)
#'  summary(M7)
#'  a <- residuals(M7)
#'  out <- M7$fit
#'  dim(out$betasp)
#' ## Plotting the spatial effects  
#' spatial_effects_plot <-function(ids=1:28){
#' boxplot(t(out$betasp[ids,]),pch=".",main="Maximum temperature",
#' xlab="Sites",ylab="Values")
#' abline(h=0, col=2)
#' }
#' a <- out$betasp
#' u <- c(t(out$betasp))
#' sn <- nrow(a)
#' itmax <- ncol(a)
#' v <- rep(1:sn, each=itmax)
#' d <- data.frame(site=as.factor(v), sp = u)
#' p <- ggplot(data=d, aes(x=site, y=sp)) + 
#' geom_boxplot(outlier.colour="black", outlier.shape=1, outlier.size=0.5) +
#' geom_abline(intercept=0, slope=0, color="blue") + 
#' labs(title= "Spatial effects of maximum temperature", x="Site", y = "Effects", size=2.5) 
#' p 
#' dim(out$betatp)
#' b <- out$betatp
#' dim(b)
#' tn <- nrow(b)
#' itmax <- ncol(b)
#' tids <- 1:tn 
#' stat <- apply(b[tids,], 1, quantile, prob=c(0.025,0.5,0.975))
#' tstat <- data.frame(tids, t(stat))
#' dimnames(tstat)[[2]] <- c("Days", "low", "median", "up")
#' head(tstat)
#' yr <- c(min(c(stat)),max(c(stat)))
#' p <- ggplot(data=tstat, aes(x=Days, y=median)) +
#' geom_point(size=3) + 
#' ylim(yr) + 
#' geom_segment(data=tstat, aes(x=Days, y=median, xend=Days, yend=low), linetype=1) +
#' geom_segment(data=tstat, aes(x=Days, y=median, xend=Days, yend=up), linetype=1) +
#' geom_abline(intercept=0, slope=0, col="blue") +
#' labs(title="Temporal effects of wind speed", x="Days", y="Temporal effects") 
#' p 
#' ## Figure: 5, Section: 4.2 ##
#' temporal_effects_plot <- function(tids=1:nrow(out$betatp)){
#' stat <- apply(out$betatp[tids,],1,quantile,prob=c(0.025,0.5,0.975))
#' plot(stat[2,],type="p",lty=3,col=1,ylim=c(min(c(stat)),max(c(stat))),
#' pch=19,ylab="",xlab="Days",axes=FALSE,main="Wind Speed",cex=0.8)
#' k <- length(tids)
#' for(i in 1:k){
#' segments(i, stat[2,i], i, stat[3,i])
#' segments(i, stat[2,i], i, stat[1,i])
#' }
#' axis(1,tids,labels=tids);axis(2)
#' # abline(v=31.5,lty=2)
#' abline(h=0, col=2)
#' #   text(15,0.32,"July");  text(45,0.32,"August");
#' }
#' # spBayes 
#' set.seed(44)
#' validt <- sort(sample(1:62, size=31))
#' valids <- c(1, 5, 10)
#' vrows <- getvalidrows(sn=28, tn=62, valids=valids, validt=validt)
#' ymat <- matrix(nysptime$y8hrmax, byrow=T, ncol=62)
#' yholdout <- ymat[valids, validt]
#' scale.transform <- "SQRT"
#' library(spTDyn)
#' M7.v <- Bsptime(package="sptDyn", model="GP", formula=f3, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#' mchoice=T,  N=5000, validrows=vrows, n.report=2)
#' M8 <- Bsptime(package="spBayes",  formula=f2, data=nysptime, 
#' prior.sigma2=c(2, 25), prior.tau2 =c(2, 25),
#' prior.sigma.eta =c(2, 0.001),
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#' mchoice=T,  N=5000,  n.report=200)
#' M8.v <- Bsptime(package="spBayes",  formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#' prior.sigma2=c(2, 25), prior.tau2 =c(2, 25),
#' prior.sigma.eta =c(2, 0.001),mchoice=T,  N=5000, validrows=vrows,  
#' n.report=200)
#' summary(M8.v)
#' a1 <- c(M7$mchoice, unlist(M7.v$stats))
#' M8$mchoice <- as.list(M8$mchoice)
#' b1 <- c(M8$mchoice$gof, M8$mchoice$penalty, M8$mchoice$pmcc,  
#' unlist(M8.v$stats))
#' table7.10 <- rbind(a1, b1)
#' round(table7.10, 2)
#' # spBayes parameter plots 
#' modfit <- M8$fit
#' N <- 5000
#' burn.in <- 1000
#' tn <- 62
#' ysamples <- modfit$p.y.samples[,burn.in:N]
#' # if (scale.transform=="SQRT") 
#' ysamples <- (ysamples)^2
#' quant95 <- function(x){
#' quantile(x, prob=c(0.025, 0.5, 0.975))
#' }
#' y.hat <- apply(ysamples, 1, quant95) ## 3 by 1736
#' dim(y.hat)
#' ## Extract the ones for hold out
#' y.hat.med <- matrix(y.hat[1, ], ncol=tn)
#' y.hat.up <- matrix(y.hat[3,], ncol=tn)
#' y.hat.low <- matrix(y.hat[2,], ncol=tn)
#' 
#' y.valids.med <- y.hat.med[valids,]
#' y.valids.low <- y.hat.low[valids,]
#' y.valids.up <- y.hat.up[valids,]
#' dim(y.valids.up)
#' yobs <- ymat[valids,]
#' p1 <- fig11.13.plot(yobs[1, ], y.valids.low[1, ], y.valids.med[1, ], y.valids.up[1, ], misst=validt)
#' p2 <- fig11.13.plot(yobs[2, ], y.valids.low[2, ], y.valids.med[2, ], y.valids.up[2, ], misst=validt)
#' p3 <- fig11.13.plot(yobs[3, ], y.valids.low[3, ], y.valids.med[3, ], y.valids.up[3, ], misst=validt)
#' p1
#' p2
#' p3
#' ## spBayes parameter plots 
#' beta <- apply(modfit$p.beta.samples[burn.in:N,], 2, quant95)
#' theta <- apply(modfit$p.theta.samples[burn.in:N,], 2, quant95)
#' sigma.sq <- theta[,grep("sigma.sq", colnames(theta))]
#' tau.sq <- theta[,grep("tau.sq", colnames(theta))]
#' phi <- theta[,grep("phi", colnames(theta))]
#' cat("\nOnly first 3 theta parameters are plotted\n")
#' par(mfrow=c(3,1))
#' plot(1:tn, sigma.sq[1,], pch=19, cex=0.5, xlab="Time", ylab="sigma.sq", ylim=range(sigma.sq))
#' arrows(1:tn, sigma.sq[1,], 1:tn, sigma.sq[3,], length=0.02, angle=90)
#' arrows(1:tn, sigma.sq[1,], 1:tn, sigma.sq[2,], length=0.02, angle=90)
#' plot(1:tn, tau.sq[1,], pch=19, cex=0.5, xlab="Time", ylab="tau.sq", ylim=range(tau.sq))
#' arrows(1:tn, tau.sq[1,], 1:tn, tau.sq[3,], length=0.02, angle=90)
#' arrows(1:tn, tau.sq[1,], 1:tn, tau.sq[2,], length=0.02, angle=90)
#' plot(1:tn, 3/phi[1,], pch=19, cex=0.5, xlab="Time", ylab="eff. range (km)", ylim=range(3/phi))
#' arrows(1:tn, 3/phi[1,], 1:tn, 3/phi[3,], length=0.02, angle=90)
#' arrows(1:tn, 3/phi[1,], 1:tn, 3/phi[2,], length=0.02, angle=90)
#' par(mfrow=c(1,1))
#' adat <- data.frame(x=1:tn, low=sigma.sq[1, ], med=sigma.sq[2, ], up=sigma.sq[3, ])
#' head(adat)
#' psigma <- ggplot() + 
#' geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
#' geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
#' theme(legend.position ="none") + 
#' labs(y = "sigma2", x = "Days") 
#' psigma
#' adat <- data.frame(x=1:tn, low=tau.sq[1, ], med=tau.sq[2, ], up=tau.sq[3, ])
#' head(adat)
#' ptau <- ggplot() + 
#' geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
#' geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
#' theme(legend.position ="none") + 
#' labs(y = "tau2", x = "Days") 
#' ptau
#' adat <- data.frame(x=1:tn, low=3/phi[3, ], med=3/phi[2, ], up=3/phi[1, ])
#' head(adat)
#' prange <- ggplot() + 
#' geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
#' geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
#' theme(legend.position ="none") + 
#' labs(y = "Range", x = "Days") 
#' prange
#' library(ggpubr)
#' ggarrange(psigma, ptau, prange, common.legend = TRUE, legend = "none", nrow = 3, ncol = 1)
#' vnames <- all.vars(f2)
#' xnames <- vnames[-1]
#' k <- 4
#' cat("\nOnly first 4 beta parameters are plotted\n")
#' beta.0 <- beta[,grep("Intercept", colnames(beta))]
#' par(mfrow=c(k,1))
#' plot(1:tn, beta.0[1,], pch=19, cex=0.5, xlab="Time", ylab="Intercept", ylim=range(beta.0))
#' arrows(1:tn, beta.0[1,], 1:tn, beta.0[3,], length=0.02, angle=90)
#' arrows(1:tn, beta.0[1,], 1:tn, beta.0[2,], length=0.02, angle=90) 
#' for (j in 2:k) {
#' betaj <- beta[,grep(xnames[j-1], colnames(beta))]
#' plot(1:tn, betaj[1,], pch=19, cex=0.5, xlab="Time", ylab=xnames[j-1], ylim=range(betaj))
#' abline(h=0, col=2)
#' arrows(1:tn, betaj[1,], 1:tn, betaj[3,], length=0.02, angle=90)
#' arrows(1:tn, betaj[1,], 1:tn, betaj[2,], length=0.02, angle=90)
#' }
#' adat <- data.frame(x=1:tn, low=beta.0[1,], med=beta.0[2,], up=beta.0[3,]
#' head(adat)
#' pint <- ggplot() + 
#' geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
#' geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
#' geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
#' theme(legend.position ="none") + labs(y = "Intercept", x = "Days") 
#' pint
#' j <- 2
#' betaj <- beta[,grep(xnames[j-1], colnames(beta))]
#' adat <- data.frame(x=1:tn, low=betaj[1,], med=betaj[2,], up=betaj[3,])
#' ptmp <- ggplot() + 
#' geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
#' geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
#' geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
#' theme(legend.position ="none") + 
#' labs(y = "Max temp", x = "Days") 
#' ptmp
#' j <- 3
#' betaj <- beta[,grep(xnames[j-1], colnames(beta))]
#' adat <- data.frame(x=1:tn, low=betaj[1,], med=betaj[2,], up=betaj[3,])
#' head(adat)
#' pwdsp <- ggplot() + 
#' geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
#' geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
#' theme(legend.position ="none") +  
#' geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
#' labs(y = "Wind speed", x = "Days") 
#' pwdsp
#' j <- 4
#' betaj <- beta[,grep(xnames[j-1], colnames(beta))]
#' adat <- data.frame(x=1:tn, low=betaj[1,], med=betaj[2,], up=betaj[3,])
#' head(adat)
#' prh <- ggplot() + 
#' geom_point(data = adat, aes(x =x, y = med, shape=19), shape=19, col="blue", size = 2) + 
#' geom_ribbon(data = adat, aes(x =x, ymin =low, ymax = up), alpha = 0.2, color = "grey50") +
#' theme(legend.position ="none") +  
#' geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)+
#' labs(y = "Rel humidity", x = "Days") 
#' prh
#' library(ggpubr)
#' ggarrange(pint, ptmp, pwdsp, prh, common.legend = TRUE, legend = "none", nrow = 2, ncol = 2)
#' ## GPP
#' M9 <-  Bsptime(package="spTimer", model="GPP", g_size=5, 
#' formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT")
#' table7.11 <- M9$params
#' table7.11
#' aresid <- residuals(M9)
#' valids <- c(8,11,12,14,18,21,24,28)
#' vrows <- getvalidrows(sn=28, tn=62, valids=valids, allt=T)
#' M9.5 <-   Bsptime(package="spTimer", model="GPP", g_size=5, 
#' formula=f2, data=nysptime, validrow=vrows, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", mchoice=T)
#' M9.4 <-   Bsptime(package="spTimer", model="GPP", g_size=4, 
#' formula=f2, data=nysptime, validrow=vrows, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", mchoice=T)
#' M9.3 <-   Bsptime(package="spTimer", model="GPP", g_size=3, 
#' formula=f2, data=nysptime, validrow=vrows, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", mchoice=T)
#' u3 <-  c(3, M9.3$mchoice, unlist(M9.3$stats))
#' u4 <-  c(4, M9.4$mchoice, unlist(M9.4$stats))
#' u5 <-  c(5, M9.5$mchoice, unlist(M9.5$stats))
#' table7.12 <- rbind(u3, u4, u5)
#' round(table7.12, 2)
#' # Validation table 
#' head(nysptime)
#' f2 <- y8hrmax~xmaxtemp+xwdsp+xrh
#' f3 <- y8hrmax~ xmaxtemp + sp(xmaxtemp)+ tp(xwdsp) + xrh
#' valids <- c(8,11,12,14,18,21,24,28)
#' vrows <- getvalidrows(sn=28, tn=62, valids=valids, allt=T)
#' M1.v <- Bsptime(model="lm", formula=f2, data=nysptime, 
#' scale.transform = "SQRT", validrows=vrows, N=5000, mchoice=T)
#' M2.v <- Bsptime(model="separable",  formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, phi.s=0.005, phi.t=0.05, 
#' scale.transform = "SQRT", validrows=vrows, N=5000, mchoice=T)
#' M3.v <- Bsptime(package="spTimer", model="GP", formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#' mchoice=T, validrows=vrows, N=5000,  n.report=2)
#' M4.v <- Bsptime(package="stan",formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#' N=5000, mchoice=T, validrows=vrows, verbose = F)
#' # 8 mins 
#' M5.v <- Bsptime(package="spTimer", model="AR", formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#' validrows=vrows, mchoice=T,  N=5000, n.report = 2)
#' M6.v <- Bsptime(package="inla", model="AR", formula=f2, data=nysptime,
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#' validrows=vrows,  N=5000, mchoice=T)
#' #6 mins 41s
#' library(spTDyn)
#' M7.v <- Bsptime(package="sptDyn", model="GP", formula=f3, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#' mchoice=T,  N=5000, validrows=vrows, n.report=2)
#' M8.v <- Bsptime(package="spBayes",  formula=f2, data=nysptime, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", 
#' prior.sigma2=c(2, 25),prior.tau2 =c(2, 25),
#' prior.sigma.eta =c(2, 0.001),mchoice=T,  N=5000, validrows=vrows,  
#' n.report=200)
#' M9.v <-   Bsptime(package="spTimer", model="GPP", g_size=5, 
#' formula=f2, data=nysptime, validrow=vrows, 
#' coordtype="utm", coords=4:5, scale.transform = "SQRT", mchoice=T, 
#' n.report=20, N=5000)
#' results <- cbind.data.frame(lm=unlist(M1.v$stats), 
#' separable=unlist(M2.v$stats), spTimerGP=unlist(M3.v$stats), 
#' stan=unlist(M4.v$stats),inla=unlist(M6.v$stats),
#' spTimerAR=unlist(M5.v$stats), spTDyn=unlist(M7.v$stats), 
#' spBayes=unlist(M8.v$stats),  sptimerGPP=unlist(M9.v$stats))
#' mcres <- cbind.data.frame(lm=unlist(M1.v$mchoice)[9:11], 
#' separable=unlist(M2.v$mchoice)[9:11], 
#' spTimerGP=unlist(M3.v$mchoice)[9:11], 
#' stan=unlist(M4.v$mchoice)[9:11],
#' inla=unlist(M6.v$mchoice)[5:7],
#' spTimerAR=unlist(M5.v$mchoice), spTDyn=unlist(M7.v$mchoice), 
#' spBayes=unlist(M8.v$mchoice),  sptimerGPP=unlist(M9.v$mchoice))                   
#' results
#' round(results, 2)
#' mcres
#' allres <- rbind.data.frame(results, mcres)
#' round(allres, 2)
#' allres <- allres[, -8]
#' round(allres, 2)
#' }
#' @export
Bsptime <- function(formula, # =y8hrmax~xmaxtemp+xwdsp+xrh, 
                   data, #=nysptime,
                   package="none", 
                   model="GP",  
                   coordtype=NULL, #="utm",  
                   coords=NULL, #=4:5, 
                   validrows=NULL, 
                   scale.transform="NONE", 
                   prior.beta0=0, prior.M=0.0001, prior.sigma2=c(2, 1),
                   prior.tau2 = c(2, 0.1),  prior.sigma.eta =c(2, 0.001),
                   phi.s = NULL,  phi.t =NULL, prior.phi = "Gamm", 
                   prior.phi.param = NULL, phi.tuning =NULL, phi.npoints=NULL,
                   prior.range= c(1, 0.5), prior.sigma = c(1, 0.005), 
                   offset = c(10, 140), max.edge=c(50, 1000),  
                   rhotp = 0,  time.data = NULL, truncation.para = list(at = 0,lambda = 2), 
                   newcoords = NULL, newdata = NULL, annual.aggrn = "NONE",
                   cov.model = "exponential",  g_size = NULL, knots.coords = NULL,
                   tol.dist = 0.005, 
                   N=2000, burn.in=1000, rseed =44, n.report = N/2, 
                   no.chains =1, ad.delta = 0.80, t.depth=15, s.size=0.01, 
                   verbose=FALSE, plotit=T, mchoice=FALSE, ...){
  
 set.seed(rseed)
 start.time<-proc.time()[3]
 data <- as.data.frame(data)
# Error checking 
  if (!is.data.frame(data)) stop("Need a data frame in the data argument")
  
  if (model !="lm") { # Will fit spatio-temporal models 
    s1 <- length(coords)
    if (s1 !=2) stop("Need 2-dimensional spatial coordinates in the coords argument\n")
    coords <-  as.matrix(unique(data[, coords])) 
    sn <- nrow(coords)
    tn <- nrow(data)/sn
    a <- abs(tn - floor(tn))
    if (a>0) stop("Unequal number of time points: check numbers of locations and times")
  } else {
    sn <- nrow(data)
    tn <- 0
  }
  

 implemented <- c("inla", "spBayes", "stan", "spTimer", "spTDyn")
 a <- grepl(package, x=implemented, ignore.case = TRUE)
 if (any(a)) { 
   package <- implemented[which(a)]
   } else { 
  imodel <- c("lm", "separable")
   b <-  grepl(model, x=imodel, ignore.case = TRUE)
   if (any(b)) { 
     model <- imodel[which(b)]
     package <- "none"
   } else { stop("Wrong package or model. Please see helpfile")}
   }
 

  if (package=="none") { 
     if (model=="lm") {
          results <- Blm_sptime(formula=formula, data=data,
                      validrows=validrows, scale.transform=scale.transform,
                      prior.beta0=prior.beta0, prior.M=prior.M, prior.sigma2=prior.sigma2,
                      N=N, plotit=plotit, rseed =rseed,
                      verbose=verbose, mchoice=mchoice)
     } else if (model=="separable") { 
      results <-  Bsp_sptime(formula=formula, data=data, validrows=validrows,
                        coordtype=coordtype, coords=coords, 
                        phi.s=phi.s, phi.t = phi.t,   
                        scale.transform=scale.transform,
                        prior.beta0=prior.beta0, prior.M=prior.M, 
                        prior.sigma2=prior.sigma2,
                        mchoice=mchoice, N=N, verbose =verbose, rseed =rseed,
                        plotit=plotit)
     } else { stop("Model not implemented.") }
  } else  if (package=="spBayes") { 
     results <- BspBayes_sptime(formula=formula, data=data,
                        validrows=validrows,  scale.transform =scale.transform,
                            coordtype=coordtype, coords=coords,
                            prior.beta0=prior.beta0, prior.M=prior.M,
                            prior.sigma2 = prior.sigma2,  
                            prior.tau2 = prior.tau2,    prior.sigma.eta =prior.sigma.eta,
                            prior.phi.param=prior.phi.param,
                            phi.tuning = phi.tuning,   
                            n.report = n.report,
                            verbose = verbose,
                            mchoice=mchoice,
                            N=N, burn.in=burn.in, rseed =rseed, plotit=plotit)
    } else if (package=="stan") { 
      results <- Bstan_sptime(formula=formula, data=data,
                    validrows=validrows,  scale.transform =scale.transform,
                           coordtype=coordtype, 
                           coords=coords,
                           # prior.beta0=0, prior.M=0.0001,
                           prior.sigma2=prior.sigma2,
                           prior.tau2=prior.tau2,
                           prior.phi = prior.phi, 
                           prior.phi.param = prior.phi.param, 
                           verbose = verbose,
                           mchoice=mchoice,
                           no.chains =no.chains,
                           rseed =rseed,  ad.delta =ad.delta, s.size=s.size, t.depth=t.depth,
                           N=N, burn.in=burn.in, plotit=plotit)
    } else if (package=="inla") { 
      results <- Binla_sptime(formula=formula, data=data,
                           validrows=validrows,  
                           scale.transform =scale.transform,
                           coordtype=coordtype, coords=coords,
                           verbose =verbose ,
                           mchoice=mchoice,
                           prior.tau2=prior.tau2,
                           prior.range= prior.range,
                           prior.sigma = prior.sigma,
                           offset=offset, 
                           max.edge = max.edge, 
                           N=N, burn.in=burn.in, rseed=rseed,  plotit=plotit)
    } else if (package=="spTimer") { 
      results <- BspTimer_sptime(data=data, formula=formula, model=model,
               coordtype=coordtype, coords=coords,
               validrows=validrows,   scale.transform =scale.transform,
               prior.beta0=prior.beta0, prior.M=prior.M, prior.sigma2 =prior.sigma2,
               prior.phi=prior.phi,
               prior.phi.param =prior.phi.param, 
               phi.tuning=phi.tuning, phi.npoints=phi.npoints,
               cov.fnc=cov.model, tol.dist = tol.dist, 
               time.data = time.data, newcoords = newcoords, newdata =newdata,
               truncation.para = truncation.para, annual.aggrn = annual.aggrn,
               # fitted.values=fitted.values, 
               N=N, burn.in=burn.in, plotit=plotit,
               mchoice=mchoice, verbose=verbose, rseed=rseed,
               g_size =g_size, knots.coords = knots.coords, n.report=n.report) 
  } else if (package=="spTDyn") { 
     results <- BspTDyn_sptime(data=data, formula=formula, model=model,
                                coordtype=coordtype, coords=coords, time.data=time.data, 
                                validrows=validrows,   scale.transform =scale.transform,
                                prior.beta0=prior.beta0, prior.M=prior.M, prior.sigma2 =prior.sigma2,
                                prior.phi=prior.phi,
                                prior.phi.param =prior.phi.param, 
                                phi.tuning=phi.tuning, phi.npoints=phi.npoints,
                                rhotp = rhotp, cov.fnc=cov.model, truncation.para=truncation.para, 
                                N=N, burn.in=burn.in, plotit=plotit,
                                mchoice=mchoice, verbose=verbose, rseed=rseed, n.report=n.report)
  } else { 
    cat("Implemented packages are none,", implemented, "\n")
    cat("\nIf package is none then the implemented models are lm and separable\n")
    stop("But, sorry, the package or model opted for  has not been implemented yet!")
     
  } 
 
 u <- getXy(formula=formula, data=data)
 y <- u$y  
 if (scale.transform == "SQRT") y <- sqrt(y)
 if (scale.transform == "LOG") y <- log(y)
 
 if  (model=="truncatedGP") { ## Transform the y's 
  truncpara <- results$fit$truncation.para
  at <- truncpara$at 
  lambda <- truncpara$lambda
  y <- truncated.fnc(y, at=at[1], lambda=lambda)$zm
 }  
 
#  k <- length(results$fitteds)
#  if (k<nrow(data)) {
#    allX <- u$X
#    p <- ncol(allX)
#    betastar <- results$params[1:p, 1]
#    fits <- as.vector(allX %*% betastar)
#    results$fitteds <- fits 
# } # else fitteds are there already as in spBayes 
 
 results$residuals <- y - results$fitteds
 results$sn <- sn 
 results$tn <- tn 
 results$formula <- formula 
 results$scale.transform <- scale.transform
 results$package <- package
 results$model <- model 
 results$call <- match.call()
 colnames(results$params) <- c("mean", "sd", "2.5%", "97.5%")
 
 end.time <- proc.time()[3]
 comp.time<-end.time-start.time
 comp.time<-fancy.time(comp.time)
 results$computation.time <- comp.time
 print(comp.time)
 class(results) <- "bmstdr"
 results 
  
}