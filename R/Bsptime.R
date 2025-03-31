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
#' organized in the order: (s1, t1), (s1, t2), ... (s1, T), ... (sn, t1), ... (sn, T). 
#' The data frame should also contain two columns giving the coordinates of the
#' locations for spatio temporal model fitting.  
#' @param package Which package is to be used in model fitting? Currently available 
#' packages are:
#' \itemize{  
#' \item "spBayes": The model implemented is the dynamic spatio-temporal 
#' model fitted using the spDynLM function in the \code{spBayes} 
#' package.  
#' \item "stan": The model implemented is the marginal independent GP model.
#' \item "inla" The only model implemented is the AR model.
#' \item "spTimer": All possible models in this package can be fitted.
#' \item "sptDyn": All possible models in this package can be fitted.
#' \item "none".: In this case case, the argument  \code{model} must be 
#' specified either as "lm" or "separable". See below.
#' }
#' Further details and more examples are provided in Chapters 7-9 of the book 
#' \insertCite{Sahubook;textual}{bmstdr}.
#' @param model The model to be fitted. This argument is passed to the fitting package. 
#' In case the package is none, then it can be either "lm" or "separable". 
#' The "lm" option is for an independent error regression model 
#' while the other option fits a separable model without any nugget effect.
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
#' which should be used for validation. When the model is "separable" this argument 
#' must include all the time points for the sites to be validated.  Otherwise, the 
#' user is allowed to select the row numbers of the data frame validation as they wish. 
#' The default NULL value instructs that validation will not be performed. 
#' \itemize{ 
#' \item lm model: If package is "none" and the model is "lm", this argument is 
#' a vector of row indices of the data frame which should be used for validation. 
#' \item separable model: If package is "none" and the model is "separable" this argument is 
#'  a vector of site indices which should be used for validation. The "separable" 
#'  model  does not allow some sites to be used for both fitting and validation.
#'  Thus it is not possible to validate at selected time points using the separable model.
#'  Further details are provided in Chapter 7 of 
#'  \insertCite{Sahubook;textual}{bmstdr}.
#' }
#' @param scale.transform Transformation of the response variable. It can take three values: SQRT, LOG or NONE.
#' Default value is "NONE". 
#' @param prior.beta0 A scalar value or a vector providing the prior mean for beta parameters.
#' @param prior.M Prior precision value (or matrix) for beta.  Defaults to a diagonal 
#' matrix with diagonal values 10^(-4).
#' @param prior.sigma2 Shape and scale parameter value for the gamma prior on 1/sigma^2, the precision. 
#' @param prior.tau2 Shape and scale parameter value for the gamma prior on tau^2, the nugget effect. 
#' @param prior.sigma.eta Shape and scale parameter value for the inverse gamma prior 
#' distribution for sigma^2 eta; only used in the spBayes package. 
#' @param phi.s Only used if the model is "separable". The value of the 
#'  fixed spatial decay parameter for the exponential covariance function.
#' If this is not provided then a value is chosen which corresponds to an effective range
#' which is the maximum distance between the data locations.
#' @param phi.t  Only used if the model is "separable". 
#' The fixed decay parameter for the exponential covariance function in the temporal domain.
#' If this is not provided then a value is chosen which corresponds to an effective temporal
#' range which is the maximum time of the data set.
#' @param prior.phi Specifies the prior distribution for \eqn{\phi} only when 
#' package is one of Stan, spTimer or spTDyn.  Distribution options 
#' uniform specified by "Unif" and gamma specified by "Gamm" have been implemented in 
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
#' Tuning parameter fo sampling phi. See the help file for spT.Gibbs
#' @param phi.npoints  Only relevant for spTimer and spTDyn models. 
#' Number of points for the discrete uniform prior distribution on phi. See the help file for spT.Gibbs
#' @param g_size   Only relevant for GPP models fitted by either spTimer or spTDyn. 
#' The grid size c(m, n) for the knots for the GPP model. A square grid is assumed 
#' if this is passed on as a scalar. This does not need to be given if knots.coords is given instead.
#' @param knots.coords Only relevant for GPP models fitted by either spTimer or spTDyn. 
#' Optional two column matrix of UTM-X and UTM-Y coordinates of the knots in kilometers.
#' It is preferable to specify the g_size parameter instead.
#' @param time.data	Defining the segments of the time-series set up using the function 
#' spT.time in the spTimer package. Only used with the spTimer package. 
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
#' @param offset Only used in INLA based modeling.  Offset parameter. See documentation for \code{inla.mesh.2d}.
#' @param max.edge Only used in INLA based modeling. See documentation for \code{inla.mesh.2d}.
#' @param cov.model Model for the covariance function. Only relevant for the spBayes, spTimer and the spTDyn packages.  Default is the exponential model. 
#' See the documentation for spLM in the package spBayes. 
#' @param tol.dist	Minimum separation distance between any two locations out of those specified by 
#' coords, knots.coords and pred.coords. The default is 0.005. The program 
#' will exit if the minimum distance is less than the non-zero specified value. 
#' This will ensure non-singularity of the covariance matrices.
#' @param N MCMC sample size. 
#' @param burn.in How many initial iterations to discard. 
#' Only relevant for MCMC based model fitting, i.e., when package is spBayes or Stan.  
#' @param rseed Random number seed that controls the starting point for the random number stream.
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
#'     those have been requested by the input argument \code{mchoice=TRUE}.
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
#'    \item validationplots - Present only if validation has been performed. 
#'    Contains three validation plots with or without segment and 
#'    an ordinary plot.  See \code{\link{obs_v_pred_plot}} for more. 
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
#' @references
#' \insertAllCited{}
#' @example inst/examples/bsptime_examples.R
#' @export
Bsptime <- function(formula, 
                   data,
                   package="none", 
                   model="GP",  
                   coordtype=NULL,
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
                   N=2000, burn.in=1000, rseed =44, n.report = 2, 
                   no.chains =1, ad.delta = 0.80, t.depth=15, s.size=0.01, 
                   verbose=FALSE, plotit=TRUE, mchoice=FALSE, ...){
  
 set.seed(rseed)
 start.time<-proc.time()[3]
 data <- as.data.frame(data)
# Error checking 
 if (!is.data.frame(data)) stop("Need a data frame in the data argument")
 if (!inherits(formula, "formula")) stop("Need a valid formula")
 if (!is.null(coordtype)) coordtype <- match.arg(coordtype, 
                                                 choices=c("utm", "lonlat", "plain"))
 if (!is.null(scale.transform)) scale.transform <- match.arg(scale.transform, 
                                                             choices=c("NONE", "SQRT", "LOG"))
  
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
 
 # implemented <- c("none", "inla", "spBayes", "stan", "spTimer", "spTDyn")
 # package <-  match.arg(arg = package, choices = implemented)
 
 if (package=="inla") {
   if (inlabru::bru_safe_inla()) {
     if (verbose) message("INLA will be used.")
   } else stop("The chosen package INLA is not available.")
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
     } else { stop("Model not implemented for package=none.") }
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
    message("Implemented packages are: none,", implemented, "\n")
    message("If package is none then the implemented models are lm and separable\n")
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
 message(comp.time)
 class(results) <- "bmstdr"
 results 
  
}