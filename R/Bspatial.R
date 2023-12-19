#' Bayesian regression model fitting for point referenced spatial data. 
#' Calculates parameter estimates, validation statistics, and 
#' estimated values of several Bayesian model choice criteria. 
#' @param formula An object of class "formula" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param data The data frame for which the model formula is to be fitted. 
#' If a spatial model is to be fitted then the data frame should contain 
#' two columns containing the locations of the coordinates. See the \code{coords} argument below. 
#' @param package Which package is to be used in model fitting? Currently implemented 
#' packages are: "spBayes", "stan",  "inla" and "none". In case, package is chosen as "none"
#' then the argument  \code{model} must be specified either as "lm" or "spat". See below. 
#' @param package Which package is to be used in model fitting? Currently available 
#' packages are:
#' \itemize{  
#' \item "spBayes": The model implemented is the marginal model with 
#' nugget effect using the \code{spLM} function.  
#' \item "stan": The model implemented is the full spatial random effect model 
#' with nugget effect where the decay parameter has been assumed to be fixed. 
#' \item "inla": The model fitted is the spatial random effect model with the 
#' nugget effect.
#' \item "none":  In this case case, the argument  \code{model} must be 
#' specified either as "lm" or "spat". See below.
#' }
#'  Further details and more examples are provided in Chapter 6 of the book 
#' \insertCite{Sahubook;textual}{bmstdr}.
#' @param model Only used when the package has been chosen to be "none". 
#' It can take one of two values: either "lm" or "spat". The "lm" option is for an independent error 
#' regression model while the "spat" option fits a  spatial model without any nugget effect.  
#' @param coordtype Type of coordinates: utm, lonlat or plain with utm 
#' (supplied in meters) as the default. Distance will be calculated in units of kilometer
#' if this argument is either utm or lonlat. Euclidean distance will be calculated 
#' if this is given as the third type plain.  If  distance in meter is to be calculated 
#' then coordtype should be passed on as plain although the coords are supplied in UTM. 
#' @param coords A vector of size two identifying the two column numbers 
#' of the data frame to take as coordinates. 
#' Or this can be given as a  matrix of number of sites by 2 providing the coordinates of all the
#' data locations. 
#' @param validrows A vector of site indices which should be used for validation. 
#' This function does not allow some sites to be used for both fitting and validation.
#' The remaining observations will be used for model fitting. The default NULL value instructs that
#' validation will not be performed.
#' @param scale.transform Transformation of the response variable. It can take three values: SQRT, LOG or NONE.
#' @param prior.beta0 A scalar value or a vector providing the prior mean for beta parameters.
#' @param prior.M Prior precision value (or matrix) for beta.  Defaults to a diagonal 
#' matrix with diagonal values 10^(-4).
#' @param prior.sigma2 Shape and scale parameter value for the gamma prior on 1/sigma^2, the precision. 
#' @param prior.tau2 Shape and scale parameter value for the gamma prior on tau^2, the nugget effect. 
#' @param phi The spatial decay parameter for the exponential covariance function. Only 
#' used if the package is  Stan or the model is a spatial model "spat" without nugget effect when the 
#' \code{package} is "none". 
#' @param prior.phi.param Lower and upper limits of the uniform prior distribution for
#' \eqn{\phi},  the spatial decay parameter when the package is \code{spBayes}. 
#' If this is not specified the default values are chosen so that the effective range is
#' uniformly distributed between 25\% and 100\% of the maximum distance between data locations.
#' @param prior.range A length 2 vector, with (range0, Prange) specifying 
#' that \eqn{P(\rho < \rho_0)=p_{\rho}}, where \eqn{\rho} is the spatial range of 
#' the random field. If Prange is NA, then range0 is used as a fixed range value. 
#' If this parameter is unspecified then range0=0.90 * maximum distance 
#' and Prange =0.95. If instead a single value is specified then the range is set at the single value.
#' @param prior.sigma A length 2 vector, with (sigma0, Psigma) specifying 
#' that \eqn{P(\sigma > \sigma_0)=p_{\sigma}}, where \eqn{\sigma} is the marginal 
#' standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.
#' @param offset Only used in INLA based modeling.  Offset parameter. See documentation for \code{inla.mesh.2d}.
#' @param max.edge Only used in INLA based modeling. See documentation for \code{inla.mesh.2d}.
#' @param cov.model Only relevant for the spBayes package.  Default is the exponential model. 
#' See the documentation for \code{\link{spLM}} in the package spBayes. 
#' @param N MCMC sample size. Default value 5000. 
#' @param burn.in How many initial iterations to discard. Default value 1000. 
#' Only relevant for MCMC based model fitting, i.e., when package is spBayes or Stan.  
#' @param rseed Random number seed that controls the starting point for the random number stream.
#' A set value is required to help reproduce the results.
#' @param n.report  How many times to report in MCMC progress. This is used only when the package is spBayes. 
#' @param no.chains Number of parallel chains to run in Stan. 
#' @param ad.delta Adaptive delta controlling the behavior of Stan during fitting. 
#' @param t.depth Maximum allowed tree depth in the fitting process of Stan. 
#' @param s.size step size in the fitting process of Stan. 
#' @param plotit  Logical scalar value: whether to plot the predictions against the observed values.
#' @param verbose Logical scalar value: whether to print various estimates and statistics.
#' @param mchoice Logical scalar value: whether model choice statistics should be calculated.
#' @param ... Any additional arguments that may be passed on to the fitting package. 
#' @return A list containing:
#'  \itemize{
#'    \item params -  A table of parameter estimates 
#'    \item  fit -  The fitted model object. This is present only if a named 
#'    package, e.g.   \code{spBayes}  has been used. 
#'    \item  max.d -  Maximum distance between data locations. 
#'    This is in unit of kilometers unless the \code{coordtype} argument 
#'    is set as \code{plain}.     
#'    \item  fitteds -  A vector of fitted values.   
#'     \item  mchoice -  Calculated model choice statistics if those have been 
#'     requested by the input argument \code{mchoice=TRUE}. Not all model fits will contain 
#'     all the model choice statistics. 
#'     \item  stats -  The four validation statistics: rmse, mae, crps and coverage. 
#'      This is present only if model validation has been performed. 
#'    \item  yobs_preds -  A data frame containing the validation rows of the model 
#'    fitting data  frame. The last five columns of this data frame contains 
#'    the validation prediction summaries: mean, sd, median, and 95\% prediction interval. 
#'    This is present only if model validation has been performed. 
#'    \item  valpreds -  A matrix containing the MCMC samples of the validation predictions. 
#'    The dimension of this matrix is the number of validations times the number of retained 
#'    MCMC samples. This is present only if model validation has been performed.  
#'    \item validationplots - Present only if validation has been performed. 
#'    Contains three validation plots with or without segment and 
#'    an ordinary plot.  See \code{\link{obs_v_pred_plot}} for more. 
#'    \item  residuals -  A vector of residual values.  
#'    \item  sn -  The number of data locations used in fitting.
#'    \item  tn  Defaults to 1. Used for plotting purposes. 
#'    \item  phi -  If present this contains the fixed value of 
#'    the spatial decay parameter \eqn{phi} used to fit the model. 
#'    \item  prior.phi.param -   If present this contains the values of the hyperparameters 
#'    of the prior distribution for the spatial decay parameter \eqn{phi}.  
#'    \item  prior.range -   Present only if the \code{INLA} package has been used 
#'    in model fitting.  This contains the values of the hyperparameters 
#'    of the prior distribution for the range.  
#'    \item  logliks -   A list containing the log-likelihood values used in calculation 
#'    of the model choice statistics if those have been requested in the first place. 
#'    \item  formula -  The input formula for the regression part of the model.  
#'     \item  scale.transform -  The transformation adopted by the input argument with the 
#'    same name.  
#'    \item  package -  The name of the package used for model fitting.  
#'    \item  model -  The name of the fitted model.   
#'    \item  call -  The command used to call the model fitting function.  
#'    \item  computation.time -  Computation time required to run the model fitting.  
#' }
#' @seealso \code{\link{Bsptime}} for Bayesian spatio-temporal model fitting.
#' @example inst/examples/bspat_examples.R
#' @export
Bspatial <- function(formula, # =yo3~xmaxtemp+xwdsp+xrh, 
                   data, # =nyspatial,
                   package="none", 
                   model="lm",  
                   coordtype=NULL, # "utm",  
                   coords=NULL, # 4:5, 
                   validrows=NULL, 
                   scale.transform="NONE",
                   prior.beta0=0, prior.M=0.0001, 
                   prior.sigma2=c(2, 1),
                   prior.tau2 = c(2, 0.1),  phi = NULL,
                   prior.phi.param = NULL, 
                   prior.range= c(1, 0.5),
                   prior.sigma = c(1, 0.005),
                   offset = c(10, 140), max.edge=c(50, 1000),  
                   cov.model = "exponential",  
                   N=5000, burn.in=1000, rseed =44, n.report = 500, 
                   no.chains =1, ad.delta = 0.99, s.size=0.01,  t.depth=15, 
                   verbose=TRUE, plotit=TRUE, mchoice=FALSE, ...){
  
 start.time<-proc.time()[3]
 set.seed(rseed)
 data <- as.data.frame(data)
 
 # Error checking 
 if (!is.data.frame(data)) stop("Need a data frame in the data argument")
 if (!inherits(formula, "formula")) stop("Need a valid formula")
 if (!is.null(coordtype)) coordtype <- match.arg(coordtype, 
                                    choices=c("utm", "lonlat", "plain"))
 if (!is.null(scale.transform)) scale.transform <- match.arg(scale.transform, 
                                                             choices=c("NONE", "SQRT", "LOG"))
 
 
  implemented <- c("inla", "spBayes", "stan")
  a <- grepl(package, x=implemented, ignore.case = TRUE)
  if (any(a)) { 
   package <- implemented[which(a)]
   } else { 
  imodel <- c("lm", "spat")
   b <-  grepl(model, x=imodel, ignore.case = TRUE)
  if (any(b)) { 
     model <- imodel[which(b)]
      package <- "none"
    } else { stop("Wrong package or model. Please see helpfile")}
   }
  # package <- c("none", "inla", "spBayes", "stan")
  # package <-  match.arg(package, several.ok = TRUE)
  # message("package=", package, "model=", model)
  
  if (model !="lm") { # Will fit spatio-temporal models 
     s1 <- length(coords)
     if (s1 !=2) stop("Need 2-dimensional spatial coordinates in the coords argument\n")
  }
 
 if (package=="inla") {
   if (inlabru::bru_safe_inla()) {
     if (verbose) message("INLA will be used.")
   } else stop("The chosen package INLA is not available.")
}

  if (package=="none") { 
     if (model=="lm") {
          results <- Blm_sp(formula=formula, data=data,
                      validrows=validrows, scale.transform=scale.transform,
                      prior.beta0=prior.beta0, prior.M=prior.M, prior.sigma2=prior.sigma2,
                      N=N, plotit=plotit, rseed =rseed,
                      verbose=verbose, mchoice=mchoice)
     } else if (model=="spat") { 
      results <-  Bsp_sp(formula=formula, data=data,validrows=validrows,
                        coordtype=coordtype, coords=coords, phi=phi, 
                        scale.transform=scale.transform,
                        prior.beta0=prior.beta0, prior.M=prior.M, 
                        prior.sigma2=prior.sigma2,
                        mchoice=mchoice, N=N, verbose =verbose, rseed =rseed,
                        plotit=plotit)
        
     } else { stop("Model not implemented.") }
  } else  if (package=="spBayes") { 
     results <- BspBayes_sp(formula=formula, data=data,
                        validrows=validrows, scale.transform =scale.transform,
                            coordtype=coordtype, coords=coords,
                            prior.beta0=prior.beta0, prior.M=prior.M,
                            prior.sigma2 = prior.sigma2,
                            prior.tau2 = prior.tau2,
                            prior.phi.param=prior.phi.param, cov.model = cov.model,
                            n.report = n.report,
                            verbose = verbose,
                            mchoice=mchoice,
                            N=N, burn.in=burn.in, rseed =rseed, plotit=plotit)
    
    } else if (package=="stan") { 
      results <- Bstan_sp(formula=formula, data=data,
                    validrows=validrows, scale.transform =scale.transform,
                           coordtype=coordtype, 
                           coords=coords,
                           # prior.beta0=0, prior.M=0.0001,
                           phi=phi, prior.sigma2=prior.sigma2,
                           prior.tau2=prior.tau2,
                           verbose = verbose,
                           mchoice=mchoice,
                           no.chains =no.chains,
                           rseed =rseed,  ad.delta =ad.delta, s.size=s.size,  t.depth=t.depth,
                           N=N, burn.in=burn.in, plotit=plotit)
    } else if (package=="inla") { 
      results <- Binla_sp (formula=formula, data=data,
                           validrows=validrows,
                           scale.transform =scale.transform,
                           coordtype=coordtype, coords=coords,
                           verbose =verbose ,
                           mchoice=mchoice,
                           prior.tau2=prior.tau2,
                           prior.range= prior.range,
                           prior.sigma = prior.sigma,
                           offset = offset, max.edge=max.edge,  
                           N=N, burn.in=burn.in, rseed=rseed,  plotit=plotit)
    } else { 
      message("Implemented packages are none,", implemented, "\n")
      message("\n If package is none then the implemented models are lm and spat\n")
      stop("But, sorry, the package or model opted for has not been implemented yet!")
    } 
 
 u <- getXy(formula=formula, data=data)
 y <- u$y  
 if (scale.transform == "SQRT") y <- sqrt(y)
 if (scale.transform == "LOG") y <- log(y)
 
#  k <- length(results$fitteds)
#  if (k<nrow(data)) {
#   allX <- u$X
#   p <- ncol(allX)
#   betastar <- results$params[1:p, 1]
#   fits <- as.vector(allX %*% betastar)
#   results$fitteds <- fits 
# } 
 
 results$residuals <- y - results$fitteds
 results$sn <- nrow(data)
 results$tn <- 1 
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
