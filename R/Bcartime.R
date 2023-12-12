#' Bayesian regression model fitting for areal and areal spatio-temporal data. 
#' Calculates parameter estimates, validation statistics, and 
#' estimated values of several Bayesian model choice criteria. 
#' @param formula An object of class "formula" (or one that can be coerced to that class):
#' a symbolic description of the regression model to be fitted.
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
#' \item "inla": INLA model fitting for areal data.
#'  See \insertCite{gomezrubio;textual}{bmstdr}.
#' \item "CARBayes": All possible models in this package can be fitted. 
#' See \insertCite{LeeCARBayes2021;textual}{bmstdr}.
#' \item "CARBayesST": All possible models in this package can be fitted.
#' See \insertCite{CarBayesST;textual}{bmstdr}.
#' Further details and more examples are provided in Chapters 10 and 11 of the book 
#' \insertCite{Sahubook;textual}{bmstdr}.
#' }
#' @param link The link function to use for INLA based model fitting. This is 
#' ignored for the CARBayes and   CARBayesST models. 
#' @param model The specific spatio temporal model to be fitted. 
#' If the package is "INLA" then the model argument should be a vector with two elements 
#' giving the spatial model as the first component. The alternatives for the 
#' spatial model are: "bym", "bym2", "besag", "besag2", "besagproper", "besagproper2", "iid" 
#' and "none". The temporal model as the second second component can  be one of 
#'  "iid", "rw1",  "rw2", ar1" or "none". 
#'  In case the model component is "none" then no spatial/temporal random effects 
#'  will be fitted. No temporal random effects will be fitted in case \code{model} is
#'  supplied as a  scalar, i.e. not a vector of two values. 
#' @param family	One of either "gaussian", "binomial","poisson" or "zip", 
#' which respectively specify a Gaussian, binomial likelihood model with the 
#' logistic link function, a Poisson likelihood model with a log link function, 
#' or a zero-inflated Poisson model with a log link function.
#' @param trials Only used if family="binomial". 
#' Either the name (character) or number of the column in 
#' the supplied data frame containing the total number of trials 
#' The program will try to access data[, trials] to get 
#' the binomial denominators. 
#' @param residtype Residual type, either "response" or "pearson",  
#' in GLM fitting with the packages CARBayes and CARBayesST.
#' Default is "response" type observed minus fitted. The other option "pearson" is for 
#' Pearson residuals in GLM. For INLA based model fitting only the default response 
#' residuals are calculated.  
#' @param W	A non-negative K by K neighborhood matrix (where K is the number of spatial units). 
#' Typically a binary specification is used, where the jkth element equals one if areas (j, k) 
#' are spatially close (e.g. share a common border) and is zero otherwise. 
#' The matrix can be non-binary, but each row must contain at least one non-zero entry.
#' This argument may not need to be specified if \code{adj.graph} is specified instead. 
#' @param adj.graph	Adjacency graph which may be specified instead of the adjacency matrix 
#' matrix. This argument is used if  \code{W} has not been supplied. The argument 
#'  \code{W} is used in case both \code{W} and \code{adj.graph} are supplied.  
#' @param validrows A vector providing the rows of the data frame which 
#' should be used for validation. 
#' The default NULL value instructs that validation will not be performed. 
#' @param scol Either the name (character) or number of the column in the supplied data frame 
#' identifying the spatial units. The program will try to access data[, scol] 
#' to identify the spatial units. If this is omitted, no spatial modeling will be performed. 
#' @param tcol Like the \code{scol} argument for the time identifier. 
#'  Either the name (character) or number of the column in the supplied data frame 
#' identifying the  time indices. The program will try to access data[, tcol] 
#' to identify the time points. If this is omitted, no temporal modeling will be performed.  
#' @inheritParams CARBayes::S.glm 
#' @inheritParams CARBayes::S.CARbym
#' @inheritParams CARBayes::S.CARdissimilarity
#' @inheritParams CARBayes::S.CARleroux
#' @inheritParams CARBayes::S.CARlocalised
#' @inheritParams CARBayes::S.CARmultilevel
#' @inheritParams CARBayesST::ST.CARlinear
#' @inheritParams CARBayesST::ST.CARanova
#' @inheritParams CARBayesST::ST.CARsepspatial
#' @inheritParams CARBayesST::ST.CARar
#' @inheritParams CARBayesST::ST.CARlocalised
#' @inheritParams CARBayesST::ST.CARadaptive
#' @inheritParams CARBayesST::ST.CARclustrends
#' @param prior.mean.delta A vector of prior means for the regression parameters delta 
#' (Gaussian priors are assumed) for the zero probability logistic regression 
#' component of the model. Defaults to a vector of zeros.
#' @param prior.var.delta	A vector of prior variances for the regression parameters 
#' delta (Gaussian priors are assumed) for the zero probability logistic 
#' regression component of the model. Defaults to a vector with values 100000.
#' @param MALA	Logical, should the function use Metropolis adjusted Langevin algorithm (MALA) 
#' updates (TRUE) or simple random walk (FALSE, default) updates for the 
#' regression parameters and random effects.
#' @param G	The maximum number of distinct intercept terms (groups) to allow in the localised model. 
#' @param rho.slo	The value in the interval [0, 1] that the spatial dependence parameter 
#' for the slope of the linear time trend, rho.slo, is fixed at if it should not be estimated. 
#' If this argument is NULL then rho.slo is estimated in the model.
#' @param rho.int	The value in the interval [0, 1] that the spatial dependence parameter for 
#' the intercept of the linear time trend, rho.int, is fixed at if it should not 
#' be estimated. If this argument is NULL then rho.int is estimated in the model.
#' @param interaction	TRUE or FALSE indicating whether the spatio-temporal interaction 
#' random effects should be included. Defaults to TRUE unless family="gaussian" in which 
#' case interactions are not allowed.
#' @param rho	The value in the interval [0, 1] that the spatial dependence parameter rho 
#' is fixed at if it should not be estimated. If this argument is NULL 
#' then rho is estimated in the model. Setting rho=1, reduces the random effects 
#' prior to the intrinsic CAR model but does require epsilon>0.
#' @param epsilon	Diagonal ridge parameter to add to the random effects prior precision matrix, 
#' only required when rho = 1, and the prior precision is improper. Defaults to 0. Only used for adaptive 
#' model fitting in CARBayesST. 
#' @param rho.S The value in the interval [0, 1] that the spatial dependence parameter rho.S is 
#' fixed at if it should not be estimated. If this argument is NULL then rho.S is 
#' estimated in the model.
#' @param rho.T The value in the interval [0, 1] that the temporal dependence parameter 
#' rho.T is fixed at if it should not be estimated. If this argument is NULL 
#' then rho.T is estimated in the model.
#' @param offsetcol Only used in INLA based modeling. The column name or number 
#' in the data frame that should be used as the offset.  
#' @param N MCMC sample size. 
#' @param thin The level of thinning to apply to the MCMC samples to reduce 
#' their temporal autocorrelation. Defaults to 1 (no thinning).
#' @param burn.in How many initial iterations to discard. 
#' Only relevant for MCMC based model fitting, i.e., when package is spBayes or Stan.  
#' @param rseed Random number seed that controls the starting point for the random number stream.
#' A set value is required to help reproduce the results.
#' @param plotit  Logical scalar value: whether to plot the predictions against the observed values.
#' @param verbose Logical scalar value: whether to print various estimates and statistics.
#' @return A list containing:
#' \itemize{
#'    \item params -   A table of parameter estimates  
#'    \item fit -   The fitted model object.    
#'    \item fitteds -   A vector of fitted values.    
#'     \item mchoice -   Calculated model choice statistics if those have been 
#'     requested by the input argument \code{mchoice=T}. Not all model fits will contain 
#'     all the model choice statistics.  
#'    \item residuals -   A vector of residual values.   
#'     \item stats -   The four validation statistics: rmse, mae, crps and coverage. 
#'      This is present only if model validation has been performed.  
#'    \item yobs_preds -   A data frame containing the validation rows of the model 
#'    fitting data  frame. The last five columns of this data frame contains 
#'    the validation prediction summaries: mean, sd, median, and 95\% prediction interval. 
#'    This is present only if model validation has been performed.  
#'    \item valpreds -   A matrix containing the MCMC samples of the validation predictions. 
#'    The dimension of this matrix is the number of validations times the number of retained 
#'    MCMC samples. This is present only if model validation has been performed. 
#'    \item validationplots - Present only if validation has been performed. 
#'    Contains three validation plots with or without segment and 
#'    an ordinary plot.  See \code{\link{obs_v_pred_plot}} for more.   
#'    \item sn -   The number of areal units used in fitting.   
#'    \item tn -  The number of time points used in fitting.  
#'    \item formula - The input formula for the regression part of the model.   
#'     \item scale.transform -   It is there for compatibility with \code{Bsptime} output.    
#'    \item package -   The name of the package used for model fitting.   
#'    \item model -   The name of the fitted model.    
#'    \item call -   The command used to call the model fitting function.   
#'    \item computation.time -   Computation time required to run the model fitting.   
#' }
#' @references
#' \insertAllCited{}
#' @example inst/examples/bcar_examples.R
#' @export
Bcartime <- function(formula, 
                   data,   
                   family, 
                   link = NULL, 
                   trials = NULL, 
                   offsetcol = NULL, 
                   formula.omega = NULL, 
                   scol = NULL, 
                   tcol = NULL,
                   package="CARBayes", 
                   model="glm",  
                   AR = 1, 
                   W = NULL, 
                   adj.graph=NULL,  
                   residtype="response", 
                   interaction =TRUE, 
                   Z = NULL, 
                   W.binary = NULL, 
                   changepoint=NULL, 
                   knots = NULL, 
                   validrows = NULL, 
                  # prior.beta0=NULL,  
                   prior.mean.delta =NULL, 
                   prior.mean.beta =NULL, 
                   prior.var.beta =NULL, 
                   prior.mean.gamma=NULL, 
                   prior.var.gamma=NULL,
                   prior.sigma2=NULL, 
                   prior.tau2 = c(2,1),  
                   prior.delta =NULL, 
                   prior.var.delta = NULL, 
                   prior.lambda =NULL, 
                   prior.nu2 =c(2,1),
                   epsilon=0,   G=NULL, 
                   ind.area =NULL, 
                   trends = NULL, 
                   rho.T =NULL, rho.S=NULL, rho = NULL, rho.slo=NULL, rho.int=NULL, 
                   MALA = FALSE, 
                   N=2000, burn.in=1000, thin=10, rseed =44,
                   Nchains=4, 
                   verbose=TRUE, plotit=TRUE){
  
 set.seed(rseed)
 start.time<-proc.time()[3]
 data <- as.data.frame(data)
 if (!is.data.frame(data)) stop("Need a data frame in the data argument")
 if (!inherits(formula, "formula")) stop("Need a valid formula")
 
 if (family=="binomial") { 
   if (length(trials)<1) stop("Need number of trials for binomial family\n")
   if (length(trials)>1) stop("Too many values to identify the number of trials for binomial family!\n")
   data$Ntrials <- data[, trials]
 } else { 
   data$Ntrials <- 1
 }
 
 
 s1 <- length(scol)
 t1 <- length(tcol)
 indep <- F
 spatialonly <- F
 sptemporal <- F 
 
 if ( (s1==0) & (t1==0) ) { 
   if (verbose) {
   message("No column has been identified as either spatial or temporal identifier columns\n")
   message("Only independent error Bayesian glm will be fitted using the CARBayes package\n")
   message("If this is a mistake, please specify at least the 
       scol argument for spatial data and the tcol argument if data are temporal too.\n") 
   }
   indep <- T
   sn <- nrow(data)
   tn <- 0
 } else if ( (s1==0) & (t1>0) ) { 
   if (verbose) {
   message("It is not possible to fit models for temporal data only.\n")
   message("Only independent error Bayesian glm will be fitted using the CARBayes package\n")
   }
   indep <- T
   sn <- nrow(data)
   tn <- 0
 } else if ( (s1>0) & (t1==0) ) { 
   spatialonly <- T
   if (verbose) {
   message("No temporal column has been supplied, hence only spatial models will be fitted. \n")
   }
   sids <- data[, scol]
   scode <- unique(sids)
   tn <- 1
   implemented <- c("inla", "CARBayes")
   a <- grepl(package, x=implemented, ignore.case = TRUE)
   if (any(a)) { 
     package <- implemented[which(a)]
   } else stop("Wrong package. Please see helpfile")
 } else { 
   sptemporal <- T
   if (verbose) {
   message("Spatio-temporal models will be fitted.\n")
   }
   sids <- data[, scol]
   timeids <- data[, tcol]
   scode <- unique(sids)
   tcode <- unique(timeids)
   tn <- length(tcode)
   implemented <- c("inla", "CARBayesST")
   a <- grepl(package, x=implemented, ignore.case = TRUE)
   if (any(a)) { 
     package <- implemented[which(a)]
   } else stop("Wrong package. Please see the helpfile.")
   
   if (!is.numeric(timeids)) {  
     message("I am not fitting any temporal model using INLA\n")
     message("since the tcol column is not numeric. \n") 
     message("To fit temporal random effects please supply  \n
          numerical values for the temporal identifier column tcol\n")
   }   
   data$timeids <- eval(data[, tcol])
 }
 
 if (package=="inla") {
   if (inlabru::bru_safe_inla()) {
     if (verbose) {
     message("INLA will be used.")
     }
   } else stop("Sorry, INLA is chosen but not available.")
 }
 ## Prepare for validation 
 u <- getXy(formula=formula, data=data)
 ynavec <- u$y
 data$ynavec <- ynavec
 orig_formula <- formula
 
 nvalid <- length(validrows)
 if (nvalid>0) { 
   yholdout <- ynavec[validrows]
   ynavec[validrows] <- NA
   formula <- update(formula, ynavec ~ . ) 
   data$ynavec <- ynavec
   #x <- sample(1:10)
   # [1]  4  5  9  3  8  1  6 10  7  2
   #match(c(4,8),x)
   #x <- sample(1:4,10,replace=TRUE)
   # [1] 3 4 3 3 2 3 1 1 2 2
   #which(x %in% c(2,4))
   # ynavec <- c(NA, 2, NA, 4, NA, 6, NA)
   # validrows <- c(5, 7)
   allmissings <- which(is.na(ynavec))
   allmissings
   validsamongallmissings  <- match(validrows, allmissings) 
   # Find the position of the validation
   # if (verbose) print(validsamongallmissings)
 }
 nmissing <- length(which(is.na(ynavec)))
 # cat("Number missing:", nmissing, "\n")
 
if (indep ==T) {
  results <- BcarBayes(formula=formula, formula.omega=formula.omega, 
                       family=family, data=data, 
                       model = "glm", nmissing = nmissing, 
                       burn.in=burn.in, N=N,
                       thin=thin, 
                       prior.mean.beta=prior.mean.beta, 
                       prior.var.beta=prior.var.beta, prior.nu2=prior.nu2,
                       prior.mean.delta=prior.mean.delta, 
                       prior.var.delta=prior.var.delta, MALA=MALA, 
                       verbose=verbose)
  
} else { 
  
  if (is.null(W)) {
    if (is.null(adj.graph) ) stop("You must specify either the W matrix or 
                                  the adjacency graph")
    if (inlabru::bru_safe_inla()) {
      a <- INLA::inla.read.graph(adj.graph)
      W <- INLA::inla.graph2matrix(a)
    } else stop("Either the adjacency matrix or the INLA package is required to convert adjacency 
                to a an adjacency matrix")
    
  } 
  
  n1 <- nrow(W)
  n2 <- ncol(W)
  if (n1 != n2) stop("W matrix not square! \n")
  sn <- length(scode)
  if (sn !=n1)  stop("Number of areal units and the dimension of W do not match \n")
  if (sn * tn !=nrow(data) ) stop("There is more than one data for each 
                                  space-time combination")
 
 if (spatialonly) { 
   if (package == "CARBayes") {
    results <- BcarBayes(formula = formula, data=data,
    family=family, model=model,  W=W, Z =Z, 
    W.binary = W.binary, G = G, 
    nmissing = nmissing, 
    formula.omega = formula.omega, 
    prior.mean.beta=prior.mean.beta, 
    prior.var.beta=prior.var.beta, prior.tau2=prior.tau2,
    prior.sigma2=prior.sigma2, prior.mean.delta=prior.mean.delta, 
   prior.var.delta=prior.var.delta, 
   prior.nu2 =prior.nu2, rho=rho, 
   ind.area=ind.area,
   MALA=MALA, N=N, burn.in=burn.in, thin=thin,  
   verbose=verbose)
 } else if (package=="inla")  { 
   inlaN <- (N-burn.in)/thin
   newresults <- Bcarinla(data=data, formula=formula,  
                      scol = scol, W=W, adj.graph = adj.graph,
                      sptemporal=FALSE, offsetcol=offsetcol, 
                      family=family, link = link, 
                      prior.nu2 =prior.nu2, prior.tau2 =prior.tau2,
                       verbose=verbose, N=inlaN) 
 } else stop("No other package implemented yet \n")
 }  # Done spatial only 
 if (sptemporal) { 
    # message("I am in sptemporal\n")
     if (package=="CARBayesST") {
       results <- BcarBayesST(formula=formula, 
                              data=data,
                              W= W, 
                              family=family, 
                              model=model,
                              AR = AR, 
                              nmissing = nmissing, 
                              interaction = interaction, 
                              prior.mean.delta =prior.mean.delta, 
                              prior.sigma2=prior.sigma2,
                              prior.tau2 = prior.tau2,  
                              prior.nu2 =prior.nu2,
                              knots=knots, changepoint = changepoint, trends=trends, 
                              Nchains = Nchains, prior.lambda=prior.lambda,
                              prior.mean.gamma=prior.mean.gamma, prior.var.gamma=prior.var.gamma,
                              prior.var.delta = prior.var.delta,  epsilon=epsilon,   G=G, 
                              prior.delta = prior.delta, 
                              rho = rho,  rho.slo=rho.slo, rho.int=rho.int, 
                              rho.S = rho.S, rho.T=rho.T, 
                              N=N, burn.in=burn.in, thin=thin, 
                              verbose=verbose)
     } else if (package=="inla") { 
       inlaN <- (N-burn.in)/thin
       newresults <- Bcarinla(data=data, formula=formula, sptemporal=TRUE, 
                           scol = scol, tcol=tcol, W=W,   offsetcol=offsetcol, 
                           adj.graph = adj.graph, 
                           model=model, family=family, link = link, 
                           prior.nu2 =prior.nu2,  prior.tau2 =prior.tau2, 
                           N=inlaN,  verbose = verbose) 
  } else stop("No other package implemented yet\n")
 } # Done sptemporal 
} # independent model or not 
 
 if (package !="inla") {
 params <- results$summary.results 
 mchoice <- results$modelfit
 newresults <- list(params=params, fit=results, fitteds=fitted(results), mchoice=mchoice)
 if (residtype=="response") { 
   newresults$residuals <- results$residuals$response
 } else if (residtype=="pearson"){ 
   newresults$residuals <- results$residuals$pearson
 } else { 
   if (verbose) {
   message(paste("Your requested residual type ", residtype, " has not been implemented."), "\n")
   message("Returning  the response residuals.\n")
   }
   newresults$residuals <- results$residuals$response
   
 }
  if (nvalid>0) { 
   yits <- t(results$samples$Y[, validsamongallmissings]) ## pickout the validation rows 
  #message("dim yits=", dim(yits), "length yholdout=", length(yholdout))
   sds <- apply(yits, 1, sd) 
   means <- apply(yits, 1, mean) 
   ipreds <- apply(yits, 1, quantile, probs=c(0.5, 0.025, 0.975))
   dim(ipreds)
   psums <- data.frame(obs=yholdout, meanpred=means,  sdpred=sds, medianpred=ipreds[1,], low=ipreds[2, ], up=ipreds[3, ])
  } # nvalid >0 loop
 } else { # Doing it for INLA
   if (nvalid>0) {  
     modfit <- newresults$fit 
     ipreds <- modfit$summary.fitted.values[validrows, ]
     # print(head(ipreds))
     psums <- data.frame(obs=yholdout, meanpred=ipreds$mean,  sdpred=ipreds$sd, medianpred=ipreds$`0.5quant`, low=ipreds$`0.025quant`, up=ipreds$`0.975quant`)
     yits <- matrix(NA, ncol=N, nrow=nvalid)
     for (i in 1:nvalid) {
       yits[i, ]  <- INLA::inla.rmarginal(N, modfit$marginals.fitted.values[[validrows[i]]]) 
     }
    if (family=="binomial")  {
      ntrialsholdout <- data$Ntrials[validrows]
      for (j in 1:ncol(yits)) yits[,j] <- yits[,j] * ntrialsholdout
      for (j in 2:ncol(psums)) psums[,j] <- psums[,j] * ntrialsholdout
    }
   } # nvalid loop done
 } # INLA loop done 
 ## Now doing it for all packages 
 if (nvalid>0) { 
   if (verbose) {
   message("Calculating validation statistics\n This may take a while. \n")
   }
   if (family=="gaussian") {
     bstat <- calculate_validation_statistics(yholdout, yits, summarystat="mean")
   } else { 
     bstat <- calculate_validation_statistics(yholdout, yits, summarystat="median")
   }
   valframe <- data[validrows, ]
   yvalids <- data.frame(valframe,  psums)
   # newresults$stats  <- c(rmse=rmse, mae=mae, cvg=cvg)
   newresults$stats <- bstat$stats
   newresults$yobs_preds <- yvalids 
   newresults$valpreds <- yits
   
   allvplots <- obs_v_pred_plot(yobs=yholdout, predsums = psums, plotit=plotit)
   newresults$validationplots <- allvplots
  
   if (verbose) print(newresults$stats)
   
 }
 
 newresults$sn <- sn
 newresults$tn <- tn
 newresults$formula <- orig_formula 
 newresults$scale.transform <- "NONE"
 newresults$package <- package
 newresults$model <- model 
 newresults$call <- match.call()
 
 class(newresults) <- "bmstdr"
 end.time <- proc.time()[3]
 comp.time<-end.time-start.time
 comp.time<-fancy.time(comp.time)
 newresults$computation.time <- comp.time
 message(comp.time)
 
 newresults 
}

 
BcarBayes <- function(formula,  data,
                      family, 
                      model="glm",  nmissing =0, 
                      W=NULL, 
                      Z = NULL, 
                      W.binary = NULL,
                      G = NULL, 
                      formula.omega=NULL,
                      valids=NULL, 
                      prior.mean.beta=NULL, 
                      prior.var.beta=NULL, 
                      prior.tau2=NULL,
                      prior.sigma2=NULL, 
                      prior.mean.delta=NULL, 
                      prior.var.delta=NULL,  
                      prior.delta = NULL, 
                      prior.nu2 = NULL,
                      rho=NULL, 
                      ind.area=NULL,
                      MALA=FALSE,
                      N=2000, burn.in=1000, thin=10, rseed =44, 
                      verbose=TRUE) {
  ###
  ###
  # message("in BcarBayes N=", N, " burn in=", burn.in, " thin =", thin, "\n")
  if (model=="glm") {
    if (family=="binomial") {
    results <-  CARBayes::S.glm(formula=formula, formula.omega=formula.omega, 
                      family=family, data=data, 
                      trials=data$Ntrials, 
                      burnin=burn.in, n.sample=N,
                      thin=thin, 
                      prior.mean.beta=prior.mean.beta, 
                      prior.var.beta=prior.var.beta, prior.nu2=prior.nu2,
                      prior.mean.delta=prior.mean.delta, 
                      prior.var.delta=prior.var.delta, MALA=MALA, 
                      verbose=verbose)
    } else { 
      results <-  CARBayes::S.glm(formula=formula, formula.omega=formula.omega, 
                                  family=family, data=data, 
                                  burnin=burn.in, n.sample=N,
                                  thin=thin, 
                                  prior.mean.beta=prior.mean.beta, 
                                  prior.var.beta=prior.var.beta, prior.nu2=prior.nu2,
                                  prior.mean.delta=prior.mean.delta, 
                                  prior.var.delta=prior.var.delta, MALA=MALA, 
                                  verbose=verbose)
    }
    
  } else if (model=="bym") { 
    if (family=="binomial") {
       results <-  CARBayes::S.CARbym(formula = formula, formula.omega=formula.omega, 
                family=family, data=data, trials=data$Ntrials, W=W, burnin=burn.in,
             n.sample=N, thin=thin, prior.mean.beta=prior.mean.beta, 
             prior.var.beta=prior.var.beta, 
             prior.tau2=prior.tau2,
             prior.sigma2=prior.sigma2, 
             prior.mean.delta=prior.mean.delta, 
             prior.var.delta=prior.var.delta, 
             MALA=MALA,
             verbose=verbose)
    } else {
      results <-  CARBayes::S.CARbym(formula = formula, formula.omega=formula.omega, 
                                     family=family, data=data,  W=W, burnin=burn.in,
                                     n.sample=N, thin=thin, prior.mean.beta=prior.mean.beta, 
                                     prior.var.beta=prior.var.beta, 
                                     prior.tau2=prior.tau2,
                                     prior.sigma2=prior.sigma2, 
                                     prior.mean.delta=prior.mean.delta, 
                                     prior.var.delta=prior.var.delta, 
                                     MALA=MALA,
                                     verbose=verbose)
    }
  } else if (model=="dissimilarity") {
    results <-    CARBayes::S.CARdissimilarity(formula=formula, family=family, data=data, 
                  trials=data$Ntrials, W=W,
                  Z=Z, W.binary=W.binary, 
                  burnin=burn.in, 
                  n.sample=N, thin=thin, 
                  prior.mean.beta=prior.mean.beta,
                  prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, 
                  prior.tau2=prior.tau2, MALA=MALA, verbose=verbose)
    } else if (model=="leroux") { 
      if (family=="binomial") {
      results <- CARBayes::S.CARleroux(formula=formula, formula.omega=formula.omega, 
                family=family, data=data, trials=data$Ntrials, W=W, burnin=burn.in,
              n.sample=N, thin=thin, prior.mean.beta=prior.mean.beta, 
              prior.var.beta=prior.var.beta,
              prior.nu2= prior.nu2, prior.tau2=prior.tau2, 
              prior.mean.delta=prior.mean.delta, 
              prior.var.delta=prior.var.delta,
              rho=rho, MALA=MALA, verbose=verbose)
      } else {
        results <- CARBayes::S.CARleroux(formula=formula, formula.omega=formula.omega, 
                                         family=family, data=data,  W=W, burnin=burn.in,
                                         n.sample=N, thin=thin, prior.mean.beta=prior.mean.beta, 
                                         prior.var.beta=prior.var.beta,
                                         prior.nu2= prior.nu2, prior.tau2=prior.tau2, 
                                         prior.mean.delta=prior.mean.delta, 
                                         prior.var.delta=prior.var.delta,
                                         rho=rho, MALA=MALA, verbose=verbose)
      }
    } else if (model=="localised") {
       if (nmissing>0) stop("Validation and/or missing response values are not allowed for this model\n")
      if (family=="binomial") {
      results <- CARBayes::S.CARlocalised(formula=formula, family=family, 
                data=data, G=G, trials=data$Ntrials, W=W,
                   burnin=burn.in, n.sample=N, thin=thin, 
                prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta,
                   prior.tau2=prior.tau2,prior.delta=prior.delta, MALA=MALA, 
                verbose=verbose)
      } else {
        results <- CARBayes::S.CARlocalised(formula=formula, family=family, 
                                            data=data, G=G,  W=W,
                                            burnin=burn.in, n.sample=N, thin=thin, 
                                            prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta,
                                            prior.tau2=prior.tau2,prior.delta=prior.delta, MALA=MALA, 
                                            verbose=verbose)
      }
    } else if (model=="multilevel") {
      results <- CARBayes::S.CARmultilevel(formula=formula, family=family, 
                    data=data, trials=data$Ntrials, W=W, ind.area=ind.area,
                     burnin=burn.in, 
                    n.sample=N, thin=thin, 
                    prior.mean.beta=prior.mean.beta, 
                    prior.var.beta=prior.var.beta,
                       prior.nu2=prior.nu2, 
                    prior.tau2=prior.tau2, 
                    prior.sigma2=prior.sigma2, 
                    rho=rho, 
                    MALA=MALA,
                    verbose=verbose)
    }  else { 
      message("Your model has not been implemented in the CARBayes package. \n")
      stop("Try some other model?")
    } 
  results
}

BcarBayesST <- function(formula, data, family, 
                      model="glm",  nmissing=0, 
                      AR =AR, 
                      W=NULL, 
                      Z = NULL, 
                      W.binary = NULL,
                      interaction=TRUE, 
                      G = NULL, 
                      formula.omega=NULL,
                      valids=NULL, 
                      prior.mean.beta=NULL, 
                      prior.var.beta=NULL, 
                      prior.tau2=NULL,
                      prior.sigma2=NULL, 
                      prior.mean.delta=NULL, 
                      prior.var.delta=NULL, 
                      prior.delta = NULL, 
                      prior.nu2 = NULL,
                      rho.slo=NULL, rho.int=NULL, 
                      rho.S=NULL, rho.T = NULL, 
                      prior.mean.alpha=NULL,  prior.var.alpha=NULL, 
                      rho=NULL, 
                      ind.area=NULL,
                      trends=NULL,  
                      changepoint =NULL, 
                      knots = NULL, 
                      MALA=FALSE, epsilon=0,
                      prior.mean.gamma=NULL, prior.var.gamma=NULL, 
                      prior.lambda=NULL,Nchains=4,
                      N=2000, burn.in=1000, thin=10, rseed, 
                      verbose=TRUE) {
  ###
  ###
  if (model=="linear") {
    if (family=="binomial") {
    results <-  CARBayesST::ST.CARlinear(formula=formula,
                      family=family, data=data, 
                      W=W, 
                      trials=data$Ntrials, 
                      burnin=burn.in, n.sample=N,
                      thin=thin, prior.mean.alpha=prior.mean.alpha, 
                      prior.mean.beta=prior.mean.beta, 
                      prior.var.beta=prior.var.beta,
                      prior.var.alpha = prior.var.alpha, 
                      prior.nu2 = prior.nu2, prior.tau2=prior.tau2, 
                      rho.slo=rho.slo, 
                      rho.int=rho.int, 
                      MALA=MALA, 
                      verbose=verbose)
    } else {
      results <-  CARBayesST::ST.CARlinear(formula=formula,
                                           family=family, data=data, 
                                           W=W, 
                                           burnin=burn.in, n.sample=N,
                                           thin=thin, prior.mean.alpha=prior.mean.alpha, 
                                           prior.mean.beta=prior.mean.beta, 
                                           prior.var.beta=prior.var.beta,
                                           prior.var.alpha = prior.var.alpha, 
                                           prior.nu2 = prior.nu2, prior.tau2=prior.tau2, 
                                           rho.slo=rho.slo, 
                                           rho.int=rho.int, 
                                           MALA=MALA, 
                                           verbose=verbose)
    }
  } else if (model=="anova") {
    if (family=="binomial") {
    results <- CARBayesST::ST.CARanova(formula = formula, 
                         family=family, data=data, trials=data$Ntrials, W=W,
                         interaction =interaction, 
                         burnin=burn.in,
                         n.sample=N, thin=thin, prior.mean.beta=prior.mean.beta, 
                         prior.var.beta=prior.var.beta, 
                         prior.tau2=prior.tau2,
                         prior.nu2 = prior.nu2, 
                         rho.S = rho.S, 
                         rho.T = rho.T, 
                         MALA=MALA,
                         verbose=verbose)
    } else {
      results <- CARBayesST::ST.CARanova(formula = formula, 
                                         family=family, data=data,  W=W,
                                         interaction =interaction, 
                                         burnin=burn.in,
                                         n.sample=N, thin=thin, prior.mean.beta=prior.mean.beta, 
                                         prior.var.beta=prior.var.beta, 
                                         prior.tau2=prior.tau2,
                                         prior.nu2 = prior.nu2, 
                                         rho.S = rho.S, 
                                         rho.T = rho.T, 
                                         MALA=MALA,
                                         verbose=verbose)
    }
  } else if (model=="sepspatial") { 
    if (nmissing>0) stop("Validation and/or missing response values are not allowed")
    if (family=="binomial") {
     results <-  CARBayesST::ST.CARsepspatial(formula=formula, family=family, data=data, 
                                    trials=data$Ntrials, W=W,
                                    burnin=burn.in, 
                                    n.sample=N, thin=thin, 
                                    prior.mean.beta=prior.mean.beta,
                                    prior.var.beta=prior.var.beta, 
                                    prior.tau2 = prior.tau2, 
                                    rho.T=rho.T, 
                                    rho.S=rho.S, 
                                    MALA=MALA, verbose=verbose)
    } else {
      results <-  CARBayesST::ST.CARsepspatial(formula=formula, family=family, data=data, 
                                                W=W,
                                               burnin=burn.in, 
                                               n.sample=N, thin=thin, 
                                               prior.mean.beta=prior.mean.beta,
                                               prior.var.beta=prior.var.beta, 
                                               prior.tau2 = prior.tau2, 
                                               rho.T=rho.T, 
                                               rho.S=rho.S, 
                                               MALA=MALA, verbose=verbose)
    }
  } else if (model=="ar") {  
    if (family=="binomial") {
    results <- CARBayesST::ST.CARar(formula=formula, 
                        family=family, data=data,  
                        AR = AR, 
                        trials=data$Ntrials, W=W, burnin=burn.in, 
                        n.sample=N, thin=thin, 
                      prior.mean.beta=prior.mean.beta, 
                      prior.var.beta=prior.var.beta, 
                      prior.nu2=prior.nu2, 
                      prior.tau2=prior.tau2,
                      rho.S=rho.S, rho.T=rho.T, MALA=MALA, verbose=verbose)
    } else {
      results <- CARBayesST::ST.CARar(formula=formula, 
                                      family=family, data=data,  
                                      AR = AR, 
                                      W=W, burnin=burn.in, 
                                      n.sample=N, thin=thin, 
                                      prior.mean.beta=prior.mean.beta, 
                                      prior.var.beta=prior.var.beta, 
                                      prior.nu2=prior.nu2, 
                                      prior.tau2=prior.tau2,
                                      rho.S=rho.S, rho.T=rho.T, MALA=MALA, verbose=verbose)
    }
  } else if (model=="localised") {  
    if (nmissing>0) stop("Validation and/or missing response values are not allowed for this model\n")
    if (family=="binomial") {
     results <- CARBayesST::ST.CARlocalised(formula=formula, family=family, 
                              data=data, G=G, trials=data$Ntrials, W=W,
                              burnin=burn.in, n.sample=N, thin=thin, 
                              prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta,
                              prior.tau2=prior.tau2, prior.delta=prior.delta, MALA=MALA, 
                              verbose=verbose)
    } else{
      results <- CARBayesST::ST.CARlocalised(formula=formula, family=family, 
                                             data=data, G=G,W=W,
                                             burnin=burn.in, n.sample=N, thin=thin, 
                                             prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta,
                                             prior.tau2=prior.tau2, prior.delta=prior.delta, MALA=MALA, 
                                             verbose=verbose)
    }
  } else if (model=="adaptive") {
    if (family=="binomial") {
    results <- CARBayesST::ST.CARadaptive(formula=formula, family=family, 
                              data=data, trials=data$Ntrials, W=W, 
                               burnin=burn.in, 
                              n.sample=N, thin=thin, 
                              prior.mean.beta=prior.mean.beta, 
                              prior.var.beta=prior.var.beta,
                              prior.nu2=prior.nu2, 
                              prior.tau2=prior.tau2, 
                              rho=rho, 
                              epsilon=epsilon, 
                              MALA=MALA,
                              verbose=verbose) 
    } else {
      results <- CARBayesST::ST.CARadaptive(formula=formula, family=family, 
                                            data=data,  W=W, 
                                            burnin=burn.in, 
                                            n.sample=N, thin=thin, 
                                            prior.mean.beta=prior.mean.beta, 
                                            prior.var.beta=prior.var.beta,
                                            prior.nu2=prior.nu2, 
                                            prior.tau2=prior.tau2, 
                                            rho=rho, 
                                            epsilon=epsilon, 
                                            MALA=MALA,
                                            verbose=verbose) 
    }
    
  } else if (model=="clustertrends") {
    if (family=="binomial") {
    results <- CARBayesST::ST.CARclustrends(formula=formula, family=family, 
                                           data=data, trials=data$Ntrials, W=W, 
                                           burnin=burn.in, 
                                           n.sample=N, thin=thin, 
                                           trends=trends,  
                                           changepoint = changepoint,
                                           knots = knots, 
                                           prior.mean.beta=prior.mean.beta, 
                                           prior.var.beta=prior.var.beta,
                                           prior.tau2=prior.tau2, 
                                           prior.mean.gamma=prior.mean.gamma, 
                                           prior.var.gamma=prior.var.gamma,
                                           prior.lambda=prior.lambda, 
                                           MALA=MALA,  Nchains = Nchains, 
                                           verbose=verbose)
    } else {
      results <- CARBayesST::ST.CARclustrends(formula=formula, family=family, 
                                              data=data,  W=W, 
                                              burnin=burn.in, 
                                              n.sample=N, thin=thin, 
                                              trends=trends,  
                                              changepoint = changepoint,
                                              knots = knots, 
                                              prior.mean.beta=prior.mean.beta, 
                                              prior.var.beta=prior.var.beta,
                                              prior.tau2=prior.tau2, 
                                              prior.mean.gamma=prior.mean.gamma, 
                                              prior.var.gamma=prior.var.gamma,
                                              prior.lambda=prior.lambda, 
                                              MALA=MALA,  Nchains = Nchains, 
                                              verbose=verbose)
    }
  }  else { 
    message("Your model has not been implemented in the CARBayes package. \n")
    stop("Try some other model?")
  } 
  results
}



Bcarinla <- function(
  formula, data, W=NULL, adj.graph=NULL,  scol ="spaceid", tcol=NULL, 
  model=c("bym", "iid"),  sptemporal = FALSE,  offsetcol=NULL, 
  link="log",  family="poisson", prior.nu2 =c(2, 1), prior.tau2 =c(2, 1),
  N=1000,  verbose = TRUE) {
  ###

  if (is.null(W)) {
     if (is.null(adj.graph) ) stop("You must specify either the W matrix or the adjacency graph")
     inla.adj <- adj.graph  
  } else {  
    a <- INLA::inla.read.graph(W)
    # INLA::inla.write.graph(a, filename ="inla.graph")
    # file.path(tempdir(), 'inla.graph') <- "inla.graph"
    INLA::inla.write.graph(a, filename =file.path(tempdir(), 'inla.graph'))
    inla.adj <- file.path(tempdir(), 'inla.graph')
  }
  
  spaceid <- data[, scol]
  data$spaceid <- data[, scol]
   # prior for random effect variance 
  prec.prior <- list(prec = list(prior = "loggamma", param = c(prior.tau2[1], prior.tau2[2])))
  if (family=="gaussian") hyper 	<- list(prec = list(prior = "loggamma", param = c(prior.nu2[1], prior.nu2[2])))
  else hyper <- NULL
  
  if (model[1]=="bym") {
    newformula <- update(formula, . ~ . + f(spaceid, model="bym", graph=file.path(tempdir(), 'inla.graph'), constr=TRUE))
  }  else  if (model[1]=="bym2") {
    newformula <- update(formula, . ~ . + f(spaceid, model="bym2", graph=file.path(tempdir(), 'inla.graph'), constr=TRUE))
} else  if (model[1]=="besag") {
  newformula <- update(formula, . ~ . + f(spaceid, model="besag", graph=file.path(tempdir(), 'inla.graph'), constr=TRUE))
} else  if (model[1]=="besag2") {
  newformula <- update(formula, . ~ . + f(spaceid, model="besag2", graph=file.path(tempdir(), 'inla.graph'), constr=TRUE))
} else  if (model[1]=="besagproper") {
  newformula <- update(formula, . ~ . + f(spaceid, model="besagproper", graph=file.path(tempdir(), 'inla.graph'), constr=TRUE))
}else  if (model[1]=="besagproper2") {
  newformula <- update(formula, . ~ . + f(spaceid, model="besagproper2", graph=file.path(tempdir(), 'inla.graph'), constr=TRUE))
  message("Your model has not been implemented in the CARBayes package. \n")
} else  if (model[1]=="iid") {
  newformula <- update(formula, . ~ . + f(spaceid, model="iid", constr=TRUE))
} else  if (model[1]=="none") {
  message("INLA is not fitting any spatial model\n")
} else {
  stop("Your spatial model has not been implemented. Try one of the alternatives?\n")
 } 


  if (sptemporal) {  
    k <- length(model)
    if (k < 2) {
      message("I am not fitting any temporal model using INLA\n")
      message("since the model argument does not contain a temporal model. \n")
    } else { 
      if (model[2] =="ar1") {
       newformula <- update(newformula, . ~ . + f(timeids, model="ar1",constr=TRUE)) 
      } else if (model[2] =="rw1") {
          newformula <- update(newformula, . ~ . + f(timeids, model="rw1",constr=TRUE)) 
      }  else if (model[2] =="rw2") {
          newformula <- update(newformula, . ~ . + f(timeids, model="rw2",constr=TRUE)) 
      } else if (model[2] =="iid") {
          newformula <- update(newformula, . ~ . + f(timeids, model="iid",constr=TRUE)) 
      }   else if (model[2] =="none") {
        message("INLA is not fitting any temporal model\n")
      } else {
          stop("Your temporal model has not been implemented. Try some other alternatives model?\n")
      }
    }
  }
    
 # print("Here are ntrials \n")
 #  ntrials <- as.numeric(data$Ntrials)
 # print(ntrials)
  
 
  if (!is.null(offsetcol)) {   
   ifit <- INLA::inla(newformula, family=family, data=data, 
               offset = data[, offsetcol], Ntrials = Ntrials, 
                 control.family=list(link=link, hyper=hyper),
                 control.predictor=list(link=1, compute=TRUE), 
                 control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, 
                                      return.marginals.predictor=TRUE))
  } else { 
    ifit <- INLA::inla(newformula,family=family,data=data, 
                Ntrials =Ntrials,  control.family=list(link=link, hyper=hyper),
                 control.predictor=list(link=1, compute=TRUE), 
                 control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, 
                return.marginals.predictor=TRUE))
    }
  
  
  # message("Finished INLA fitting \n")
 
  # Fixed effects betas
  fixed.out <- round(ifit$summary.fixed,3)
  if (verbose) print(fixed.out)
  
  p <- nrow(fixed.out)
  beta.samp <- matrix(NA, nrow=N, ncol=p)
  for (i in 1:p) {
    beta.samp[, i] <-  as.vector(INLA::inla.rmarginal(N, ifit$marginals.fixed[[i]]))
  } 
  # colnames(beta.samp) <- rownames(fixed.out)
  samps <- data.frame(beta.samp)
  colnames(samps) <-  rownames(fixed.out)
  head(samps)
  # Hyperparameters sigma2eps and AR(1) a
  rnames <- rownames(ifit$summary.hyperpar)
  rnames
  no_h <- length(rnames) # Number of hyper-parameters
  

  
  if (model[1]=="iid") { 
    a <- grepl("Precision for spaceid", x=rnames, ignore.case = TRUE)
    k <- which(a)
    if (any(a)) { 
      prec.samp 		<- INLA::inla.rmarginal(N, ifit$marginals.hyperpar[[k]])
      sigmasqS.samp <-  1/prec.samp
      summary(sigmasqS.samp)
      samps$tau2 <- sigmasqS.samp 
    } 
  }

  a <- grepl("spatial", x=rnames, ignore.case = TRUE)
  k <- which(a)
  if (any(a)) { 
    sd.samp 		<- INLA::inla.rmarginal(N, ifit$marginals.hyperpar[[k]])
    sigmasq.samp 		<- sd.samp^2
    summary(sigmasq.samp)
    samps$tau2 <- sigmasq.samp
  } 
  # tau2 for spatial variance 
  
  a <- grepl("iid", x=rnames, ignore.case = TRUE)
  k <- which(a)
  if (any(a)) { 
    prec.samp 		<- INLA::inla.rmarginal(N, ifit$marginals.hyperpar[[k]])
    tausq.samp 		<- 1/prec.samp
    summary(tausq.samp)
    samps$sigma2 <- tausq.samp
  } 

  # sigma2 for iid error variance 
  
  a <- grepl("Rho", x=rnames, ignore.case = TRUE)
  k <- which(a)
  if (any(a)) { 
    rho.samp <-  INLA::inla.rmarginal(N, ifit$marginals.hyperpar[[k]])
    summary(rho.samp)
    samps$rho <- rho.samp
  }
  
  a <- grepl("Precision for timeids", x=rnames, ignore.case = TRUE)
  k <- which(a)
  if (any(a)) { 
    prec.samp 		<- INLA::inla.rmarginal(N, ifit$marginals.hyperpar[[k]])
    sigmasqT.samp <-  1/prec.samp
    summary(sigmasqT.samp)
    samps$tau2T <- sigmasqT.samp
  }
  
  a <- grepl("Precision for the Gaussian observation", x=rnames, ignore.case = TRUE)
  k <- which(a)
  if (any(a)) { 
    prec.samp 		<- INLA::inla.rmarginal(N, ifit$marginals.hyperpar[[k]])
    sigmasqGauss.samp <-  1/prec.samp
    summary(sigmasqGauss.samp)
    samps$nu2 <- sigmasqGauss.samp
  }
  params <- get_parameter_estimates(samps)
  allres <- list(params=params, fit=ifit)
  if (verbose) print(round(allres$params, 3))
  n <- nrow(data)
  allres$fitteds  <- ifit$summary.fitted.values$mean[1:n]
  if (family=="binomial")  allres$fitteds <- allres$fitteds * data$Ntrials
  u <- getXy(formula=formula, data=data)
  y <- u$y
  allres$residuals <- y - allres$fitteds  
  # if (mchoice)  {
    means <- ifit$summary.fitted.values$mean[1:n]
    vars <- ifit$summary.fitted.values$sd[1:n]
    gof <- sum((y-means)^2, na.rm=T)
    penalty <- sum(vars[!is.na(y)])
    pmcc <- gof+penalty
    umod <- c(unlist(ifit$dic$p.eff), unlist(ifit$dic$dic), unlist(ifit$waic$p.eff), unlist(ifit$waic$waic),
              gof, penalty, pmcc)
    names(umod) <- c("pdic", "dic", "pwaic", "waic", "gof", "penalty", "pmcc")
    
    allres$mchoice <- umod
    
    if (verbose) print(round(allres$mchoice, 2))
  # }
    
 # message("All done here\n")  
  allres

}