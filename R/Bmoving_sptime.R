#' Model fitting and validation for spatio-temporal data from moving sensors in time. 
#' @inheritParams Bsptime
#' @param coords A vector of size 2 giving the column numbers of the data 
#' frame which contain the coordinates of the data locations.  
#' Here the supplied data frame must contain a column named `time` which 
#' should indicate the time index of the data row. The values in the column `time` 
#' should be positive integers starting from 1.  
#' @param validrows Either a number of randomly selected data rows to validate
#' or a vector giving the row numbers of the data set for validation.
#' @param predspace A 0-1 flag indicating whether spatial predictions are to be made.
#' @param newdata A new data frame with the same column structure as the model fitting data set. 
#' @seealso \code{\link{Bsptime}} for spatio-temporal  model fitting.
#' @return A list containing:
#'   \itemize{
#'    \item params - A table of parameter estimates 
#'    \item fit  -   The fitted model object. 
#'    \item datatostan - A list containing all the information sent to the 
#'    rstan package. 
#'    \item prior.phi.param  -    This contains the values of 
#'    the hyperparameters of the prior distribution for the spatial 
#'    decay parameter \eqn{phi}.   
#'    \item prior.phi  -    This contains the name of  of the prior 
#'    distribution for the spatial decay parameter \eqn{phi}.  
#'    \item validationplots - Present only if validation has been performed. 
#'    Contains three validation plots with or without segment and 
#'    an ordinary plot.  See \code{\link{obs_v_pred_plot}} for more. 
#'    \item fitteds  -   A vector of fitted values.   
#'    \item residuals   -   A vector of residual values.  
#'    \item package   -   The name of the package used for model fitting.  
#'    This is  always stan for this function. 
#'    \item model   -   The name of the fitted model.   
#'    \item call  -   The command used to call the model fitting function.  
#'    \item formula   -   The input formula for the regression part of 
#'    the model.  
#'    \item scale.transform  -   The transformation adopted by the 
#'     input argument with the same name.  
#'    \item sn   -   The number of data locations used in fitting.  
#'    \item tn   -   The number of time points used in fitting for each location.  
#'    \item computation.time  -   Computation time required 
#'    to run the model fitting.     
#' }
#' @examples
#' \donttest{
#' deep <- argo_floats_atlantic_2003[argo_floats_atlantic_2003$depth==3, ]
#' deep$x2inter <- deep$xinter*deep$xinter
#' deep$month <- factor(deep$month)
#' deep$lat2 <- (deep$lat)^2
#' deep$sin1 <- round(sin(deep$time*2*pi/365), 4)
#' deep$cos1 <- round(cos(deep$time*2*pi/365), 4)
#' deep$sin2 <- round(sin(deep$time*4*pi/365), 4)
#' deep$cos2 <- round(cos(deep$time*4*pi/365), 4)
#' deep[, c( "xlat2", "xsin1", "xcos1", "xsin2", "xcos2")] <- 
#' scale(deep[,c("lat2", "sin1", "cos1", "sin2", "cos2")])
#' f2 <- temp ~ xlon + xlat + xlat2+ xinter + x2inter 
#' M2 <- Bmoving_sptime(formula=f2, data = deep, coordtype="lonlat", 
#' coords = 1:2, N=11, burn.in=6, validrows =NULL, mchoice = FALSE)
#' summary(M2)
#' plot(M2)
#' names(M2)
#' }
#' @export
#'
Bmoving_sptime <-  function(formula,  data, coordtype,  coords,  
                        prior.sigma2=c(2, 1),  
                        prior.tau2 = c(2, 1), 
                        prior.phi=NULL, 
                        prior.phi.param = NULL,
                        scale.transform ="NONE",
                        ad.delta = 0.80,  t.depth=12,   s.size=0.01,  
                        N=2500,   burn.in=1000,    no.chains=1,  
                        validrows = 10, predspace =FALSE, newdata=NULL, 
                        mchoice=TRUE,   plotit=FALSE,   rseed=44,   verbose=TRUE, 
                        knots.coords = NULL,  g_size = 5)
{
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
  s1 <- length(coords)
  if (s1 !=2) stop("Need 2-dimensional spatial coordinates in the coords argument\n")
  
  dcoords <- data[, coords]
  ## For predictions at new locations
  data$mod <- 1 ## Flag for modelling data
  data$val <- 0
    
  nvalid <- length(validrows)
  if (length(validrows)==1) {   
    nvalid <- validrows
    n <- nrow(data)
    validrows <- sample(n, nvalid)
  }
  data$val[validrows] <- 1
  vdat <-  data[validrows, ]
  
  if (predspace) { 
    message("For spatial predictions make sure the new data for prediction is in 
        the exact same form as fitting data. Otherwise this method will not work.")
    newdata$mod <- 0
    newdata$val <- 0
    r1 <- ncol(data)
    r2 <- ncol(newdata)
    if (r1 !=r2) stop("new data for predictions should have the 
                      same number of columns as the modelling data")
    newdata$mod <- 0 ## these data won't be modeled 
    newdata$val <- 0
    
    a <- rbind(data, newdata)
    data <- a[order(a$time), ]
    
    dcoords <- data[, coords]
  }
  
 dim(dcoords)
 dim(data)
 #summary(data$temp)
 #summary(a$temp)
  
  ## For knots 
  kmn <- length(knots.coords[,1]) ## length of supplied knots.coords
  gmn <- length(g_size) ## length of given  grid size
    #  message("kmn= ", kmn, " gmn =", gmn, "\n")
    if (gmn ==1) g_size <- rep(g_size, 2)
    if ( (kmn ==0) & (gmn ==0))
      stop("Need either the knots (knots.coords) or the grid size (g_size) for the GPP model")
    
    if ( (kmn > 0) & (gmn >0)) { # both given
      if (kmn != (g_size[1] * g_size[2])) stop("Conflict in knots.coords and grid size.
                                            Specify only one of those two")
    }
    # if ( (kmn == 0) & (gmn >0)) { # only grid size given
    if (gmn >0) { # only grid size given
      xcoord <- c(max(dcoords[,1]), min(dcoords[,1]))
      ycoord <- c(max(dcoords[,2]), min(dcoords[,2]))
      knots.coords <- spTimer::spT.grid.coords(xcoord, ycoord, by=g_size)
    }
    knots.coords <- as.matrix(knots.coords)
    
    b <- c(t(knots.coords)) # b is a vector of longs/lats
    n <- nrow(dcoords)
    B <- matrix(rep(b, each=n), nrow=n) # Replimessagees B as rows n times 
    a <- cbind(dcoords, B) 
    rd <- apply(a, 1, row_dists, coordtype=coordtype)
    Cdist <- t(rd)
    m <- nrow(knots.coords)
    dmat <- dist_mat(knots.coords, coordtype)
    # max(dmat)
    max.d <- max(Cdist)

  # head(data)
  # dim(data)
  
  nvec <- data$time 
  a <- table(nvec)
  length(a)

  t1 <- min(nvec) 
  if (t1 != 1) stop("Time should start at 1") 
  
  tn <- max(nvec)
  k <- sort(unique(nvec))
  
  alltime <- 1:tn
  ontime <- length(k)
  
  if (ontime<tn) { 
    message("There are missing times\n")
    misst <- alltime[-k]
    } else misst <- NULL
  
  start_row <- rep(NA, ontime)
  fin_row <- rep(NA, ontime)
  
  for (i in 1:ontime) { 
   rowids <- which(data$time==k[i])  
   start_row[i] <- min(rowids)
   fin_row[i] <- max(rowids)
  }
    
  u <- as.data.frame(table(nvec)) ## How many observations at each time point 
  u$srow <- start_row
  u$frow <- fin_row
  ots <- as.numeric(as.character(u[,1]))
  nts <- as.numeric(u[,2])
  v <- sort(c(ots, misst)) 
  summary(v-alltime)
 #  u
  
  sn <- nrow(unique(coords))
  n <- nrow(data)
  
  vnames <- all.vars(formula)
 
  u <- getXy(formula=formula, data=data)
  X <- u$X 
  xnames <- colnames(X)
  y <- u$y  
  summary(y)
  yorig <- y
  y[data$mod<1] <- NA ## setting the prediction data to NA
  summary(y)
  vdaty <- y[data$val>0]  ## Saving the y's for validations 
  y[data$val>0] <- NA ## setting the validations to NA
  # print(summary(y))
 
 ## Figure out which values of y are missing
 
  missing_flag <- rep(0, n)
  missing_flag[is.na(y)] <- 1
  ntmiss <- sum(missing_flag)
  ntobs <- n - ntmiss
  
  data_miss_idx <- which(is.na(y))
  data_obs_idx <- which(!is.na(y))
  yobs <- y[data_obs_idx]
  
  data$newy <- y 
  
  a <- data[data$mod>0, ]
  omiss <- which(is.na(y) & (data$mod>0) & (data$val<1))
  vmiss <- which(is.na(y) & (data$val>0))
  pmiss <- which(is.na(y) & (data$mod<1))
  
  a <- 1:ntmiss
  b <- sort(c(omiss, vmiss, pmiss))
  
  valindex <- match(vmiss, data_miss_idx)
  pindex <- match(pmiss, data_miss_idx)
  oindex <- match(omiss, data_miss_idx)
  
  # b <- sort(c(valindex, oindex, pindex))
  # summary(b-a)
  
  ynavec <- y
  if (scale.transform == "SQRT") { 
    if (min(yobs, na.rm=T) < 0) stop("Can't use the square root transformation.  
            \n Negative observations are there in the response. ") 
    yobs <- sqrt(yobs)
    ynavec <- sqrt(ynavec)  ## keeps the modelling y's
  } 
  
  
  if (scale.transform == "LOG") {
    if (min(yobs, na.rm=T) < 0) stop("Can't use the log transformation. 
    \n Negative observations are there in the response.") 
    yobs <- log(yobs)
    ynavec <- log(ynavec) ## keeps the modelling y's
  } 
  
    k <- length(prior.phi) 
    if (k==0) prior.phi <- "Unif"
    if (k>1) stop("Too many prior distributions for phi")
    if (k==1) {
        u <- match(prior.phi, c("Unif", "Gamm", "Cauchy"))
        if (is.na(u)) stop("Sorry, can't handle that prior distribution for phi.\n
        Please specify it as one of the three: Unif, Gamm or Cauchy")
    }
    k <- length(prior.phi.param)
    if (k>2) stop("Too many prior hyper parameters for phi")
    if (prior.phi=="Unif") { 
        if (k<2) prior.phi.param <- 3 * c(1, 100)/max.d
        phidist <- 0
     }
     if (prior.phi=="Gamm") { 
         if (k<2) prior.phi.param <- c(2, 1)
         phidist <- 1
     }
     if (prior.phi=="Cauchy") { 
         if (k<2) prior.phi.param <- c(0, 1)
         phidist <- 2 
     }
  
  p <- ncol(X)
  missing <- 0
  if (ntmiss>0) missing <- 1
  # print(prior.phi.param)
  
  datatostan <- list(n=n, tn=tn, m2=nrow(knots.coords), p=p, 
                     missing=missing,  ntmiss=ntmiss, ntobs = ntobs, 
                     data_miss_idx=as.vector(data_miss_idx),  data_obs_idx =  as.vector(data_obs_idx), 
                     time =data$time, nots=length(ots),  ots = ots, nts=nts, start_row=start_row, fin_row=fin_row,  
                     n_misst=length(misst), misst=misst,  
                     Cdist=Cdist, dmat = dmat,  
                     yobs=yobs,  X=X,
                     sigma2_prior =prior.sigma2, 
                     tau2_prior = prior.tau2, phidist = phidist,
                     prior_phi_param =prior.phi.param)
  
  initfun <- function() {
    # starting values near the lm estimates
    # variations will work as well
    list(sigma_sq = 1, tau_sq=1, beta=rep(0, p), phi=mean(prior.phi.param))
  }
  message("ATTENTION: this run can be computationally intensive!\n")
 
  mfit <- rstan::sampling(stanmodels$ind_gpp_marginal, data=datatostan,  seed =rseed, init=initfun,
                           chains = no.chains, iter = N, warmup = burn.in,
                           control = list(adapt_delta = ad.delta, stepsize=s.size, max_treedepth=t.depth))
  
  # params <- rstan::summary(mfit, pars =c("beta", "tau_sq", "sigma_sq", "phi"), probs = c(.025, .975))
  # b <- round(params$summary, 5)
  # b
  listofdraws <- rstan::extract(mfit)
 
  beta <- listofdraws$beta # N by p 
  tau_sq <- listofdraws$tau_sq
  sigma_sq <- listofdraws$sigma_sq
  phi <- listofdraws$phi
  
  samps <- cbind(beta, tau_sq, sigma_sq, phi)
  params <- get_parameter_estimates(samps)
  rownames(params) <- c(xnames, "tausq", "sigmasq", "phi")
  
  allres <- list(params=params, fit=mfit, datatostan=datatostan)
  allres$prior.phi.param <- prior.phi.param
  allres$prior.phi <- prior.phi
  
  if (verbose)  print(allres$params$summary)
  
  if (nvalid>0) {
    message("validating ", length(vdaty), " space time observations", "\n")
  
    # a <- cbind(missing_flag, val_flag)
    # b <- a[a[,1]>0, ]
    # valindex <- which(b[,2]>0)
    
    ypreds <- listofdraws$z_miss[, valindex]
    
    if (scale.transform == "SQRT") ypreds <-  ypreds^2
    if (scale.transform == "LOG")  ypreds <-  exp(ypreds)
    
    a <- calculate_validation_statistics(yval=vdaty, yits=t(ypreds))
    predsums <- get_validation_summaries(ypreds)
    
    yvalids <- data.frame(vdat, predsums)
    allres$stats <- a$stats
    allres$yobs_preds <- yvalids
    allres$valpreds <- ypreds
    # Added May 17 2022
    allvplots <- obs_v_pred_plot(vdaty, predsums)
    allres$validationplots <- allvplots
    
    if (plotit)  plot(allvplots$pwithseg)
    if (verbose) print(allres$stats)
    
    
  }
  if (predspace) {
    message("Retrieving the requested predictions \n")
    
    ypreds <- listofdraws$z_miss[, pindex]
    if (scale.transform == "SQRT") ypreds <-  ypreds^2
    if (scale.transform == "LOG")  ypreds <-  exp(ypreds)
    predsums <- get_validation_summaries(ypreds)
    preds <-  data.frame(newdata, predsums)
    allres$preds <- preds
    allres$mcmcpred <- ypreds
  }
  
  if (mchoice) {
    ## logliks <- loo::extract_log_lik(gp_fit_stan)
    ## allres$mchoice <- loo::waic(logliks)
  
    message("Calculating model choice statistics\n")
    v <- logliks_from_moving_gpp_marginal_stanfit (y=ynavec, d=datatostan, stanfit=mfit)
    waic_results <- calculate_waic(v$loglik)
    dic_results <- calculate_dic(v$log_full_like_at_thetahat, v$log_full_like_vec)
    pmcc_results <- v$pmcc
    stanmod <- c(unlist(dic_results), unlist(waic_results), unlist(pmcc_results))
    allres$mchoice <-  stanmod
    if (verbose) print(allres$mchoice)
    allres$logliks <- v
    
  } 
  u <- getXy(formula=formula, data=data)
  k <- length(allres$fitteds)
  if (k<nrow(data)) {
    allX <- u$X
    p <- ncol(allX)
    betastar <- allres$params[1:p, 1]
    fits <- as.vector(allX %*% betastar)
    allres$fitteds <- fits 
  } # else fitteds are there already as in spBayes 
  
  y <- u$y  
  if (scale.transform == "SQRT") y <- sqrt(y)
  if (scale.transform == "LOG") y <- log(y)
  allres$residuals <- y - allres$fitteds
  allres$package <- "stan"
  allres$model <- "indep_gpp_marginal" 
  allres$call <- match.call()
  allres$formula <- formula
  allres$scale.transform <- scale.transform
  allres$sn <- nrow(data) 
  allres$tn <- 1 # This is not in standard format 
  colnames(allres$params) <- c("mean", "sd", "2.5%", "97.5%")
  class(allres) <- "bmstdr"
  
  end.time <- proc.time()[3]
  comp.time<-end.time-start.time
  comp.time<-fancy.time(comp.time)
  allres$computation.time<-comp.time
  message(comp.time)
  
  allres
}




## Provides the likelihood evaluations for calculating DIC and WAIC
## Provides the conditional log-likelihoods from the marginal
## model. These are used to calculate WAIC.
## Inputs are: y (on the modelling scale), sn, tn, and the stanfitted object
## The output is a list of
## (i)   log full likelihood at theta hat
## (ii)  N (MCMC) dimensional vector of log full likelihood at N theta samples
## (iii) matrix of N (MCMC) by n (observed data points)
# 
logliks_from_moving_gpp_marginal_stanfit <- function(y, d, stanfit) {
  
 #  y <- ynavec
 # d <- datatostan 
  #  stanfit <- mfit 
  
  listofdraws <- rstan::extract(stanfit)
  
  phi <- listofdraws$phi
  itmax <- length(phi)
  beta <- listofdraws$beta # N by p
  X <- d$X
#  print(dim(beta))
#  print(dim(t(X)))
  xbeta <- beta %*% t(X) # N by n 
  tau_sq <- listofdraws$tau_sq
  sigma_sq <- listofdraws$sigma_sq
  zmiss <- listofdraws$z_miss
  
  ntobs <- d$ntobs
  loglik <- matrix(NA, nrow=itmax, ncol=ntobs)
  yrep <- matrix(NA, nrow=itmax, ncol=ntobs)
  log_full_like_vec <- numeric()
  m2 <- d$m2
  
  nna <- length(y[is.na(y)])
  
  for (it in 1:itmax) {
   # it <- 1
    sigma2 <- sigma_sq[it]
    tau2   <- tau_sq[it]
    phi_it <- phi[it]
    yimputed  <- y
    if (nna >0 ) yimputed[is.na(y)] <- zmiss[it, ]
    
    Sigma <- exp(-phi_it * d$dmat)
    diag(Sigma) <- 1
    Swinv <-  solve(Sigma) 
    Cmat <-  exp(-phi_it * d$Cdist)
   
    cond_mean_vec <- numeric()
    cond_var_vec <- numeric()
    
    u <- 0.0
    for (i in 1:d$nots) { 
     
      zt <-  yimputed[d$start_row[i]:d$fin_row[i]] 
      mut <- xbeta[it, d$start_row[i]:d$fin_row[i]]

      if (d$nts[i] >1) {  #do conditioning 
        Ct <- Cmat[d$start_row[i]:d$fin_row[i] ,]
        St <-  sigma2 * Ct %*% Swinv %*% t(Ct)
        St <- St + tau2 * diag(1, nrow=d$nts[i], ncol=d$nts[i])
        St <- 0.5* (St + t(St))
        Qmat <- solve(St)
        meanmult <- diag(1/diag(Qmat), nrow=d$nts[i], ncol=d$nts[i]) %*% Qmat
        condmean <- zt - meanmult %*% (zt - mut)
        condvar <- 1/diag(Qmat)
        logden_contr <- mnormt::dmnorm(zt, mean=mut, varcov= St, log=T)
        } else   { 
          condmean <- mut 
          condvar <- sigma2+tau2
          logden_contr <- dnorm(zt[1], mean=mut[1], sd =sqrt(sigma2+tau2), log=T)  
        }
        u <- u +  logden_contr
        # message("i=", i, " logden= ", logden_contr, "\n")
     #   message("i=", i, " mean= ", condmean, "\n")
      #  message("i=", i, " var= ", condvar, "\n")
        
        cond_mean_vec[d$start_row[i]:d$fin_row[i]] <- condmean
        cond_var_vec[d$start_row[i]:d$fin_row[i]] <- condvar
     } # i loop 

    log_full_like_vec[it] <- u
    yobs <- y[!is.na(y)]
    ymean <-   cond_mean_vec[!is.na(y)]
    yvar <- cond_var_vec[!is.na(y)]
    loglik[it, ] <- dnorm(yobs, mean=ymean, sd=sqrt(yvar), log=T)
    yrep[it, ] <- ymean + rnorm(ntobs) * sqrt(yvar)
  }
  # print(calculate_waic(loglik))
  
  ## to calculate log full likelihood at theta hat
  ## calculate theta hat first
  
  sigma2 <- mean(sigma_sq)
  tau2   <- mean(tau_sq)
  phi_mean <- mean(phi)
  yimputed  <- y
  
  if (nna >0) {
  zmissmean <- apply(zmiss, 2, mean)
  yimputed[is.na(y)] <- zmissmean
  }
  
  Sigma <- exp(-phi_it * d$dmat)
  diag(Sigma) <- 1
  Swinv <- solve(Sigma)
  meanxbeta <-  apply(xbeta, 2, mean)
  Cmat <-  exp(-phi_mean * d$Cdist)
  
  u <- 0.0
  for (i in 1:d$nots) { 
    zt <-  yimputed[d$start_row[i]:d$fin_row[i]] 
    mut <- meanxbeta[d$start_row[i]:d$fin_row[i]]
    
    if (d$nts[i] >1) { 
      Ct <- Cmat[d$start_row[i]:d$fin_row[i] ,] 
      St <-  sigma2 * Ct %*% Swinv %*% t(Ct)
      St <- St + tau2 * diag(1, nrow=d$nts[i], ncol=d$nts[i])
      St <- 0.5*(St + t(St))
      logden_contr <- mnormt::dmnorm(zt, mean=mut, varcov =St, log=T)
    } else { 
      logden_contr <- dnorm(zt, mean=mut, sd =sqrt(sigma2+tau2), log=T)  
    }
    u <- u +  logden_contr
  } # i loop 
  log_full_like_at_thetahat <- u
  ##
  
  yrepmeans <- as.vector(apply(yrep, 2, mean))
  yrepvars <- as.vector(apply(yrep, 2, var))
  yobs <- y[!is.na(y)]
  gof <-   sum((yobs-yrepmeans)^2)
  penalty <- sum(yrepvars)
  pmcc <- list(gof=gof, penalty=penalty, pmcc=gof+penalty)
  
  list(log_full_like_at_thetahat=log_full_like_at_thetahat,  log_full_like_vec=log_full_like_vec,
       loglik=loglik, pmcc=pmcc)
}

