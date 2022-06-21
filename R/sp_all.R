Blm_sp <- function(formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial,
                   validrows=NULL, scale.transform="NONE",
                   prior.beta0=0, prior.M=0.0001, prior.sigma2=c(2, 1),
                   N=5000, plotit=TRUE, rseed =44,
                   verbose=TRUE, mchoice=TRUE){
  

  r <- length(validrows)
  if (r>0) {
  fdat <- data[(-validrows), ]  ## Fitting data set
  vdat <- data[validrows, ]     ## validation data set
  u <- getXy(formula=formula, data=vdat)
  xpreds <- u$X
  vdaty <- u$y
  } else fdat <- data

  
  fdat <- na.omit(fdat)  ## Removes the NA's
  u <- getXy(formula=formula, data=fdat)
  X <- u$X
  y <- u$y

  n <- length(y)
  p <- ncol(X)
  xnames <- colnames(X)

  if (!is.vector(prior.beta0)) stop("beta0 must be a vector or scalar")
  if (length(prior.beta0) != p ) prior.beta0 <- rep(prior.beta0[1], p)
  if (!is.matrix(prior.M)) prior.M <- diag(prior.M, ncol=p, nrow=p)

  if (scale.transform == "SQRT") y <- sqrt(y)
  if (scale.transform == "LOG")  y <- log(y)


  Mstar <- prior.M + t(X) %*% X
  Mstar_inv <- solve(Mstar)
  # Mstar_inv
  betastar <- Mstar_inv %*% (prior.M %*% prior.beta0 + t(X) %*% y)
  if (verbose) {
  mod1 <- lm(formula=formula, data=fdat)
  round(cbind(mod1$coefficients,  betastar), 2) ## Estimates change slightly
  }

  two_bstar <- 2 * prior.sigma2[2] + t(prior.beta0) %*% prior.M %*% prior.beta0  + sum(y^2) - t(betastar) %*% Mstar %*% betastar
  two_bstar <- as.numeric(two_bstar)
  Var_beta_given_y <- two_bstar * Mstar_inv / (n+ 2* prior.sigma2[1] -2)
  round(sqrt(diag(Var_beta_given_y)), 2)
  gammas <- sqrt(diag(Mstar_inv)) * sqrt(two_bstar/(n+2*prior.sigma2[1]))
  crs_low <- betastar - qt(0.975, df=n+2*prior.sigma2[1]) * gammas
  crs_up <-  betastar + qt(0.975, df=n+2*prior.sigma2[1]) * gammas

  ## Estimation of sigma^2

  sparameter <- 0.5*n+prior.sigma2[1] ## Shape parameter for the posterior distribution of 1/sigma^2
  rparameter <- 0.5* two_bstar ## Rate parameter for the posterior distribution of 1/sigma^2
  sigma2_mean <-  rparameter /(sparameter -1)  ## Estimate of sigma^2
  sigma2_var <- rparameter^2 /((sparameter-1)^2 * (sparameter-2))
  sigma2_low <- 1/qgamma(0.975, shape=sparameter, rate=rparameter)
  sigma2_up <- 1/qgamma(0.025, shape=sparameter, rate=rparameter)

  sigma2_stats <- c(sigma2_mean, sqrt(sigma2_var), sigma2_low, sigma2_up)
  a <- cbind(betastar, sqrt(diag(Var_beta_given_y)),  crs_low, crs_up)
  params_table_lm  <- data.frame(rbind(a, sigma2_stats))
  pnames <- c(paste("beta", 0:(p-1), sep=""), "sigma2")
  dimnames(params_table_lm) <- list(pnames, c("mean", "sd", "low", "up"))
  round(params_table_lm, 2)

  pnames <- c(xnames,  "sigma2")
  dimnames(params_table_lm) <- list(pnames, c("mean", "sd", "low", "up"))
  # round(params_table_sp, 2)
 
  u <- getXy(formula=formula, data=data)
  allX <- u$X
  fits <- as.vector(allX %*% betastar)
  allres <- list(params=params_table_lm, fit=NULL, max.d=NULL, fitteds=fits)
  
  if (mchoice) { # Should we calculate model choice statistics
    ## Going to sampling ...
    set.seed(rseed)
    ##
    sigma2.samples <- 1.0/rgamma(N, shape=sparameter, rate=rparameter)
    summary(sigma2.samples)
    quant(sigma2.samples)
    
    v <- matrix(rnorm(p*N), nrow=p, ncol=N)
    dim(v)
    sqrtmat <- chol(Mstar_inv)
    sqrtmat
    v <- t(sqrtmat) %*% v
    dim(v)
    sigmas <- sqrt(sigma2.samples)
    sigmamat <-  matrix(rep(sigmas, p), nrow=p, ncol=N, byrow=TRUE)
    # sigmamat[, 1:5]
    # sigmas[1:4]
    betamat <- matrix(rep(as.vector(betastar), N), nrow=p, ncol=N)
    # betamat[, 1:5]
    betasamples  <- v * sigmamat + betamat
    ###
    
    lmsamps <- t(rbind(betasamples, sigma2.samples))
    # mc.lmsamps <- as.mcmc(lmsamps)
    # summary(mc.lmsamps)
    # round(params_table_lm, 2)
    
    ## Now calculate the Model choice criteria using sampl
    
    ## DIC calculation for the independent error linear model
    # pars <- params_table_lm[,1] ## Exact pars
    pars <- as.vector(apply(lmsamps, 2, mean)) ## Alternative sampling estimated parameters
    # Does not matter which one we choose
    
    log_lik_at_theta_hat <- log_full_likelihood(pars, y=y, Xmat=X)
    log_liks <- apply(lmsamps, 1, log_full_likelihood, y=y, Xmat=X)
    # log_liks
    # summary(log_liks)  # var(log_liks)
    dic_results_lm <- calculate_dic(log_lik_at_theta_hat, log_liks)
    dic_results_lm
    ## WAIC calculation for the independent regression model
    ## assumes we already have the samples lmsamps
    pars <-  as.vector(apply(lmsamps, 2, mean))
    fitmeans <- X %*% as.vector(pars[1:p])
    logdens <- log_likelihoods_lm(pars=pars, y=y, Xmat=X)
    # dens
    v <- apply(lmsamps, 1, log_likelihoods_lm, y=y, Xmat=X) ## n by N
    waic_results_lm <- calculate_waic(t(v))
    waic_results_lm
    ## waic(t(v))  ## matches with the results from the WAIC function in loo package
    
    ## PMCC for the lm
    set.seed(rseed)
    v <- apply(lmsamps, 1, pred_samples_lm, Xmat=X) ## n by N
    expected_ys <- apply(v, 1, mean)
    var_ys <- apply(v, 1, var)
    gof <- sum((y-expected_ys)^2)
    penalty <- sum(var_ys)
    pmcc <- gof + penalty
    pmcc_results_lm <- list(gof=gof, penalty=penalty, pmcc=pmcc)
    pmcc_results_lm
    #####
    lmmod <- c(unlist(dic_results_lm), unlist(waic_results_lm), unlist(pmcc_results_lm))
    round(lmmod, 2)
    allres$mchoice <-  lmmod
    allres$logliks <- list(log_full_like_at_thetahat=log_lik_at_theta_hat ,  log_full_like_vec=log_liks, loglik=t(v))
    if (verbose) print(round(allres$mchoice, 2))
    
  } # Model choice

  ## validation statistics
  if (r>0) { # Perform validation if there are sites set up for validation


  deltasquare <- diag(1, r, r)
  meanpred <- xpreds %*% betastar
  # meanpred
  vpred <- deltasquare + xpreds %*% Mstar_inv %*% t(xpreds)
  vpred <- vpred * two_bstar /(n+2*prior.sigma2[1] - 2)
  sdpred <- sqrt(diag(vpred))

  ## To calculate CRPS we need to sample from the posterior predictive distribution
  ## We just use the t distribution obtained in the book
  ypreds <- matrix(rt(n=r*N, df=n+2*prior.sigma2[1]), nrow=r, ncol=N)
  sdmat <- matrix(rep(sdpred, N), nrow=r, ncol=N)
  means <- matrix(rep(meanpred, N), nrow=r, ncol=N)
  ypreds <- ypreds*sdmat + means


  if (scale.transform == "NONE") {
      low <- meanpred - qt(0.975, df=n+2*prior.sigma2[1]) * sdpred
      upr <- meanpred + qt(0.975, df=n+2*prior.sigma2[1]) * sdpred
      predsums <- data.frame(meanpred=meanpred, sdpred=sdpred, medianpred=meanpred, low=low, up=upr)

      rmseBlm  <- sqrt(mean((vdaty-meanpred)^2, na.rm=TRUE))
      maeBlm <- mean(abs(vdaty-meanpred), na.rm=TRUE)
      cvgBlm <- spT.pCOVER(vdaty, zup=upr, zlow=low)
      tmp <- cbind(vdaty,ypreds)
      tmp <- na.omit(tmp)
      crpslm <- crpscpp(tmp) # Use C++ to calculate CRPS
      #crpslm <- crps(vdaty, ypreds)
      a <- list(rmse=rmseBlm, mae=maeBlm, crps =crpslm, cvg=cvgBlm)
      results <- list(stats=a)
  }

  if (scale.transform == "SQRT") {
      ypreds <-  ypreds^2
      results  <- calculate_validation_statistics(yval=vdaty, yits=ypreds)
      predsums <- get_validation_summaries(t(ypreds))
    }
    if (scale.transform == "LOG")  {
      ypreds <-  exp(ypreds)
      results  <- calculate_validation_statistics(yval=vdaty, yits=ypreds)
      predsums <- get_validation_summaries(t(ypreds))
    }

    yvalidrows <- cbind(vdat, predsums)
    allres$stats <- results$stats
    allres$yobs_preds <- yvalidrows
    allres$valpreds <- t(ypreds)
    
    allvplots <- obs_v_pred_plot(vdaty, predsums)
    allres$validationplots <- allvplots
    if (plotit)  plot(allvplots$pwithseg)
    
    # if (plotit)  obs_v_pred_plot(vdaty, predsums)
    if (verbose) print(round(unlist(allres$stats), 3))

  } # Validation complete

 
  ####
  allres
}


Bsp_sp <- function(formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial,validrows=NULL,
                            coordtype="utm", coords=4:5, phi=NULL, scale.transform="NONE",
                            prior.beta0=0, prior.M=0.0001, prior.sigma2=c(2, 1),
                           mchoice=TRUE,N=5000, verbose =TRUE,rseed =44,
                           plotit=TRUE){

  set.seed(rseed)
  if (length(coords)==2) coords <-  as.matrix(unique(data[, coords]))  
     
  if (length(validrows)>0) {
    fdat <- data[(-validrows), ]  ## Fitting data set
    vdat <- data[validrows, ]     ## validation data set
    r <- nrow(vdat)
    n <- nrow(fdat)
    
   
    u <- getXy(formula=formula, data=vdat)
    xpreds <- u$X
    vdaty <- u$y
    
    fcoords <- coords[-validrows, ]
    vcoords <- coords[validrows, ]
    allcoords <- rbind(vcoords, fcoords)
    alldistmat <- dist_mat(allcoords, coordtype) # distances in kilometers
    
    d12 <- alldistmat[1:r, (r+1):(r+n)] # is r by n
    d11 <-  alldistmat[1:r, 1:r] # is r by r
  } else { ## We are not doing validation
   fdat <- data
   fcoords <- coords
   n <- nrow(fdat)
  }
  dist_matrix <- dist_mat(fcoords,  coordtype) # distances in kilometers
  summary(as.vector(dist_matrix))
  max.d <- max(dist_matrix)


  u <- getXy(formula=formula, data=fdat)
  X <- u$X
  y <- u$y
  k <- length(y[is.na(y)])
  if (k>0) stop("Can't handle NAs in fitting this model. Please remove the NA data rows and try again")
  
  xnames <- colnames(X)
  p <- ncol(X)

  if (scale.transform == "SQRT") y <- sqrt(y)
  if (scale.transform == "LOG")  y <- log(y)



  if (!is.vector(prior.beta0)) stop("prior.beta0 must be a vector or scalar")
  if (length(prior.beta0) != p ) prior.beta0 <- rep(prior.beta0[1], p)
  if (!is.matrix(prior.M)) prior.M <- diag(prior.M, ncol=p, nrow=p)
  if (length(phi) ==0) phi <- 3/max.d

  ## Calculate correlation matrix

  H <- exp(-phi*dist_matrix)
  Hinv <- solve(H) ## H inverse
  summary(as.vector(H)) ## Looks good values

  Mstar <- prior.M + t(X) %*% Hinv %*% X
  Mstar_inv <- solve(Mstar)
  # Mstar_inv


  betastar <- Mstar_inv %*% (prior.M %*% prior.beta0 + t(X) %*% Hinv %*% y)
  ##
  two_bstar <- 2 * prior.sigma2[2] + t(prior.beta0) %*% prior.M %*% prior.beta0  + quadform(y, Hinv) - t(betastar) %*% Mstar %*% betastar
  two_bstar <- as.numeric(two_bstar)
  Var_beta_given_y <- two_bstar * Mstar_inv / (n+ 2* prior.sigma2[1] -2)
  gammas <- sqrt(diag(Mstar_inv)) * sqrt(two_bstar/(n+2*prior.sigma2[1]))
  gammas
  crs_low <- betastar - qt(0.975, df=n+2*prior.sigma2[1]) * gammas
  crs_up <-  betastar + qt(0.975, df=n+2*prior.sigma2[1]) * gammas

  sparameter <- 0.5 * n + prior.sigma2[1]
  rparameter <- 0.5 * two_bstar
  sigma2_mean <-  rparameter / (sparameter -1)
  # sigma2_mean
  sigma2_var <- rparameter^2 /((sparameter-1)^2 * (sparameter-2))
  sigma2_low <- 1/qgamma(0.975, shape=sparameter, rate=rparameter)
  sigma2_up <- 1/qgamma(0.025, shape=sparameter, rate=rparameter)
  # sigma2_low
  # sigma2_up

  sigma2_stats <- c(sigma2_mean, sqrt(sigma2_var), sigma2_low, sigma2_up)
  a <- cbind(betastar, sqrt(diag(Var_beta_given_y)),  crs_low, crs_up)
  params_table_sp  <- data.frame(rbind(a, sigma2_stats))
  pnames <- c(xnames,  "sigma2")
  dimnames(params_table_sp) <- list(pnames, c("mean", "sd", "low", "up"))
  
  u <- getXy(formula=formula, data=data)
  allX <- u$X
  fits <- as.vector(allX %*% betastar)
  allres <- list(params=params_table_sp, fit=NULL, max.d=max.d, fitteds=fits)
  
  if (verbose) {
    print(round(params_table_sp, 2))
  }
  if (mchoice) { # Perform model choice by sampling
    ## First draw samples from the posterior distribution
    sigma2.samples <- 1.0/rgamma(N, shape=sparameter, rate=rparameter)
    summary(sigma2.samples)
    # quant(sigma2.samples)
    
    v <- matrix(rnorm(p*N), nrow=p, ncol=N)
    dim(v)
    sqrtmat <- chol(Mstar_inv)
    # sqrtmat
    v <- t(sqrtmat) %*% v
    # dim(v)
    sigmas <- sqrt(sigma2.samples)
    sigmamat <-  matrix(rep(sigmas, p), nrow=p, ncol=N, byrow=TRUE)
    betamat <- matrix(rep(as.vector(betastar), N), nrow=p, ncol=N)
    betasamples  <- v * sigmamat + betamat
    ###
    
    spsamps <- t(rbind(betasamples, sigma2.samples))
    # mc.spsamps <- coda::as.mcmc(spsamps)
    # summary(mc.spsamps)
    round(params_table_sp, 2)
    ###
    ## Have to treat y as a single observation from the multivariate normal distribution
    pars <- as.vector(apply(spsamps, 2, mean)) ## this is thetahat Bayes
    
    loglik_at_thetahat <- log_full_likelihood_sp(pars, y=y, Xmat=X, H=H)
    log_liks <- apply(spsamps, 1, log_full_likelihood_sp, y=y, Xmat=X, H=H)
    #log_liks
    # summary(log_liks)
    # var(log_liks)
    # dim(log_liks)
    ## Has two arguments: (1) log full likelihood at thetahat and (2) vector of log-likelihood at the theta samples
    ##
    dic_results_sp <- calculate_dic(loglik_at_thetahat, log_liks)
    dic_results_sp
    
    ## WAIC calculation for the spatial model
    ## assumes we already have the samples
    # dens <- likelihoods_sp(pars=pars)
    # dens
    v <- as.matrix(apply(spsamps, 1, log_likelihoods_sp, y=y, Xmat=X, Hinv=Hinv)) ## n by N
    waic_results_sp <- calculate_waic(t(v))
    waic_results_sp
    # waic(t(v))  ## matches Needs N by n input
    ####
    
    # pars <-  as.vector(apply(spsamps, 2, mean))
    # u <- pred_samples_sp(pars)
    v <- apply(spsamps, 1, pred_samples_sp, Xmat=X, H=H) ## n by N
    expected_ys <- apply(v, 1, mean)
    var_ys <- apply(v, 1, var)
    gof <- sum( (y-expected_ys)^2)
    penalty <- sum(var_ys)
    pmcc <- gof + penalty
    pmcc_results_sp <- list(gof=gof, penalty=penalty, pmcc=pmcc)
    pmcc_results_sp
    spatmod <- c(unlist(dic_results_sp), unlist(waic_results_sp), unlist(pmcc_results_sp))
    allres$mchoice <-  spatmod
    allres$logliks <- list(log_full_like_at_thetahat=loglik_at_thetahat,  log_full_like_vec=log_liks, loglik=t(v))
    if (verbose) print(round(allres$mchoice, 2))
  }
  
  if (length(validrows)>0) { ## We are performing validation
  s12 <-  exp(-phi*d12)
  s11 <- exp(-phi * d11)
#  round(s12, 2)
#  round(s11, 1)

  deltasquare <- diag(1, r, r) - s12 %*% Hinv %*% t(s12)
  meanfac <- s12 %*% Hinv %*% X
  s12hinvy <- s12 %*% Hinv %*% y
  gprime <- xpreds - meanfac


  meanpred <- xpreds %*% betastar + s12hinvy - meanfac %*% betastar
 # meanpred
  vpred <- deltasquare + gprime %*% Mstar_inv %*% t(gprime)
  vpred <- vpred * two_bstar /(n+2*prior.sigma2[1] - 2)
  sdpred <- sqrt(diag(vpred))
  ypreds <- matrix(rt(n=r*N, df=n+2*prior.sigma2[1]), nrow=r, ncol=N)
  sdmat <- matrix(rep(sdpred, N), nrow=r, ncol=N)
  means <- matrix(rep(meanpred, N), nrow=r, ncol=N)
  ypreds <- ypreds*sdmat + means
  dim(ypreds)

   if (scale.transform =="NONE") {
    low <- meanpred - qt(0.975, df=n+2*prior.sigma2[1]) * sdpred
    upr <- meanpred + qt(0.975, df=n+2*prior.sigma2[1]) * sdpred
    predsums <- data.frame(meanpred=meanpred, sdpred=sdpred, medianpred=meanpred, low=low, up=upr)

    rmseBsp  <- sqrt(mean((vdaty-meanpred)^2, na.rm=TRUE))
    maeBsp <- mean(abs(vdaty-meanpred), na.rm=TRUE)
    cvgBsp <- cal_cvg(vdaty, yup=upr, ylow=low)
    tmp <- cbind(vdaty,ypreds)
    tmp <- na.omit(tmp)
    crpsBsp <- crpscpp(tmp) # Use C++ to calculate CRPS
    #crpsBsp <- crps(vdaty, ypreds)
    a <- list(rmse=rmseBsp, mae=maeBsp, crps =crpsBsp, cvg=cvgBsp)
    results <- list(stats=a)
    }

   if (scale.transform == "SQRT") {
      ypreds <-  ypreds^2
      results  <- calculate_validation_statistics(yval=vdaty, yits=ypreds)
      predsums <- get_validation_summaries(t(ypreds))
    }
    if (scale.transform == "LOG")  {
      ypreds <-  exp(ypreds)
      results  <- calculate_validation_statistics(yval=vdaty, yits=ypreds)
      predsums <- get_validation_summaries(t(ypreds))
    }

  yvalidrows <- data.frame(vdat, predsums)
  allres$stats <- results$stats
  allres$yobs_preds <- yvalidrows
  allres$valpreds <- t(ypreds)
  allvplots <- obs_v_pred_plot(vdaty, predsums)
  allres$validationplots <- allvplots
  if (plotit)  plot(allvplots$pwithseg)
  
  # if (plotit)  obs_v_pred_plot(vdaty, predsums)
  if (verbose) print(round(unlist(allres$stats), 3))
  }

  allres$phi <- phi  
 
  ####

  allres
}

BspBayes_sp <- function(formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial,
                                  validrows=NULL,
                                   scale.transform ="NONE",
                                   coordtype="utm", 
                                   coords=4:5,
                        prior.beta0=0, prior.M=0.0001,
                        prior.sigma2 = c(2, 2),
                        prior.tau2 = c(2, 0.1),
                        prior.phi.param=NULL,
                                cov.model = "exponential",
                                n.report = 500,
                                verbose = FALSE,
                                    mchoice=TRUE,
                                N=5000, burn.in=1000, rseed =44, plotit=TRUE,...){

  set.seed(rseed)
  if (length(coords)==2) coords <-  as.matrix(unique(data[, coords]))  
 
    if (length(validrows)>0) {
    fdat <- data[(-validrows), ]  ## Fitting data set
    vdat <- data[validrows, ]     ## validation data set
    r <- nrow(vdat)
    n <- nrow(fdat)
  
    u <- getXy(formula=formula, data=vdat)
    xpreds<- u$X
    vdaty <- u$y

    vcoords <- coords[validrows, ]
    fcoords <- coords[-validrows,]
  } else { ## We are not doing validation
    fdat <- data
    n <- nrow(fdat)
    fcoords <- coords
  }
  alldistmat <- dist_mat(coords, coordtype)
  # summary(as.vector(dist_matrix))
  max.d <- max(alldistmat)

  distmat <-  dist_mat(fcoords, coordtype)
  
  u <- getXy(formula=formula, data=fdat)
  X <- u$X
  y <- u$y

  p <- ncol(X)
  xnames <- colnames(X)[-1]
  if (coordtype=="utm") coords <- coords/1000
  
  ## prior distributions and tuning
  if (length(prior.phi.param) <2) prior.phi.param <- c(0.25, 1) * 3/max.d  #  max(0.25 * 3/max.d, 0.01)

  priors <- list("beta.Norm"=list(rep(prior.beta0, p), diag(prior.M^(-1), p)),
                 "phi.Unif"= prior.phi.param, "sigma.sq.IG"=prior.sigma2,
                 "tau.sq.IG"=prior.tau2)

  phistart <- mean(prior.phi.param)
  starting <- list("phi"=phistart, "sigma.sq"=1, "tau.sq"=1)
  tuning <- list("phi"=0.3, "sigma.sq"=0.3, "tau.sq"=0.3)
  ##

 ## if we transform
  newy <- y
  if (scale.transform == "SQRT") newy <- sqrt(y)
  if (scale.transform == "LOG")  newy <- log(y)
  fdat$newy <- newy

  #b <- "newy ~ "
  #for (j in 1:length(xnames)) {
  #  if (j < length(xnames))
  #    b <- paste(b, xnames[j],  " + ", sep="")
  #  else b <- paste(b, xnames[j], sep="")
  #}
  #newformula <- as.formula(b)
  ## finish transforming
  
  newformula <- update(formula, newy ~ .)
  
  modfit <- spBayes::spLM(formula=newformula, data=fdat, coords=fcoords, starting=starting,
               tuning=tuning, priors=priors, cov.model=cov.model,
               n.samples=N, verbose=verbose, n.report=N/n.report)

  modrecover <- spBayes::spRecover(modfit, start=burn.in+1, verbose=FALSE)
  samps <- data.frame(modrecover$p.beta.recover.samples, modrecover$p.theta.recover.samples)
  params <- get_parameter_estimates(samps)
  
  u <- getXy(formula=formula, data=data)
  allX <- u$X
  fits <- as.vector(allX %*% params[1:p, 1])
  allres <- list(params=params, fit=modfit, max.d=max.d, fitteds=fits)
  
  # if (mchoice)   allres$mchoice <-  spDiag(modrecover)
 
  if (mchoice) {
    # samps <- data.frame(modrecover$p.beta.recover.samples, modrecover$p.theta.recover.samples)
    sn <- n
    betas <- samps[, 1:p]
    tau_sq <- samps$tau.sq
    sigma_sq <- samps$sigma.sq
    itmax <- length(tau_sq)
    yrep <- matrix(NA, nrow=itmax, ncol=n)
    xbeta <- t(X %*% t(betas) )
    loglik <- matrix(NA, nrow=itmax, ncol=n)
    
    log_full_like_vec <- numeric()
    
    for (it in 1:itmax) {
      sigma2 <- sigma_sq[it]
      tau2   <- tau_sq[it]
      phi <- samps$phi[it]
      Sigma <- diag(tau2, nrow=sn, ncol=sn) + sigma2 * exp(-phi * distmat)
      Qmat <- solve(Sigma)
      meanvec <- xbeta[it, ]
      meanmult <- diag(1/diag(Qmat), nrow=sn, ncol=sn) %*% Qmat
      
      condmean <- newy - meanmult %*% (newy - meanvec)
      condvar <- 1/diag(Qmat)
      logden <- mnormt::dmnorm(newy, mean=meanvec, varcov =Sigma, log=TRUE)
      log_full_like_vec[it] <- logden
      
      loglik[it, ] <- dnorm(newy, mean=condmean, sd=sqrt(condvar), log=TRUE)
      yrep[it, ] <- condmean + rnorm(n) * sqrt(condvar)
    }
    # print(calculate_waic(loglik))
    
    ## to calculate log full likelihood at theta hat
    ## calculate theta hat first
    
    sigma2 <- mean(sigma_sq)
    tau2   <- mean(tau_sq)
    Sigma <- diag(tau2, nrow=sn, ncol=sn) + sigma2 * exp(-phi * distmat)
    meanxbeta <-  apply(xbeta, 2, mean)
    logden <- mnormt::dmnorm(newy, mean=meanxbeta, varcov =Sigma, log=TRUE)
    log_full_like_at_thetahat <- logden
    
    yrepmeans <- as.vector(apply(yrep, 2, mean))
    yrepvars <- as.vector(apply(yrep, 2, var))
    gof <-   sum((newy-yrepmeans)^2)
    penalty <- sum(yrepvars)
    pmcc_results <- list(gof=gof, penalty=penalty, pmcc=gof+penalty)
    ###
    waic_results <- calculate_waic(loglik)
    dic_results <- calculate_dic(log_full_like_at_thetahat, log_full_like_vec)
    
    sbBayesmod <- c(unlist(dic_results), unlist(waic_results), unlist(pmcc_results))
    allres$mchoice <-  sbBayesmod
    allres$logliks <- list(log_full_like_at_thetahat=log_full_like_at_thetahat,  log_full_like_vec=log_full_like_vec, loglik=loglik)
    
    allres$spDiag.mchoice <- spBayes::spDiag(modrecover)
    if (verbose) {
      print(allres$mchoice)
      print(allres$spDiag.mchoice)
    }
  }

  if (length(validrows)>0) {
  modpred <- spBayes::spPredict(modfit, pred.covars=xpreds, pred.coords=vcoords,
                        start=burn.in +1, verbose = FALSE)
  ypreds <- modpred$p.y.predictive.samples

  if (scale.transform == "SQRT") ypreds <-  ypreds^2
  if (scale.transform == "LOG")  ypreds <-  exp(ypreds)

  predsums  <- get_validation_summaries(t(ypreds))
  yvalidrows <- data.frame(vdat, predsums)

  b <- calculate_validation_statistics(vdaty, ypreds)
  
 
  allres$stats <- b$stats
  allres$yobs_preds <- yvalidrows
  allres$valpreds <- t(ypreds)
  allvplots <- obs_v_pred_plot(vdaty, predsums)
  allres$validationplots <- allvplots
  if (plotit)  plot(allvplots$pwithseg)
  
  # if (plotit)   obs_v_pred_plot(vdaty, predsums)
  }
  
allres$prior.phi.param <- prior.phi.param

allres
}

# ## @export
Bstan_sp <- function(formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial,
                            validrows=NULL,
                            scale.transform ="NONE",
                            coordtype="utm", 
                            coords=4:5,
                           # prior.beta0=0, prior.M=0.0001,
                            phi=NULL,
                            prior.sigma2=c(2, 2),
                            prior.tau2=c(2, 0.1),
                            verbose = TRUE,
                            mchoice=TRUE,
                            no.chains =1,
                            rseed =44,  ad.delta = 0.99, s.size=0.01,  t.depth=15,
                            N=1500, burn.in=500, plotit=TRUE,...){

 
  if (length(coords)==2) coords <-  as.matrix(unique(data[, coords]))
  
  r <- length(validrows)
  if (r>0) {
    fdat <- data[(-validrows), ]  ## Fitting data set
    vdat <- data[validrows, ]     ## validation data set
    #yXpred <- Formula.matrix(formula=formula, data=vdat)
    #vdaty <- as.vector(yXpred[[1]])
    # X0 <- as.matrix(yXpred[[2]])
    
    u <- getXy(formula=formula, data=vdat)
    X0 <- u$X
    vdaty <- u$y

    vcoords <- coords[validrows, ]
    fcoords <- coords[-validrows, ]
    n <- nrow(fdat)
    allcoords <- rbind(vcoords, fcoords)
    alldistmat <- dist_mat(allcoords, coordtype)
    d12 <- alldistmat[1:r, (r+1):(r+n)] # is r by n
    d11 <-  alldistmat[1:r, 1:r] # is r by r
    distmat <- alldistmat[(r+1):(r+n), (r+1):(r+n)]
  } else { ## We are not doing validation
    fdat <- data
    fcoords <- coords
    n <- nrow(fdat)
    distmat <-  dist_mat(fcoords, coordtype)
  }
  max.d <- max(distmat)

  #yX <- Formula.matrix(formula=formula, data=fdat)
  #y <- as.vector(yX[[1]])
  #X <- as.matrix(yX[[2]])
  
  u <- getXy(formula=formula, data=fdat)
  X <- as.matrix(u$X)
  y <- as.numeric(u$y)

  xnames <- colnames(X)
  p <- ncol(X)
  ## prior distributions and tuning
  if (length(phi)==0) phi <- 3/max.d # max(3/max.d, 0.01)

  ## if we transform
  newy <- y
  if (scale.transform == "SQRT") newy <- sqrt(y)
  if (scale.transform == "LOG")  newy <- log(y)

  datatostan <- list(n=n, p=p, y = newy,  X=X,
                     priorsigma2 = prior.sigma2,
                     priortau2 = prior.tau2,
                     phi=phi, dist=distmat)


  initfun <- function() {
    # starting values near the lm estimates
    # variations will work as well
    mod1 <- lm(formula=newy~-1+X, data=fdat)
    list(sigma_sq = 1, tau_sq=1, beta=coefficients(mod1))
  }

 # message("You must keep the supplied file spatial_model.stan in the sub-folder stanfiles \n")
 #  message("below the current working directory, getwd(). It will give an error if the file is not found.\n")
  if (verbose) message("ATTENTION: this run is likely to be computationally intensive!\n")
  #message("The run with supplied default arguments takes about 15 minutes to run in a fast PC\n")
    
#  stanfit <- stan(data=datatostan, file = "stanfiles/spatial_model.stan", seed = rseed,
#                  chains = no_chains, iter = N, warmup = burn.in, init=initfun,
#                  control = list(adapt_delta = ad.delta, stepsize=0.01,  max_treedepth=t.depth), ...)

  stanfit <-  rstan::sampling(stanmodels$spatial_model,  data=datatostan, seed = rseed,
                  chains = no.chains, iter = N, warmup = burn.in, init=initfun,
                  control = list(adapt_delta = ad.delta, stepsize=s.size,  max_treedepth=t.depth))

  stanestimates  <- rstan::summary(stanfit, pars =c("beta", "tau_sq", "sigma_sq"), probs = c(.025, .975))
  names(stanestimates)
  params <- data.frame(stanestimates$summary)

  listofdraws <- rstan::extract(stanfit)
  beta <- listofdraws$beta # N by p 
  tau_sq <- listofdraws$tau_sq
  sigma_sq <- listofdraws$sigma_sq
  # phi <- listofdraws$phi
  
  samps <- cbind(beta, tau_sq, sigma_sq)
  params <- get_parameter_estimates(samps)
  rownames(params) <- c(xnames, "tausq", "sigmasq")
  
  u <- getXy(formula=formula, data=data)
  allX <- u$X
  fits <- as.vector(allX %*% params[1:p, 1])
  allres <- list(params=params, fit=stanfit, max.d=max.d, fitteds=fits)
  
  
  if (mchoice)   {
    sn <- n
    listofdraws <- rstan::extract(stanfit)
    names(listofdraws)
    betas <- listofdraws$beta
    tau_sq <- listofdraws$tau_sq
    sigma_sq <- listofdraws$sigma_sq
    itmax <- length(tau_sq)
    yrep <- matrix(NA, nrow=itmax, ncol=n)
    xbeta <- betas %*% t(X) # N by n 
    loglik <- matrix(NA, nrow=itmax, ncol=n)

    log_full_like_vec <- numeric()

    for (it in 1:itmax) {
      sigma2 <- sigma_sq[it]
      tau2   <- tau_sq[it]
      Sigma <- diag(tau2, nrow=sn, ncol=sn) + sigma2 * exp(-phi * distmat)
      Qmat <- solve(Sigma)
      meanvec <- as.numeric(xbeta[it, ])
      meanmult <- diag(1/diag(Qmat), nrow=sn, ncol=sn) %*% Qmat
      
      condmean <- newy - meanmult %*% (newy - meanvec)
      condvar <- 1/diag(Qmat)
      logden <- mnormt::dmnorm(newy, mean=meanvec, varcov =Sigma, log=TRUE)
      log_full_like_vec[it] <- logden

      loglik[it, ] <- dnorm(newy, mean=condmean, sd=sqrt(condvar), log=TRUE)
      yrep[it, ] <- condmean + rnorm(n) * sqrt(condvar)
    }
    # print(calculate_waic(loglik))

    ## to calculate log full likelihood at theta hat
    ## calculate theta hat first

    sigma2 <- mean(sigma_sq)
    tau2   <- mean(tau_sq)
    Sigma <- diag(tau2, nrow=sn, ncol=sn) + sigma2 * exp(-phi * distmat)
    meanxbeta <-  apply(xbeta, 2, mean)
    logden <- mnormt::dmnorm(newy, mean=meanxbeta, varcov =Sigma, log=TRUE)
    log_full_like_at_thetahat <- logden

    yrepmeans <- as.vector(apply(yrep, 2, mean))
    yrepvars <- as.vector(apply(yrep, 2, var))
    gof <-   sum((newy-yrepmeans)^2)
    penalty <- sum(yrepvars)
    pmcc_results <- list(gof=gof, penalty=penalty, pmcc=gof+penalty)
    ###
    waic_results <- calculate_waic(loglik)
    dic_results <- calculate_dic(log_full_like_at_thetahat, log_full_like_vec)

    stanspmod <- c(unlist(dic_results), unlist(waic_results), unlist(pmcc_results))
    allres$mchoice <-  stanspmod
    allres$logliks <- list(log_full_like_at_thetahat=log_full_like_at_thetahat,  log_full_like_vec=log_full_like_vec, loglik=loglik)
    if (verbose) print(round(allres$mchoice, 2))
  }

  if (r>0) { ## perform validation

    listofdraws <- rstan::extract(stanfit)
    names(listofdraws)
    betas <- listofdraws$beta
    tau_sq <- listofdraws$tau_sq
    sigma_sq <- listofdraws$sigma_sq
    itmax <- length(tau_sq)
    ypreds <- matrix(NA, nrow=itmax, ncol=r)

    S12 <-  exp(-phi * d12)
    S11 <-  exp(-phi * d11)
    S22 <-  exp(-phi * distmat)
    S22_inv <- solve(S22)
    omega12 <- S12 %*% S22_inv
    omega11 <-  S11 - omega12 %*% t(S12)

    for (it in 1:itmax) {
      #   it <- 1
      sigma2 <- sigma_sq[it]
      tau2   <- tau_sq[it]
      etait <- as.vector(listofdraws$eta[it,])
      mustar <- omega12 %*% etait
      a <- MASS::mvrnorm(n = 1, mu=mustar, Sigma= sigma2 * omega11, tol=10^(-8))
      eta_pred <- as.vector(a)
      betait <- betas[it, ]
      mu <- X0 %*% betait + eta_pred
      v <- rnorm(n=r, 0,  sd=sqrt(tau2) )
      ypreds[it, ] <- (v+mu)
    }

    if (scale.transform == "SQRT") ypreds <-  ypreds^2
    if (scale.transform == "LOG")  ypreds <-  exp(ypreds)

    predsums  <- get_validation_summaries(ypreds)
    yvalidrows <- data.frame(vdat,  predsums)
    b <- calculate_validation_statistics(vdaty, t(ypreds))
   
    allres$stats <- b$stats
    allres$yobs_preds <- yvalidrows
    allres$valpreds <- ypreds
    
    allvplots <- obs_v_pred_plot(vdaty, predsums)
    allres$validationplots <- allvplots
    if (plotit)  plot(allvplots$pwithseg)
    
    # if (plotit)   obs_v_pred_plot(vdaty, predsums)
  } # Else returning the fitted model

  
  allres
}

Binla_sp <- function(formula=yo3~xmaxtemp+xwdsp+xrh, data=nyspatial,
                            validrows=NULL,
                            scale.transform ="NONE",
                            coordtype="utm", 
                            coords=4:5,
                            verbose = TRUE,
                            mchoice=TRUE,
                            prior.tau2=c(2, 1),
                            prior.range= c(1, 0.5),
                            prior.sigma = c(1, 0.005),
                            offset = c(10, 140), 
                            max.edge=c(50, 1000),  
                            N=2000, burn.in =1000,  rseed=44,  plotit=TRUE){


  set.seed(rseed)
  r <- length(validrows)
  Ns <- N-burn.in 
  
  if (length(coords)==2) coords <-  as.matrix(unique(data[, coords]))
  if (coordtype=="lonlat")  stop("Please either supply the coordinates in UTM meters \n 
                                or set the coordtype = plain and re-run. We are not sure if inla can calculate 
                                 geodetic distance from points with latitude and longitude pairs\n. ")
  distmat <-  dist_mat(coords, coordtype=coordtype)
  max.d <- max(distmat)
  
  if (coordtype=="utm") coords <- as.matrix(coords)/1000 ## distance will be in kilometers  
  

  
  if (r>0) {
    fdat <- data[(-validrows), ]  ## Fitting data set
    vdat <- data[validrows, ]     ## validation data set

   
    u <- getXy(formula=formula, data=vdat)
    X0 <- u$X
    vdaty <- u$y

    vcoords <- coords[validrows, ]
    fcoords <- coords[-validrows, ]
    allcoords <- rbind(vcoords, fcoords)
  } else { ## We are not doing validation
    fdat <- data
    fcoords <- coords
  }
  n <- nrow(fdat)
  u <- getXy(formula=formula, data=fdat)
  X <- u$X
  y <- u$y
  xnames <- colnames(X)
  p <- ncol(X)
  ## if we transform
  newy <- y
  if (scale.transform == "SQRT") newy <- sqrt(y)
  if (scale.transform == "LOG")  newy <- log(y)

  mesh <- INLA::inla.mesh.2d(loc=fcoords, offset=offset, max.edge=max.edge)

  spde		<- INLA::inla.spde2.pcmatern(mesh = mesh, alpha = 1.5,
                               prior.range = prior.range, prior.sigma = prior.sigma)

  A 		<- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(fcoords))
  Xcov		<- as.matrix(data.frame(intercept = 1, X[,-1]))

  if (r>0) {
    A.val		<- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(vcoords))
    Xcov.val	<- as.matrix(data.frame(intercept = 1, X0[,-1]))
  }

  stack		<- INLA::inla.stack(tag = 'est', data = list(y = newy), A = list(A, 1),
                       effects = list(se = 1:spde$n.spde, Xcov = Xcov))

  hyper 	<- list(prec = list(prior = "loggamma", param = c(prior.tau2[1], prior.tau2[2])))
  newformula	<- y ~ -1 + Xcov + f(se, model = spde)

  ifit		<- INLA::inla(formula=newformula, data = INLA::inla.stack.data(stack),
                family = "gaussian", control.family = list(hyper = hyper),
                control.predictor = list(A = INLA::inla.stack.A(stack), compute = T),
                control.compute = list(config = TRUE, dic = mchoice, waic = mchoice),
                verbose = F)
  #summary(ifit)

  prec.marg.samp 		<- INLA::inla.rmarginal(Ns, ifit$marginals.hyperpar[[1]])
  marg.tausq.samp 		<- 1/prec.marg.samp
  # summary(marg.tausq.samp)

  range.marg.samp 		<- INLA::inla.rmarginal(Ns, ifit$marginals.hyperpar[[2]])
  marg.phi.samp 		<- 3/range.marg.samp
  # summary(marg.phi.samp)

  sd.marg.samp 		<- INLA::inla.rmarginal(Ns, ifit$marginals.hyperpar[[3]])
  marg.sigmasq.samp 		<- sd.marg.samp^2
  # summary(marg.sigmasq.samp)

  beta.samp <- matrix(NA, nrow=Ns, ncol=p)

  for (i in 1:p) {
    beta.samp[, i] <-  INLA::inla.rmarginal(Ns, ifit$marginals.fixed[[i]])
  }

  dimnames(beta.samp)[[2]] <- xnames
  samps <- data.frame(beta.samp, phi=marg.phi.samp, sigmasq=marg.sigmasq.samp, tausq=marg.tausq.samp)

  params <- get_parameter_estimates(samps)
  
  u <- getXy(formula=formula, data=data)
  allX <- u$X
  fits <- as.vector(allX %*% params[1:p, 1])
  allres <- list(params=params, fit=ifit, max.d=max.d, fitteds=fits)
  
  
  if (mchoice)  {
    n <- length(newy)
    means <- ifit$summary.fitted.values$mean[1:n]
    vars <- (ifit$summary.fitted.values$sd[1:n])^2
    gof <- sum((newy-means)^2, na.rm=TRUE)
    penalty <- sum(vars[!is.na(newy)])
    pmcc <- gof+penalty
    umod <- c(unlist(ifit$dic$p.eff), unlist(ifit$dic$dic), unlist(ifit$waic$p.eff), unlist(ifit$waic$waic),
              gof, penalty, pmcc)
    names(umod) <- c("pdic", "dic", "pwaic", "waic", "gof", "penalty", "pmcc")
    
    allres$mchoice <- umod
    
  }

  if (r>0) {

    ps 		<- INLA::inla.posterior.sample(Ns, ifit) 	## posterior sampling
    contents 	<- ifit$misc$configs$contents

    idX 		<- contents$start[which(contents$tag == "Xcov1")]-1 + (1:p)
    idSpace 	<- contents$start[which(contents$tag == "se")]-1 +
      (1:contents$length[which(contents$tag == "se")])


    xLatent 	<- matrix(0, nrow = length(ps[[1]]$latent), ncol = Ns)
    xHyper  	<- matrix(0, nrow = length(ps[[1]]$hyperpar), ncol = Ns)

    for(j in 1:Ns){
      xLatent[,j] <- ps[[j]]$latent
      xHyper[,j]  <- ps[[j]]$hyperpar
    }

    xSpace 	<- xLatent[idSpace,]
    xX 	 	<- xLatent[idX,]

    sample.IIDval <- matrix(0, r, Ns)

    ID.precision 	   	<- xHyper[1, ] ## precision samples or 3? Nope 1
    tau_sample <- 1/sqrt(ID.precision)
    errs <- matrix(rnorm(r*Ns), nrow=r, ncol=Ns)
    tau_mat  <- matrix(rep(tau_sample, each=r), ncol=Ns)
    err_samp <- errs  * tau_mat

    ypreds	<- as.matrix(A.val %*% xSpace) + as.matrix(Xcov.val) %*% xX + err_samp

    if (scale.transform == "SQRT") ypreds <-  ypreds^2
    if (scale.transform == "LOG")  ypreds <-  exp(ypreds)

    predsums  <- get_validation_summaries(t(ypreds))
    yvalidrows <- data.frame(vdat, predsums)

    b <- calculate_validation_statistics(vdaty, ypreds)
    
    allres$stats <- b$stats
    allres$yobs_preds <- yvalidrows
    allres$valpreds <- t(ypreds)
    allvplots <- obs_v_pred_plot(vdaty, predsums)
    allres$validationplots <- allvplots
    if (plotit)  plot(allvplots$pwithseg)
    # if (plotit)   obs_v_pred_plot(vdaty, predsums)
    if (verbose) print(allres$stats)
  } # Else returning the fitted model


  allres
}
