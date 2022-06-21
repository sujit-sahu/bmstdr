Blm_sptime <- function(formula=y8hrmax~xmaxtemp+xwdsp+xrh, data=nysptime,
                       scale.transform="SQRT",
                       prior.beta0=0, prior.M=0.0001, prior.sigma2=c(2,1),
                       verbose =TRUE, plotit=TRUE,
                       N=1000, validrows=NULL, mchoice=TRUE, rseed=44){
  ###

  nvalid <- length(validrows) 
  n <- nrow(data)

  
  if (nvalid>n) stop("Can't validate after the data rows")
  if (nvalid > 0) { ## perform validation
    fdat <- data[-validrows, ]
    vdat <- data[validrows, ]
    u <- getXy(formula=formula, data=vdat)
    xpreds <- u$X
    vdaty <- u$y
  } else { ## We are not doing validation
    fdat <- data
  }
  n <- nrow(fdat)
  u <- getXy(formula=formula, data=fdat)
  X <- u$X
  y <- u$y
  
  xnames <- colnames(X)
  p <- ncol(X)

  miss <- which(is.na(y))
  nmiss <- length(miss)
  if (length(miss)>0) { # Impute the missing observations
    omean <- mean(y, na.rm=TRUE)
    y[miss] <- omean
    if (verbose) message("Replaced ", nmiss,  " missing observations by the grand mean of the data\n")
  }
  if (scale.transform == "SQRT") { 
    if (min(y, na.rm=TRUE) < 0) stop("Can't use the square root transformation. \n 
                                There are negative observations.") 
    else y <- sqrt(y)
  } 

  if (scale.transform == "LOG") {
    if (min(y, na.rm=TRUE) < 0) stop("Can't use the log transformation. \n 
                                  There are negative observations. \n") 
    else y <- log(y)
  } 


  if (!is.vector(prior.beta0)) stop("prior.beta0 must be a vector or scalar")
  if (length(prior.beta0) != p ) prior.beta0 <- rep(prior.beta0[1], p)
  if (!is.matrix(prior.M)) prior.M <- diag(prior.M, ncol=p, nrow=p)

  Mstar <- prior.M + t(X) %*% X
  Mstar_inv <- solve(Mstar)
  # Mstar_inv
  betastar <- Mstar_inv %*% (prior.M %*% prior.beta0 + t(X) %*% y)


  two_bstar <- 2 * prior.sigma2[2] + t(prior.beta0) %*% prior.M %*% prior.beta0  + sum(y^2) - t(betastar) %*% Mstar %*% betastar
  two_bstar <- as.numeric(two_bstar)
  Var_beta_given_y <- two_bstar * Mstar_inv / (n+ 2* prior.sigma2[1] -2)
  # Var_beta_given_y
  ## Sds of parameter estimates
  round(sqrt(diag(Var_beta_given_y)), 3)


  gammas <- sqrt(diag(Mstar_inv)) * sqrt(two_bstar/(n+2*prior.sigma2[1]))
  ## gammas
  ## 95% Credible interval for betas
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
  pnames <- c(xnames, "sigma2")
  dimnames(params_table_lm) <- list(pnames, c("mean", "sd", "low", "up"))
  round(params_table_lm, 2)
  mod1 <- lm(formula=y~-1+X)
  summary(mod1)
 
  u <- getXy(formula=formula, data=data)
  allX <- u$X
  fitmeans <- as.vector(allX %*% betastar)
  allres <- list(params=params_table_lm, fit=NULL, max.d=NULL, fitteds=fitmeans)
 
  ## Variance to calculate later


  if (verbose)  print(round(allres$params, 3))
  if (mchoice) { # Should we calculate model choice statistics
    ## Going to sampling ...
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
    #  mc.lmsamps <- as.mcmc(lmsamps)
    #  summary(mc.lmsamps)
    #  round(params_table_lm, 2)
    
    ## Now calculate the Model choice criteria using sampl
    
    ## DIC calculation for the independent error linear model
    # pars <- params_table_lm[,1] ## Exact pars
    pars <- as.vector(apply(lmsamps, 2, mean)) ## Alternative sampling estimated parameters
    # Does not matter which one we choose
    
    log_lik_at_theta_hat <- log_full_likelihood(pars, y=y, Xmat=X)
    log_liks <- apply(lmsamps, 1, log_full_likelihood, y=y, Xmat=X)
    # log_liks
    # summary(log_liks)
    # var(log_liks)
    
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
    
    v <- apply(lmsamps, 1, pred_samples_lm, Xmat=X) ## n by N
    expected_ys <- apply(v, 1, mean)
    var_ys <- apply(v, 1, var)
    gof <- sum((y-expected_ys)^2)
    penalty <- sum(var_ys)
    pmcc <- gof + penalty
    pmcc_results_lm <- list(gof=gof, penalty=penalty, pmcc=pmcc)
    pmcc_results_lm
    #####
    lmmod_hand <- c(unlist(dic_results_lm), unlist(waic_results_lm), unlist(pmcc_results_lm))
    round(lmmod_hand, 2)
    allres$mchoice <-  lmmod_hand
    allres$logliks <- list(log_full_like_at_thetahat=log_lik_at_theta_hat,  log_full_like_vec=log_liks, loglik=t(v))
    
    if (verbose) print(round(allres$mchoice, 2))
  } # Model choice
  
  ## validation statistics
  if (nvalid>0) { # Perform validation if there are sites set up for validation
    if (verbose) message("validating ", length(vdaty), " space time observations", "\n")
    k <- nrow(vdat)
    deltasquare <- diag(1, k, k)
    meanpred <- xpreds %*% betastar
    # meanpred
    vpred <- deltasquare + xpreds %*% Mstar_inv %*% t(xpreds)
    vpred <- vpred * two_bstar /(n+2*prior.sigma2[1] - 2)
    sdpred <- sqrt(diag(vpred))

    ypreds <- matrix(rt(n=k*N, df=n+2*prior.sigma2[1]), nrow=k, ncol=N)
    sdmat <- matrix(rep(sdpred, N), nrow=k, ncol=N)
    means <- matrix(rep(meanpred, N), nrow=k, ncol=N)
    ypreds <- ypreds * sdmat + means

    if (scale.transform == "NONE") {
      low <- meanpred - qt(0.975, df=n+2*prior.sigma2[1]) * sdpred
      upr <- meanpred + qt(0.975, df=n+2*prior.sigma2[1]) * sdpred
      # predsums <- data.frame(meanpred=meanpred,  sdpred=sdpred, medianpred=meanpred, low=low, up=upr)
      predsums <- get_validation_summaries(t(ypreds))
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
    # Added May 17 2022
    allvplots <- obs_v_pred_plot(vdaty, predsums)
    allres$validationplots <- allvplots
    if (plotit)  plot(allvplots$pwithseg)
    
    #if (plotit)  obs_v_pred_plot(vdaty, predsums)
    if (verbose) print(round(unlist(allres$stats), 3))

  } # Validation complete

  
  
  allres
}

Bsp_sptime <- function(formula=y8hrmax~xmaxtemp+xwdsp+xrh, data=nysptime, coordtype="utm", coords=4:5,
                       scale.transform ="SQRT", phi.s=NULL, phi.t =NULL,
                       prior.beta0=0, prior.M=0.0001, prior.sigma2=c(2, 1),
                       verbose =TRUE,  plotit=TRUE,
                       N=1000, validrows=NULL, mchoice=TRUE, rseed=44){
    ##
  
 
  nvalid <- length(validrows)
  if (length(coords)==2) coords <-  as.matrix(unique(data[, coords])) 

    
  if (nvalid>0) { ## perform validation
  
    sn <- nrow(coords)
    tn <- nrow(data)/sn
    a <- abs(tn - floor(tn))
    if (a>0) stop("Unequal number of time points: check numbers of locations and times")
    k <- length(data$s.index)
    if (k==0) data$s.index <- rep(1:sn, each=tn)
    
    valids <- unique(data$s.index[validrows])
    r <- nvalid/tn                                
    a <- abs(length(valids) - r)
    if (a>0) stop("Unequal number of validation time points. 
    Please include all the time points in the selected validation sites\n")
   
    fdat <- spT.subset(data=data, var.name=c("s.index"), s=valids, reverse=TRUE)
    vdat <- spT.subset(data=data, var.name=c("s.index"), s=valids)
    
    u <- getXy(formula=formula, data=vdat)
    xpreds <- u$X
    vdaty <- u$y

    fcoords <- coords[-valids, ]
    vcoords <- coords[valids, ]
    allcoords <- rbind(vcoords, fcoords)
  } else { ## We are not doing validation
   fdat <- data
   fcoords <- coords
   allcoords <- coords
   r <- 0
  }

  sn <- nrow(fcoords)
  n <- nrow(fdat)
  tn <- n/sn
  a <- abs(tn - floor(tn))
  if (a>0) stop("Unequal number of time points: check numbers of locations and times")
  # alldistmat <- as.matrix(dist(allcoords)) # distances in kilometers
  alldistmat <- dist_mat(allcoords, coordtype)
  dist_matrix <- alldistmat[(r+1):(r+sn), (r+1):(r+sn)]
  summary(as.vector(alldistmat))
  max.d <- max(alldistmat)

  if (length(phi.s) == 0) phi.s <- 3/max.d
  if (length(phi.t) == 0) phi.t <- 3/tn

  
  u <- getXy(formula=formula, data=fdat)
  X <- u$X
  y <- u$y

  a <- rep(1:tn, each=tn)
  b <- rep(1:tn, tn)
  a <- abs(a-b)
  tcov <- matrix(exp(-phi.t * a), ncol=tn)
  scov <- exp(-phi.s * dist_matrix)
  invscov <- solve(scov)
  invtcov <- solve(tcov)


  xnames <- colnames(X)
  p <- ncol(X)

  miss <- which(is.na(y))
  nmiss <- length(miss)
  if (length(miss)>0) { # Impute the missing observations
    if (verbose) message("This model fitting method cannot handle missing data\n")
    omean <- mean(y, na.rm=TRUE)
    y[miss] <- omean 
    if (verbose) message("Replaced ", nmiss,  " missing observations by the grand mean of the data.
        Because this function cannot handle missing values.\n")
  }
  if (scale.transform == "SQRT") { 
    if (min(y, na.rm=TRUE) < 0) stop("Can't use the square root transformation. \n Negative values in response.") 
    else y <- sqrt(y)
  } 
  
  
  if (scale.transform == "LOG") {
    if (min(y, na.rm=TRUE) < 0) stop("Can't use the log transformation. \n Negative values in response.") 
    else y <- log(y)
  } 
  


  if (!is.vector(prior.beta0)) stop("prior.beta0 must be a vector or scalar")
  if (length(prior.beta0) != p ) prior.beta0 <- rep(prior.beta0[1], p)
  if (!is.matrix(prior.M)) prior.M <- diag(prior.M, ncol=p, nrow=p)


  XHinv <- matrix(0, p, sn*tn)
  for (k in 1:p) {
    Xs <- matrix(X[,k], byrow=TRUE, ncol=tn)
    # print(dim(scov))
    # print(dim(Xs))
    temp1 <- invscov %*% Xs
    temp <- temp1 %*% invtcov
    XHinv[k, ]<- as.vector(t(temp))
  }

  XHinvX <- XHinv %*% X
  invXHX <- solve(XHinvX)


  Mstar <- prior.M +  XHinvX
  Mstar_inv <- solve(Mstar)
  # Mstar_inv
  betastar <- Mstar_inv %*% (prior.M %*% prior.beta0 + XHinv %*% y)
  # mod1 <- lm(formula=y~-1+X)
  #round(cbind(mod1$coefficients,  betastar), 3) ## Estimates change
  ## Error here
  ymat <- matrix(y, byrow=TRUE, ncol=tn)
  temp1 <- invscov %*% ymat
  temp <- temp1 %*% invtcov
  ytHinv <- as.vector(t(temp))

  two_bstar <- 2 * prior.sigma2[2] + t(prior.beta0) %*% prior.M %*% prior.beta0  + sum(y*ytHinv) - t(betastar) %*% Mstar %*% betastar
  two_bstar <- as.numeric(two_bstar)
  Var_beta_given_y <- two_bstar * Mstar_inv / (n+ 2* prior.sigma2[1] -2)
  # Check_beta_given_y
  # round(sqrt(diag(Var_beta_given_y)), 2)
  #summary(mod1)


  gammas <- sqrt(diag(Mstar_inv)) * sqrt(two_bstar/(n+2*prior.sigma2[1]))
  gammas
  crs_low <- betastar - qt(0.975, df=n+2*prior.sigma2[1]) * gammas
  crs_up <-  betastar + qt(0.975, df=n+2*prior.sigma2[1]) * gammas

  # Parameter Estimates

  sparameter <- 0.5*n+prior.sigma2[1]
  rparameter <- 0.5* two_bstar
  sigma2_mean <-  rparameter /(sparameter -1)
  # sigma2_mean
  sigma2_var <- rparameter^2 /((sparameter-1)^2 * (sparameter-2))
  sigma2_low <- 1/qgamma(0.975, shape=sparameter, rate=rparameter)
  sigma2_up <- 1/qgamma(0.025, shape=sparameter, rate=rparameter)
  # sigma2_low
  # sigma2_up

  sigma2_stats <- c(sigma2_mean, sqrt(sigma2_var), sigma2_low, sigma2_up)
  a <- cbind(betastar, sqrt(diag(Var_beta_given_y)),  crs_low, crs_up)
  params_table_sp  <- data.frame(rbind(a, sigma2_stats))
  pnames <- c(xnames, "sigma2")
  dimnames(params_table_sp) <- list(pnames, c("mean", "sd", "low", "up"))

  u <- getXy(formula=formula, data=data)
  allX <- u$X
  fitmeans <- as.vector(allX %*% betastar)
  
  allres <- list(params=params_table_sp, fit=NULL, max.d=max.d, fitteds=fitmeans)
 
  # allres$fit <- fitmeans

  if (verbose)  print(round(allres$params, 3))

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
    #  mc.spsamps <- as.mcmc(spsamps)
    #  summary(mc.spsamps)
    # round(params_table_sp, 2)
    ###
    
    ## Have to treat y as a single observation from the multivariate normal distribution
    
    pars <- as.vector(apply(spsamps, 2, mean)) ## this is thetahat Bayes
    
    ### start here ...
    logdetSinv <- log(det(invscov))
    logdetTinv <- log(det(invtcov))
    
    loglik_at_thetahat <- log_full_likelihood_sptime(pars, y=y, Xmat=X,
                                                     Sinv=invscov, Tinv = invtcov, log_detSinv=logdetSinv, log_detTinv=logdetTinv)
    
    
    log_liks <- apply(spsamps, 1, log_full_likelihood_sptime, y=y, Xmat=X,
                      Sinv=invscov, Tinv = invtcov, log_detSinv=logdetSinv, log_detTinv=logdetTinv)
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
    
    #  v <- as.matrix(apply(spsamps, 1, log_likelihoods_sp, y=y, Xmat=X, H=H)) ## n by N
    v <- as.matrix(apply(spsamps, 1, log_likelihoods_sptime, y=y, Xmat=X,
                         Sinv=invscov, Tinv = invtcov)) ## n by N
    
    waic_results_sp <- calculate_waic(t(v))
    waic_results_sp
    # waic(t(v)  ## matches
    ####
  
    Ls <- chol(scov)
    Lt <- chol(tcov)
    
    pars <-  as.vector(apply(spsamps, 2, mean))
    u <- pred_samples_sptime(pars, y=y, Xmat=X,  Sinv=invscov, Tinv = invtcov)
    
    v <- apply(spsamps, 1, pred_samples_sptime, y=y, Xmat=X,  Sinv=invscov, Tinv = invtcov) ## n by N
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
  
  
  if (r>0) { ## We are performing validation
    if (verbose) message("validating ", length(vdaty), " space time observations", "\n")

    S <- exp(-phi.s * alldistmat)
    S12 <-  as.matrix(S[1:r, (r+1):(r+sn)]) # is r by n
    if (r==1) S12 <- t(S12)
    S11 <-  as.matrix(S[1:r, 1:r]) # is r by r
    S22 <-  as.matrix(S[(r+1):(r+sn), (r+1):(r+sn)]) # sn by sn

    fitmeans <- X %*% as.vector(betastar)
    errbs <- y - fitmeans
    materrbs <- matrix(errbs, byrow=TRUE, ncol=tn) ## This is sn by tn

    # Here is the additional contribution to the mean from Kriging
    # print(dim(S12))
    # print(dim(S22))
    # print(dim(S))
    val.mean.mult <- S12 %*% solve(S22)

    valmean.add <-  val.mean.mult %*% materrbs # this is nvalid by tn
    vmeanadd.vec <- as.vector(t(valmean.add)) # This is r*tn by 1
    meanpred <- xpreds %*% betastar + vmeanadd.vec
    summary(meanpred)
  
    
    u <- S11- S12 %*% solve(S22) %*% t(S12)
    valas <- diag(u)  # this as(s') in paper and nvalid by 1
    csptp <- matrix(rep(valas, tn), ncol=tn) # this is nvalid by tn

    temp <- matrix(0, nrow=r*tn, ncol=p)
    m <- 1
    for (j in 1:r) {
      for (k in 1:tn) {
        for (i in 1:p) {
          a <- as.matrix(X[,i], ncol=sn)
          temp[m, i] <- sum(a[k, ] * val.mean.mult[j, ])
        }
        m <- m+1
      }
    }
    gmat <- xpreds - temp

    gpMinvg <- matrix(NA, nrow=r, ncol=tn)
    m <- 1
    for (j in 1:r) {
      for (k in 1:tn) {
        gpMinvg[j,k] <- quadform(gmat[m, ], Mstar_inv)
        m <- m+1
      }
    }
    vpred <- (csptp + gpMinvg) * two_bstar /(n+2*prior.sigma2[1] - 2)
    sdpred <- sqrt(as.vector(t(vpred)))

    ypreds <- matrix(rt(n=r*tn*N, df=n+2*prior.sigma2[1]), nrow=r*tn, ncol=N)
    sdmat <- matrix(rep(sdpred, N), nrow=r*tn, ncol=N)
    means <- matrix(rep(meanpred, N), nrow=r*tn, ncol=N)
    ypreds <- ypreds * sdmat + means
    #dim(ypreds)
    # range(ypreds)
    # a <- apply(ypreds, 1, mean)

  if (scale.transform =="NONE") {
    low <- meanpred - qt(0.975, df=n+2*prior.sigma2[1]) * sdpred
    upr <- meanpred + qt(0.975, df=n+2*prior.sigma2[1]) * sdpred
    predsums <- data.frame(meanpred=meanpred, sdpred=sdpred, medianpred=meanpred,  low=low, up=upr)
    # predsums <- get_validation_summaries(t(ypreds))  # Alternative using sampling  
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
  # Added May 17 2022
  allvplots <- obs_v_pred_plot(vdaty, predsums)
  allres$validationplots <- allvplots
  if (plotit)  plot(allvplots$pwithseg)
  
  # if (plotit)  obs_v_pred_plot(vdaty, predsums)
  if (verbose) print(round(unlist(allres$stats), 3))


  } ## validation complete

  allres$phi.s <- phi.s
  allres$phi.t <- phi.t
  
  allres
}
##

BspBayes_sptime <- function(formula=y8hrmax~xmaxtemp+xwdsp+xrh, data=nysptime,
                           coordtype="utm", coords=4:5,  scale.transform ="SQRT",
                           prior.beta0=0, prior.M = 0.00001,
                           prior.sigma2=c(2, 25),
                           prior.tau2 =c(2, 25),
                           prior.sigma.eta =c(2, 0.001),
                           prior.phi.param = NULL, phi.tuning =NULL, 
                           cov.model="exponential", 
                           verbose =TRUE, plotit=FALSE,        
                           N=2000,  burn.in=N-999, n.report=N/2, 
                           validrows=NULL, mchoice=TRUE,  rseed=44) {
  ###
 
  nvalid <- length(validrows)
  if (length(coords)==2) coords <-  as.matrix(unique(data[, coords]))
  
 #  message("I am here\n")
  sn <- nrow(coords)
  n <- nrow(data)
  tn <- n/sn
  a <- abs(tn - floor(tn))
  if (a>0) stop("Unequal number of time points: check numbers of locations and times")

  max.d <- max(dist_mat(coords, coordtype)) ## Distance in kilometers
  
  if (coordtype=="utm") coords <- coords/1000 ## Needed to pass to spBayes
  
  u <- getXy(formula=formula, data=data)
  X <- u$X
  y <- u$y
  vnames <- all.vars(formula)
  xnames <- colnames(X)[-1]
  
 
  p <- ncol(X)
  xmat <- matrix(X[,2], byrow=TRUE, ncol=tn)
  if (p > 2) {
    for (j in 2:length(xnames)) {
        a <- matrix(X[,j+1], byrow=TRUE, ncol=tn)
        xmat <- cbind(xmat, a)
    }
  }

  a <- NULL
  for (j in 1:length(xnames)) {
    a <- c(a, paste(xnames[j], 1:tn, sep="."))
  }
  dimnames(xmat)[[2]] <- a

  u <- NULL
  for (i in 1:tn) {
    b <- paste(vnames[1], ".", i, " ~ ", sep="")
    for (j in 1:length(xnames)) {
      if (j < length(xnames))
       b <- paste(b, xnames[j],".", i, " + ", sep="")
      else b <- paste(b, xnames[j],".", i, sep="")
    }
    u[i] <- b
  }

  # mods <- lapply(paste("y8hrmax.", 1:tn, "~xmaxtemp.", 1:tn, "+xwdsp.", 1:tn, "+xrh.", 1:tn, sep=""), as.formula)

  nmods <- lapply(u, as.formula)

  if (nvalid >0) {
    vdaty <- y[validrows]   
    y[validrows] <- NA
    val_flag <- rep(0, n)
    val_flag[validrows] <- 1
    vdat <- data[val_flag>0, ]
  }
  
  ymat <- matrix(y, byrow=TRUE, ncol=tn)
  dimnames(ymat)[[2]] <- paste(vnames[1], 1:tn, sep=".")
  

  if (scale.transform == "SQRT") { 
    if (min(c(ymat), na.rm=TRUE) < 0) stop("Can't use the square root transformation.  
            \n Negative observations are there in the response. ") 
    else ymat <- sqrt(ymat)
  } 
  
  
  if (scale.transform == "LOG") {
    if (min(c(ymat), na.rm=TRUE) < 0) stop("Can't use the log transformation. 
    \n Negative observations are there in the response.") 
    else ymat <- log(ymat)
  } 
  
  fdat <- cbind.data.frame(ymat, xmat)
  head(fdat)

  k <- length(prior.phi.param)
  if (k < 2) { # prior.phi.param <- 3 * c(1/max.d, 100/max.d)
    prior.phi.param <- 3 * c(1/(0.9*max.d), 1/(0.05*max.d))
  } else  stop("Too many prior hyper parameters for phi")
    
  phi.mid <- mean(prior.phi.param[1:2])
    
  k <- length(phi.tuning)
  if (k==1)  tuning <- list("phi"=rep(phi.tuning, tn))
  if (k==0)  tuning <- list("phi"=rep(0.1*phi.mid, tn))
  if (k==tn) tuning <- phi.tuning
  if (k>2 & k!=tn) stop("Correctly specify the tuning parameter for phi")
    
   
  starting <- list("beta"=rep(0, tn*p), # "phi"=rep(phi.mid, tn),
                   "phi"=rep(3/(0.5*max.d), tn),
                   "sigma.sq"=rep(2,tn), "tau.sq"=rep(1, tn),
                   "sigma.eta"=diag(rep(0.01, p)))

  priors <- list("beta.0.Norm"=list(rep(prior.beta0, p), diag(prior.M^(-1),p)),
                 "phi.Unif"=list(rep(prior.phi.param[1], tn), rep(prior.phi.param[2], tn)),
                 "sigma.sq.IG"=list(rep(prior.sigma2[1],tn), rep(prior.sigma2[2],tn)),
                 "tau.sq.IG"=list(rep(prior.tau2[1], tn), rep(prior.tau2[2],tn)),
                 "sigma.eta.IW"=list(prior.sigma.eta[1], diag(prior.sigma.eta[2], p)))

  m.1 <- spBayes::spDynLM(nmods, data=fdat, coords=coords,
                 starting=starting, tuning=tuning, priors=priors, get.fitted =TRUE,
                 cov.model=cov.model, n.samples=N, n.report=N/n.report)

  allres <- list(fit=m.1, max.d=max.d)
  a1 <- m.1$p.beta.samples[burn.in:N,]
  a2 <- m.1$p.theta.samples[burn.in:N,]
  v1 <- get_parameter_estimates(a1)
  v2 <- get_parameter_estimates(a2)
  beta <- apply(m.1$p.beta.samples[burn.in:N,], 2, quant)
  theta <- apply(m.1$p.theta.samples[burn.in:N,], 2, quant)
  allres$params <- rbind(v1, v2)
  # Work on parameters 
  ## plotting removed. Went into sptime R commands in Rfiles 
  
  y <- as.vector(t(ymat))
  ysamples <- m.1$p.y.samples
  means <- apply(ysamples, 1, mean)
  allres$fitteds <- as.vector(means)
  
  if (mchoice==TRUE) {
  vars <- apply(ysamples, 1, var)
  gof <- sum((y-means)^2, na.rm=TRUE)
  penalty <- sum(vars[!is.na(y)]) 
  pmcc <- gof+penalty
  pmcc_results <- list(gof=gof, penalty=penalty, pmcc=pmcc)
  umod <- c(unlist(pmcc_results))
  allres$mchoice <-  umod
  if (verbose) print(round(unlist(allres$mchoice), 2))
  
}


  if (nvalid>0) {
    dim(m.1$p.y.samples)
    a <- m.1$p.y.samples
    ## This is in spTimer arrangement
    ## (s1, t1), (s1, t2), (s1, t3), ... (sn, T)
    dim(a)
    ypreds <- a[val_flag>0, ]
    dim(ypreds)

    if (scale.transform == "SQRT") ypreds <-  ypreds^2
    if (scale.transform == "LOG")  ypreds <-  exp(ypreds)


    if (verbose) message("validating ", length(vdaty), " space time observations", "\n")
    a <- calculate_validation_statistics(yval=vdaty, yits=ypreds)
    predsums <- get_validation_summaries(t(ypreds))

    yvalidrows <- data.frame(vdat, predsums)
    allres$stats  <-  a$stats
    allres$yobs_preds <- yvalidrows
    allres$valpreds <- t(ypreds)
    
    # Added May 17 2022
    allvplots <- obs_v_pred_plot(vdaty, predsums)
    allres$validationplots <- allvplots
    if (plotit)  plot(allvplots$pwithseg)
    
    # if (plotit)  obs_v_pred_plot(vdaty, predsums)
    if (verbose) print(round(unlist(allres$stats), 3))

  } # Finished validation
  
  allres$prior.phi.param <- prior.phi.param
  allres$phi.tuning <- phi.tuning 
  
 allres
}

Bstan_sptime <-  function(data=nysptime, formula=y8hrmax~xmaxtemp+xwdsp+xrh,
                          coordtype="utm", coords=4:5,
                          validrows=NULL, 
                         # prior.beta0=0, beta_prior_var=10^4,              
                          prior.sigma2=c(2, 10), prior.tau2 = c(2, 5),
                          prior.phi = "Unif", prior.phi.param = NULL, 
                          scale.transform ="SQRT",
                          ad.delta = 0.80, t.depth=15, s.size=0.01,
                          N=1100, burn.in=100,  no.chains=1,
                          mchoice=TRUE, plotit=FALSE, rseed=44, verbose=F) {
###
  
 
  nvalid <- length(validrows)
  
  if (length(coords)==2) coords <-  as.matrix(unique(data[, coords]))
  
  sn <- nrow(coords)
  n <- nrow(data)
  tn <- n/sn
  a <- abs(tn - floor(tn))
  if (a>0) stop("Unequal number of time points: check numbers of locations and times")
  nT <- sn*tn

  alldistmat <- as.matrix(dist_mat(coords, coordtype)) # distances in kilometers
  max.d <- max(alldistmat)
  
  u <- getXy(formula=formula, data=data)
  X <- as.matrix(u$X)
  y <- as.numeric(u$y)
  p <- ncol(X)
  vnames <- all.vars(formula)
  xnames <- colnames(X)
  ynavec <- y

  if (nvalid >0) {
    ynavec[validrows] <- NA
    val_flag <- rep(0, nT)
    val_flag[validrows] <- 1
    vdaty <- y[val_flag>0]
    vdat <- data[val_flag>0, ]
  }

## Figure out which values of y are missing
missing_flag <- rep(0, nT)
missing_flag[is.na(ynavec)] <- 1
ntmiss <- sum(missing_flag)
ntobs <- nT - ntmiss
if (ntmiss >0 ) { missing <- 1
} else missing <- 0

data_miss_idx <- which(is.na(ynavec))
data_obs_idx <- which(!is.na(ynavec))
yobs <- y[data_obs_idx]


if (scale.transform == "SQRT") { 
  if (min(yobs, na.rm=TRUE) < 0) stop("Can't use the square root transformation.  
            \n Negative observations are there in the response. ") 
  yobs <- sqrt(yobs)
  ynavec <- sqrt(ynavec)  ## keeps the modelling y's
} 


if (scale.transform == "LOG") {
  if (min(yobs, na.rm=TRUE) < 0) stop("Can't use the log transformation. 
    \n Negative observations are there in the response.") 
  yobs <- log(yobs)
  ynavec <- log(ynavec) ## keeps the modelling y's
} 

    k <- length(prior.phi)
    if (k>1) stop("Too many prior distributions for phi")
    if (k==0) prior.phi <- "Unif"
    if (k==1) {
        u <- match(prior.phi, c("Unif", "Gamm", "Cauchy"))
        if (is.na(u)) stop("Sorry, can't handle that prior distribution for phi.\n
        Please specify it as one of: Unif, Gamma or Cauchy for this model fitting")
    }
    k <- length(prior.phi.param)
    if (k>2) stop("Too many prior hyper parameters for phi")
    if (prior.phi=="Unif") { 
        if (k<2) prior.phi.param <- 3 * c(1/max.d, 100/max.d)
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
        
datatostan <- list(sn=sn, tn=tn, nT=nT, p=p, ntmiss=ntmiss, ntobs = ntobs, missing=missing,
             data_miss_idx=data_miss_idx,  data_obs_idx =  data_obs_idx,
             yobs=yobs,  X=X,
             sigma2_prior=prior.sigma2,
             tau2_prior = prior.tau2,
             phidist = phidist,
             prior_phi_param =prior.phi.param,
             dist = alldistmat, verbose=as.numeric(verbose))

initfun <- function() {
  list(sigma_sq = 1, tau_sq=1, beta=rep(0, p), phi =mean(prior.phi.param))
}
# message("You must keep the supplied file gp_marginal.stan in the sub-folder stanfiles\n")
# message("below the current working directory, getwd(). It will give an error if the file is not found.\n")
if (verbose) message("ATTENTION: the run is likely to be computationally intensive!\n")
# message("The run with supplied default arguments takes about an hour and 10 minutes to run in a fast PC\n")
# gp_fit_stan <- rstan::stan(data=datatostan, file = "gp_marginal.stan", seed =rseed, init=initfun,
#                      chains = no_chains, iter = N, warmup = burn.in,
# 		     control = list(adapt_delta = ad.delta, stepsize=s.size, max_treedepth=t.depth))

gp_fit_stan <- rstan::sampling(stanmodels$gp_marginal, data=datatostan,  seed =rseed, init=initfun,
                           chains = no.chains, iter = N, warmup = burn.in,
                           control = list(adapt_delta = ad.delta, stepsize=s.size, max_treedepth=t.depth))

# u <- rstan::summary(gp_fit_stan, pars =c("beta", "tau_sq", "sigma_sq", "phi"), probs = c(.025, .975))

# u <- a$fit

listofdraws <- rstan::extract(gp_fit_stan)

beta <- listofdraws$beta # N by p 
tau_sq <- listofdraws$tau_sq
sigma_sq <- listofdraws$sigma_sq
phi <- listofdraws$phi

samps <- cbind(beta, tau_sq, sigma_sq, phi)
params <- get_parameter_estimates(samps)
rownames(params) <- c(xnames, "tausq", "sigmasq", "phi")
params
# summary(post.gp)$summary
# summary(lm(yX$Y~-1+yX$X))

betastar <- params[1:p, 1]
fits <- as.vector(X %*% betastar)

allres <- list(params=params, fit=gp_fit_stan, max.d=max.d, fitteds=fits)

if (verbose)  print(allres$params)

  if (mchoice) {
      ## logliks <- loo::extract_log_lik(gp_fit_stan)
      ## allres$mchoice <- loo::waic(logliks)
      if (verbose) message("Calculating model choice statistics\n")
      v <- logliks_from_gp_marginal_stanfit(y=ynavec, X=X, sn=sn, tn=tn,  distmat=alldistmat, stanfit=gp_fit_stan)
      waic_results <- calculate_waic(v$loglik)
      dic_results <- calculate_dic(v$log_full_like_at_thetahat, v$log_full_like_vec)
      pmcc_results <- v$pmcc
      
      stanmod <- c(unlist(dic_results), unlist(waic_results), unlist(pmcc_results))
      allres$mchoice <-  stanmod
      if (verbose) print(allres$mchoice)
      allres$logliks <- v
  }
    
if (nvalid>0) {
  if (verbose) message("validating ", length(vdaty), " space time observations", "\n")
  listofdraws <- rstan::extract(gp_fit_stan)

  a <- cbind(missing_flag, val_flag)
  b <- a[a[,1]>0, ]
  valindex <- which(b[,2]>0) # alternative
 

  preds <- listofdraws$z_miss
  ypreds <- preds[, valindex] # N by r


  if (scale.transform == "SQRT") ypreds <-  ypreds^2
  if (scale.transform == "LOG")  ypreds <-  exp(ypreds)

  a <- calculate_validation_statistics(yval=vdaty, yits=t(ypreds))
  predsums <- get_validation_summaries(ypreds)

  yvalidrows <- data.frame(vdat, predsums)
  allres$stats <- a$stats
  allres$yobs_preds <- yvalidrows
  allres$valpreds <- t(ypreds)
  # Added May 17 2022
  allvplots <- obs_v_pred_plot(vdaty, predsums)
  allres$validationplots <- allvplots
  if (plotit)  plot(allvplots$pwithseg)
  
  # if (plotit)  obs_v_pred_plot(vdaty, predsums)
  if (verbose) print(allres$stats)

}
 
 allres$prior.phi.param <- prior.phi.param
 allres$prior.phi <- prior.phi
allres
}

Binla_sptime <- function(data=nysptime, formula=y8hrmax~xmaxtemp+xwdsp+xrh,
                                         coordtype="utm", coords=4:5,
                                         scale.transform ="SQRT",
                                         validrows=NULL,
                                         prior.tau2 =c(2, 1),
                                        prior.range= c(1, 0.5),
                                         prior.sigma = c(1, 0.005),
                                         offset = c(10, 140), 
                                        max.edge=c(50, 1000),  
                                         N=1000, burn.in=500, mchoice=TRUE,  plotit=TRUE,
                                         verbose = TRUE, rseed=44) {
 ###

 
 nvalid <- length(validrows)
 Ns <- N -burn.in
 
 if (length(coords)==2) coords <-  as.matrix(unique(data[, coords])) 
 if (coordtype=="lonlat")  stop("Please either supply the coordinates in UTM meters \n 
                                or set the coordtype = plain and re-run.")
 if (coordtype=="utm")coords <- as.matrix(coords)/1000 ## distance will be in kilometers
 
  sn <- nrow(coords)
  n <- nrow(data)
  tn <- n/sn
  a <- abs(tn - floor(tn))
  if (a>0) stop("Unequal number of time points: check numbers of locations and times")

  u <- getXy(formula=formula, data=data)
  X <- u$X
  y <- u$y
  
  xnames <- colnames(X)
  vnames <- all.vars(formula)
  
  p <- ncol(X)
  times <- rep(1:tn, each=sn)
  all.locs <- matrix(rep(coords, each=tn), byrow=F, ncol=2) ## can come from data columns too

  if (nvalid >0) {
  # message("Will perform validation\n")
  val_flag <- rep(0, sn*tn)
  val_flag[validrows] <- 1
  vdaty <- y[val_flag>0]
  valframe <- data[val_flag>0, ] 
  y[validrows] <- NA
  }
  
  if (scale.transform == "SQRT") { 
    if (min(y, na.rm=TRUE) < 0) stop("Can't use the square root transformation.  
            \n Negative observations are there in the response. ") 
    else y <- sqrt(y)
  } 
  
  if (scale.transform == "LOG") {
    if (min(y, na.rm=TRUE) < 0) stop("Can't use the log transformation. 
    \n Negative observations are there in the response.") 
    else y <- log(y)
  } 
  
  alldistmat <- as.matrix(dist(coords)) # distances in kilometers
  max.d <- max(alldistmat)
  k <- length(prior.range)
  if (k<2) {
      if (k==0) prior.range = c(0.50*max.d,  0.95)
      if (k==1) prior.range = c(0.50*max.d, NA)
  }

  # max.edge   <- diff(range(coords[,1]))/15
  # bound.outer 	<- diff(range(coords[,2]))/3

 #  mesh <- INLA::inla.mesh.2d(  loc = coords, max.edge = c(1,5)*max.edge, offset = c(max.edge, bound.outer), cutoff = max.edge/5)
  
 mesh <- INLA::inla.mesh.2d(loc=coords, offset=offset, max.edge=max.edge)
 old.par <- par(no.readonly = TRUE)
 on.exit(par(old.par))
 
 if (plotit)  { 
   par(mfrow=c(1, 1))
   plot(mesh)
   points(coords[,1], coords[,2], pch=20, cex=2)
 }

  spde		<- INLA::inla.spde2.pcmatern(mesh = mesh, alpha = 1.5,
                               prior.range = prior.range, prior.sigma = prior.sigma)

  hyper 	<- list(prec = list(prior = "loggamma", param = c(prior.tau2[1], prior.tau2[2])))

  A_est <- INLA::inla.spde.make.A(mesh=mesh, loc=all.locs, group=times, n.group=tn)
  dim(A_est)

  s_index <- INLA::inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde, n.group=tn)
  names(s_index)


  u <- sum(abs(X[,1]-1))
  if (u==0) { 
     Xframe <- data.frame(X[, -1])
  } else Xframe <- data.frame(X)  
  
  
  stack_est <- INLA::inla.stack(data=list(y=y), A=list(A_est, 1),
                         effects=list(c(s_index,list(Intercept=1)), list(Xframe)), tag="est")
            

  stack <- INLA::inla.stack(stack_est)

  xnames <- colnames(Xframe)

  b <- paste0(" y  ~ ")
    for (j in 1:length(xnames)) {
      if (j < length(xnames)) b <- paste0(b, xnames[j], " + ")
      else b <- paste0(b, xnames[j])
    }
  b <- as.formula(b)
  newformula <- update(b, .  ~ . -1 + f(spatial.field, model = spde, group=spatial.field.group, control.group=list(model="ar1")))
  
  # newformula	<- y ~ -1 + Xcov + f(spatial.field, model = spde, group=spatial.field.group, control.group=list(model="ar1"))


  if (verbose) message("ATTENTION: this INLA run is likely to be computationally intensive!\n")
 
  ifit <- INLA::inla(newformula, data=INLA::inla.stack.data(stack, spde=spde), family="gaussian",
               control.family = list(hyper = hyper),
               control.predictor=list(A=INLA::inla.stack.A(stack), compute=TRUE),
               control.compute = list(config = TRUE, dic = mchoice, waic = mchoice))
  if (verbose) message("Finished INLA fitting. \n")

  # Fixed effects betas
  fixed.out <- round(ifit$summary.fixed,3)
  if (verbose) print(fixed.out)
  
  p <- nrow(fixed.out)
  beta.samp <- matrix(NA, nrow=Ns, ncol=p)
  for (i in 1:p) {
    beta.samp[, i] <-  as.vector(INLA::inla.rmarginal(Ns, ifit$marginals.fixed[[i]]))
  } 
  colnames(beta.samp) <- rownames(fixed.out)
  samps <- data.frame(beta.samp)
  
  # Hyperparameters sigma2eps and AR(1) a
  rnames <- rownames(ifit$summary.hyperpar)
  rnames
  no_h <- length(rownames(ifit$summary.hyperpar)) +1 # number of hyper parameters plus 1 
  
  a <- grepl("Rho", x=rnames, ignore.case = TRUE)
  k <- which(a)
  if (any(a)) { 
    rho.samp <-  INLA::inla.rmarginal(Ns, ifit$marginals.hyperpar[[k]])
    summary(rho.samp)
    samps$rho <- rho.samp
  }
  ###
  prec.samp 		<- INLA::inla.rmarginal(Ns, ifit$marginals.hyperpar[[1]])
  tausq.samp 		<- 1/prec.samp
  summary(tausq.samp)
  samps$sigma2eps <- tausq.samp
  
  a <- grepl("Stdev", x=rnames, ignore.case = TRUE)
  k <- which(a)
  if (any(a)) { 
    sd.samp 		<- INLA::inla.rmarginal(Ns, ifit$marginals.hyperpar[[k]])
    sigmasq.samp 		<- sd.samp^2
    summary(sigmasq.samp)
    samps$sig2eta <- sigmasq.samp
  } else { 
    samps$sig2eta <- (prior.sigma[1])^2
  }
  
  a <- grepl("Range", x=rnames, ignore.case = TRUE)
  k <- which(a)
  if ( any(a) ) { 
    #message("I am here")
    range.samp 		<- INLA::inla.rmarginal(Ns, ifit$marginals.hyperpar[[k]])
    phi.samp 		<- 3/range.samp
    summary(phi.samp)
    samps$phi <- phi.samp
  } else { 
    samps$phi <- 1/prior.range[1]
  }
  


  params <- get_parameter_estimates(samps)

  allres <- list(params=params, fit=ifit, max.d=max.d)
  if (verbose) print(round(allres$params, 3))
  
  n <- length(y)
  allres$fitteds  <- ifit$summary.fitted.values$mean[1:n]
  
  if (mchoice)  {
    n <- length(y)
    means <- ifit$summary.fitted.values$mean[1:n]
    vars <- (ifit$summary.fitted.values$sd[1:n])^2
    gof <- sum((y-means)^2, na.rm=TRUE)
    penalty <- sum(vars[!is.na(y)])
    
    #allres$mchoice <- list(pdic=unlist(ifit$dic$p.eff), dic=unlist(ifit$dic$dic), 
    #                      pwaic=unlist(ifit$waic$p.eff), waic=unlist(ifit$waic$waic),
    #                      gof=gof, penalty=penalty, pmcc = gof+penalty)
    pmcc <- gof+penalty
    umod <- c(unlist(ifit$dic$p.eff), unlist(ifit$dic$dic), unlist(ifit$waic$p.eff), unlist(ifit$waic$waic),
              gof, penalty, pmcc)
    names(umod) <- c("pdic", "dic", "pwaic", "waic", "gof", "penalty", "pmcc")
    
    allres$mchoice <- as.data.frame(t(umod))
    
    if (verbose) print(allres$mchoice)
  }
  
  
  ###
  if (nvalid>0) {  
   if (verbose) message("validating ", length(vdaty), " space time observations", "\n")
  
    ypreds <- matrix(NA,  nrow=nvalid, ncol=Ns)
    for (i in 1:nvalid) {
      isamples <- INLA::inla.rmarginal(Ns, ifit$marginals.fitted.values[[validrows[i]]]) 
      ypreds[i, ] <- isamples
    }
  
    if (scale.transform == "SQRT") {
      ypreds <-  (ypreds)^2
    }
     if (scale.transform == "LOG")  {
      ypreds <-  exp(ypreds)
    }
    predsums <- get_validation_summaries(t(ypreds))
    b <- calculate_validation_statistics(vdaty, ypreds)
   
    yvalidrows <- data.frame(valframe,  predsums)
    allres$stats  <- b$stats
    allres$yobs_preds <- yvalidrows
    allres$valpreds <- t(ypreds)
    
    allvplots <- obs_v_pred_plot(vdaty, predsums)
    allres$validationplots <- allvplots
    if (plotit)  plot(allvplots$pwithseg)
    
    # if (plotit)  obs_v_pred_plot(vdaty, predsums)
    if (verbose) print(round(unlist(allres$stats), 3))

  }
 
  allres$prior.range <- prior.range
 
  allres
}
##

BspTDyn_sptime <- function(data=nysptime, formula=y8hrmax~xmaxtemp+sp(xmaxtemp)+tp(xwdsp)+xrh, model="GP",
                            coordtype="utm", coords=4:5, time.data=NULL, 
                            validrows=NULL,   scale.transform ="SQRT",
                            prior.beta0=0, prior.M=0.001, prior.sigma2 =c(2, 1),
                           prior.phi="Gamm", prior.phi.param =NULL,
                           phi.tuning=NULL, phi.npoints=NULL,
                            rhotp = 0, cov.fnc="exponential", truncation.para = NULL, 
                            N=5000, burn.in=1000, plotit=TRUE, n.report=10, 
                            mchoice=TRUE, verbose=TRUE, rseed=44) {
  ###
  ###



  nvalid <- length(validrows)

  if (length(coords)==2) coords <-  unique(data[, coords]) 
  
  sn <- nrow(coords)
  n <- nrow(data)
  tn <- n/sn
  a <- abs(tn - floor(tn))
  if (a>0) stop("Unequal number of time points: check numbers of locations and times")
  nT <- sn*tn


  alldistmat <- dist_mat(coords, coordtype)
  max.d <- max(alldistmat)
  priors <- spTDyn::priors(inv.var.prior=Gamm(prior.sigma2[1], prior.sigma2[2]),
                                beta.prior=Norm(prior.beta0, prior.M^(-1)), rho.prior =Norm(0,10^10))

   k <- length(prior.phi)
    if (k>1) stop("Too many prior distributions for phi")
    if (k==0) prior.phi <- "Gamm"

     if (k==1) {
        u <- match(prior.phi, c("Unif", "Gamm", "FIXED"))
        if (is.na(u)) stop("Sorry, can't handle that prior distribution for phi.\n
        Please specify it as one of the three: Unif, Gamm or FIXED")
    }
    k <- length(prior.phi.param)
    if (k>2) stop("Too many prior hyper parameters for phi")
    if (prior.phi=="Unif") { 
        if (k<2) prior.phi.param <- 3 * c(1/max.d, 100/max.d)
        if (length(phi.npoints)==0) phi.npoints <- 10
     }
     if (prior.phi=="Gamm") { 
         if (k<2) prior.phi.param <- c(2, 1)
         if (length(phi.tuning)==0) phi.tuning <- mean(prior.phi.param) * 0.02
         spatial.decay <- spTimer::spT.decay(distribution=Gamm(prior.phi.param[1], prior.phi.param[2]), tuning=phi.tuning)
     }
    if (prior.phi == "FIXED") {
        if (k==1) phi <- prior.phi.param[1]
        else  phi <- 3/max.d 
        spatial.decay <- spTimer::spT.decay(distribution="FIXED", value= phi)
    }   
    if (prior.phi == "Unif") {
        if (k<2) prior.phi.param <- 3 * c(1, 100)/max.d 
        if (length(phi.npoints)==0) phi.npoints <- 10
        spatial.decay <- spTimer::spT.decay(distribution=Unif(prior.phi.param[1], prior.phi.param[2]),
                                             npoints=phi.npoints)
   }
  
  initials <- spTDyn::initials(rhotp=0, rho=0, sig2eps=0.01, sig2eta=0.1, sig2beta=0.1, sig2delta=0.1, phi=mean(prior.phi.param))

  distance.method <- "euclidean"
  if (coordtype=="utm") { 
  coords <- coords/1000
  } 
  if (coordtype=="lonlat") distance.method <- "geodetic:km"
  
 
  u <- getXy(formula=formula, data=data)
  X <- u$X
  y <- u$y
  vnames <- all.vars(formula)
  xnames <- colnames(X)
  yname <- vnames[1]
  p <- ncol(X)
  ynavec <- y

  if (nvalid >0) {
    ynavec[validrows] <- NA
    val_flag <- rep(0, nT)
    val_flag[validrows] <- 1
    vdaty <- y[val_flag>0]
    vdat <- data[val_flag>0, ]
  } 

  data$ynavec <- ynavec
 

  newformula <- update(formula, ynavec ~ .)

  fit <- spTDyn::GibbsDyn(formula=newformula, data=data,
                               model=model, coords=coords, time.data=time.data, 
                               nItr =N, nBurn=burn.in, distance.method=distance.method,
                               priors=priors, spatial.decay=spatial.decay,
                               initials = initials, cov.fnc =cov.fnc, 
                               truncation.para = truncation.para, 
                               scale.transform=scale.transform, report=n.report)

  allres <- list(params=fit$parameter[,-2], fit=fit, max.d=max.d)
  fits <- fit$fitted
  allres$fitteds <- fits[,1] 
  if (verbose) print(round(allres$params, 3))

  if (mchoice) {
    if (verbose) message("Calculating model choice statistics\n")
    pmcc_results <- list(gof=fit$PMCC[1], penalty=fit$PMCC[2], pmcc=fit$PMCC[3])
    
    sptimermod <- c(unlist(pmcc_results))
    allres$mchoice <-  sptimermod
    if (verbose) print(round(allres$mchoice, 2))
  }
  
  if (nvalid>0) {
    
    itmax <- fit$iterations - fit$nBurn
    if (verbose) message("validating ", nvalid, " space time observations", "\n")
    
    a <- matrix(rnorm(nvalid*itmax), nrow=nvalid, ncol=itmax)
    if  (model=="truncated") { 
      ovalues <- fit$op[val_flag>0, ]
      truncpara <- fit$truncation.para
      at <- truncpara$at 
      lambda <- truncpara$lambda
      ypreds <- reverse.truncated.fnc(ovalues,at=at[1],lambda=lambda,at2=at[2])
    } else { 
      ovalues <- fit$op
      sig2eps <-  fit$sig2ep
      meanmat <- ovalues[val_flag>0, ]
      dim(meanmat)
      sige <- sqrt(sig2eps)
      # a <- 1:3
      # matrix(rep(a, each=4), byrow=F, ncol=3)
      sigemat <- matrix(rep(sige, each=nvalid), byrow=F, ncol=itmax)
      ypreds <- meanmat + a * sigemat
      ##
    }
    ####
    #ovalues <- fit$op
    #sig2eps <-  fit$sig2ep
    #ovalidation <- ovalues[val_flag>0, ]
    #itmax <- ncol(ovalidation)
    #nvalidrows <- nrow(ovalidation)
    #dim(ovalidation)
    #sige <- sqrt(sig2eps)
    # a <- 1:3
    # matrix(rep(a, each=4), byrow=F, ncol=3)
    #sigemat <- matrix(rep(sige, each=nvalidrows), byrow=F, ncol=itmax)
    #a <- matrix(rnorm(nvalidrows*itmax), nrow=nvalidrows, ncol=itmax)
    #ypreds <- ovalidation + a * sigemat
    ##

    if (scale.transform == "SQRT")  ypreds <-  (ypreds)^2
    if (scale.transform == "LOG")  ypreds <-  exp(ypreds)

    predsums <- get_validation_summaries(t(ypreds))
    b <- calculate_validation_statistics(vdaty, ypreds)

    ##
    yvalidrows <- data.frame(vdat, predsums)
    allres$stats  <- b$stats
    allres$yobs_preds <- yvalidrows
    allres$valpreds <- t(ypreds)
    
    allvplots <- obs_v_pred_plot(vdaty, predsums)
    allres$validationplots <- allvplots
    if (plotit)  plot(allvplots$pwithseg)
    
    # if (plotit)  obs_v_pred_plot(vdaty, predsums)
    if (verbose) print(round(unlist(allres$stats), 3))
  }
 
  allres$prior.phi.param <- prior.phi.param
  allres$prior.phi <- prior.phi
  allres

}


##

BspTimer_sptime <- function(data=nysptime, formula=y8hrmax~xmaxtemp+xwdsp+xrh, model="GP",
                            coordtype="utm", coords=4:5,
                            validrows=NULL,   scale.transform ="SQRT",
                            prior.beta0=0, prior.M=0.001, prior.sigma2 =c(2, 1),
                            prior.phi="Gamm",
                            prior.phi.param =NULL, 
                            phi.tuning=NULL, phi.npoints=NULL,
                            cov.fnc="exponential",  tol.dist=0.005, 
                            time.data = NULL, newcoords = NULL, newdata =NULL,
                            truncation.para = NULL, annual.aggrn = "NONE",
                            N=5000, burn.in=1000, plotit=TRUE,
                            mchoice=TRUE, verbose=TRUE, rseed=44,
                            g_size = NULL, knots.coords = NULL, n.report=10) {
  ###
  ###

  
  nvalid <- length(validrows)
  
  if (length(coords)==2) coords <-  as.matrix(unique(data[, coords]))
  
  sn <- nrow(coords)
  n <- nrow(data)
  tn <- n/sn
  a <- abs(tn - floor(tn))
  if (a>0) stop("Unequal number of time points: check numbers of locations and times")
  nT <- sn*tn
  
  alldistmat <- dist_mat(coords, coordtype) # distances in kilometers
  max.d <- max(alldistmat)
  
  
  priors <- spTimer::spT.priors(model=model, inv.var.prior=Gamm(prior.sigma2[1], prior.sigma2[2]),
                                beta.prior=Norm(prior.beta0, prior.M^(-1)))
  
  k <- length(prior.phi)
  if (k>1) stop("Too many prior distributions for phi")
  if (k==0) prior.phi <- "Gamm"
  
  if (k==1) {
    u <- match(prior.phi, c("Unif", "Gamm", "FIXED"))
    if (is.na(u)) stop("Sorry, can't handle that prior distribution for phi.\n
        Please specify it as one of the three: Unif, Gamm or FIXED")
  }
  k <- length(prior.phi.param)
  if (k>2) stop("Too many prior hyper parameters for phi")
  if (prior.phi=="Unif") { 
    if (k<2) prior.phi.param <- 3 * c(1/max.d, 100/max.d)
    if (length(phi.npoints)==0) phi.npoints <- 10
  }
  if (prior.phi=="Gamm") { 
    if (k<2) prior.phi.param <- c(2, 1)
    if (length(phi.tuning)==0) phi.tuning <- mean(prior.phi.param) * 0.02
    spatial.decay <- spTimer::spT.decay(distribution=Gamm(prior.phi.param[1], prior.phi.param[2]), tuning=phi.tuning)
  }
  if (prior.phi == "FIXED") {
    if (k==1) phi <- prior.phi.param[1]
    else  phi <- 3/max.d 
    spatial.decay <- spTimer::spT.decay(distribution="FIXED", value= phi)
  }   
  if (prior.phi == "Unif") {
    if (k<2) prior.phi.param <- 3 * c(1, 100)/max.d 
    if (length(phi.npoints)==0) phi.npoints <- 5
    spatial.decay <- spTimer::spT.decay(distribution=Unif(prior.phi.param[1], prior.phi.param[2]),
                                        npoints=phi.npoints)
  }
  
  u <- getXy(formula=formula, data=data)
  X <- u$X
  y <- u$y
  vnames <- all.vars(formula)
  xnames <- colnames(X)
  
  p <- ncol(X)
  ynavec <- y
  
  ## Additional code for GPP
  ## Checks and sets up the knots.coords
  if ( (model == "GPP") || (model=="truncatedGPP") )  {
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
      xcoord <- c(max(coords[,1]), min(coords[,1]))
      ycoord <- c(max(coords[,2]), min(coords[,2]))
      knots.coords <- spTimer::spT.grid.coords(xcoord, ycoord, by=g_size)
    }
    knots.coords <- as.matrix(knots.coords)
    
    if (coordtype=="utm") knots.coords <-  knots.coords/1000
  }
  
  distance.method <- "euclidean"
  if (coordtype=="utm") { 
    coords <- coords/1000
  } 
  if (coordtype=="lonlat") distance.method <- "geodetic:km"
  
  if (nvalid >0) {
    ynavec[validrows] <- NA
    val_flag <- rep(0, nT)
    val_flag[validrows] <- 1
    vdaty <- y[val_flag>0]
    vdat <- data[val_flag>0, ]
  }
  
  data$ynavec <- ynavec
  newformula <- update(formula, ynavec ~ . )
  
  if ( (model=="GPP") || (model=="truncatedGPP") ) { 
    gp_fit <- spTimer::spT.Gibbs(formula=newformula, data=data,
                                 model=model, coords=coords,
                                 nItr =N, nBurn=burn.in, distance.method=distance.method,
                                 priors=priors, spatial.decay=spatial.decay,
                                 scale.transform=scale.transform, knots.coords=knots.coords, 
                                 cov.fnc=cov.fnc,  tol.dist = tol.dist, 
                                 time.data = time.data, newcoords = newcoords, newdata =newdata,
                                 truncation.para = truncation.para, annual.aggrn = annual.aggrn,
                                 report=n.report)
  } else { 
    gp_fit <- spTimer::spT.Gibbs(formula=newformula, data=data,
                                 model=model, coords=coords,
                                 nItr =N, nBurn=burn.in, distance.method=distance.method,
                                 priors=priors, spatial.decay=spatial.decay, tol.dist = tol.dist, 
                                 scale.transform=scale.transform,   cov.fnc=cov.fnc, 
                                 time.data = time.data, newcoords = newcoords, newdata =newdata,
                                 truncation.para = truncation.para, annual.aggrn = annual.aggrn,
                                 report=n.report)
  }
  
  allres <- list(params=gp_fit$parameter[,-2], fit=gp_fit, max.d=max.d)
  if ( (model=="GPP") || (model=="truncatedGPP") ) allres$knots.coords <- knots.coords
  if (verbose) print(round(allres$params, 3))
  
  
  if  (model=="truncatedGP") { ## Fitteds at the transformed scale
    ovalues <- gp_fit$op
    fits <- get_parameter_estimates(t(ovalues))
  } else { 
    fits <- gp_fit$fitted 
  }
  allres$fitteds <- fits[,1] 
  
  if (mchoice) {
    if (verbose) message("Calculating model choice statistics\n")
    pmcc_results <- list(gof=gp_fit$PMCC[1], penalty=gp_fit$PMCC[2], pmcc=gp_fit$PMCC[3])
    
    if (model=="GP")  {
      v <- logliks_from_full_gp_spTimer(gpfit=gp_fit)
      allres$logliks <- v
      waic_results <- calculate_waic(v$loglik)
      dic_results <- calculate_dic(v$log_full_like_at_thetahat, v$log_full_like_vec)
      sptimermod <- c(unlist(dic_results), unlist(waic_results), unlist(pmcc_results))
    } else sptimermod <- c(unlist(pmcc_results))
    allres$mchoice <-  sptimermod
    if (verbose) print(round(allres$mchoice, 2))
  }
  
  
  if (nvalid>0) {
 
    itmax <- gp_fit$iterations-gp_fit$nBurn
    if (verbose) message("validating ", nvalid, " space time observations", "\n")
    
    
    a <- matrix(rnorm(nvalid*itmax), nrow=nvalid, ncol=itmax)
    if  (model=="truncatedGP") { 
      ovalues <- gp_fit$op[val_flag>0, ]
      truncpara <- gp_fit$truncation.para
      at <- truncpara$at 
      lambda <- truncpara$lambda
      ypreds <- reverse.truncated.fnc(ovalues,at=at[1],lambda=lambda,at2=at[2])
    } else if (model=="GPP") {
      ## Generating the ypreds by approximation from the fitteds
      v <- fitted(gp_fit)[val_flag>0,]
      meanmat <- matrix(rep(v$Mean, each=itmax), byrow=TRUE, ncol=itmax)
      sigemat <- matrix(rep(v$SD, each=itmax), byrow=TRUE, ncol=itmax) 
      ypreds <- meanmat + a * sigemat
    } else {
      ovalues <- gp_fit$op
      sig2eps <-  gp_fit$sig2ep
      meanmat <- ovalues[val_flag>0, ]
      dim(meanmat)
      sige <- sqrt(sig2eps)
      # a <- 1:3
      # matrix(rep(a, each=4), byrow=F, ncol=3)
      sigemat <- matrix(rep(sige, each=nvalid), byrow=F, ncol=itmax)
      ypreds <- meanmat + a * sigemat
      ##
    }
    
    if (scale.transform == "SQRT")  ypreds <-  (ypreds)^2
    if (scale.transform == "LOG")  ypreds <-  exp(ypreds)
    
    predsums <- get_validation_summaries(t(ypreds))
    b <- calculate_validation_statistics(vdaty, ypreds)
    ##
    yvalidrows <- data.frame(vdat, predsums)
    allres$stats  <- b$stats
    allres$yobs_preds <- yvalidrows
    allres$valpreds <- t(ypreds)
    allvplots <- obs_v_pred_plot(vdaty, predsums)
    allres$validationplots <- allvplots
    if (plotit)  plot(allvplots$pwithseg)
    
    # if (plotit)  obs_v_pred_plot(vdaty, predsums)
    if (verbose) print(round(unlist(allres$stats), 3))
  }
 
  allres$prior.phi.param <- prior.phi.param
  allres$prior.phi <- prior.phi
  allres
  
}

