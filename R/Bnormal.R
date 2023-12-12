#' N(theta, sigma2):  Using different methods. 
#' @param y A vector of data values. Default is 28 ydata values from the package bmstdr
#' @param mu0 The value of the prior mean if kprior=0. Default is the data mean.
#' @param kprior A scalar providing how many data standard deviation the prior
#' mean is from the data mean. Default value is 0.
#' @param prior.M Prior sample size, defaults to 10^(-4).#' 
#' @param prior.sigma2 Shape and scale parameter value for the gamma prior on 1/sigma^2, the precision.
#' @param package Which package (or method) to use. Possibilities are: 
#' \itemize{  
#' \item "exact": Use exact theoretical calculation.  
#' \item "RGibbs": Use Gibbs sampler using R code.  
#' \item "stan": Use HMC by implementing in Stan.
#' \item "inla": Use the INLA package.
#' }
#' @param N is the number of Gibbs sampling iterations
#' @param burn.in is the number of initial iterations to discard before 
#' making inference. 
#' @param rseed is the random number seed defaults to 44. 
#' @return A list containing the exact posterior means and variances of theta and sigma2
#' @example inst/examples/bnormal_examples.R
#' @export
Bnormal <- function(package="exact", y=ydata, mu0=mean(y), kprior=0, prior.M=0.0001, 
                     prior.sigma2=c(0, 0),  N=2000, burn.in=1000, rseed=44){
  if (!is.vector(y)) {
    stop("The argument must be a vector\n")
  }
  set.seed(rseed)
  implemented <- c("exact",  "RGibbs", "stan", "inla")
  a <- grepl(package,  x=implemented, ignore.case = TRUE)
  if (any(a)) { 
    package <- implemented[which(a)]
  } else { stop("Wrong package or model. Please see helpfile")}
 
  if (package=="exact") { 
     results <- normal_theta_sigma2_exact(y=y, mu0=mean(y), kprior=kprior,prior.M=prior.M, prior.sigma2=prior.sigma2)
     message("Results from exact methods.\n")
  } else if (package=="RGibbs") { 
    results <- normal_theta_sigma2_gibbs(y=y, mu0=mean(y), kprior=kprior,  prior.M=prior.M, prior.sigma2=prior.sigma2, N=N)   
    message("Results from Gibbs sampler coded in R.\n")
  } else if (package=="stan") { 
    results <- normal_theta_sigma2_by_stan(y, mu0=mean(y), kprior=kprior, prior.M=prior.M, prior.sigma2=prior.sigma2, N=N, burn.in=burn.in) 
    message("Results from Hamiltonian Monte Carlo in Stan.\n")
  } else if (package=="inla") {   
    if (inlabru::bru_safe_inla()) {
    results <- normal_theta_sigma2_by_inla(y, mu0=mean(y), kprior=kprior, prior.M=prior.M, prior.sigma2=prior.sigma2, N=N)  
    message("Results from INLA.\n")
    } else {
      stop("The chosen package INLA is not available.")
    }
  } else { stop("Either the method or package has not been implemented.")}
  results
}
## N(theta, sigma^2) example computed using several methods
## Assume the hierarchical prior distributions:
## pi(theta | sigma^2) and pi(lambda^2) where lambda^2=1/sigma^2
## Use the 28 air pollution values from New York to illustrate

## #' N(theta, sigma^2):  Exact estimates of theta and sigma2
## #'
## #' @inheritParams Bnormal
## # #' @export
normal_theta_sigma2_exact <- function(y=ydata, mu0=mean(y), kprior=0, prior.M=0.0001, prior.sigma2=c(0,0)){
  if (!is.vector(y)) {
    stop("The argument must be a vector\n")
  }
  y <- y[!is.na(y)] ## Removes the NA's in y
  n <- length(y)
  ybar <- mean(y)
  s2y <- var(y)
  mu <- mu0 + kprior * sqrt(s2y/n)
  mup <- (n*ybar + prior.M * mu)/(n+prior.M)
  astar <- prior.sigma2[1] + 0.5 * n
  bstar <- prior.sigma2[2] + 0.5 * (n-1) * s2y + 0.5 * n * prior.M *(ybar - mu)^2/(n+prior.M)
  var_theta <- 2 * bstar /((n+2*prior.sigma2[1] -2) * (n+prior.M))
  mean_sigma2 <- bstar / (astar -1)
  var_sigma2 <- (mean_sigma2)^2 /(astar-2)

  list(mean_theta = mup, var_theta = var_theta, mean_sigma2=mean_sigma2, var_sigma2=var_sigma2)
}

## #' N(theta, sigma^2):  Gibbs sampling for estimation of theta and sigma2
## #'
##@inheritParams Bnormal
##@return Prints a list containing the estimated posterior means and variances of theta, sigma2 and
##the co-efficient of variation sigma/theta. A N by 2 matrix containing the sampled values are returned.
##@examples
##psamps <- normal_theta_sigma2_gibbs(y=ydata, kprior=1, prior.M=1, prior.sigma2=c(2,1))
##psamps <- normal_theta_sigma2_gibbs(mu0=0) ## non-informative prior
##psamps <- normal_theta_sigma2_gibbs()
##psamps <- normal_theta_sigma2_gibbs(kprior=10, prior.M=1, prior.sigma2=c(2, 1))
# ##@export
normal_theta_sigma2_gibbs <- function(y=ydata, mu0=mean(y), kprior=0, prior.M=0.0001, prior.sigma2=c(0, 0), N=10000){
  if (!is.vector(y)) {
    stop("The argument must be a vector.\n")
  }
  y <- y[!is.na(y)] ## Removes the NA's in y
  n <- length(y)
  ybar <- mean(y)
  s2y <- var(y)
  mu <- mu0 + kprior * sqrt(s2y/n)
  mup <- (n*ybar + prior.M * mu)/(n+prior.M)

  ## Setup initial values
  theta <- numeric() ## to hold the sampled values of theta
  sigma2 <- numeric() ## to hold the sampled value of precision =1/sigma^2

  theta[1] <- ybar ## Any reasonable starting value will do
  sigma2[1] <- 1  ## Any reasonable starting value will do
  # set.seed(44) ## Set the random number seed to reproduce results
  for (j in 1:(N-1)) {
    theta[j+1] <- rnorm(n=1, mean=mup, sd=sqrt(sigma2[j]/(n+prior.M)))
    u <-  prior.sigma2[2] + 0.5  * ( sum((y-theta[j+1])^2) + prior.M * (theta[j+1] - mu)^2)  ## Keep this
    sigma2[j+1] <- 1.0/rgamma(n=1, shape=0.5*(n+1) + prior.sigma2[1], rate=u)
  }
  print(list(mean_theta = mean(theta), var_theta_estimate = var(theta),  mean_sigma2=mean(sigma2), var_sigma2_estimate=var(sigma2)))
  # Return the samples as a matrix
 cbind(theta, sigma2)
}

##N(theta, sigma^2):  Using Stan to estimate theta and sigma2
##@inheritParams Bnormal
##@return Returns a table containing the exact (first row) and the Stan estimated posterior means and variances of
##theta, sigma2 and the co-efficient of variation sigma/theta.
##@examples
##normal_theta_sigma2_by_stan(N=1000, kprior=1, prior.M=1, prior.sigma2=c(2,1))
# ##@export
normal_theta_sigma2_by_stan <- function(y=ydata, mu0=mean(y), kprior=1, prior.M=1, prior.sigma2=c(0, 0), N=1000, burn.in=500){
  if (!is.vector(y)) {
    stop("The argument must be a vector.\n")
  }
  y <- y[!is.na(y)] ## Removes the NA's in y
  n <- length(y)
  ybar <- mean(y)
  s2y <- var(y)
  mu <- mu0 + kprior * sqrt(s2y/n) # prior mean

 data_to_send_to_stan <- list(n=n, y=y, mu=mu, mprior=prior.M, aprior=prior.sigma2[1], bprior = prior.sigma2[2])

  init_fun <- function() {
    list(theta = mean(y), sigma2 = 1)
  }

  results_from_stan <- rstan::sampling(stanmodels$normal, data=data_to_send_to_stan, 
                                       seed = 44, chains =1, iter = N, 
                                       warmup = burn.in, init=init_fun)
  
  #results_from_stan <- stan(data=data_to_send_to_stan, file = "stanfiles/normal.stan", seed = 44,
   #                                chains =1, iter = N, warmup = Nburn, init=init_fun)

  # normal_stan_estimates <- summary(normal_results_from_stan,  pars =c("theta", "sigma2"), probs = c(.025, .975))
  # normal_stan_estimates

  sres <- rstan::extract(results_from_stan, pars=c("theta", "sigma2"), permuted=TRUE)
  samps <- cbind(sres$theta, sres$sigma2)

  eresults <- normal_theta_sigma2_exact(y=y, mu0=mu0, kprior=kprior, prior.M=prior.M, prior.sigma2=prior.sigma2)
  ## Run Gibbs sampling
  stanres <- list(mean_theta = mean(samps[,1]), var_theta = var(samps[,1]), mean_sigma2=mean(samps[,2]), var_sigma2=var(samps[,2]))
  stanlow <- list(theta_low=quantile(samps[,1], probs=0.025), var_theta=NA, sigma2_low =quantile(samps[,2], probs=0.025), var_sigma2=NA)
  stanupp <- list(theta_low=quantile(samps[,1], probs=0.975), var_theta=NA, sigma2_low =quantile(samps[,2], probs=0.975), var_sigma2=NA)
  a <- rbind(unlist(eresults), unlist(stanres), unlist(stanlow), unlist(stanupp))
  cvsamp <- sqrt(samps[,2])/samps[,1]
  cv <- c(NA, mean(cvsamp), quantile(cvsamp, probs=c(0.025, 0.975)))
  u <- cbind(a, cv)  
  rownames(u) <- c("Exact", "Estimate", "2.5%", "97.5%")
  u
}
## b <- normal_theta_sigma2_by_stan(N=1000, kprior=1, prior.M=1, prior.sigma2[1]=2, prior.sigma2[2]=1)
## xtable(b, digits = 3)

##N(theta, sigma^2):  Using INLA to estimate theta and sigma2
##It is not possible to implement this simple hierarchical prior model in INLA.
##Instead  an independent prior model for theta and sigma2 is implemented
##@inheritParams Bnormal
##@return Returns a table containing the exact (first row) and the Stan estimated posterior means and variances of
##theta, sigma2 and the co-efficient of variation sigma/theta.
##@examples
##normal_theta_sigma2_by_inla(N=1000, kprior=1, prior.M=1, prior.sigma2=c(2,1))
# ##@export
normal_theta_sigma2_by_inla <- function(y=ydata, mu0=mean(y), kprior=1, prior.M=1, prior.sigma2=c(2,1), N=1000){
  if (!is.vector(y)) {
    stop("The argument must be a vector.\n")
  }
  y <- y[!is.na(y)] ## Removes the NA's in y
  n <- length(y)
  ybar <- mean(y)
  s2y <- var(y)
  mu <- mu0 + kprior * sqrt(s2y/n) # prior mean

  # set.seed(44)
  ydf <- data.frame(y=y)
  formula <- y~1
  hyper <- list(prec=list(prior="loggamma",param=prior.sigma2))
  mod1 <- INLA::inla(formula, data=ydf,  family="gaussian", control.family = list(hyper=hyper),
             control.fixed = list(mean.intercept=mu, prec.intercept=s2y/prior.M))
  names(mod1)
  summary(mod1)
  prec.marg.samp <- INLA::inla.rmarginal(N, mod1$marginals.hyperpar[[1]])
  sigma2 <- 1/prec.marg.samp #variance
  theta <- INLA::inla.rmarginal(N, mod1$marginals.fixed[[1]])
  #cv
  cv <- sqrt(sigma2)/theta
  samps <- cbind(theta, sigma2)
  eresults <- normal_theta_sigma2_exact(y=y, mu0=mu0, kprior=kprior, prior.M=prior.M, prior.sigma2=prior.sigma2)
  ## Run Gibbs sampling
  stanres <- list(mean_theta = mean(samps[,1]), var_theta = var(samps[,1]), mean_sigma2=mean(samps[,2]), var_sigma2=var(samps[,2]))
  stanlow <- list(theta_low=quantile(samps[,1], probs=0.025), var_theta=NA, sigma2_low =quantile(samps[,2], probs=0.025), var_sigma2=NA)
  stanupp <- list(theta_low=quantile(samps[,1], probs=0.975), var_theta=NA, sigma2_low =quantile(samps[,2], probs=0.975), var_sigma2=NA)
  a <- rbind(unlist(eresults), unlist(stanres), unlist(stanlow), unlist(stanupp))
  cvsamp <- sqrt(samps[,2])/samps[,1]
  cv <- c(NA, mean(cvsamp), quantile(cvsamp, probs=c(0.025, 0.975)))
  u <- cbind(a, cv)  
  rownames(u) <- c("Exact", "Estimate", "2.5%", "97.5%")
  u
}



