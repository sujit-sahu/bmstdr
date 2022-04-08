#' Model choice criteria calculation for univariate 
#' normal model for both known and unknown sigma^2
#' @param case One of the three cases: 
#' \itemize{  
#' \item{"Exact.sigma2.known" }{Use exact theoretical calculation.}  
#' \item{"MC.sigma2.known" }{Use Monte Carlo methods for drawing samples from the 
#' posterior assuming known sigma2.}  
#' \item{"MC.sigma2.unknown" }{Use the Gibbs sampler to generate samples 
#' from the joint posterior distribution of theta and sigma^2.}
#' }
#' @param y A vector of data values. Default is 28 ydata values from the package bmstdr
#' @param mu0 The value of the prior mean if kprior=0. Default is the data mean.
#' @param sigma2 Value of the known data variance; defaults to sample variance of the data. This is ignored 
#' in the third case when sigma2 is assumed to be unknown.
#' @param kprior A scalar providing how many data standard deviation the prior
#' mean is from the data mean. Default value is 0.
#' @param prior.M Prior sample size, defaults to 10^(-4).
#' @param prior.sigma2 Shape and scale parameter value for the gamma prior on 1/sigma^2, the precision.
#' @param N The number of samples to generate. 
#' @param rseed The random number seed. Defaults to 44 to reproduce the results
#'  in the book  \insertCite{Sahubook;textual}{bmstdr}.  
#' @return A list containing the exact values of pdic, dic, pdicalt, dicalt,
#'  pwaic1, waic1, pwaic2, waic2, gof, penalty and pmcc.
#'  Also prints out the posterior mean and variance.
#'  @references
#' \insertAllCited{}
#' @examples
#' Bmchoice()
#' b1 <- Bmchoice(case="Exact.sigma2.known")
#' b2 <- Bmchoice(case="MC.sigma2.known")
#' d1 <- Bmchoice(case="MC.sigma2.unknown")
#' d2 <- Bmchoice(y=rt(100, df=8),  kprior=1, prior.M=1)
#' 
#' @export
Bmchoice <- function(case="Exact.sigma2.known", y=ydata, mu0=mean(y), sigma2=22, 
                    kprior=1, prior.M=1, prior.sigma2=c(2, 1),  N=10000, rseed=44){
  if (!is.vector(y)) {
    stop("The argument must be a vector\n")
  }
  set.seed(rseed)
  implemented <- c("Exact.sigma2.known",  "MC.sigma2.known", "MC.sigma2.unknown")
  a <- grepl(case,  x=implemented, ignore.case = TRUE)
  if (any(a)) { 
    case <- implemented[which(a)]
  } else { stop("Wrong case. Please see helpfile")}
  
  if (case=="Exact.sigma2.known") { 
    results <- exact_normal_theta_model_choice_values(sigma2=sigma2,  kprior=kprior, prior.M=prior.M)
    message("Results by using  exact theoretical calculations assuming known sigma2.\n")
  } else if (case=="MC.sigma2.known") { 
    results <- normal_theta_unknown_mc_values(sigma2=sigma2,  kprior=kprior, prior.M=prior.M, N=N)
    message("Results by  using Monte Carlo Samples from the posterior distribution of theta for known sigma2.\n")
  } else if (case=="MC.sigma2.unknown") { 
    results <- normal_both_unknown_mc_values(kprior=kprior, prior.M=prior.M, prior.sigma2=prior.sigma2, N=N)
    message("Results by  using Monte Carlo Samples from the joint posterior distribution of theta and sigma2.\n")
  } else { stop("The supplied case has not been implemented.")}
  results
}


# Model choice criteria calculation
## For known sigma^2
## Calculates the exact values of the different model choice criteria
##
## @param y A vector of data values. Default is 28 ydata values from the package bmstdr
## @param mu0 The value of the prior mean if k=0. Default is the data mean.
## @param kprior A scalar providing how many data standard deviation the prior
## mean is from the data mean. Default value is 0.
## @param prior.M Prior sample size, defaults to 10^(-4).
## @param sigma2 Value of the known data variance; defaults to sample variance of the data
## @return A list containing the exact values of aic, bic, p_dic, dic, p_dic alt, dic alt,
##  p_waic1, waic1, p_waic2, waic2, gof, penalty and pmcc.
##  Also prints out the posterior mean and variance.
## @examples
## exact_normal_theta_model_choice_values()
## a <- exact_normal_theta_model_choice_values(sigma2=22,  kprior=1, prior.M=1)
## exact_normal_theta_model_choice_values(rt(100, df=8),  kprior=1, prior.M=1)
# ## @export
exact_normal_theta_model_choice_values <- function(y=ydata, mu0=mean(y), kprior=1, prior.M=10^(-4), sigma2=var(y)) {
  if (!is.vector(y)) {
    stop("The argument must be a vector\n")
    }
y <- y[!is.na(y)] ## Removes the NA's in y
ybar <- mean(y)
s2y <- var(y)
n <- length(y)
sumsqsy <- (n-1) * s2y
# Prior mean
mu <- mu0 +  kprior * sqrt(sigma2/n)
# Prior variance
tau2 <- sigma2/prior.M
# Posterior variance
sigma2p <- 1/(n/sigma2 + 1/tau2)
# Posterior mean
mup <- sigma2p * (n*mean(y)/sigma2 + mu/tau2)
print(list(posterior.mean=mup, posterior.var=sigma2p, data.mean=mean(y),
           sigma2byn=sigma2/n))
aic <- n * log(2 * pi * sigma2) + sumsqsy/sigma2 + 2 ## AIC
bic <- n * log(2 * pi * sigma2) + sumsqsy/sigma2 + log(n) ## BIC

pdic <- n *  sigma2p/sigma2
dic <-  n * log(2 * pi * sigma2) + sum((y-mup)^2)/sigma2 + 2 * pdic
# list(pdic=pdic, dic=dic)

pdicalt <- n^2 * sigma2p * (sigma2p + 2 * (mup-ybar)^2)/sigma2^2
dicalt <-  n * log(2 * pi * sigma2) + sum((y-mup)^2)/sigma2 + 2 * pdicalt
# list(pdicalt=pdicalt, dicalt=dicalt)

# dic_exact <- list(pdic=pdic, pdicalt=pdicalt, dic1=dic, dicalt=dicalt)

sigma2star <-  sigma2 + sigma2p

pwaic1 <-  n * log(sigma2/sigma2star) + n * sigma2p/sigma2   + sigma2p * sum((y-mup)^2)/(sigma2 * sigma2star) # 0.9167
pwaic2 <- n * sigma2p^2 / sigma2^2 + sigma2p/sigma2^2 * sum((y-mup)^2) ##  0.9647962

waic_first_term <- n * log(2 * pi * sigma2star) + sum((y-mup)^2)/sigma2star # 165.0721

waic1 <- waic_first_term + 2 * pwaic1 # 166.9056
waic2 <- waic_first_term + 2 * pwaic2 # 167.0017
waic_exact <- list(pwaic1=pwaic1, pwaic2=pwaic2, waic1=waic1, waic2=waic2)

# pmcc

gof <- sum((y-mup)^2)
penalty <- n * sigma2star

pmcc <- gof + penalty
pmcc_exact <- list(gof=gof, penalty=penalty, pmcc=pmcc)
print(list(aic=aic, bic=bic))
list(pdic=pdic,  pdicalt=pdicalt, dic=dic, dicalt=dicalt, pwaic1=pwaic1, pwaic2=pwaic2, waic1=waic1, waic2=waic2,
     gof=gof, penalty=penalty, pmcc=pmcc)

}
## Model choice criteria calculation using sampling
## For known sigma^2
## Calculates approximate values of the different model choice criteria
##
## @param y A vector of data values. Default is 28 ydata values from the package bmstdr
## @param mu0 The value of the prior mean if k=0. Default is the data mean.
## @param kprior A scalar providing how many data standard deviation the prior
## mean is from the data mean. Default value is 0.
## @param prior.M Prior sample size, defaults to 10^(-4).
## @param sigma2 Value of the known data variance; defaults to sample variance of the data
## @param N The number of samples to generate
## @return A list containing the exact values of aic, bic, p_dic, dic, p_dic alt, dic alt,
##  p_waic1, waic1, p_waic2, waic2, gof, penalty and pmcc.
##  Also prints out the posterior mean and variance.
## @examples
## a2 <- normal_theta_unknown_mc_values(sigma2=22,  kprior=1, prior.M=1)
# ## @export
normal_theta_unknown_mc_values <- function(y=ydata, mu0=mean(y), kprior=1, prior.M=10^(-4), sigma2=var(y), N=10000) {
  if (!is.vector(y)) {
    stop("The argument must be a vector\n")
  }
  y <- y[!is.na(y)] ## Removes the NA's in y
  ybar <- mean(y)
  s2y <- var(y)
  n <- length(y)
  sumsqsy <- (n-1) * s2y
  # Prior mean
  mu <- mu0 +  kprior * sqrt(sigma2/n)
  # Prior variance
  tau2 <- sigma2/prior.M
  # Posterior variance
  sigma2p <- 1/(n/sigma2 + 1/tau2)
  # Posterior mean
  mup <- sigma2p * (n*mean(y)/sigma2 + mu/tau2)
  print(list(posterior.mean=mup, posterior.var=sigma2p, data.mean=mean(y),
             sigma2byn=sigma2/n))

  ### Generate theta from the posterior distribution
  theta <- rnorm(n=N, mean=mup, sd=sqrt(sigma2p))  ## This is the exact sampling case

  ## dic calculation for normal with known variance but using sampling
  samps <- cbind(theta, rep(sigma2, N))
  a2 <- normal_theta_sigma2_model_choice_values(samps, y=y)
  a2
}
## Model choice criteria calculation using sampling
## when both parameters are unknown
## Calculates approximate values of the different model choice criteria
##
## @param y A vector of data values. Default is 28 ydata values from the package bmstdr
## @param mu0 The value of the prior mean if k=0. Default is the data mean.
## @param kprior A scalar providing how many data standard deviation the prior
## mean is from the data mean. Default value is 0.
## @param prior.M Prior sample size, defaults to 10^(-4).
## @param prior.sigma2 Shape and scale parameter value for the gamma prior on 1/sigma^2, the precision.
## @param sigma2 Value of the known data variance; defaults to sample variance of the data
## @param N The number of samples to generate
## @return A list containing the exact values of p_dic, p_dic alt, dic, dic alt,
##  p_waic1, waic1, p_waic2, waic2, gof, penalty and pmcc.
## @examples
## a3 <- normal_both_unknown_mc_values(kprior=1, prior.M=1)
# ## @export
normal_both_unknown_mc_values <- function(y=ydata, mu0=mean(y), kprior=1, prior.M=10^(-4), 
                                          prior.sigma2=c(2,1), N=10000) {
  if (!is.vector(y)) {
    stop("The argument must be a vector\n")
  }
  y <- y[!is.na(y)] ## Removes the NA's in y
  samps <- normal_theta_sigma2_gibbs(y=y, mu0=mu0, kprior=kprior, prior.M=prior.M, prior.sigma2=prior.sigma2, N=N)
  a3 <- normal_theta_sigma2_model_choice_values(samps, y=y)
  a3
}

## Model choice criteria calculation using sampling for both theta and sigma2
## Calculates estimated values of the different model choice criteria
##
## @param samps A N by 2 matrix of sampled theta and sigma2 values
## @param y A vector of data values
normal_theta_sigma2_model_choice_values <- function(samps, y=ydata) {
  thetahat <- mean(samps[,1])
  sigma2hat <- mean(samps[,2])
  pars <- c(thetahat, sigma2hat)

  log_full_likelihood_N_theta <- function(pars, y=y) { ## log likelihood value for a given theta sigma2 as first and second components of pars
    logden <- unlist(lapply(y, dnorm, mean=pars[1], sd=sqrt(pars[2]), log=TRUE))
    sum(logden)
  }
  loglikatthetahat <- log_full_likelihood_N_theta(pars, y=y)  ## just check everything works
  log_liks <- apply(samps, 1, log_full_likelihood_N_theta, y=y)
  dic <- calculate_dic(loglikatthetahat, log_liks)

  ## WAIC calculations

  likelihoods_N_theta <- function(pars, y=y) { ## likelihood value for a given theta sigma2 as first and second components of pars
    unlist(lapply(y, dnorm, mean=pars[1], sd=sqrt(pars[2]), log=TRUE))
  }
  tv <- apply(samps, 1, likelihoods_N_theta, y=y) ## n by N
  v <- t(tv)
  waic <- calculate_waic(v)

  ## PMCC calculation for known sigma2 but sampling the thetas

  n <- length(y)
  pred_samples_N_theta <- function(pars, n) { ## Draw n samples from the posterior predictive distribution
    rnorm(n, mean=pars[1], sd=sqrt(pars[2]))
  }
  v <- apply(samps, 1, pred_samples_N_theta, n=n) ## n by N

  expected_ys <- apply(v, 1, mean)
  var_ys <- apply(v, 1, var)

  gof <- sum( (y-expected_ys)^2)
  penalty <- sum(var_ys)

  pmcc <- gof + penalty

  list(pdic=dic$pdic,  pdicalt=dic$pdicalt, dic=dic$dicorig,
       dicalt=dic$dicalt, pwaic1=waic$pwaic1, pwaic2=waic$pwaic2, waic1=waic$waic1, waic2=waic$waic2,
       gof=gof, penalty=penalty, pmcc=pmcc)
}





