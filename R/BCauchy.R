#' Cauchy prior simulation example. 
#' @param method Which method or package to use. Possibilities are: 
#' \itemize{  
#' \item{"exact" }{Use exact numerical integration.}  
#' \item{"importance" }{Use importance sampling with the prior distribution as the importance sampling distribution.}
#' \item{"rejection" }{Use rejection sampling with the prior distribution as the importance sampling distribution.}
#' \item{"independence" }{Use the Metropolis-Hastings independence sampler with the prior distribution as the proposal distribution.}
#' \item{"randomwalk" }{Use the Metropolis-Hastings random-walk sampler with normal distribution 
#' with mean 0 and variance (tuning.sd)^2 as the increment distribution.}
#' }
#' @param true.theta True value of theta with a default value of 5.
#' @param n Data sample size; defaults to 100. 
#' @param N is the number of Monte Carlo samples.  
#' @param rseed is the random number seed for drawing data samples.  
#' @param tuning.sd is the standard deviation of the proposal increment distribution for the random walk sampler.  
#' @return A list containing the estimated posterior mean, ybar (the data mean)
#' and the values of the numerator and the denominator integrals
#' The routine simulates n observations from N(theta, 1). 
#' Mean of the simulated data values are returned as ybar.
#' @examples
#' \donttest{
#' BCauchy(true.theta = 1, n=25) 
#' BCauchy(true.theta = 5, n=100) 
#' BCauchy(method="importance", true.theta = 1, n=25) 
#' BCauchy(method="importance", true.theta = 1, n=25, N=20000) 
#' BCauchy(method="rejection", true.theta = 1, n=25) 
#' BCauchy(method="independence", true.theta = 1, n=25) 
#' BCauchy(method="randomwalk", true.theta = 1, n=25, tuning.sd =1) 
#'}
#' @export
BCauchy <- function(method="exact", true.theta=1, n=25, N=10000, rseed=44, tuning.sd=1) {
  implemented <- c("exact", "importance", "rejection", "independence", "randomwalk")
  a <- grepl(method, x=implemented, ignore.case = TRUE)
  if (any(a)) { 
    method <- implemented[which(a)]
  } else { stop("Wrong method for the Cauchy prior example. Please see helpfile")}
  set.seed(rseed)
  if (method=="exact") { 
    results <- Cauchy_prior_exact_estimate(true.theta=true.theta, n=n, rseed=rseed)
    message("Results from exact integration.\n")
  } else if (method=="importance") { 
    results <- Cauchy_prior_IS_estimate(true.theta=true.theta, n=n,  N=N)
    message("Results from using importance sampling.\n")
  } else if (method=="rejection") { 
    results <- Cauchy_prior_rejection_estimate(true.theta=true.theta, n=n,  N=N)
    message("Results from using rejection sampling.\n")
  } else if (method=="independence") {   
    results <- Cauchy_prior_independence_sampler(true.theta=true.theta, n=n,  N=N) 
    message("Results from independence sampler.\n")
  }  else if (method=="randomwalk") {   
    results <- Cauchy_prior_random_walk_sampler(true.theta=true.theta, n=n,  N=N, tuning.sd=tuning.sd) 
    message("Results from the random walk sampler.\n")
  } else { stop("The supplied method has not been implemented.")}
  results
}
## Estimate the exact value of theta by numerical integration for the
## normal(theta, 1) data distribution and a Cauchy prior on theta
## @inheritParams BCauchy
## @return A list containing the estimated posterior mean, ybar (the data mean)
## and the values of the numerator and the denominator integrals
## The routine simulates n observations from N(theta, 1) where the default
## seed is 44. Simulated values are not returned.
## @family Cauchy_prior_IS_estimate Cauchy_prior_rejection_estimate
## @examples
## Cauchy_prior_exact_estimate(true.theta = 1, n=25) 
## Cauchy_prior_exact_estimate(true.theta = 5, n=100) 
## Cauchy_prior_exact_estimate(true.theta = -1, n=100)
## @export
Cauchy_prior_exact_estimate <- function(true.theta=5, n=100, rseed=44) {
set.seed(rseed)
y <- rnorm(n, mean=true.theta, sd=1)
ybar <- mean(y)
Cauchy_post_den <- function(theta, ybar=ybar, n=n) {
  u <- n * 0.5 * (theta-ybar)^2
  exp(-u) /(1+theta^2)
}

Cauchy_theta_star_postden<- function(theta, ybar=ybar, n=n) {
  u <- n * 0.5 * (theta-ybar)^2
  theta * exp(-u) /(1+theta^2)
}
I1 <- integrate(Cauchy_post_den, ybar=ybar, n=n, lower=-Inf, upper=Inf)
I2 <- integrate(Cauchy_theta_star_postden, ybar=ybar, n=n, lower=-Inf, upper=Inf)
# print(I2$value/I1$value)
post.mean <- as.numeric(I2$value)/as.numeric(I1$value)
list(exact.post.mean=post.mean, ybar=ybar, true.theta=true.theta, I2=I2$value, I1=I1$value)
}

## Illustrate results from  importance sampling for the
## normal(theta, 1) data distribution and a Cauchy prior on theta
## @inheritParams BCauchy
## @return A list containing the estimated posterior mean using importance sampling and the
## estimated value using numerical integration.
## Cauchy_prior_IS_estimate(true.theta = 1, n=25) \cr
## Cauchy_prior_IS_estimate(true.theta = 5, n=100) \cr
## Cauchy_prior_IS_estimate(true.theta = -1, n=100)
# ## @export
Cauchy_prior_IS_estimate <- function(true.theta=5, n=100, N=10000) {
 eresults <- Cauchy_prior_exact_estimate(true.theta =true.theta, n=n)
 ybar <- eresults$ybar
 ## Simulate N observations from the prior
 theta <- rcauchy(N)
 w <- exp(-n * 0.5 * (theta-ybar)^2) ## calculate the weights
 top <- mean(theta*w)  ## Estimate of the top integral
 bot <- mean(w)  ## Estimate of the bottom integral
 list(is_estimate =top/bot, exact_estimate=eresults$exact.post.mean, ybar=eresults$ybar)
}
## Illustrate results from  rejection sampling for the
## normal(theta, 1) data distribution and a Cauchy prior on theta
## @inheritParams BCauchy
## @return A list containing the estimated posterior mean using importance sampling and the
## estimated value using numerical integration.
##
## Cauchy_prior_rejection_estimate(true.theta = 1, n=25) \cr
## Cauchy_prior_rejection_estimate(true.theta = 5, n=100) \cr
## Cauchy_prior_rejection_estimate(true.theta = -1, n=100)
# ## @export
Cauchy_prior_rejection_estimate <- function(true.theta=5, n=100, N=10000) {
  eresults <- Cauchy_prior_exact_estimate(true.theta =true.theta, n=n)
  ybar <- eresults$ybar
  ## Generate samples from the importance density
  theta <- rcauchy(N)
  w <- exp(-n * 0.5 * (theta-ybar)^2)
  u <-  runif(N)
  test <- ifelse(u<w, 1, 0)
  retained <- theta[test==1]
  list(rejection_estimate =mean(retained), exact_estimate=eresults$exact.post.mean,
          acceptance_rate=100*length(retained)/N, 
       interval = quantile(retained, probs=c(0.025, 0.975)), ybar=eresults$ybar)
}
## Illustrate results from  the Metropolis-Hastings independence sampler
## normal(theta, 1) data distribution and a Cauchy prior on theta
## @inheritParams BCauchy
## @return A list containing the estimated posterior mean using importance sampling and the
## estimated value using numerical integration.
##
## Cauchy_prior_independence_sampler(true.theta =1, n=25) \cr
## Cauchy_prior_independence_sampler(true.theta = 5, n=100) \cr
# ## Cauchy_prior_independence_sampler(true.theta = -1, n=100)
# ## @export
Cauchy_prior_independence_sampler <- function(true.theta=5, n=100, N=10000) {
  eresults <- Cauchy_prior_exact_estimate(true.theta =true.theta, n=n)
  ybar <- eresults$ybar
  theta <- numeric()
  theta[1] <- ybar ## sensible starting value
  # Note that the likelihood function is density of the normal
  # distribution with mean ybar and variance 1/n.
  count <- 1 ## counts how many we have accepted

  for (j in 1:(N-1)) {
  # j indexes iteration
  # the code below assumes we are at iteration j and the value is theta[j]
  # we would like to perform one iteration to go to iteration j+1
  phi <- rcauchy(n=1) ## Draw one proposal

  top <- dnorm(phi, mean=ybar, sd=1/sqrt(n)) # Numerator at the proposal phi
  bot <- dnorm(theta[j], mean=ybar, sd=1/sqrt(n)) # Denominator at the current theta

  alpha <- min(1, top/bot) # Acceptance probability
  if (runif(n=1) < alpha) { # Test whether to accept
    theta[j+1] <- phi  ## Accepted the proposal
    count <- count +1
  }
  else theta[j+1] <- theta[j] ## Rejected, stay where we are
  }
 list(MH_ind_samp_estimate =mean(theta), exact_estimate=eresults$exact.post.mean,
       acceptance_rate=100*count/N, interval = quantile(theta, probs=c(0.025, 0.975)), ybar=ybar)
}

## Illustrate results from  the Metropolis-Hastings independence sampler
## normal(theta, 1) data distribution and a Cauchy prior on theta
## @inheritParams BCauchy
## @return A list containing the estimated posterior mean using importance sampling and the
## estimated value using numerical integration.
## Cauchy_prior_random_walk_sampler(true.theta = 1, n=25) \cr
## Cauchy_prior_random_walk_sampler(true.theta = 5, n=100) \cr
## Cauchy_prior_random_walk_sampler(true.theta = -1, n=100)
# ## @export
Cauchy_prior_random_walk_sampler <- function(true.theta=5, n=100, N=10000,  tuning.sd=0.5) {
  eresults <- Cauchy_prior_exact_estimate(true.theta =true.theta, n=n)
  ybar <- eresults$ybar

  theta <- numeric()
  theta[1] <- ybar ## reasonable starting value
  count <- 1  # keeps track of how many acceptances

  logpost <- function(theta) { ## evaluates the log posterior density
    u <- dnorm(theta, mean=ybar, sd=sqrt(1/n), log=TRUE) ## log likelihood
    v <- dcauchy(theta, log=TRUE)   ## log of the prior density
    logden <- u+v  ## log of the posterior
    logden
  }

  for (j in 1:(N-1)) {
    phi <- rnorm(n=1, mean=theta[j], sd=tuning.sd) # generate proposal
    top <- logpost(phi) # Numerator in the ratio
    bot <- logpost(theta[j]) # Denominator in the ratio
    dratio <- exp(top-bot)
    alpha <- min(1, dratio) # Acceptance probability
    if (runif(n=1) < alpha) {  # Test whether to accept
      theta[j+1] <- phi
      count <- count + 1
    } else
      theta[j+1] <- theta[j]
  }
  # End j loop: Metropolis algorithm has finished
  list(MH_rand_walk_estimate =mean(theta), exact_estimate=eresults$exact.post.mean,
       acceptance_rate=100*count/N, interval = quantile(theta, probs=c(0.025, 0.975)),
       ybar=ybar, tuning.sd=tuning.sd)
}
