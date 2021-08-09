# Find the posterior mean and variance  using `exact' methods - no sampling 
# or approximation  
Bnormal(kprior = 1, prior.M = 1, prior.sigma2 = c(2, 1))
# Use default non-informative prior
Bnormal(mu0 = 0) 
# Start creating table
y <-  ydata
mu0 <-  mean(y)
kprior <-  1
prior.M  <-  1
prior.sigma2 <- c(2, 1)
N  <-  10000
eresults <- Bnormal(package = "exact", y = y, mu0 = mu0, kprior = kprior,
    prior.M = prior.M, prior.sigma2 = prior.sigma2)
eresults 
# Run Gibbs sampling
samps <- Bnormal(package = "RGibbs", y = y, mu0 = mu0, kprior = kprior,
    prior.M = prior.M, prior.sigma2 = prior.sigma2, N = N)
gres <- list(mean_theta = mean(samps[, 1]), var_theta = var(samps[, 1]),
    mean_sigma2 = mean(samps[, 2]), var_sigma2 = var(samps[, 2]))
glow <- list(theta_low = quantile(samps[, 1], probs = 0.025), var_theta = NA,
    sigma2_low = quantile(samps[, 2], probs = 0.025), var_sigma2 = NA)
gupp <- list(theta_low = quantile(samps[, 1], probs = 0.975), var_theta = NA,
    sigma2_low = quantile(samps[, 2], probs = 0.975), var_sigma2 = NA)
a <- rbind(unlist(eresults), unlist(gres), unlist(glow), unlist(gupp))
cvsamp <- sqrt(samps[, 2])/samps[, 1]
cv <- c(NA, mean(cvsamp), quantile(cvsamp, probs = c(0.025, 0.975)))
u <- data.frame(a, cv)
rownames(u) <- c("Exact", "Estimate", "2.5%", "97.5%")
print(u)
# End create table 
##
# Compute using the model fitted by Stan 
u <- Bnormal(package = "stan", kprior = 1, prior.M = 1, prior.sigma = c(2, 1),
    N = 2000, burn.in = 1000)
print(u)
###
# Compute using INLA 
if (require(INLA)) {
    u <- Bnormal(package = "inla", kprior = 1, prior.M = 1, prior.sigma = c(2,
        1), N = 1000)
    print(u)
}
