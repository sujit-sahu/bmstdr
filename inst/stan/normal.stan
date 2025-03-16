data {
int<lower=0> n; // number of sites
array[n] real y;
real mu;   // for the prior
real mprior;  // for the prior
real aprior; // for the prior
real bprior;  // for the prior
}

parameters {
  real theta;
  real<lower=0> sigma2;
}

transformed parameters {
real tau;
tau = sqrt(sigma2/mprior);
}

model {
sigma2 ~ inv_gamma(aprior, bprior);
 theta ~ normal(mu, tau);  // parameterised with mean and sd
 y ~ normal(theta, sqrt(sigma2));
}
