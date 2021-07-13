// data block. This must contain the same variables as data list
// that will come from R

data {
int<lower=0> n; // number of sites
int<lower=0> p; // number of sites
vector[n] y;
matrix[n, p] X;
matrix[n, n] dist; // to hold n by n distance matrix
real<lower=0.00001> phi;
vector<lower=0>[2] priorsigma2;
vector<lower=0>[2] priortau2;
}

// These will never change during the MCMC computation. 
transformed data {
  real delta=1e-5;
  vector[n] mu_0 = rep_vector(0, n);
}

// Declare all the parameters to be sampled here

parameters {
  vector[p] beta;
  real<lower=0> sigma_sq;
  vector[n] eta;
  real<lower=0> tau_sq;
}


// Model specification

model {
vector[n] xbmodel;
vector[n] dats;
matrix[n, n] L;
matrix[n, n] Sigma;
real u;

// print(beta)
// print(y)

xbmodel = X * beta; 
 for (i in 1:n) {
   for (j in 1:n) {
     Sigma[i, j] = sigma_sq * exp((-1)*phi*dist[i,j]);
   }
    Sigma[i, i] = Sigma[i, i] + delta;
  }
 L = cholesky_decompose(Sigma);
 
 eta ~   multi_normal_cholesky(mu_0, L);
 sigma_sq ~ inv_gamma(priorsigma2[1], priorsigma2[2]);
 tau_sq ~ inv_gamma(priortau2[1], priortau2[2]);
 // dats ~ normal(xbmodel, sqrt(tau_sq));
 y ~ normal(xbmodel+eta, sqrt(tau_sq));
 // dats = y - eta; 
 
}


