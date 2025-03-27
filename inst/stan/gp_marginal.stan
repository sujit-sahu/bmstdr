
/*
datatostan <- list(sn=sn, tn=tn, nT=nT, p=p, ntmiss=ntmiss, ntobs = ntobs, missing=missing,
             data_miss_idx=data_miss_idx,  data_obs_idx =  data_obs_idx,
               yobs=yobs,  X=X,
             sigma2_prior=prior.sigma2,
             tau2_prior = prior.tau2,
             phidist = phidist,
             prior_phi_param =prior.phi.param,
             dist = alldistmat, verbose)
	     */
data {
int<lower=0> sn; // number of sites
int<lower=0> tn; // number of times
int<lower=0> nT;
int<lower=0> p;
int<lower=0> ntmiss;
int<lower=0> ntobs;
int <lower=0, upper=1> missing;
array[ntmiss] int data_miss_idx;
array[ntobs] int data_obs_idx;
array[ntobs] real yobs;
matrix[nT, p] X;
array[2] real<lower=0> sigma2_prior;
array[2] real<lower=0> tau2_prior;
int<lower=0> phidist; 
array[2] real <lower=0> prior_phi_param;
matrix[sn, sn] dist; // to hold sn by sn distance matrix
int <lower=0> verbose;
}

transformed data {
  vector[sn] mu_0 = rep_vector(0, sn);
}

parameters {
  vector[p] beta;
  real<lower=0> phi;
  real<lower=0> sigma_sq;
  real<lower=0> tau_sq;
  array[missing ? ntmiss : 0] real z_miss;  // Only define z_miss if there is missing data 
}


model {
vector[nT] xbmodel;
matrix[tn, sn] mus;
matrix[tn, sn] dats;
array[nT] real z1;
matrix [sn, sn] L;
matrix [sn, sn] Sigma;
real u;


 sigma_sq ~ inv_gamma(sigma2_prior[1], sigma2_prior[2]);
 tau_sq ~ inv_gamma(tau2_prior[1], tau2_prior[2]);
 if (phidist ==0)   phi ~ uniform(prior_phi_param[1], prior_phi_param[2]); 
 if (phidist ==1)   phi ~ gamma(prior_phi_param[1], prior_phi_param[2]); 
 if (phidist ==2)   phi ~ cauchy(prior_phi_param[1], prior_phi_param[2]); 
 if (phidist>2) reject("Wrong prior distribution for phi; found phidist=", phidist);


for(i in 1:(sn-1)){
for(j in (i+1):sn){
     Sigma[i,j] = sigma_sq * exp((-1)*phi*dist[i,j]);
   //  Sigma[i,j] = 0.0;
     Sigma[j,i] = Sigma[i,j];
   }
 }
for(i in 1:sn) Sigma[i, i] = sigma_sq + tau_sq;
L = cholesky_decompose(Sigma);

xbmodel = X * beta;

for (i in 1:ntobs)
    z1[data_obs_idx[i]] = yobs[i];
    
if (missing>0)  { 
 for (k in 1:ntmiss)
    z1[data_miss_idx[k]] = z_miss[k];
}

for (i in 1:tn) {
  for (j in 1:sn) {
     dats[i, j] = z1[i + (j-1) * tn];
     mus[i, j]=  xbmodel[i + (j-1) * tn];
   }
}


 for (t  in 1:tn) {
 dats[t] ~ multi_normal_cholesky(mus[t], L);

 }

if (verbose >0){
 print("beta= ", beta); 
 print("sigma sq = ", sigma_sq); 
 print("tau sq = ", tau_sq); 
 print("range = ", 3.0/phi);
}


}
