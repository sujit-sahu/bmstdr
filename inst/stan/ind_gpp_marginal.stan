// Edited 14th September 2022 to solve missing time problems 
/* 
  datatostan <- list(n=n, tn=tn, m2=nrow(knots.coords), p=p, 
                     missing=missing,  ntmiss=ntmiss, ntobs = ntobs, 
                     data_miss_idx=as.vector(data_miss_idx),  data_obs_idx =  as.vector(data_obs_idx), 
                     time =data$time, nots=length(ots),  ots = ots, nts=nts, start_row=start_row, fin_row=fin_row,  
                     n_misst=n_misst, 
                     Cdist=Cdist, dmat = dmat,  
                     yobs=yobs,  X=X,
                     sigma2_prior =prior.sigma2, 
                     tau2_prior = prior.tau2, phidist = phidist,
                     prior_phi_param =prior.phi.param)
		     */
data {
int<lower=0> n; // number of rows
int<lower=0> tn; // number of times
int<lower=0> m2; // how many knot locations 
int<lower=0> p;
int<lower = 0, upper = 1> missing;  // Is there any missing observation? 
int<lower=0> ntmiss;
int<lower=0> ntobs;
int data_miss_idx[ntmiss];
int data_obs_idx[ntobs];
int time[n]; // gives the time t for a data row 
int nots;  // How many observed times 
int ots[nots]; // Observed times 
int nts[nots]; // Number of observations in time ots 
int start_row[nots];  // Starting row minus 1 for the observations at t th time 
int fin_row[nots];  // Finish row  for the observations at t th time 
// blank line
int n_misst;  // Number of time points without any observations  
// int misst[n_misst];  // Which times have missing observations 
// Do not need this here. 
matrix[n, m2] Cdist; // to hold n by m2 distance matrix
matrix[m2, m2] dmat; // to hold m2 by m2 distance matrix
real yobs[ntobs];
matrix[n, p] X;
real<lower=0> sigma2_prior[2];
real<lower=0> tau2_prior[2];
int<lower=0> phidist; 
real <lower=0> prior_phi_param[2];
}

transformed data { 
  int sumnt2=0; 
  
  for (i in 1:nots) {
    sumnt2 = sumnt2 + nts[i] * nts[i]; 
  }
} 

parameters {
  vector[p] beta;
  real<lower=0> phi;
  real<lower=0> sigma_sq;
  real z_miss[missing ? ntmiss : 0];  // Only define z_miss if there is missing data 
  // vector[include_alpha ? N : 0] alpha;
  real<lower=0> tau_sq;
}

transformed parameters {
vector[n] xbmodel;
real bigS[sumnt2]; 
{ 
matrix [m2, m2] Sigma;
matrix [m2, m2] Swinv;
matrix [n, m2] Cmat; 


xbmodel = X * beta;

for(i in 1:(m2-1)){
for(j in (i+1):m2){
     Sigma[i,j] = exp((-1)*phi*dmat[i,j]);
   //  Sigma[i,j] = 0.0;
     Sigma[j,i] = Sigma[i,j];
   }
 }
for(i in 1:m2) Sigma[i, i] = 1.0;

Swinv = inverse_spd(Sigma); 

for (i in 1:n) { 
  for (j in 1:m2) { 
    Cmat[i, j] = exp((-1)*phi*Cdist[i,j]);
    }
  }
  
 { 
 int m=1; 
 for (i in 1:nots) { 
   // nts[i] dimensional normal distribution 
   matrix [nts[i], nts[i]] St; 
   matrix [nts[i], m2] Ct; 
   
   Ct = Cmat[start_row[i]:fin_row[i], ];
   St = sigma_sq * Ct * Swinv * Ct'; 
   for (j in 1:nts[i]) { 
     St[j, j] += tau_sq; 
   } /* j loop */
   
    for (k in 1:nts[i]) { 
      for (j in 1:nts[i]) { 
      bigS[m] = St[k, j]; 
      m = m +1; 
      }
   }

}
}
}
 /*
 print("beta= ", beta); 
 print("sigma sq = ", sigma_sq); 
 print("tau sq = ", tau_sq); 
 print("range = ", 3.0/phi);
 */ 
 }

model {
vector[n] z1;
int m=1; 


 sigma_sq ~ inv_gamma(sigma2_prior[1], sigma2_prior[2]);
 tau_sq ~ inv_gamma(tau2_prior[1], tau2_prior[2]);
 if (phidist ==0)   phi ~ uniform(prior_phi_param[1], prior_phi_param[2]); 
 if (phidist ==1)   phi ~ gamma(prior_phi_param[1], prior_phi_param[2]); 
 if (phidist ==2)   phi ~ cauchy(prior_phi_param[1], prior_phi_param[2]); 
 if (phidist>2) reject("Wrong prior distribution for phi; found phidist=", phidist);

  for (i in 1:ntobs)
    z1[data_obs_idx[i]] = yobs[i];
 
 if (missing>0 ) { 
  for (k in 1:ntmiss)
    z1[data_miss_idx[k]] = z_miss[k];
 }

 for (i in 1:nots) { 
   // nts[i] dimensional normal distribution 
   vector[nts[i]] zt;
   vector[nts[i]] mut; 
   matrix [nts[i], nts[i]] St; 
   
   
   zt = z1[start_row[i]:fin_row[i]]; 
   mut = xbmodel[start_row[i]:fin_row[i]];
  
     for (k in 1:nts[i]) { 
      for (j in 1:nts[i]) { 
      St[k, j] = bigS[m]; 
      m = m +1; 
      }
     }
    
   zt ~ multi_normal(mut, St); 
} // i loop 

 

}

