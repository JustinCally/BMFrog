
data {
  int<lower=0> nfiles; // number of files to process (observation periods)
  int<lower=0> nsites; // number of sites
  int<lower=0> nspec; // number of species
  array[nfiles, nspec] real<lower=0,upper=1> score; // occupancy detection score
  array[nsites, nspec] int site_idx_start;
  array[nsites, nspec] int site_idx_end;
  int start_idx_0[nsites, nspec];
  int end_idx_0[nsites, nspec];
  int start_idx_1[nsites, nspec];
  int end_idx_1[nsites, nspec];
  int start_idx_2[nsites, nspec];
  int end_idx_2[nsites, nspec];
  int any_seen[nsites, nspec];
  int<lower=0> theta_nc;
  array[nspec] matrix[nfiles, theta_nc] theta_mm; // model matrix
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> sigma[2, nspec];
  real alpha[nspec];
  array[nsites, nspec] real eps_raw;
  real<lower=0> eps_sd[nspec];
  array[nspec] vector[theta_nc] beta_theta; //call rate
// mu
  real<lower=0> mu[2, nspec];
}

transformed parameters {
  array[nsites, nspec] real<lower=0,upper=1> psi;
  array[nspec] vector<lower=0,upper=1>[nfiles] theta;

for(j in 1:nspec) {
for(i in 1:nsites) {
  psi[i,j] = inv_logit(alpha[j] + eps_sd[j] * eps_raw[i,j]);
}

  theta[j] = inv_logit(theta_mm[j] * beta_theta[j]);
}

}


model {
  // priors

for(j in 1:nspec) {
  alpha[j] ~ normal(0,2);
  eps_raw[,j] ~ normal(0,1);
  eps_sd[j] ~ normal(0,1);

for(i in 1:2) {
  mu[i,j] ~ normal(0,10);
}
  sigma[,j] ~ normal(0,10);

  beta_theta[j] ~ normal(0,2);

}

for(j in 1:nspec) {
for(i in 1:nsites) {

  if(start_idx_0[i,j] != 0) {
    if(any_seen[i,j] == 0) {
      target += log1m(psi[i,j]) + beta_lpdf(score[start_idx_0[i,j]:end_idx_0[i,j],j] | mu[1, j], sigma[1,j]);
    }
      target += log(psi[i,j]) + log1m(theta[j, start_idx_0[i,j]:end_idx_0[i,j]]) + beta_lpdf(score[start_idx_0[i,j]:end_idx_0[i,j],j] | mu[1, j], sigma[1,j]);
}
  if(start_idx_1[i,j] != 0) {
    target += log(psi[i,j]) + log(theta[j, start_idx_1[i,j]:end_idx_1[i,j]]) + beta_lpdf(score[start_idx_1[i,j]:end_idx_1[i,j],j] | mu[2, j], sigma[2,j]);
}
  if(start_idx_2[i,j] != 0) {
        if(any_seen[i,j] == 0) {
    target += log1m(psi[i,j]) + beta_lpdf(score[start_idx_2[i,j]:end_idx_2[i,j],j] | mu[1, j], sigma[1,j]);
        }
    target += log_sum_exp(log(psi[i,j]) + log_sum_exp(log1m(theta[j, start_idx_2[i,j]:end_idx_2[i,j]])) + beta_lpdf(score[start_idx_2[i,j]:end_idx_2[i,j],j] | mu[1, j], sigma[1,j]),
    log(psi[i,j]) + log_sum_exp(log(theta[j, start_idx_2[i,j]:end_idx_2[i,j]])) + beta_lpdf(score[start_idx_2[i,j]:end_idx_2[i,j],j] | mu[2,j], sigma[2,j]));
}
}
}

}

generated quantities {
    int<lower=0,upper=1> z[nsites, nspec];
    real mu_pred[2,nspec];
    real theta_phen [365,nspec];
    array[nspec] vector<lower=0,upper=1>[nfiles] det_pred;
    array[nfiles, nspec] int<lower=0,upper=1> det_z;
    array[nfiles, nspec] real<lower=0,upper=1> score_pred;

for(j in 1:nspec) {
    z[,j] = bernoulli_rng(psi[,j]);
    mu_pred[,j] = beta_rng(mu[,j], sigma[,j]);
    for(i in 1:365) {
      // last two predictors must be cos and sin days
    theta_phen[i,j] = inv_logit(beta_theta[j,1] + beta_theta[j,theta_nc-1]*cos(2*pi()*i/365) +  beta_theta[j,theta_nc]*sin(2*pi()*i/365));
    }
    for(n in 1:nsites) {
      det_pred[j, site_idx_start[n,j]:site_idx_end[n,j]] = psi[n,j] * theta[j,site_idx_start[n,j]:site_idx_end[n,j]];
      det_z[site_idx_start[n,j]:site_idx_end[n,j],j] = bernoulli_rng(det_pred[j,site_idx_start[n,j]:site_idx_end[n,j]]);
    }
    for(k in 1:nfiles) {
      if(det_z[k,j] == 0) {
        score_pred[k,j] = beta_rng(mu[1,j], sigma[1,j]);
      } else if(det_z[k,j] == 1) {
        score_pred[k,j] = beta_rng(mu[2,j], sigma[2,j]);
      }
    }
}
}

