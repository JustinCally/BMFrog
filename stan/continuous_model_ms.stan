
data {
  int<lower=0> nfiles; // number of files to process (observation periods)
  int<lower=0> nsites; // number of sites
  int<lower=0> nspec; // number of species
  array[nfiles, nspec] real score; // occupancy detection score
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
// re for mu
  real mu_alpha[2, nspec];
  array[nsites, 2, nspec] real eps_raw_mu;
  real<lower=0> eps_sd_mu[2, nspec];
}

transformed parameters {
  array[nsites, nspec] real<lower=0,upper=1> psi;
  array[nspec] vector<lower=0,upper=1>[nfiles] theta;
  array[nsites, 2, nspec] real mu;
for(j in 1:nspec) {
for(i in 1:nsites) {
  psi[i,j] = inv_logit(alpha[j] + eps_sd[j] * eps_raw[i,j]);
  mu[i,1,j] = mu_alpha[1,j] + eps_sd_mu[1,j] * eps_raw_mu[i, 1,j];
  mu[i,2,j] = mu_alpha[2,j] + eps_sd_mu[2,j] * eps_raw_mu[i, 2,j];
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
  mu_alpha[i,j] ~ normal(0,10);
  eps_sd_mu[i,j] ~ normal(0,1);
  eps_raw_mu[,i,j] ~ normal(0,1);
}
  sigma[,j] ~ normal(0,10000);

  beta_theta[j] ~ normal(0,2);

}

for(j in 1:nspec) {
for(i in 1:nsites) {

  if(start_idx_0[i,j] != 0) {
    if(any_seen[i,j] == 0) {
      target += log1m(psi[i,j]) + normal_lpdf(score[start_idx_0[i,j]:end_idx_0[i,j],j] | mu[i, 1, j], sigma[1,j]);
    }
      target += log(psi[i,j]) + log1m(theta[j, start_idx_0[i,j]:end_idx_0[i,j]]) + normal_lpdf(score[start_idx_0[i,j]:end_idx_0[i,j],j] | mu[i, 1, j], sigma[1,j]);
}
  if(start_idx_1[i,j] != 0) {
    target += log(psi[i,j]) + log(theta[j, start_idx_1[i,j]:end_idx_1[i,j]]) + normal_lpdf(score[start_idx_1[i,j]:end_idx_1[i,j],j] | mu[i, 2, j], sigma[2,j]);
}
  if(start_idx_2[i,j] != 0) {
        if(any_seen[i,j] == 0) {
    target += log1m(psi[i,j]) + normal_lpdf(score[start_idx_2[i,j]:end_idx_2[i,j],j] | mu[i, 1, j], sigma[1,j]);
        }
    target += log_sum_exp(log(psi[i,j]) + log_sum_exp(log1m(theta[j, start_idx_2[i,j]:end_idx_2[i,j]])) + normal_lpdf(score[start_idx_2[i,j]:end_idx_2[i,j],j] | mu[i, 1, j], sigma[1,j]),
    log(psi[i,j]) + log_sum_exp(log(theta[j, start_idx_2[i,j]:end_idx_2[i,j]])) + normal_lpdf(score[start_idx_2[i,j]:end_idx_2[i,j],j] | mu[i, 2,j], sigma[2,j]));
}
}
}

}

generated quantities {
    int<lower=0,upper=1> z[nsites, nspec];
    real mu_pred[2,nspec];
for(j in 1:nspec) {
    z[,j] = bernoulli_rng(psi[,j]);
    mu_pred[,j] = normal_rng(mu_alpha[,j], sigma[,j]);
}
}

