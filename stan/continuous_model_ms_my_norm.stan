
data {
  int<lower=0> nfiles; // number of files to process (observation periods)
  int<lower=0> nsites; // number of site-year combos
  int<lower=0> nspec; // number of species
  int<lower=0> nyear; // number of years
  int<lower=0> nclusters; // number of clusters (pair of recorders)
  int<lower=0> nloc; // number of site locations (independent of years)
  array[nfiles, nspec] real score; // occupancy detection score
  array[nsites, nspec] int loc_code;
  array[nsites, nspec] int cluster_code;
  array[nsites, nspec] int year_code;
  array[nyear, nspec] real year_counts;
  real year_counts_file[nyear];
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
  int<lower=0> psi_nc;
  array[nspec] matrix[nsites, psi_nc] psi_mm; // occ model matrix
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // mu
  real mu_int[2];
  real<lower=0> sigma[2];
  // re for mu
  array[2, nspec] real mu_raw;
  real<lower=0> mu_sd;
  // real alpha[nspec];
  // location random effect
  array[nloc, nspec] real loc_raw;
  real<lower=0> loc_sd[nspec];
  // year random effect
  array[nyear, nspec] real year_raw;
  real<lower=0> year_sd[nspec];
  // site random effect
  array[nsites, nspec] real site_raw;
  real<lower=0> site_sd[nspec];
  // files random effect
  array[nfiles, nspec] real file_raw;
  real<lower=0> file_sd[nspec];
  // call rate coefficient
  array[nspec] vector[theta_nc] beta_theta; //call rate
  // occ coef
  array[nspec] vector[psi_nc] beta_psi; //call rate

}

transformed parameters {
  array[nsites, nspec] real<lower=0,upper=1> psi;
  array[nspec] vector<lower=0,upper=1>[nfiles] theta;
  // random effects
  array[nloc, nspec] real eps_loc;
  array[nyear, nspec] real eps_year;
  array[nsites, nspec] real eps_site;
  array[2, nspec] real mu;

for(j in 1:nspec) {
  // mu
  mu[1,j] = mu_int[1] + mu_sd * mu_raw[1,j];
  mu[2,j] = mu_int[2] + mu_sd * mu_raw[2,j];
  // locations
  for(l in 1:nloc) {
    eps_loc[l,j] = loc_sd[j] * loc_raw[l,j];
  }
  // years
    for(y in 1:nyear) {
    eps_year[y,j] = year_sd[j] * year_raw[y,j];
  }

for(i in 1:nsites) {
  eps_site[i,j] = site_sd[j] * site_raw[i,j];
  psi[i,j] = inv_logit(psi_mm[j, i] * beta_psi[j] + eps_loc[loc_code[i,j], j] + eps_year[year_code[i,j], j] + eps_site[i,j]);
}

  theta[j] = inv_logit(theta_mm[j] * beta_theta[j]); // + (file_sd[j] * to_vector(file_raw[,j])));
}

}


model {
  // priors

for(j in 1:nspec) {
  // alpha[j] ~ normal(0,2);
  // re priors
  loc_raw[,j] ~ normal(0,1);
  loc_sd[j] ~ normal(0,1);
  year_raw[,j] ~ normal(0,1);
  year_sd[j] ~ normal(0,1);
  site_raw[,j] ~ normal(0,1);
  site_sd[j] ~ normal(0,1);
  file_raw[,j] ~ normal(0,1);
  file_sd[j] ~ normal(0,1);

  beta_theta[j] ~ normal(0,2);
  beta_psi[j] ~ normal(0,2);

}
  sigma ~ normal(0,1);
  mu_int[1] ~ normal(-3,1);
  mu_int[2] ~ normal(0,1);
  mu_raw[1,] ~ normal(0,1);
  mu_raw[2,] ~ normal(0,1);

  mu_sd ~ normal(0,1);

for(j in 1:nspec) {
for(i in 1:nsites) {

  if(start_idx_0[i,j] != 0) {
    if(any_seen[i,j] == 0) {
      target += log1m(psi[i,j]) + normal_lpdf(score[start_idx_0[i,j]:end_idx_0[i,j],j] | mu[1, j], sigma[1]);
    }
      target += log(psi[i,j]) + log1m(theta[j, start_idx_0[i,j]:end_idx_0[i,j]]) + normal_lpdf(score[start_idx_0[i,j]:end_idx_0[i,j],j] | mu[1, j], sigma[1]);
}
  if(start_idx_1[i,j] != 0) {
    target += log(psi[i,j]) + log(theta[j, start_idx_1[i,j]:end_idx_1[i,j]]) + normal_lpdf(score[start_idx_1[i,j]:end_idx_1[i,j],j] | mu[2, j], sigma[2]);
}
  if(start_idx_2[i,j] != 0) {
        if(any_seen[i,j] == 0) {
    target += log1m(psi[i,j]) + normal_lpdf(score[start_idx_2[i,j]:end_idx_2[i,j],j] | mu[1, j], sigma[1]);
        }
    target += log_sum_exp(log(psi[i,j]) + log_sum_exp(log1m(theta[j, start_idx_2[i,j]:end_idx_2[i,j]])) + normal_lpdf(score[start_idx_2[i,j]:end_idx_2[i,j],j] | mu[1, j], sigma[1]),
    log(psi[i,j]) + log_sum_exp(log(theta[j, start_idx_2[i,j]:end_idx_2[i,j]])) + normal_lpdf(score[start_idx_2[i,j]:end_idx_2[i,j],j] | mu[2, j], sigma[2]));
}
}
}

}

generated quantities {
    int<lower=0,upper=1> zm[nsites, nspec];
    int<lower=0,upper=1> z[nsites, nspec];
    int<lower=0> z_sum[nyear, nspec];
    array[nyear, nspec] real z_hat;
    int<lower=0> theta_sum[nyear, nspec];
    array[nyear, nspec] real theta_hat;
    array[nsites, nspec] real theta_site;
    array[nsites, nspec] real rel_ab;
    array[nyear, nspec] real rel_ab_sum;
    array[nyear, nspec] real rel_ab_mean;
    real mu_pred[2,nspec];
    real theta_phen [365,nspec];
    array[nspec] vector<lower=0,upper=1>[nfiles] det_pred;
    array[nfiles, nspec] int<lower=0,upper=1> det_z;
    array[nfiles, nspec] real score_pred;
    array[nsites] int<lower=0> species_richness;

for(j in 1:nspec) {
    zm[,j] = bernoulli_rng(psi[,j]);

    mu_pred[, j] = normal_rng(mu[,j], sigma);
    // for(i in 1:365) {
    //   // last two predictors must be cos and sin days
    // theta_phen[i,j] = inv_logit(beta_theta[j,1] + beta_theta[j,theta_nc-1]*cos(2*pi()*i/365) +  beta_theta[j,theta_nc]*sin(2*pi()*i/365));
    // }
    for(y in 1:nyear) z_sum[y,j] = 0;
    for(y in 1:nyear) theta_sum[y,j] = 0;
    for(y in 1:nyear) rel_ab_sum[y,j] = 0;
    for(n in 1:nsites) {
      det_pred[j, site_idx_start[n,j]:site_idx_end[n,j]] = psi[n,j] * theta[j,site_idx_start[n,j]:site_idx_end[n,j]];
      det_z[site_idx_start[n,j]:site_idx_end[n,j],j] = bernoulli_rng(det_pred[j,site_idx_start[n,j]:site_idx_end[n,j]]);
      theta_site[n,j] = mean(det_z[site_idx_start[n,j]:site_idx_end[n,j],j]);
      theta_sum[year_code[n,j],j] += sum(det_z[site_idx_start[n,j]:site_idx_end[n,j],j]);
      z[n,j] = max(zm[n,j], any_seen[n,j]);
      rel_ab[n,j] = z[n,j]*theta_site[n,j];
      rel_ab_sum[year_code[n,j],j] += rel_ab[n,j];
      z_sum[year_code[n,j],j] += z[n,j];
    }
    for(k in 1:nfiles) {
      if(det_z[k,j] == 0) {
        score_pred[k,j] = normal_rng(mu[1, j], sigma[1]);
      } else if(det_z[k,j] == 1) {
        score_pred[k,j] = normal_rng(mu[2, j], sigma[2]);
      }
    }
}
for(j in 1:nspec) {
    for(y in 1:nyear) {
      z_hat[y,j] = z_sum[y,j]/year_counts[y,j];
      rel_ab_mean[y,j] = rel_ab_sum[y,j]/year_counts[y,j];
      theta_hat[y,j] = theta_sum[y,j]/year_counts_file[y];
    }
}
// site-level community score
for(n in 1:nsites) {
  species_richness[n] = sum(z[n,]);
}
}

