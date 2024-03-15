
data {
  int<lower=0> nfiles; // number of files to process (observation periods)
  int<lower=0> nsites; // number of sites
  array[nfiles] real score; // occupancy detection score
  int site_idx[nfiles]; // site
  array[nfiles] int d_idx; //detection index
  int start_idx_0[nsites];
  int end_idx_0[nsites];
  int start_idx_1[nsites];
  int end_idx_1[nsites];
  int start_idx_2[nsites];
  int end_idx_2[nsites];
  int any_seen[nsites];
  int<lower=0> theta_nc;
  matrix[nfiles, theta_nc] theta_mm; // model matrix
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu[2];
  real<lower=0> sigma[2];
  real alpha;
  vector[nsites] eps_raw;
  real<lower=0> eps_sd;
  vector[theta_nc] beta_theta; //call rate
}

transformed parameters {
  vector<lower=0,upper=1>[nsites] psi;
  vector<lower=0,upper=1>[nfiles] theta;

for(i in 1:nsites) {
  psi[i] = inv_logit(alpha + eps_sd * eps_raw[i]);
}

  theta = inv_logit(theta_mm * beta_theta);

}


model {
  // priors
  alpha ~ normal(0,2);
  eps_raw ~ normal(0,1);
  eps_sd ~ normal(0,1);

  mu ~ normal(0,5);
  sigma ~ normal(0,5);

  beta_theta ~ normal(0,2);


for(i in 1:nsites) {

  if(start_idx_0[i] != 0) {
    if(any_seen[i] == 0) {
      target += log1m(psi[i]) + normal_lpdf(score[start_idx_0[i]:end_idx_0[i]] | mu[1], sigma[1]);
    }
      target += log(psi[i]) + log1m(theta[start_idx_0[i]:end_idx_0[i]]) + normal_lpdf(score[start_idx_0[i]:end_idx_0[i]] | mu[1], sigma[1]);
}
  if(start_idx_1[i] != 0) {
    target += log(psi[i]) + log(theta[start_idx_1[i]:end_idx_1[i]]) + normal_lpdf(score[start_idx_1[i]:end_idx_1[i]] | mu[2], sigma[2]);
}
  if(start_idx_2[i] != 0) {
        if(any_seen[i] == 0) {
    target += log1m(psi[i]) + normal_lpdf(score[start_idx_2[i]:end_idx_2[i]] | mu[1], sigma[1]);
        }
    target += log_sum_exp(log(psi[i]) + log_sum_exp(log1m(theta[start_idx_2[i]:end_idx_2[i]])) + normal_lpdf(score[start_idx_2[i]:end_idx_2[i]] | mu[1], sigma[1]),
    log(psi[i]) + log_sum_exp(log(theta[start_idx_2[i]:end_idx_2[i]])) + normal_lpdf(score[start_idx_2[i]:end_idx_2[i]] | mu[2], sigma[2]));
}
}

}

generated quantities {
    int<lower=0,upper=1> z[nsites];
    z = bernoulli_rng(psi);
}

