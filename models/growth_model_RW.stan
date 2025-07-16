// RW ONLY MODEL 
data {
  int<lower=0> N_Tr;            // number of trees with rw measurements
  int<lower=0> N_years;         // number of years of data 
  int<lower=0> N_X;             // number of increment values to estimate
  int<lower=0> N_Xobs;          // total number of increments
  int<lower=0> N_taxa;          // total number of taxa
  int<lower=0> Tr2X[N_Tr];      // maps DBH measurements to estimated increment
  int<lower=0> X2Tr[N_X];       // maps estimated increment values to tree ID
  int<lower=0> X2year[N_X];     // maps estimated increment values to year ID
  real logTr[N_Tr];             // log of DBH measurements
  real logXobs[N_Xobs];         // log of increment measurements
  real<lower=0> sig_d_obs;      // mean diameter measurement error
  int<lower=0> idx_Tr[N_Tr, 3]; // indicates first and last value index for each tree ID
  int<lower=0> Xobs2X[N_Xobs];  // maps increment measurements to estimated values

  int<lower=0> X2taxon[N_X];
  //int<lower=0> Tr2taxon[N_Tr];
}

parameters {
  real beta0;
  real beta[N_Tr];
  real<lower=0> beta_sd;
  real<lower=0> beta_t_sd;
  // real<lower=1e-6> sig_x;
  real<lower=1e-6> sig_x[N_taxa];
  real<lower=1e-6> sig_x_obs;
  real<lower=-30, upper=80> D0[N_Tr];
  vector<lower=1e-6> [N_X] X;

  matrix[N_years, N_taxa] beta_t;
}
transformed parameters {
    // process evolution
  vector<lower=-30> [N_X] D;

  for (tree in 1:N_Tr){
    D[idx_Tr[tree,2]] = D0[tree] + 2.0 * X[idx_Tr[tree,2]] / 10.0;
    for (val in (idx_Tr[tree,2]+1):(idx_Tr[tree,3])){
      D[val] = D[val-1] + 2.0 * X[val] / 10.0;
    }
  }
}

model{

  beta0     ~ normal(0, 1.0/0.00001);
  sig_x_obs ~ uniform(1e-6, 2.0);
  //sig_d_obs ~ uniform(1e-6, 1000);
  
  // sig_x     ~ uniform(1e-6, 1000);
  beta_sd   ~ uniform(1e-6, 1000);
  beta_t_sd ~ uniform(1e-6, 1000);
    
  for(tree in 1:N_Tr) {
    D0[tree] ~ uniform(-30, 80);
    beta[tree] ~ normal(beta0, beta_sd);
  }
  
  // for(year in 1:N_years) {
  //   beta_t[year] ~ normal(0, beta_t_sd);
  // }

    // temporal species effect prior 
  for(year in 1:N_years) {
    beta_t[year,] ~ normal(0, beta_t_sd);
  }

  for (taxon in 1:N_taxa) {
    sig_x[taxon] ~ uniform(1e-6, 1000);
  }

  
  // increment likelihood
  
  // estimates 
  for (val in 1:N_X){
    X[val] ~ lognormal(beta[X2Tr[val]] + beta_t[X2year[val], X2taxon[val]], sig_x[X2taxon[val]]);
   }
   
  // RW estimates against RW increments
  for (inc in 1:N_Xobs){
   logXobs[inc] ~ normal(log(X[Xobs2X[inc]]), sig_x_obs);
  }
  
  // diameter estimates against diameter measurements
  for (tree in 1:N_Tr){
    if (logTr[tree] == -999){
    } else {
      logTr[tree] ~ student_t(3, log(D[Tr2X[tree]]), sig_d_obs);
    }
  }
}
