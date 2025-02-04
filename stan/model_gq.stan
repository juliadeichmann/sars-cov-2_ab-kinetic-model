// Stan code that accompanies the following manuscript:
//
// Title:   Predicting antibody kinetics and duration of protection against 
//          SARS-CoV-2 following vaccination from sparse serological data
// Authors: Julia Deichmann, Noam Barda, Michal Canetti, Yovel Peretz, Yael Ottolenghi,
//          Yaniv Lustig, Gili Regev-Yochay, Marc Lipsitch
//
// Date:    February 4, 2025
// Author:  Julia Deichmann <jdeichmann@hsph.harvard.edu>

functions {
  vector AbKin(vector t, vector tv, real A0, vector gB, vector delta, real rho, real c1, real r, real c2) {   

    vector[size(t)] y;
    vector[size(t)] delta_t;

    for (v in 1:size(tv)){

      for (i in 1:size(t)){
        delta_t[i] = max([t[i]-tv[v]-delta[v], 0.]);
      }

      if (v==1){
        y = A0 + gB[v] * (rho * (exp(-c2*delta_t) - exp(-r*delta_t)) / (r - c2) + (1-rho) * (exp(-c1*delta_t) - exp(-r*delta_t)) / (r - c1));
      }
      else {
        y = y + gB[v] * (rho * (exp(-c2*delta_t) - exp(-r*delta_t)) / (r - c2) + (1-rho) * (exp(-c1*delta_t) - exp(-r*delta_t)) / (r - c1));
      }
    }
        
    return y;
  }
  
  real transformation(real y_old){
    real y_BAU;
    y_BAU = exp(4.506+0.6634*log(y_old)-0.0852*(log(y_old))^2 + 0.0403 * (log(y_old))^3);
    return y_BAU;
  }
}


// The input data.
data {
  int<lower=0> N;          // total number of observations
  int<lower=0> np;         // number of individuals
  int<lower=0> nv;         // number of vaccination events
  
  array[np] int nd;        // number of observations person [np]
  array[np] int nb;        // number of doses per person    [np]
  
  int nc;                  // number of covariates
  matrix[np,nc] X;         // covariate matrix
  
  vector[N] t;             // observation time points
  vector[N] obs;           // observed antibody concentration
  
  vector[nv] tvacc;        // vaccination time points
  
  real sigma_prior;        // prior of sigma parameters, random effects
  
  int nt;                  // number of prediction time points
  vector[nt] tpred;        // prediction time points
}


// The parameters to be estimated.
parameters {
  real log_A0_pop;                    // pre-vacc antibody concentration
  vector[2] log_gB_pop;               // antibody boost per dose
  real<lower=-6,upper=6> rho_pop;     // fraction of short-lived APCs
  real log_r;                         // antibody decay rate
  ordered[2] log_c_pop;               // decay rate of short- and long-lived APCs
  
  vector[nc] beta_gB1;                // covariate parameters
  vector[nc] beta_gB2;
  vector[nc] beta_c1;
  vector[nc] beta_c2;
  
  real<lower=0> A0_sigma;             // standard deviation of random effects
  real<lower=0> gB1_sigma;
  real<lower=0> gB2_sigma;
  real<lower=0> rho_sigma;
  real<lower=0> c1_sigma;
  real<lower=0> c2_sigma;
  
  vector[np] u_A0;                    // individual deviation, random effects 
  vector[np] u_gB1;
  vector[np] u_gB2;
  vector[np] u_rho;
  vector[np] u_c1;
  vector[np] u_c2;
  
  real<lower=0> sigma;                // observational error
}


transformed parameters {
  
  vector[N] ypred;
  
  {
    vector[np] log_A0;
    vector[np] log_gB_1;
    vector[np] log_gB_2;
    vector[np] log_c1;
    vector[np] log_c2;
    
    log_A0 = log_A0_pop + u_A0 * A0_sigma;
    log_gB_1 = log_gB_pop[1] + X * beta_gB1 + u_gB1 * gB1_sigma;
    log_gB_2 = log_gB_pop[2] + X * beta_gB2 + u_gB2 * gB2_sigma;                       
    log_c1 = log_c_pop[1] + X * beta_c1 + u_c1 * c1_sigma;
    log_c2 = log_c_pop[2] + X * beta_c2 + u_c2 * c2_sigma;
    
    vector[np] rho;
    rho = inv_logit(rho_pop + u_rho * rho_sigma);
    
    vector[np] A0;
    vector[np] gB_1;
    vector[np] gB_2;
    real r;
    vector[np] c1;
    vector[np] c2;
      
    A0 = exp(log_A0);
    gB_1 = exp(log_gB_1);
    gB_2 = exp(log_gB_2);
    r = exp(log_r);
    c1 = exp(log_c1);
    c2 = exp(log_c2);
    
    real delta_1 = 10;
    real delta_2 = 4;

    int ndd = 1;
    int nbb = 1;
    
    for (i in 1:np){
      int steps = nd[i];
      int doses = nb[i];

      ypred[ndd:(ndd+steps-1)] = AbKin(t[ndd:(ndd+steps-1)], tvacc[nbb:(nbb+doses-1)], A0[i], [gB_1[i], gB_2[i]]', [delta_1, delta_2]', rho[i], c1[i], r, c2[i]);
      ndd += steps;
      nbb += doses;
    }
  } 
}


// The model to be estimated.
model {
  // Priors
  log_A0_pop ~ normal(log(0.05), 0.5);
  u_A0 ~ normal(0, 1);
  A0_sigma ~ normal(0, sigma_prior);
  
  log_gB_pop[1] ~ normal(log(0.8), 0.2);    
  u_gB1 ~ normal(0, 1);
  gB1_sigma ~ normal(0, sigma_prior);
  
  log_gB_pop[2] ~ normal(log(10), 0.2);
  u_gB2 ~ normal(0, 1);
  gB2_sigma ~ normal(0, sigma_prior);

  rho_pop ~ normal(3, 0.1);
  u_rho ~ normal(0, 1);
  rho_sigma ~ normal(0, sigma_prior);
  
  log_r ~ normal(log(0.033), 0.1);
  
  log_c_pop[1] ~ normal(log(0.002), 0.5);
  u_c1 ~ normal(0, 1);
  c1_sigma ~ normal(0, sigma_prior);
  
  log_c_pop[2] ~ normal(log(0.14), 0.1);
  u_c2 ~ normal(0, 1);
  c2_sigma ~ normal(0, sigma_prior);
  
  for (i in 1:nc){
    beta_gB1[i] ~ normal(0, .2);
    beta_gB2[i] ~ normal(0, .2);
    beta_c1[i] ~ normal(0, .2);
    beta_c2[i] ~ normal(0, .2);
  }
  
  // Prior - measurement error
  sigma ~ normal(0, 0.1);

  // observation model                          
  obs ~ lognormal(log(ypred), sigma);  
}


// Posterior predictive simulation
generated quantities {
  vector[np*nt] pp_yhat;
  vector[np] pp_ymax;
  vector[np] t_max;
  vector[np] t_prot;
  int ntt = 1;
  int nbb = 1;
  
  vector[1] tmp;
  
  {
    vector[np] log_A0;
    vector[np] log_gB_1;
    vector[np] log_gB_2;
    vector[np] log_c1;
    vector[np] log_c2;
    
    log_A0 = log_A0_pop + u_A0 * A0_sigma;
    log_gB_1 = log_gB_pop[1] + X * beta_gB1 + u_gB1 * gB1_sigma;  
    log_gB_2 = log_gB_pop[2] + X * beta_gB2 + u_gB2 * gB2_sigma;  
    log_c1 = log_c_pop[1] + X * beta_c1 + u_c1 * c1_sigma;   
    log_c2 = log_c_pop[2] + X * beta_c2 + u_c2 * c2_sigma;   
    
    vector[np] rho;
    rho = inv_logit(rho_pop + u_rho * rho_sigma);
    
    vector[np] A0;
    vector[np] gB_1;
    vector[np] gB_2;
    real r;
    vector[np] c1;
    vector[np] c2;
      
    A0 = exp(log_A0);
    gB_1 = exp(log_gB_1);
    gB_2 = exp(log_gB_2);
    r = exp(log_r);
    c1 = exp(log_c1);
    c2 = exp(log_c2);
    
    real delta_1 = 10;
    real delta_2 = 4;
    
    
    for (i in 1:np){          
      int doses = nb[i];
      
      // antibody trajectory
      pp_yhat[ntt:(ntt+nt-1)] = AbKin(tpred, tvacc[nbb:(nbb+doses-1)], A0[i], [gB_1[i], gB_2[i]]', [delta_1, delta_2]', rho[i], c1[i], r, c2[i]);
      
      // peak antibody level
      pp_ymax[i] = max(pp_yhat[ntt:(ntt+nt-1)]);

      // time of peak
      int j = 1;
      while (pp_yhat[ntt+j-1] < pp_ymax[i]){
        j += 1;
      }
      t_max[i] = tpred[j];
      
      // duration of protection (time to y<500)
      if (transformation(pp_yhat[ntt+j-1]) > 500){
        while (transformation(pp_yhat[ntt+j-1]) > 500 && j < nt){
          j += 1;
          }
        if (transformation(pp_yhat[ntt+j-1]) <= 500){
          t_prot[i] = tpred[j];
        }
        else {
          // extend prediction until y<500
          int k = 1;
          tmp = AbKin(to_vector([tpred[j] + k*5]), tvacc[nbb:(nbb+doses-1)], A0[i], [gB_1[i], gB_2[i]]', [delta_1, delta_2]', rho[i], c1[i], r, c2[i]);
          while (transformation(tmp[1]) > 500 && k < 80){
            k += 1;
            tmp = AbKin(to_vector([tpred[j] + k*5]), tvacc[nbb:(nbb+doses-1)], A0[i], [gB_1[i], gB_2[i]]', [delta_1, delta_2]', rho[i], c1[i], r, c2[i]);
          }
          t_prot[i] = tpred[j] + k*5;
        }
      }
      else {
        t_prot[i] = 0;
      }
      
      ntt += nt;
      nbb += doses;
    }
  }
}
