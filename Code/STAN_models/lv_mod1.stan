//Simple factor analysis model 
// BD 11/4/2024
// The following model is based off of https://github.com/cbrown5/ecological-condition-latent-model/blob/main/indicator-model.stan and should cite the associated publication...


data{
  //observed variables 
  int<lower=1> N; 
  vector[N] fric;
  vector[N] feve;
  vector[N] fdiv;
}

parameters{
  //LV
  vector[N] nu;
  
  //Indicators parameters for regression on the latent
  real beta_fric; //regression on latent 
  real beta_feve;
  real beta_fdiv;
  
  real a_fric; //intercepts
  real a_feve; //intercepts
  real a_fdiv; //intercepts
  
  // SDs for observation variables
  real<lower=0> sigma_fric;
  real<lower=0> sigma_feve;
  real<lower=0> sigma_fdiv;
  
}

transformed parameters{

 //Obs var means
 vector[N] fric_hat;
 vector[N] feve_hat;
 vector[N] fdiv_hat;
 
 //
 // Latent regressions
 //
 fric_hat = beta_fric * nu + a_fric;
 feve_hat = beta_feve * nu + a_feve;
 fdiv_hat = beta_fdiv * nu + a_fdiv;

}

model{ 

  //
  //Ecological condition 
  //
  nu ~ std_normal(); 
  //Note variance fixed to 1 by using standard normal here
  
  //Priors for indiator params 
  a_fric ~ normal(0,10);
  a_feve~ normal(0,10);
  a_fdiv ~ normal(0,10);
  
  beta_fric ~ exponential(1.8);
  beta_feve ~ normal(0,1);
  beta_fdiv ~ normal(0,1);
  //Note use of exponential() for beta_fric to ensure positive
  // only values. Helps convergence. 
  
  //Observation errors
  sigma_fric ~  exponential(0.1);
  sigma_feve ~  exponential(0.1);
  sigma_fdiv ~  exponential(0.1);
  
  // Observations
  fric ~ normal(fric_hat, sigma_fric);
  feve ~ normal(feve_hat, sigma_feve);
  fdiv ~ normal(fdiv_hat, sigma_fdiv);
}
