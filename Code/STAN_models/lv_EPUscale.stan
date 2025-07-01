data {
    // Observed variables
    int<lower=1> N;                       // Number of observations
    vector[N] fric;                // Richness
    vector[N] feve;              // Evenness
    vector[N] fdiv;              // Redundancy
    int<lower=1> T;                       // Number of years
    int<lower=1, upper=T> year[N];        // Year indices for each observation
    int<lower=1> M;                       // Number of EPUs
    int<lower=1, upper=M> epu[N];         // EPU indices for each site (grouped by EPU)
}

parameters {
    // Latent variables for each year and EPU
    matrix[M, T] nu;                      // Latent variables (one for each EPU and year)
    
    // Regression coefficients for each site
    real<lower=0> beta_fric;                // Regression coefficient for richness
    real beta_feve;              // Regression coefficient for evenness
    real beta_fdiv;              // Regression coefficient for redundancy
    
    // Intercepts for each site
    real a_fric[M];               // Intercept for richness 
    real a_feve[M];             // Intercept for evenness 
    real a_fdiv[M];             // Intercept for redundancy
    
    
    // SDs for the random intercepts and coefficients
    real<lower=0> sigma_a_fric;   // SD for richness intercept
    real<lower=0> sigma_a_feve; // SD for evenness intercept
    real<lower=0> sigma_a_fdiv; // SD for redundancy intercept
    //real<lower=0> sigma_beta_feve; // SD for evenness slope
    //real<lower=0> sigma_beta_fdiv; // SD for redundancy slope

    // SDs for the observations
    real<lower=0> sigma_fric;     // SD for richness
    real<lower=0> sigma_feve;   // SD for evenness
    real<lower=0> sigma_fdiv;   // SD for redundancy
}

transformed parameters {
    // Predicted values for each observation
    vector[N] fric_hat;
    vector[N] feve_hat;
    vector[N] fdiv_hat;
    

    // Latent variable regression predictions (for each observation)
    for (n in 1:N) {
        fric_hat[n] = beta_fric * nu[epu[n], year[n]] + a_fric[epu[n]];
        feve_hat[n] = beta_feve * nu[epu[n], year[n]] + a_feve[epu[n]];
        fdiv_hat[n] = beta_fdiv * nu[epu[n], year[n]] + a_fdiv[epu[n]];
    }
}

model {
    // Priors for latent variable nu (standard normal prior for each EPU and year)
    for (m in 1:M) {
        for (t in 1:T) {
            nu[m, t] ~ std_normal();  // Latent variable for each EPU and year
        }
    }

    // Priors for raw intercepts (before applying non-centered parameterization)
    a_fric ~ normal(0, sigma_a_fric);   // Mean intercept for richness
    a_feve ~ normal(0, sigma_a_feve); // Mean intercept for evenness
    a_fdiv ~ normal(0, sigma_a_fdiv); // Mean intercept for redundancy

    // Priors for regression coefficients (betas)
    beta_fric ~ exponential(1.8);   // Positive-only prior for richness coefficient
    beta_feve ~ normal(0,1);     // Normal prior for evenness coefficient
    beta_fdiv ~ normal(0,1);     // Normal prior for redundancy coefficient

    // Priors for the SDs of the random intercepts and regression coefficients
    sigma_a_fric ~ exponential(1);
    sigma_a_feve ~ exponential(1);
    sigma_a_fdiv ~ exponential(1);
    //sigma_beta_feve ~ exponential(1);
    //sigma_beta_fdiv ~ exponential(1);

    // Priors for the observation errors
    sigma_fric ~ exponential(1);
    sigma_feve ~ exponential(1);
    sigma_fdiv ~ exponential(1);

    // Likelihood for the observations
    fric ~ normal(fric_hat, sigma_fric);
    feve ~ normal(feve_hat, sigma_feve);
    fdiv ~ normal(fdiv_hat, sigma_fdiv);
}


