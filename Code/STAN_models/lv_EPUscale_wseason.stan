data {
    // Observed variables
    int<lower=1> N;                       // Number of observations
    vector[N] fric;                       // Richness
    vector[N] feve;                       // Evenness
    vector[N] fdiv;                       // Redundancy
    int<lower=1> T;                       // Number of years
    int<lower=1, upper=T> year[N];        // Year indices for each observation
    int<lower=1> M;                       // Number of EPUs
    int<lower=1, upper=M> epu[N];         // EPU indices for each site (grouped by EPU)
    int<lower=1> S;                       // Number of seasons
    int<lower=1, upper=S> season[N];      // Season indices for each observation
}

parameters {
    // Latent variable for each observation (indexed by observation number)
    vector[N] nu;                        // Latent variables for each observation (EPU-specific latent state)
    
    // Regression coefficients for each site
    real<lower=0> beta_fric;              // Regression coefficient for richness
    real beta_feve;                       // Regression coefficient for evenness
    real beta_fdiv;                       // Regression coefficient for redundancy
    
    // Intercepts for each site and season
    real a_fric[M];                       // EPU-specific intercept for richness
    real a_feve[M];                       // EPU-specific intercept for evenness
    real a_fdiv[M];                       // EPU-specific intercept for redundancy
    
    real a_fric_season[S];                // Season-specific intercept for richness
    real a_feve_season[S];                // Season-specific intercept for evenness
    real a_fdiv_season[S];                // Season-specific intercept for redundancy
    
    // SDs for the random intercepts and coefficients
    real<lower=0> sigma_a_fric;           // SD for richness EPU intercept
    real<lower=0> sigma_a_feve;           // SD for evenness EPU intercept
    real<lower=0> sigma_a_fdiv;           // SD for redundancy EPU intercept
    real<lower=0> sigma_a_fric_season;    // SD for richness season intercept
    real<lower=0> sigma_a_feve_season;    // SD for evenness season intercept
    real<lower=0> sigma_a_fdiv_season;    // SD for redundancy season intercept
    
    // SDs for the observations
    real<lower=0> sigma_fric;             // SD for richness
    real<lower=0> sigma_feve;             // SD for evenness
    real<lower=0> sigma_fdiv;             // SD for redundancy
}

transformed parameters {
    // Predicted values for each observation
    vector[N] fric_hat;
    vector[N] feve_hat;
    vector[N] fdiv_hat;
    
    // Predicted values based on the latent state, intercepts, and coefficients
    for (n in 1:N) {
        fric_hat[n] = beta_fric * nu[n] + a_fric[epu[n]] + a_fric_season[season[n]];
        feve_hat[n] = beta_feve * nu[n] + a_feve[epu[n]] + a_feve_season[season[n]];
        fdiv_hat[n] = beta_fdiv * nu[n] + a_fdiv[epu[n]] + a_fdiv_season[season[n]];
    }
}

model {
    // Priors for the latent variables (standard normal prior for each observation)
    nu ~ std_normal();  // Latent variable for each observation, indexed by the observation number
    
    // Priors for raw intercepts (before applying non-centered parameterization)
    a_fric ~ normal(0, sigma_a_fric);   // EPU mean intercept for richness
    a_feve ~ normal(0, sigma_a_feve);   // EPU mean intercept for evenness
    a_fdiv ~ normal(0, sigma_a_fdiv);   // EPU mean intercept for redundancy
    
    // Priors for season-specific intercepts
    // a_fric_season ~ normal(0, sigma_a_fric_season);  // Season-specific intercept for richness
    // a_feve_season ~ normal(0, sigma_a_feve_season);  // Season-specific intercept for evenness
    // a_fdiv_season ~ normal(0, sigma_a_fdiv_season);  // Season-specific intercept for redundancy
    
    a_fric_season ~ normal(0, sigma_a_fric_season);  // Season-specific intercept for richness
    a_feve_season ~ normal(0, sigma_a_feve_season);  // Season-specific intercept for evenness
    a_fdiv_season ~ normal(0, sigma_a_fdiv_season);  // Season-specific intercept for redundancy
    
    // Priors for regression coefficients (betas)
    beta_fric ~ exponential(1.8);        // Positive-only prior for richness coefficient
    beta_feve ~ normal(0, 1);            // Normal prior for evenness coefficient
    beta_fdiv ~ normal(0, 1);            // Normal prior for redundancy coefficient
    
    // Priors for the SDs of the random intercepts and regression coefficients
    sigma_a_fric ~ exponential(1.8);
    sigma_a_feve ~ exponential(1.8);
    sigma_a_fdiv ~ exponential(1.8);
    sigma_a_fric_season ~ exponential(1.8);
    sigma_a_feve_season ~ exponential(1.8);
    sigma_a_fdiv_season ~ exponential(1.8);
    
    // Priors for the observation errors
    sigma_fric ~ exponential(1);
    sigma_feve ~ exponential(1);
    sigma_fdiv ~ exponential(1);
    
    // Likelihood for the observations
    fric ~ normal(fric_hat, sigma_fric);
    feve ~ normal(feve_hat, sigma_feve);
    fdiv ~ normal(fdiv_hat, sigma_fdiv);
}

