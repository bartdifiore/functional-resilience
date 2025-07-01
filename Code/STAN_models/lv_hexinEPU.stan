data {
    // Observed variables
    int<lower=1> N;                       // Number of observations
    vector[N] fg_richness;                // Richness
    vector[N] fg_evennessSh;              // Evenness
    vector[N] fg_redundancy;              // Redundancy
    int<lower=1> K;                       // Number of sites
    int<lower=1, upper=K> group[N];       // Group indices (site-level)
    int<lower=1> T;                       // Number of years
    int<lower=1, upper=T> year[N];        // Year indices for each observation
    int<lower=1> M;                       // Number of EPUs
    int<lower=1, upper=M> epu[N];         // EPU indices for each site (grouped by EPU)
}

parameters {
    // Latent variables for each year and EPU
    matrix[M, T] nu;                      // Latent variables (one for each EPU and year)
    
    // Regression coefficients for each site
    real<lower=0> beta_fg_richness;                // Regression coefficient for richness
    real beta_fg_evennessSh;              // Regression coefficient for evenness
    real<lower=0> beta_fg_redundancy;              // Regression coefficient for redundancy
    
    // Intercepts for each site
    real a_fg_richness_raw;               // Intercept for richness (raw scale, for non-centered parameterization)
    real a_fg_evennessSh_raw;             // Intercept for evenness (raw scale)
    real a_fg_redundancy_raw;             // Intercept for redundancy (raw scale)
    
    // Non-centered random effects (raw values for intercepts)
    vector[K] z_a_fg_richness;           // Non-centered random intercepts for richness
    vector[K] z_a_fg_evennessSh;         // Non-centered random intercepts for evenness
    vector[K] z_a_fg_redundancy;         // Non-centered random intercepts for redundancy
    
    // SDs for the random intercepts and coefficients
    real<lower=0> sigma_a_fg_richness;   // SD for richness intercept
    real<lower=0> sigma_a_fg_evennessSh; // SD for evenness intercept
    real<lower=0> sigma_a_fg_redundancy; // SD for redundancy intercept
    real<lower=0> sigma_z_a_fg_richness;   // SD for richness intercept random effects
    real<lower=0> sigma_z_a_fg_evennessSh; // SD for evenness intercept random effects
    real<lower=0> sigma_z_a_fg_redundancy; // SD for redundancy intercept random effects

    // SDs for the observations
    real<lower=0> sigma_fg_richness;     // SD for richness
    real<lower=0> sigma_fg_evennessSh;   // SD for evenness
    real<lower=0> sigma_fg_redundancy;   // SD for redundancy
}

transformed parameters {
    // Predicted values for each observation
    vector[N] fg_richness_hat;
    vector[N] fg_evennessSh_hat;
    vector[N] fg_redundancy_hat;
    
    // Transformed intercepts using non-centered random effects
    vector[K] a_fg_richness = a_fg_richness_raw + sigma_z_a_fg_richness * z_a_fg_richness;
    vector[K] a_fg_evennessSh = a_fg_evennessSh_raw + sigma_z_a_fg_evennessSh * z_a_fg_evennessSh;
    vector[K] a_fg_redundancy = a_fg_redundancy_raw + sigma_z_a_fg_redundancy * z_a_fg_redundancy;

    // Latent variable regression predictions (for each observation)
    for (n in 1:N) {
        fg_richness_hat[n] = beta_fg_richness * nu[epu[n], year[n]] + a_fg_richness[group[n]];
        fg_evennessSh_hat[n] = beta_fg_evennessSh * nu[epu[n], year[n]] + a_fg_evennessSh[group[n]];
        fg_redundancy_hat[n] = beta_fg_redundancy * nu[epu[n], year[n]] + a_fg_redundancy[group[n]];
    }
}

model {
    // Priors for latent variable nu (standard normal prior for each EPU and year)
    for (m in 1:M) {
        for (t in 1:T) {
            nu[m, t] ~ normal(0, 1);  // Latent variable for each EPU and year
        }
    }

    // Non-centered random intercepts for each site
    z_a_fg_richness ~ normal(0, sigma_z_a_fg_richness); // Prior for non-centered random intercepts for richness
    z_a_fg_evennessSh ~ normal(0, sigma_z_a_fg_evennessSh); // Prior for non-centered random intercepts for evenness
    z_a_fg_redundancy ~ normal(0, sigma_z_a_fg_redundancy); // Prior for non-centered random intercepts for redundancy

    // Priors for raw intercepts (before applying non-centered parameterization)
    a_fg_richness_raw ~ normal(0, sigma_a_fg_richness);   // Mean intercept for richness
    a_fg_evennessSh_raw ~ normal(0, sigma_a_fg_evennessSh); // Mean intercept for evenness
    a_fg_redundancy_raw ~ normal(0, sigma_a_fg_redundancy); // Mean intercept for redundancy

    // Priors for regression coefficients (betas)
    beta_fg_richness ~ gamma(2,1);   // Positive-only prior for richness coefficient
    beta_fg_evennessSh ~ normal(0,1);     // Normal prior for evenness coefficient
    beta_fg_redundancy ~ gamma(2, 1);     // Normal prior for redundancy coefficient

    // Priors for the SDs of the random intercepts and regression coefficients
    sigma_a_fg_richness ~ exponential(1);
    sigma_a_fg_evennessSh ~ exponential(1);
    sigma_a_fg_redundancy ~ exponential(1);

    // Priors for the observation errors
    sigma_fg_richness ~ exponential(1);
    sigma_fg_evennessSh ~ exponential(1);
    sigma_fg_redundancy ~ exponential(1);

    // Likelihood for the observations
    fg_richness ~ normal(fg_richness_hat, sigma_fg_richness);
    fg_evennessSh ~ normal(fg_evennessSh_hat, sigma_fg_evennessSh);
    fg_redundancy ~ normal(fg_redundancy_hat, sigma_fg_redundancy);
}

