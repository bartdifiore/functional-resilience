# Use the mFD package. 
 
## The basic workflow: 
   
   1. Include categorical traits in the trait matrix, specifically trophic groups
   2. estimate the PCA of the trait x species matix using Gowers distance
   3. Estimate CWM's for each trait and for each axis
   4. Model CWM traits (e.g. length-at-maturity) and PCA axes using sdmTMB
   5. Use those same PCA axes to compute functional diversity indices
   6. Include the functional diversity indices in a latent state SEM model. Include a more complicated SEM where Temperature ----> LS -----> Total ecosystem biomass | Total commercial biomass (e.g. metrics of ecosytem services)
