
T_mat_test = matrix(nrow = 3, ncol = 1, data = c(0.5,0.9,0.2))
S_mat_test = matrix(nrow = 5, ncol = 3, data = rpois(5*3, lambda = 1))

CWM <- function(T_mat, S_mat){
  # Data are collected at the level of the trawl. Here, I estimate the CWM for each trait, k, as the sum of the value for each trait multiplied by the biomass the species, such that
  # T_mat is the trait matrix with traits in columns and species in rows
  # S_mat is the species matrix with abundance/biomass/presence-absence for each species in columns and rows are sites (here individual tows). 
  
  W = apply(S_mat, 1, sum)
  
  (S_mat %*% T_mat) / W

    
}


CWM(T_mat_test, S_mat_test)



#5x3 * 3x1 = 5*1