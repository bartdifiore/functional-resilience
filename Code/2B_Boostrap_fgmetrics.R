
df.test <- df %>% 
  mutate(id = paste(hex_id, five_year_period, sep = "-")) %>%
  select(id, EPU, abundance, svspp, haul_id) %>%
  pivot_wider(names_from = "svspp", values_from = "abundance", values_fill = 0)

id <- unique(df.test$id)


# bootstrap_richness <- function(id){
#   # Step 1. Find a SAC, and M-M fit for each hex-year combination. Save output to a dataframe. Filter out hex-year combinations where the number of observed species did NOT reach 50% of the asymptotic species richness (defined by Vm). 
  out <- list()
  for(i in 1:length(id)){
    temp = df.test[df.test$id == id[963], ]
    mat = as.matrix(temp[,-c(1:3)])
    spec_random = specaccum(mat, method = "random", permutations = 1000)
    fit_models = fitspecaccum(spec_random, model = "michaelis-menten")
    out[[i]] = c(id = id[i], 
                 Vm = median(fit_models$coefficients[1,]), 
                 K = median(fit_models$coefficients[2,])) # This value represents the number of hauls needed to observe 50% of the asymptotic species richness.
  }
#   
#   filtered_survey_mat <- survey_mat %>%
#     left_join(out) %>%
#     filter()
#   
#   # Step 2. Sample the survey matrix according to a random number of hauls in each remaining hex-year. In other words for each run, pull a random number from 1 to the maximum number of hauls in a hex-year, then sample that number of hauls from the survey matrix. 
#   
#   # Step 3. Run the calculations of the 4 resilience indicators for each run. Save the output. 
#   
#   # Step 4. Iterate the processes at least 100 times. 
# }
  
  
  
  
  
  out <- list()
  
  # Define a vector of models to try
  models_to_try <- c("michaelis-menten","arrhenius", "gleason", "gitay", "lomolino", "asymp", "gompertz", "logis", "weibull")  # Add more models as needed
  
  for (i in 1:length(id)) {
    # Use tryCatch to handle errors
    result <- tryCatch({
      # Subset the data
      temp <- df.test[df.test$id == id[i], ]
      mat <- as.matrix(temp[,-c(1:3)])  # Remove first 3 columns (id, and possibly other columns)
      
      # Run specaccum and attempt fitting with different models
      spec_random <- specaccum(mat, method = "random", permutations = 1000)
      
      # Try different models until one works
      fit_models <- NULL
      selected_model <- NULL  # Store the selected model name
      for (model in models_to_try) {
        fit_models <- tryCatch({
          fitspecaccum(spec_random, model = model, control = nls.control(maxiter = 100, minFactor = 1e-9))
        }, error = function(e) {
          message(paste("Error in model", model, "for iteration", i, ":", e$message))
          return(NULL)  # If model fails, try the next one
        })
        
        # If fitting was successful, store the model and break out of the loop
        if (!is.null(fit_models)) {
          selected_model <- model  # Store the name of the selected model
          break
        }
      }
      
      # If no model was successfully fitted, stop the process for this iteration
      if (is.null(fit_models)) {
        stop("All models failed for iteration ", i)  # Stop if all models fail
      }
      
      # Store results including the selected model
      c(id = id[i],
        p1 = as.numeric(median(fit_models$coefficients[1,])),
        p2 = as.numeric(median(fit_models$coefficients[2,])),
        selected_model = selected_model, 
        n_hauls = as.numeric(dim(temp)[1]))  # Add the selected model
      
    }, error = function(e) {
      # Handle error in the tryCatch block: print message and return NULL or NA for the iteration
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  # You can return NA or any other value if you prefer
    })
    
    # If no error occurred, save the result
    if (!is.null(result)) {
      out[[i]] <- result
    } else {
      out[[i]] <- NA  # Optional: store NA for failed iterations
    }
  }
  
  
  