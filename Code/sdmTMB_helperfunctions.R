#sdmTMB helper functions

predict_mod <- function(sdmTMB_mod, type = "response"){
  polys <- read_sf("Data/Shapefiles/prediction_grid.shp") %>%
    rename(fid = FID)
  
  pred_grid <- read_sf("Data/Spatial/predictiongrid_wcovairates.shp") %>%
    rename(bt_seasonal = men_tmp, 
           total_landings_mt = fd_lnd_) %>% 
    mutate(X = as.numeric(st_coordinates(.)[,1])/1000, 
           Y = as.numeric(st_coordinates(.)[,2])/1000, 
           depth = mean(df$depth), 
           year = as.numeric(year)) %>% 
    st_drop_geometry() %>%
    filter(year %in% unique(df$year)) %>%
    drop_na(total_landings_mt) %>%
    mutate(season = toupper(season), 
           season_temp = ifelse(season == "SPRING", 25, 75), 
           year_season = as.numeric(paste(year, season_temp, sep = ".")))
  
  pred_grid$scaled_year <- (pred_grid$year - mean(pred_grid$year)) / 10
  pred_grid$scaled_year_season <- (pred_grid$year_season - mean(pred_grid$year_season)) / 10
  pred_grid$temp_scaled <- scale(pred_grid$bt_seasonal)
  pred_grid$log_landings_scaled <- scale(log(pred_grid$total_landings_mt))
  
  if (type == "response") {
    predictions <- predict(sdmTMB_mod, newdata = pred_grid, type = "response") %>% 
      left_join(polys) %>%
      st_as_sf()
  }
  
  if (type == "link") {
    predictions <- predict(sdmTMB_mod, newdata = pred_grid, type = "link") %>% 
      left_join(polys) %>%
      st_as_sf()
  }

  
  return(predictions)
}

plot_spatiotemporal <- function(sdmTMB_mod, years = c(1975, 1985, 1995, 2005, 2010, 2019), label, seasonal_model = F){
  
  if(seasonal_model == T){
  predict_mod(sdmTMB_mod) %>%
    filter(year %in% years) %>%
    ggplot()+
    geom_sf(data = ecodata::epu_sf, aes(geometry = geometry), fill = "transparent", color = "darkred", lwd = 0.25)+
    geom_sf(aes(fill = est))+
    scale_fill_viridis_c()+
    facet_grid(season~year)+
    labs(title = label, fill = label)+
    theme_bd()
  }
  
  if(seasonal_model == F){
    predict_mod(sdmTMB_mod) %>%
      filter(year %in% years) %>%
      ggplot()+
      geom_sf(data = ecodata::epu_sf, aes(geometry = geometry), fill = "transparent", color = "darkred", lwd = 0.25)+
      geom_sf(aes(fill = est))+
      scale_fill_viridis_c()+
      facet_wrap(~year, nrow = 1)+
      labs(title = label, fill = label)+
      theme_bd()
  }
}









# Function to compute slope, SE, and CI
get_slope_ci <- function(epu, season, coefs, vc) {
  terms <- "scaled_year"
  
  if (epu != "GB") {
    terms <- c(terms, paste0("scaled_year:EPU", epu))
  }
  
  if (season == "SPRING") {
    terms <- c(terms, "scaled_year:seasonSPRING")
    if (epu != "GB") {
      terms <- c(terms, paste0("scaled_year:EPU", epu, ":seasonSPRING"))
    }
  }
  
  idx <- match(terms, coefs$term)
  beta <- coefs$estimate[idx]
  beta[is.na(beta)] <- 0
  
  # Total estimate
  est <- sum(beta)
  
  # Variance
  var_total <- 0
  for (i in seq_along(idx)) {
    for (j in seq_along(idx)) {
      if (!is.na(idx[i]) && !is.na(idx[j])) {
        var_total <- var_total + vc[idx[i], idx[j]]
      }
    }
  }
  
  se <- sqrt(var_total)
  ci_low <- est - 1.96 * se
  ci_high <- est + 1.96 * se
  
  tibble(EPU = epu, season = season, slope = est, se = se,
         ci_low = ci_low, ci_high = ci_high)
}


plot_model_outputs <- function(model, return_plot = c("p1", "p2"), .label = "", family.log = F) {
  return_plot <- match.arg(return_plot)
  # .var <- rlang::enquo(.var)  
  
  # Prep coefficients and slope table
  coefs <- tidy(model, effects = "fixed", conf.int = TRUE)
  vc <- as.matrix(vcov(model))
  
  epus <- c("GOM", "GB", "MAB", "SS")
  seasons <- c("FALL", "SPRING")
  
  slope_table <- expand_grid(EPU = epus, season = seasons) %>%
    pmap_dfr(~ get_slope_ci(..1, ..2, coefs, vc))
  
  if (return_plot == "p1") {
    # Predict over grid
    
    if(family.log == F){
      predictions <- predict(model, newdata = pred_grid, type = "response") %>%
        left_join(polys) %>%
        st_as_sf()
      
      p1 <- predictions %>%
        filter(year %in% c(1970, 1990, 2010, 2024)) %>%
        ggplot() +
        geom_sf(aes(fill = est)) +
        scale_fill_viridis_c() +
        geom_sf(data = ecodata::epu_sf, aes(geometry = geometry),
                fill = "transparent", color = "darkred", lwd = 0.25) +
        facet_grid(season ~ year) +
        coord_sf(expand = FALSE) +
        scale_x_continuous(breaks = c(-68, -72, -76)) +
        labs(fill = .label) +
        theme_bd() +
        theme(
          # Places legend cleanly above the plotting frame next to the title
          legend.position = "top", 
          legend.direction = "horizontal",
          legend.justification = "left",
          legend.margin = margin(b = -10, t = 0), # Pulls the bar up closer to title
          strip.background.y = element_blank(),
          strip.text.y = element_blank()
        )
    }
    if(family.log == T){
      predictions <- predict(model, newdata = pred_grid, type = "link") %>%
        left_join(polys) %>%
        st_as_sf()
      
      p1 <- predictions %>%
        filter(year %in% c(1970, 1990, 2010, 2024)) %>%
        ggplot() +
        geom_sf(aes(fill = est)) +
        scale_fill_viridis_c(
          # Hardcodes the log breakpoints directly within your active 11.28 to 15.52 window
          breaks = log(c(1e5, 5e5, 3e6)),
          
          # Displays them cleanly on the response scale in scientific format
          labels = c("100K", "500K", "3M")
        )+
        geom_sf(data = ecodata::epu_sf, aes(geometry = geometry),
                fill = "transparent", color = "darkred", lwd = 0.25) +
        facet_grid(season ~ year) +
        coord_sf(expand = FALSE) +
        scale_x_continuous(breaks = c(-68, -72, -76)) +
        labs(fill = .label) +
        theme_bd() +
        theme(
          # Places legend cleanly above the plotting frame next to the title
          legend.position = "top", 
          legend.direction = "horizontal",
          legend.justification = "left",
          legend.margin = margin(b = -10, t = 0), # Pulls the bar up closer to title
          strip.background.y = element_blank(),
          strip.text.y = element_blank()
        )
    }
    
    return(p1)
  }
  
  if (return_plot == "p2") {
    # Filter slope table to valid slopes only
    valid_slopes <- slope_table %>%
      filter(ci_low > 0 | ci_high < 0) %>%
      mutate(include = TRUE)
    
    # Create EPU-level prediction grid
    ndata <- expand.grid(
      year = unique(df$year),
      season = unique(df$season),
      EPU = c("GOM", "GB", "MAB")
    ) %>%
      mutate(
        season_temp = ifelse(season == "SPRING", 25, 75),
        year_season = as.numeric(paste(year, season_temp, sep = "."))
      ) %>%
      left_join(valid_slopes %>% select(EPU, season, include), by = c("EPU", "season")) %>%
      filter(include == TRUE)
    
    ndata$scaled_year <- (ndata$year - mean(ndata$year)) / 10
    
    # Point estimates and uncertainty
    est <- predict(model, newdata = ndata, re_form = NA, se_fit = TRUE)
    
    # Simulated time series predictions
    predictions2 <- predict(model, newdata = pred_grid, nsim = 1000, type = "response")
    
    pred2 <- data.frame(predictions2) %>%
      bind_cols(year = pred_grid$year, season = pred_grid$season, EPU = pred_grid$EPU) %>%
      group_by(year, season, EPU) %>%
      pivot_longer(cols = starts_with("X")) %>%
      select(-name) %>%
      mutate(value = value) %>%
      tidybayes::mean_qi(value) %>%
      filter(EPU != "SS")
    
    if(family.log == F){
      p2 <- ggplot() +
        # geom_line(data = shelf, aes(x = year, y = !!.var), color = "black")+
        geom_line(data = pred2, aes(x = year, y = value, color = EPU), linewidth = 1) +
        geom_line(data = est, aes(x = year, y = est, color = EPU), linetype = 4, linewidth = 1) +
        geom_ribbon(data = est, aes(
          x = year,
          y = est,
          ymin = est - 1.96 * est_se,
          ymax = est + 1.96 * est_se,
          group = EPU
        ), alpha = 0.1) +
        scale_color_manual(values = c("GOM" = "#0072B2", "GB" = "#E69F00", "MAB" = "#009E73")) +
        scale_y_continuous(position = "left") +
        facet_wrap(~season, nrow = 2, strip.position = "right") +
        labs(y = .label, x = "", color = "") +
        theme_bd() +
        theme(
          legend.position = "top",
          # legend.position.inside = c(0.5, 0.03),
          # legend.direction = "horizontal"
        )
    }
    
    if(family.log == T){
      p2 <- ggplot() +
        # geom_line(data = shelf, aes(x = year, y = !!.var), color = "black")+
        geom_line(data = pred2, aes(x = year, y = value, color = EPU), linewidth = 1) +
        geom_line(data = est, aes(x = year, y = exp(est), color = EPU), linetype = 4, linewidth = 1) +
        geom_ribbon(data = est, aes(
          x = year,
          y = exp(est),
          ymin = exp(est - 1.96 * est_se),
          ymax = exp(est + 1.96 * est_se),
          group = EPU
        ), alpha = 0.1) +
        scale_color_manual(values = c("GOM" = "#0072B2", "GB" = "#E69F00", "MAB" = "#009E73")) +
        scale_y_continuous(position = "left") +
        facet_wrap(~season, nrow = 2, strip.position = "right") +
        labs(y = .label, x = "", color = "") +
        theme_bd() +
        theme(
          legend.position = "top",
          # legend.position.inside = c(0.5, 0.03),
          # legend.direction = "horizontal"
        )
    }
    
    return(p2)
  }
}

estimate_decadal_slopes <- function(model, ref_epu = "GB", ref_season = "FALL") {
  library(dplyr)
  
  # 1. Extract fixed effects and terms
  tidy_coefs <- sdmTMB::tidy(model, effects = "fixed")
  beta <- tidy_coefs$estimate
  names(beta) <- tidy_coefs$term
  
  # Extract unique levels from model terms to reconstruct the grid
  all_terms <- tidy_coefs$term
  
  # Find unique EPUs present in model terms (e.g., "EPUGOM" -> "GOM")
  epu_terms <- all_terms[grepl("^EPU", all_terms) & !grepl(":", all_terms)]
  epu_levels <- c(ref_epu, gsub("^EPU", "", epu_terms))
  
  # Find unique seasons present in model terms (e.g., "seasonSPRING" -> "SPRING")
  season_terms <- all_terms[grepl("^season", all_terms) & !grepl(":", all_terms)]
  season_levels <- c(ref_season, gsub("^season", "", season_terms))
  
  # 2. Extract variance-covariance matrix
  V <- vcov(model)
  
  # 3. Create minimal grid for contrasts (scaled_year at 1 vs 0)
  grid <- expand.grid(EPU = epu_levels, season = season_levels)
  grid_1 <- grid %>% mutate(scaled_year = 1)
  grid_0 <- grid %>% mutate(scaled_year = 0)
  
  # 4. Generate model matrices based on fixed formula
  form <- ~ scaled_year * EPU * season
  X1 <- model.matrix(form, data = grid_1)
  X0 <- model.matrix(form, data = grid_0)
  
  # 5. Build contrast matrix and align columns with beta
  L <- X1 - X0
  L <- L[, names(beta), drop = FALSE]
  
  # 6. Compute log-scale slopes and standard errors
  slope_estimates <- L %*% beta
  slope_variances <- diag(L %*% V %*% t(L))
  slope_se        <- sqrt(slope_variances)
  
  # 7. Compile results and calculate % change per decade
  results <- grid %>%
    mutate(
      log_slope      = as.vector(slope_estimates),
      log_std.error  = slope_se,
      log_conf.low   = log_slope - 1.96 * slope_se,
      log_conf.high  = log_slope + 1.96 * slope_se,
      
      # Transform to Percent Change per Decade
      pct_change      = (exp(log_slope) - 1) * 100,
      pct_conf_low    = (exp(log_conf.low) - 1) * 100,
      pct_conf_high   = (exp(log_conf.high) - 1) * 100
    )
  
  return(results)
}


