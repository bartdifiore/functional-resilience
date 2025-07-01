#sdmTMB helper functions

predict_mod <- function(sdmTMB_mod){
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
  
  predictions <- predict(sdmTMB_mod, newdata = pred_grid, type = "response") %>% 
    left_join(polys) %>%
    st_as_sf()
  
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
        labs(fill = .label, title = .label) +
        theme_bd() +
        theme(
          legend.position = "inside",
          legend.position.inside = c(0.925, 0.1),
          legend.direction = "horizontal",
          legend.title = element_blank(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          legend.background = element_blank()
        )
    
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






