library(tidyverse)
library(sdmTMB) # This is the workhorse to do the spatio-temporal modeling
library(sf)
source("Code/theme.R")
source("Code/CWM_functions.R")
source("Code/sdmTMB_helperfunctions.R")

#-----------------------------
## Read in the data
#-----------------------------

year_filter <- c(seq(1963, 1969, by = 1), 2017, 2020)

meta <- read.csv("Data/Derived/trawlmeta_wBT.csv")

epu_strata <- read.csv("Data/Derived/strata_to_epu.csv")

df <- read.csv("Data/Derived/CWM_trawl_dataset.csv") %>% 
  as_tibble() %>% 
  left_join(meta) %>%
  drop_na(lat, lon) %>%
  add_utm_columns(ll_names = c("lon", "lat"), units = "km") %>%
  mutate(season = as.factor(season)) %>% 
  ungroup() %>%
  mutate(season_temp = ifelse(season == "SPRING", 25, 75), 
         year_season = as.numeric(paste(year, season_temp, sep = "."))) %>%
  left_join(epu_strata) %>%
  drop_na(EPU) %>%
  mutate(EPU = as.factor(EPU), 
         season = as.factor(season))

sf_use_s2(FALSE)

pj_crs <- st_crs(x = "EPSG:32619")

df$scaled_year <- (df$year - mean(df$year))/10 
df$scaled_year_season <- (df$year_season - mean(df$year_season))/10

df %>% 
  group_by(year, season, EPU) %>%
  summarize(mean_bt = mean(bottom_temp, na.rm = T)) %>%
  ggplot(aes(x = year, y = mean_bt))+
  geom_line(aes(color = season)) +
  facet_wrap(~EPU)

#----------------------------------------------------
## Model
#----------------------------------------------------

# Make the mesh
mesh <- sdmTMB::make_mesh(df, xy_cols = c("X", "Y"), cutoff = 30)

plot(mesh)

m_playtime <- sdmTMB(
  data = df,
  formula = length_maturity ~ bottom_temp*EPU*season,
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year_season",
  time_varying = ~ 1,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017.25, 2017.75, 2020.25, 2020.75),
  silent = F
)

m_playtime
summary(m_playtime)
sanity(m_playtime)

dat <- ggeffects::ggeffect(m_playtime, terms = c("bottom_temp", "season", "EPU[GB, MAB, GOM]"))
plot(dat)

m_playtime2 <- sdmTMB(
  data = df,
  formula = length_maturity ~ 0 + scale(bottom_temp),
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year",
  time_varying = ~ 1 + EPU*season,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)

write_rds(m_playtime2, "Data/Derived/Model_output/m_playtime2.rds", compress = "gz")
m_playtime2 <- readRDS("Data/Derived/Model_output/m_playtime2.rds")

m_playtime2
sanity(m_playtime2)

dat <- ggeffects::ggeffect(m_playtime2, terms = c("bottom_temp"))
plot(dat)

bt_grid <- readRDS("Data/Derived/prediction_grid_wBT.rds")
  

epus <- ecodata::epu_sf %>% 
  select(EPU, geometry)%>%
  st_transform(crs = pj_crs)

pred_grid <- read_sf("Data/Shapefiles/prediction_grid.shp") %>%
  st_centroid() %>% 
  st_join(epus) %>%
  mutate(X = as.numeric(st_coordinates(.)[,1])/1000, 
         Y = as.numeric(st_coordinates(.)[,2])/1000) %>%
  st_drop_geometry() %>%
  drop_na(EPU) %>%
  replicate_df(time_name = "year", time_values = 1970:2024) %>% 
  expand_grid(season = c("FALL", "SPRING")) %>%
  # filter(!year %in% year_filter) %>%
  rename(fid = FID) %>% 
  left_join(bt_grid) %>%
  mutate(EPU = as.factor(EPU), 
         season = as.factor(season), 
         bottom_temp = replace_na(bottom_temp, mean(bottom_temp, na.rm = T)))


# Extract the index TMB object
preds_spatial <- predict(m_playtime2, newdata = pred_grid)


preds_spatial %>% 
  group_by(year, season, EPU) %>%
  summarize(mean_non_rf = mean(est_non_rf), 
            mean_bt = mean(bottom_temp)) %>%
  ggplot(aes(x=mean_bt,y= exp(mean_non_rf)))+
  geom_line(aes(color = season))+
  facet_wrap(~EPU)


index_df <- sdmTMB::get_index_split(m_playtime2, newdata = pred_grid, bias_correct = F, area = 1, nsplit = 2)

# Divide the output by the total number of cells per year/group to get the true mean trait value
cells_per_year <- pred_grid %>% 
  group_by(year, EPU, season) %>% 
  summarise(n_cells = n(), .groups = "drop")

mean_trait_df <- index_df %>%
  mutate(
    mean_est = est / 140,
    mean_lwr = lwr / 140,
    mean_upr = upr / 140
  )

head(mean_trait_df)

mean_trait_df %>% 
  ggplot(aes(x = year, y = mean_est))+
  geom_line()
















# 1. Create a minimal grid (no spatial coordinates needed for just the TV effects)
nd <- expand.grid(
  year = unique(c(df$year, 2017, 2020)),
  EPU = unique(df$EPU),
  season = unique(df$season),
  bottom_temp = 0 # Held at mean
)

# 2. Predict link-space values
# est_tv isolates the time_varying component specifically
preds <- predict(m_playtime2, newdata = nd, re_form = NA) 

# 3. Plot est_tv directly
ggplot(preds, aes(x = year, y = exp(est), color = season)) +
  geom_line() + 
  facet_wrap(~EPU) +
  labs(y = "Time-Varying Intercept")


nd <- df %>%
  group_by(season, EPU) %>%
  reframe(
    expand_grid(
      bottom_temp = seq(min(bottom_temp, na.rm = T), max(bottom_temp, na.rm = T), length.out = 100), 
      year = seq(1970, 2024, by = 1)
    )
  )

# 2. Predict link-space values
# est_tv isolates the time_varying component specifically
preds <- predict(m_playtime2, newdata = nd, re_form = NA) 

# 3. Plot est_tv directly

preds %>%
  filter(year %in% c(1970, 1980, 1990, 2000, 2010, 2020, 2024)) %>%
ggplot(aes(x = bottom_temp, y = exp(est))) +
  geom_line(aes(color=as.factor(year))) +
  facet_grid(season~EPU) +
  labs(y = "Time-Varying Intercept")


nd <- df %>% 
  group_by(year, EPU, season) %>% 
  summarize(bottom_temp = mean(bottom_temp, na.rm = T)) %>%
  drop_na(bottom_temp)

preds <- predict(m_playtime2, newdata = nd, re_form = NA) 

ggplot(preds, aes(x = year, y = exp(est)))+
  geom_line(aes(color = season))+
  facet_wrap(season~EPU, scales = "free")



m_playtime3 <- sdmTMB(
  data = df,
  formula = length_maturity ~ 1,
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year",
  time_varying = ~ 0 + EPU*season,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)

m_playtime3
sanity(m_playtime3)

# 1. Create a minimal grid (no spatial coordinates needed for just the TV effects)
nd <- expand.grid(
  year = unique(c(df$year, 2017, 2020)),
  EPU = unique(df$EPU),
  season = unique(df$season),
  bottom_temp = 0 # Held at mean
)

# 2. Predict link-space values
# est_tv isolates the time_varying component specifically
preds <- predict(m_playtime3, newdata = nd, re_form = NA) 

# 3. Plot est_tv directly
ggplot(preds, aes(x = year, y = exp(est), color = season)) +
  geom_line() + 
  facet_wrap(~EPU) +
  labs(y = "Time-Varying Intercept")




m_playtime4 <- sdmTMB(
  data = df,
  formula = length_maturity ~ scaled_year*EPU*season + scale(bottom_temp),
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year",
  time_varying = ~ 0 + EPU*season,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)
sanity(m_playtime4)




dat <- ggeffects::ggeffect(m_playtime4, terms = c("bottom_temp", "scaled_year"))
plot(dat)


nd <- df %>% 
  group_by(year, EPU, season) %>% 
  summarize(bottom_temp = mean(bottom_temp, na.rm = T)) %>%
  drop_na(bottom_temp)

nd$scaled_year <- (nd$year - mean(nd$year))/10

preds <- predict(m_playtime4, newdata = nd, re_form = NA) 

preds %>%
  filter(EPU != "SS") %>%
  ggplot(aes(x = year, y = exp(est)))+
  geom_line(aes(color = EPU))+
  facet_wrap(~season)



#-------------------------------------------
## Two temp scales
#-------------------------------------------

df_rolling <- df %>%
  group_by(EPU, season, year) %>%
  summarize(annual_mean_bt = mean(bottom_temp, na.rm = T)) %>%
  ungroup() %>%
  arrange(EPU, season, year) %>%
  group_by(EPU, season) %>%
  mutate(climate_trend = zoo::rollapply(annual_mean_bt, width = 7, FUN = mean, na.rm = T, align = "center", fill = "extend")) #%>%
  # ungroup() %>%
  # mutate(temp_anomaly = bottom_temp - climate_trend)

df2 <- df %>% 
  left_join(df_rolling) %>%
  ungroup() %>%
  mutate(temp_anomaly = bottom_temp - climate_trend)

df_rolling %>%
  select(year, EPU, season, climate_trend) %>%
  distinct(year, EPU, season, climate_trend) %>%
  ggplot()+
  geom_line(aes(y = climate_trend, x = year, color = EPU))+
  facet_wrap(~season)


m_two_scales <- sdmTMB(
  data = df2,
  formula = length_maturity ~ climate_trend + temp_anomaly,
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year",
  time_varying = ~ EPU*season,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)
sanity(m_two_scales)

tidy(m_two_scales, effects = "fixed")

dat <- ggeffects::ggeffect(m_two_scales, terms = c("climate_trend", "temp_anomaly"))
plot(dat)


nd <- df2 %>% 
  select(year, season, EPU, climate_trend, temp_anomaly) %>%
  drop_na(temp_anomaly) %>%
  filter(EPU != "SS")
  
preds <- predict(m_two_scales, newdata = nd, re_form = NA) 

preds %>%
  select(year, season, EPU, climate_trend, est) %>%
  group_by(year, season, EPU) %>%
  summarize(mean_est = mean(est)) %>%
  ggplot()+
  geom_line(aes(x = year, y = exp(mean_est), color = EPU))+
  facet_wrap(~season)

ggplot(preds)+
  geom_line(aes(x = temp_anomaly, y = exp(est), color = EPU))+
  facet_grid(EPU~season)

nd <- expand.grid(year = unique(df2$year),
                  season = unique(df2$season), 
                  EPU = sort(unique(df2$EPU))[1:3])
  


m_two_scales2 <- sdmTMB(
  data = df2,
  formula = length_maturity ~ EPU * season + EPU * season : (climate_trend + temp_anomaly),
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year",
  time_varying = ~ EPU*season,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)
sanity(m_two_scales2)
tidy(m_two_scales2, effects = "fixed")


library(tidyverse)

# 1. Generate the tidy dataframe directly from your model object
# (Replace 'm_two_scales2' with your actual sdmTMB model variable name)
coef_raw <- tidy(m_two_scales2, effects = "fixed", conf.int = TRUE)

# 2. Automatically parse and clean the terms for the heatmap
df_heatmap <- coef_raw %>%
  # Filter to look only at the temperature slope rows
  filter(str_detect(term, "climate_trend|temp_anomaly")) %>%
  mutate(
    # Extract EPU: Find GOM, MAB, or SS. If none are found, it's the reference baseline (GB)
    EPU = case_when(
      str_detect(term, "GOM") ~ "GOM",
      str_detect(term, "MAB") ~ "MAB",
      str_detect(term, "SS")  ~ "SS",
      TRUE                    ~ "GB"
    ),
    # Extract Season: Look for SPRING or FALL
    Season = case_when(
      str_detect(term, "SPRING") ~ "Spring",
      str_detect(term, "FALL")   ~ "Fall"
    ),
    # Extract Temperature Scale Type
    Scale = case_when(
      str_detect(term, "climate_trend") ~ "Chronic (Climate Trend)",
      str_detect(term, "temp_anomaly")  ~ "Acute (Temp Anomaly)"
    ),
    # Identify significance using the model's computed confidence intervals
    is_significant = if_else(conf.low > 0 | conf.high < 0, TRUE, FALSE),
    # Format text label for display on the tiles (Percentage change)
    label_text = paste0(round(estimate * 100, 1), "%", if_else(is_significant, "*", ""))
  )

# 3. Generate the Heatmap Plot
ggplot(df_heatmap, aes(x = Season, y = EPU, fill = estimate)) +
  geom_tile(color = "white", lwd = 0.8) +
  geom_text(aes(label = label_text), color = "black", fontface = "bold", size = 4) +
  scale_fill_gradient2(
    high = "#2166ac", mid = "#f7f7f7", low = "#b2182b", 
    midpoint = 0, name = "Log Effect Slope\n(Per 1°C)"
  ) +
  facet_wrap(~Scale) +
  labs(
    title = "Community Trait Sensitivity to Temperature Scales",
    subtitle = "Values display % change in length-at-maturity CWM per 1°C (* indicates p < 0.05)",
    x = "Season",
    y = "Ecosystem Production Unit (EPU)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "transparent"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "right"
  )

library(tidyverse)

# 1. Generate predictions from your sdmTMB model
# We isolate the fixed effects + the AR1 time-varying process components,
# while deliberately dropping the localized spatiotemporal random field noise.
bt_grid <- readRDS("Data/Derived/prediction_grid_wBT.rds")


epus <- ecodata::epu_sf %>% 
  select(EPU, geometry)%>%
  st_transform(crs = pj_crs)

pred_grid <- read_sf("Data/Shapefiles/prediction_grid.shp") %>%
  st_centroid() %>% 
  st_join(epus) %>%
  mutate(X = as.numeric(st_coordinates(.)[,1])/1000, 
         Y = as.numeric(st_coordinates(.)[,2])/1000) %>%
  st_drop_geometry() %>%
  drop_na(EPU) %>%
  replicate_df(time_name = "year", time_values = 1970:2024) %>% 
  expand_grid(season = c("FALL", "SPRING")) %>%
  # filter(!year %in% year_filter) %>%
  rename(fid = FID) %>% 
  left_join(bt_grid) %>%
  mutate(EPU = as.factor(EPU), 
         season = as.factor(season), 
         bottom_temp = replace_na(bottom_temp, mean(bottom_temp, na.rm = T)))

pred_grid_rolling <- pred_grid %>%
  group_by(EPU, season, year) %>%
  summarize(annual_mean_bt = mean(bottom_temp, na.rm = T)) %>%
  ungroup() %>%
  arrange(EPU, season, year) %>%
  group_by(EPU, season) %>%
  mutate(climate_trend = zoo::rollapply(annual_mean_bt, width = 7, FUN = mean, na.rm = T, align = "center", fill = "extend"))
  
pred_grid2 <- pred_grid %>% 
  left_join(pred_grid_rolling) %>%
  ungroup() %>%
  mutate(temp_anomaly = bottom_temp - climate_trend)

# Extract the index TMB object
preds_spatial <- predict(m_two_scales2, newdata = pred_grid2)






library(tidyverse)

# 1. Build the clean, regional-level prediction dataframe
# We take the aggregated grid rolling averages so 'climate_trend' is accurate for each EPU/season/year,
# and we explicitly set 'temp_anomaly' to 0 to eliminate short-term temperature variance.
regional_baseline_frame <- pred_grid_rolling %>%
  select(year, EPU, season, climate_trend) %>%
  distinct() %>%
  mutate(temp_anomaly = 0)

# 2. Predict from the sdmTMB model
# We set re_form = NA to instruct sdmTMB to ignore spatial and spatiotemporal random fields,
# allowing it to run smoothly on a non-spatial dataframe.
preds_baseline <- predict(
  m_two_scales2, 
  newdata = regional_baseline_frame, 
  re_form = NA
)

# 3. Exponentiate the link-scale estimate and plot the counterfactual time series
preds_baseline %>%
  mutate(predicted_cwm = exp(est)) %>%
  ggplot(aes(x = year, y = predicted_cwm, color = season)) +
  # Reconstructed baseline trendlines
  geom_line(lwd = 1.3) +
  geom_point(size = 2) +
  # Face-wrap by EPU to handle varying community scales cleanly
  facet_wrap(~EPU, scales = "free_y") +
  scale_color_manual(
    values = c("SPRING" = "#1f78b4", "FALL" = "#ff7f00"),
    labels = c("SPRING" = "Spring", "FALL" = "Fall"),
    name = "Survey Season"
  ) +
  labs(
    title = "Model-Predicted Trait Dynamics (Counterfactual Baseline Trend)",
    subtitle = "Evaluated at the chronic climate trend line with acute temp anomalies set to 0",
    x = "Year",
    y = "Predicted CWM Length-at-Maturity"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "transparent"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    legend.position = "bottom"
  )




library(tidyverse)

# 1. Create a dense grid of climate trends for every single actual year
# This satisfies the time requirement while matching your observed data ranges.
df_marginal_with_time <- df2 %>%
  group_by(EPU, season) %>%
  summarize(
    min_trend = min(climate_trend, na.rm = TRUE),
    max_trend = max(climate_trend, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(EPU, season) %>%
  reframe(
    climate_trend = seq(min_trend, max_trend, length.out = 100)
  ) %>%
  # Cross-join with every unique year from your dataset to satisfy the AR1 time parameter
  expand_grid(year = unique(df2$year)) %>%
  mutate(temp_anomaly = 0) # Lock anomaly at 0

# 2. Predict from your sdmTMB model (ignoring spatial fields)
preds_raw <- predict(
  m_two_scales2, 
  newdata = df_marginal_with_time, 
  re_form = NA
)

# 3. Average out the year effect (AR1 drift) to isolate the pure marginal fixed slopes
preds_marginal <- preds_raw %>%
  group_by(EPU, season, climate_trend) %>%
  summarize(
    # Average on the log scale first, then transform to response scale
    predicted_cwm = exp(mean(est, na.rm = TRUE)),
    .groups = "drop"
  )

# 4. Plot the conditional effects
ggplot(preds_marginal, aes(x = climate_trend, y = predicted_cwm, color = EPU)) +
  geom_line(lwd = 1.5) +
  facet_grid(~season, scales = "free") + 
  # scale_color_manual(
  #   values = c("SPRING" = "#1f78b4", "FALL" = "#ff7f00"),
  #   labels = c("SPRING" = "Spring", "FALL" = "Fall"),
  #   name = "Survey Season"
  # ) +
  labs(
    title = "Conditional Effect of Chronic Climate Warming on Community Traits",
    subtitle = "Averaged across AR1 year intercepts with acute temp anomalies set to 0",
    x = "Chronic Climate Trend (7-Year Bottom Temp Average, °C)",
    y = "Predicted CWM Length-at-Maturity"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "transparent"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    legend.position = "bottom"
  )


dat <- ggeffects::ggeffect(m_two_scales2, terms = c("climate_trend", "EPU[GOM, MAB, GB]", "season"))
plot(dat)

#----------------------------------
## Fecundity
#----------------------------------


m_two_scales_fecundity <- sdmTMB(
  data = df2,
  formula = fecundity ~ EPU * season + EPU * season : (climate_trend + temp_anomaly),
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year",
  time_varying = ~ 0 + EPU*season,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)
sanity(m_two_scales_fecundity)
tidy(m_two_scales_fecundity, effects = "fixed")


dat <- ggeffects::ggeffect(m_two_scales_fecundity, terms = c("climate_trend", "EPU[GOM, MAB, GB]", "season"))
plot(dat)

# 1. Generate the tidy dataframe directly from your model object
# (Replace 'm_two_scales2' with your actual sdmTMB model variable name)
coef_raw <- tidy(m_two_scales_fecundity, effects = "fixed", conf.int = TRUE)

# 2. Automatically parse and clean the terms for the heatmap
df_heatmap <- coef_raw %>%
  # Filter to look only at the temperature slope rows
  filter(str_detect(term, "climate_trend|temp_anomaly")) %>%
  mutate(
    # Extract EPU: Find GOM, MAB, or SS. If none are found, it's the reference baseline (GB)
    EPU = case_when(
      str_detect(term, "GOM") ~ "GOM",
      str_detect(term, "MAB") ~ "MAB",
      str_detect(term, "SS")  ~ "SS",
      TRUE                    ~ "GB"
    ),
    # Extract Season: Look for SPRING or FALL
    Season = case_when(
      str_detect(term, "SPRING") ~ "Spring",
      str_detect(term, "FALL")   ~ "Fall"
    ),
    # Extract Temperature Scale Type
    Scale = case_when(
      str_detect(term, "climate_trend") ~ "Chronic (Climate Trend)",
      str_detect(term, "temp_anomaly")  ~ "Acute (Temp Anomaly)"
    ),
    # Identify significance using the model's computed confidence intervals
    is_significant = if_else(conf.low > 0 | conf.high < 0, TRUE, FALSE),
    # Format text label for display on the tiles (Percentage change)
    label_text = paste0(round(estimate * 100, 1), "%", if_else(is_significant, "*", ""))
  )

# 3. Generate the Heatmap Plot
ggplot(df_heatmap, aes(x = Season, y = EPU, fill = estimate)) +
  geom_tile(color = "white", lwd = 0.8) +
  geom_text(aes(label = label_text), color = "black", fontface = "bold", size = 4) +
  scale_fill_gradient2(
    high = "#2166ac", mid = "#f7f7f7", low = "#b2182b", 
    midpoint = 0, name = "Log Effect Slope\n(Per 1°C)"
  ) +
  facet_wrap(~Scale) +
  labs(
    title = "Community Trait Sensitivity to Temperature Scales",
    subtitle = "Values display % change in length-at-maturity CWM per 1°C (* indicates p < 0.05)",
    x = "Season",
    y = "Ecosystem Production Unit (EPU)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "transparent"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "right"
  )


preds_baseline <- predict(
  m_two_scales_fecundity, 
  newdata = regional_baseline_frame, 
  re_form = NA
)

# 3. Exponentiate the link-scale estimate and plot the counterfactual time series
preds_baseline %>%
  mutate(predicted_cwm = exp(est)) %>%
  ggplot(aes(x = year, y = predicted_cwm, color = season)) +
  # Reconstructed baseline trendlines
  geom_line(lwd = 1.3) +
  geom_point(size = 2) +
  # Face-wrap by EPU to handle varying community scales cleanly
  facet_wrap(~EPU, scales = "free_y") +
  scale_color_manual(
    values = c("SPRING" = "#1f78b4", "FALL" = "#ff7f00"),
    labels = c("SPRING" = "Spring", "FALL" = "Fall"),
    name = "Survey Season"
  ) +
  labs(
    title = "Model-Predicted Trait Dynamics (Counterfactual Baseline Trend)",
    subtitle = "Evaluated at the chronic climate trend line with acute temp anomalies set to 0",
    x = "Year",
    y = "Predicted CWM Length-at-Maturity"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "transparent"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    legend.position = "bottom"
  )






#----------------------------------
## L infinity
#----------------------------------

m_two_scales_linf <- sdmTMB(
  data = df2,
  formula = l_inf ~ EPU * season + EPU * season : (climate_trend + temp_anomaly),
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year",
  time_varying = ~ 0 + EPU*season,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)
sanity(m_two_scales_linf)
tidy(m_two_scales_linf, effects = "fixed")


dat <- ggeffects::ggeffect(m_two_scales_linf, terms = c("climate_trend", "EPU[GOM, MAB, GB]", "season"))
plot(dat)


#----------------------------------
## k
#----------------------------------

m_two_scales_k <- sdmTMB(
  data = df2,
  formula = k ~ EPU * season + EPU * season : (climate_trend + temp_anomaly),
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year",
  time_varying = ~ 0 + EPU*season,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)
sanity(m_two_scales_k)
tidy(m_two_scales_k, effects = "fixed")

library(emmeans)

# 1. Slopes for climate_trend across EPU and season
climate_slopes <- emtrends(m_two_scales_k, 
                           specs = ~ EPU * season, 
                           var = "climate_trend")

# 2. Slopes for temp_anomaly across EPU and season
temp_slopes <- emtrends(m_two_scales_k, 
                        specs = ~ EPU * season, 
                        var = "temp_anomaly")

# View the estimates
summary(climate_slopes)
summary(temp_slopes)

dat <- ggeffects::ggeffect(m_two_scales_k, terms = c("climate_trend", "EPU[GOM, MAB, GB]", "season"))
plot(dat)

temp <- ggeffects::predict_response(m_two_scales_k, terms = c("EPU[GOM, MAB, GB]", "season"))

plot(temp)
