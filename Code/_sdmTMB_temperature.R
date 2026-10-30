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
         season = as.factor(season)) %>%
  ungroup() %>%
  group_by(EPU, season, year) %>%
  mutate(annual_mean_bt = mean(bottom_temp, na.rm = T)) %>%
  ungroup() %>%
  mutate(temp_anomaly = bottom_temp - annual_mean_bt)

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

# Make the mesh
mesh <- sdmTMB::make_mesh(df, xy_cols = c("X", "Y"), cutoff = 30)

plot(mesh)


df_rolling <- df %>%
  group_by(EPU, season, year) %>%
  summarize(annual_mean_bt = mean(bottom_temp, na.rm = T)) %>%
  ungroup() %>%
  arrange(EPU, season, year) %>%
  group_by(EPU, season) %>%
  mutate(
    climate_trend = zoo::rollapply(annual_mean_bt, width = 5, FUN = mean, na.rm = T, align = "center", fill = "extend"),
    # Calculate the anomaly purely at the YEARLY scale
    annual_temp_anomaly = annual_mean_bt - climate_trend 
  )

df2 <- df %>% 
  left_join(df_rolling, by = c("EPU", "season", "year")) %>%
  ungroup()
# Do NOT subtract individual bottom_temp here. 
# Use annual_temp_anomaly in your sdmTMB formula instead.

df_rolling %>%
  select(year, EPU, season, climate_trend) %>%
  distinct(year, EPU, season, climate_trend) %>%
  ggplot()+
  geom_line(aes(y = climate_trend, x = year, color = EPU))+
  facet_wrap(~season)

df %>%
  group_by(EPU, season, year) %>%
  mutate(annual_mean_bt = mean(bottom_temp, na.rm = T)) %>%
  ggplot()+
  geom_line(aes(x = year, y = annual_mean_bt, color = EPU))+
  facet_wrap(~season)


#----------------------------------------------------
## Model: length-at-maturity
#----------------------------------------------------



plot(mesh)

df_lm <- df2 %>%
  drop_na(length_maturity, climate_trend, annual_temp_anomaly)

mesh <- sdmTMB::make_mesh(df_lm, xy_cols = c("X", "Y"), cutoff = 30)

m_two_scales_lm <- sdmTMB(
  data = df_lm,
  formula = length_maturity ~ EPU * season *(climate_trend + annual_temp_anomaly),
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year",
  time_varying = ~  1,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)
sanity(m_two_scales_lm)
tidy(m_two_scales_lm, effects = "fixed") %>% View()

dat <- ggeffects::ggeffect(m_two_scales_lm, terms = c("temp_anomaly", "EPU[GOM, MAB, GB]", "season"))
plot(dat)

# 1. Extract coefficients and the variance-covariance matrix
# Works for sdmTMB, lm, glm, etc.
coefs <- tibble::enframe(coef(m_two_scales_lm), name = "term", value = "estimate")
vc <- as.matrix(vcov(m_two_scales_lm))

# 2. Define your categories (update matching your exact level names)
epus <- c("GOM", "GB", "MAB")       # Example EPUs
seasons <- c("FALL", "SPRING")            # Example seasons
covariates <- c("climate_trend", "temp_anomaly")

# 3. Create a combination grid for calculation
slope_grid <- expand_grid(EPU = epus, season = seasons, covariate = covariates)

# 5. Map over the grid to build your final results table
slopes_df <- slope_grid %>%
  purrr::pmap_dfr(~ bind_cols(tibble(EPU = ..1, season = ..2, covariate = ..3), 
                              get_slope_ci(..1, ..2, ..3, coefs, vc)))

#----------------------------------------------------
## Model: fecundity
#----------------------------------------------------

m_two_scales_fecundity <- sdmTMB(
  data = df2,
  formula = fecundity ~ EPU * season *(climate_trend + temp_anomaly),
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
tidy(m_two_scales_fecundity, effects = "fixed") %>% View()

dat <- ggeffects::ggeffect(m_two_scales_fecundity, terms = c("temp_anomaly","EPU[GOM, MAB, GB]", "season"))
plot(dat)

# 1. Extract coefficients and the variance-covariance matrix
# Works for sdmTMB, lm, glm, etc.
coefs <- tibble::enframe(coef(m_two_scales_fecundity), name = "term", value = "estimate")
vc <- as.matrix(vcov(m_two_scales_fecundity))

# 2. Define your categories (update matching your exact level names)
epus <- c("GOM", "GB", "MAB")       # Example EPUs
seasons <- c("FALL", "SPRING")            # Example seasons
covariates <- c("climate_trend", "temp_anomaly")

# 3. Create a combination grid for calculation
slope_grid <- expand_grid(EPU = epus, season = seasons, covariate = covariates)

# 5. Map over the grid to build your final results table
slopes_df <- slope_grid %>%
  purrr::pmap_dfr(~ bind_cols(tibble(EPU = ..1, season = ..2, covariate = ..3), 
                              get_slope_ci(..1, ..2, ..3, coefs, vc)))
slopes_df


#----------------------------------------------------
## Model: l_inf
#----------------------------------------------------

m_two_scales_linf <- sdmTMB(
  data = df2,
  formula = l_inf ~ EPU * season *(climate_trend + temp_anomaly),
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
tidy(m_two_scales_linf, effects = "fixed") %>% View()

dat <- ggeffects::ggeffect(m_two_scales_linf, terms = c("temp_anomaly","EPU[GOM, MAB, GB]", "season"))
plot(dat)

# 1. Extract coefficients and the variance-covariance matrix
# Works for sdmTMB, lm, glm, etc.
coefs <- tibble::enframe(coef(m_two_scales_linf), name = "term", value = "estimate")
vc <- as.matrix(vcov(m_two_scales_linf))

# 2. Define your categories (update matching your exact level names)
epus <- c("GOM", "GB", "MAB")       # Example EPUs
seasons <- c("FALL", "SPRING")            # Example seasons
covariates <- c("climate_trend", "temp_anomaly")

# 3. Create a combination grid for calculation
slope_grid <- expand_grid(EPU = epus, season = seasons, covariate = covariates)

# 5. Map over the grid to build your final results table
slopes_df <- slope_grid %>%
  purrr::pmap_dfr(~ bind_cols(tibble(EPU = ..1, season = ..2, covariate = ..3), 
                              get_slope_ci(..1, ..2, ..3, coefs, vc)))
slopes_df

#----------------------------------------------------
## Model: k
#----------------------------------------------------

m_two_scales_k <- sdmTMB(
  data = df2,
  formula = k ~ EPU * season *(climate_trend + temp_anomaly),
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
tidy(m_two_scales_k, effects = "fixed") %>% View()

dat <- ggeffects::ggeffect(m_two_scales_k, terms = c("temp_anomaly","EPU[GOM, MAB, GB]", "season"))
plot(dat)

# 1. Extract coefficients and the variance-covariance matrix
# Works for sdmTMB, lm, glm, etc.
coefs <- tibble::enframe(coef(m_two_scales_k), name = "term", value = "estimate")
vc <- as.matrix(vcov(m_two_scales_k))

# 2. Define your categories (update matching your exact level names)
epus <- c("GOM", "GB", "MAB")       # Example EPUs
seasons <- c("FALL", "SPRING")            # Example seasons
covariates <- c("climate_trend", "temp_anomaly")

# 3. Create a combination grid for calculation
slope_grid <- expand_grid(EPU = epus, season = seasons, covariate = covariates)

# 5. Map over the grid to build your final results table
slopes_df <- slope_grid %>%
  purrr::pmap_dfr(~ bind_cols(tibble(EPU = ..1, season = ..2, covariate = ..3), 
                              get_slope_ci(..1, ..2, ..3, coefs, vc)))
slopes_df

#----------------------------------------------------
## Model: offspring_size
#----------------------------------------------------

m_two_scales_os <- sdmTMB(
  data = df2,
  formula = offspring_size ~ EPU * season *(climate_trend + temp_anomaly),
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
sanity(m_two_scales_os)
tidy(m_two_scales_os, effects = "fixed") %>% View()

dat <- ggeffects::ggeffect(m_two_scales_os, terms = c("temp_anomaly","EPU[GOM, MAB, GB]", "season"))
plot(dat)

# 1. Extract coefficients and the variance-covariance matrix
# Works for sdmTMB, lm, glm, etc.
coefs <- tibble::enframe(coef(m_two_scales_os), name = "term", value = "estimate")
vc <- as.matrix(vcov(m_two_scales_os))

# 2. Define your categories (update matching your exact level names)
epus <- c("GOM", "GB", "MAB")       # Example EPUs
seasons <- c("FALL", "SPRING")            # Example seasons
covariates <- c("climate_trend", "temp_anomaly")

# 3. Create a combination grid for calculation
slope_grid <- expand_grid(EPU = epus, season = seasons, covariate = covariates)

# 5. Map over the grid to build your final results table
slopes_df <- slope_grid %>%
  purrr::pmap_dfr(~ bind_cols(tibble(EPU = ..1, season = ..2, covariate = ..3), 
                              get_slope_ci(..1, ..2, ..3, coefs, vc)))
slopes_df

#----------------------------------------------------
## Model: trophic level
#----------------------------------------------------

m_two_scales_tl <- sdmTMB(
  data = df2,
  formula = trophic_level ~ EPU * season *(climate_trend + temp_anomaly),
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
sanity(m_two_scales_tl)
tidy(m_two_scales_tl, effects = "fixed") %>% View()

dat <- ggeffects::ggeffect(m_two_scales_tl, terms = c("temp_anomaly","EPU[GOM, MAB, GB]", "season"))
plot(dat)

# 1. Extract coefficients and the variance-covariance matrix
# Works for sdmTMB, lm, glm, etc.
coefs <- tibble::enframe(coef(m_two_scales_tl), name = "term", value = "estimate")
vc <- as.matrix(vcov(m_two_scales_tl))

# 2. Define your categories (update matching your exact level names)
epus <- c("GOM", "GB", "MAB")       # Example EPUs
seasons <- c("FALL", "SPRING")            # Example seasons
covariates <- c("climate_trend", "temp_anomaly")

# 3. Create a combination grid for calculation
slope_grid <- expand_grid(EPU = epus, season = seasons, covariate = covariates)

# 5. Map over the grid to build your final results table
slopes_df <- slope_grid %>%
  purrr::pmap_dfr(~ bind_cols(tibble(EPU = ..1, season = ..2, covariate = ..3), 
                              get_slope_ci(..1, ..2, ..3, coefs, vc)))
slopes_df


#----------------------------------------------------
## Model: max length
#----------------------------------------------------

m_two_scales_ml <- sdmTMB(
  data = df2,
  formula = max_length ~ EPU * season *(climate_trend + temp_anomaly),
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
sanity(m_two_scales_ml)
tidy(m_two_scales_ml, effects = "fixed") %>% View()

dat <- ggeffects::ggeffect(m_two_scales_ml, terms = c("temp_anomaly","EPU[GOM, MAB, GB]", "season"))
plot(dat)

# 1. Extract coefficients and the variance-covariance matrix
# Works for sdmTMB, lm, glm, etc.
coefs <- tibble::enframe(coef(m_two_scales_ml), name = "term", value = "estimate")
vc <- as.matrix(vcov(m_two_scales_ml))

# 2. Define your categories (update matching your exact level names)
epus <- c("GOM", "GB", "MAB")       # Example EPUs
seasons <- c("FALL", "SPRING")            # Example seasons
covariates <- c("climate_trend", "temp_anomaly")

# 3. Create a combination grid for calculation
slope_grid <- expand_grid(EPU = epus, season = seasons, covariate = covariates)

# 5. Map over the grid to build your final results table
slopes_df <- slope_grid %>%
  purrr::pmap_dfr(~ bind_cols(tibble(EPU = ..1, season = ..2, covariate = ..3), 
                              get_slope_ci(..1, ..2, ..3, coefs, vc)))
slopes_df
