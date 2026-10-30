library(ecodata)
library(sf)

epu <- ecodata::epu_sf %>%
  select(EPU)

bg <- ecodata::bottom_temp_model_gridded %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(epu)) %>%
  st_join(epu)

bt_means <- bg %>%
  filter(Var %in% c("fall", "spring"), EPU != "SS") %>%
  rename(year = Time, bottom_temp = Value, season = Var) %>%
  st_drop_geometry() %>%
  mutate(season = toupper(season)) %>%
  group_by(EPU, season, year) %>%
  summarize(annual_mean_bt = mean(bottom_temp, na.rm = T)) %>%
  ungroup() %>%
  arrange(EPU, season, year) %>%
  group_by(EPU, season) %>%
  mutate(
    climate_trend = zoo::rollapply(annual_mean_bt, width = 7, FUN = mean, na.rm = T, align = "center", fill = "extend"),
    # Calculate the anomaly purely at the YEARLY scale
    annual_temp_anomaly = annual_mean_bt - climate_trend 
  )

write_rds(bt_means, "Data/Derived/bottom_temp_means.rds")
 

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
  select(-bottom_temp) %>%
  filter(EPU != "SS") %>%
  left_join(bt_means) %>%
  mutate(id = as.factor(paste(EPU, year, season, sep = "-")))

sf_use_s2(FALSE)

pj_crs <- st_crs(x = "EPSG:32619")

df$scaled_year <- (df$year - mean(df$year))/10 
df$scaled_year_season <- (df$year_season - mean(df$year_season))/10

summary(df)

df_lm <- df%>%
  drop_na(length_maturity, climate_trend, annual_temp_anomaly)

mesh <- sdmTMB::make_mesh(df_lm, xy_cols = c("X", "Y"), cutoff = 30)

m_two_scales_lm <- sdmTMB(
  data = df_lm,
  formula = length_maturity ~ EPU * season *(climate_trend + annual_temp_anomaly) + (1|id),
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "on",
  time = "year",
  # time_varying = ~  1,
  # time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017, 2020),
  silent = F
)
sanity(m_two_scales_lm)
tidy(m_two_scales_lm, effects = "fixed") %>% View()

dat <- ggeffects::ggeffect(m_two_scales_lm, terms = c("climate_trend", "EPU[GOM, MAB, GB]", "season"))
plot(dat)

# 1. Extract coefficients and the variance-covariance matrix
# Works for sdmTMB, lm, glm, etc.
coefs <- tibble::enframe(coef(m_two_scales_lm), name = "term", value = "estimate")
vc <- as.matrix(vcov(m_two_scales_lm))

# 2. Define your categories (update matching your exact level names)
epus <- c("GOM", "GB", "MAB")       # Example EPUs
seasons <- c("FALL", "SPRING")            # Example seasons
covariates <- c("climate_trend", "annual_temp_anomaly")

# 3. Create a combination grid for calculation
slope_grid <- expand_grid(EPU = epus, season = seasons, covariate = covariates)

# 5. Map over the grid to build your final results table
slopes_df <- slope_grid %>%
  purrr::pmap_dfr(~ bind_cols(tibble(EPU = ..1, season = ..2, covariate = ..3), 
                              get_slope_ci(..1, ..2, ..3, coefs, vc)))


