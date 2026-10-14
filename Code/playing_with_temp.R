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

meta <- readRDS("Data/NMFS_trawl/survdat_pull_20250206.rds")$survdat %>%
  janitor::clean_names() %>%
  mutate(haul_id = paste(cruise6, station, sep = "_")) %>%
  select(cruise6, station, haul_id, stratum, tow, year, season, lat, lon) %>%
  distinct()

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
  drop_na(EPU)

sf_use_s2(FALSE)

pj_crs <- st_crs(x = "EPSG:32619")

df$scaled_year <- (df$year - mean(df$year))/10 
df$scaled_year_season <- (df$year_season - mean(df$year_season))/10

hist(df$length_maturity)



# Other stuff I will need

polys <- read_sf("Data/Shapefiles/prediction_grid.shp") %>%
  rename(fid = FID)

epus <- ecodata::epu_sf %>% 
  select(EPU, geometry) %>%
  st_transform(crs = pj_crs)

btg <- ecodata::bottom_temp_model_gridded %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"))

bt <- ecodata::bottom_temp_model_anom %>% 
  pivot_wider(names_from = Var, values_from = Value) %>% 
  janitor::clean_names() %>%
  rename(FALL = fall_bottom_temp_anomaly, SPRING = spring_bottom_temp_anomaly) %>% 
  select(time, epu, source, FALL, SPRING) %>% 
  group_by(time, epu, source) %>%
  pivot_longer(cols = c(FALL, SPRING)) %>%
  filter(source != "MOM6") %>% 
  rename(year = time, EPU = epu, bottom_temp = value, season = name)

ggplot(bt)+
  geom_line(aes(x = time, y = value, color = source))+
  facet_wrap(~epu)

df_wtemp <- df %>% 
  left_join(bt)

pred_grid <- read_sf("Data/Shapefiles/predictiongrid_wcovairates.shp") %>%
  st_centroid() %>% 
  st_join(epus) %>%
  mutate(X = as.numeric(st_coordinates(.)[,1])/1000, 
         Y = as.numeric(st_coordinates(.)[,2])/1000) %>%
  st_drop_geometry() %>%
  drop_na(EPU) %>%
  replicate_df(time_name = "year", time_values = 1970:2024) %>% 
  expand_grid(season = c("FALL", "SPRING")) %>%
  mutate(season_temp = ifelse(season == "SPRING", 25, 75), 
         year_season = as.numeric(paste(year, season_temp, sep = "."))) %>% 
  filter(!year %in% year_filter) %>%
  rename(fid = FID)

pred_grid$scaled_year <- (pred_grid$year - mean(pred_grid$year)) / 10

shelf <- read.csv("Data/Derived/CWM_trawl_dataset_wholeshelf.csv")

#----------------------------------------------------
## Model PCA1
#----------------------------------------------------

# Make the mesh
mesh <- sdmTMB::make_mesh(df_wtemp, xy_cols = c("X", "Y"), cutoff = 30)

plot(mesh)

m_playtime <- sdmTMB(
  data = df_wtemp,
  formula = length_maturity ~ bottom_temp*as.factor(EPU)*as.factor(season),
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year_season",
  time_varying = ~1,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017.25, 2017.75, 2020.25, 2020.75),
  silent = F
)

m_playtime
summary(m_playtime)
sanity(m_playtime)

ggeffects::ggemmeans(m_playtime, terms = c("bottom_temp", "season", "EPU"))
dat <- ggeffects::ggeffect(m_playtime, terms = c("bottom_temp", "season", "EPU"))
plot(dat)

estimate_decadal_slopes(m_playtime)


df_wtemp <- df_wtemp %>%
  mutate(EPU = as.factor(EPU), 
         season = as.factor(season))


m_playtime <- sdmTMB(
  data = df_wtemp,
  formula = length_maturity ~ bottom_temp*EPU*season,
  # formula = length_maturity ~ 0,
  mesh = mesh,
  family = gengamma(link = "log"), 
  spatial = "off",
  time = "year_season",
  time_varying = ~0 + EPU,
  time_varying_type = "ar1",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  extra_time = c(2017.25, 2017.75, 2020.25, 2020.75),
  silent = F
)

m_playtime
summary(m_playtime)
sanity(m_playtime)

tidy(m_playtime, effects = "fixed")

dat <- ggeffects::ggeffect(m_playtime, terms = c("bottom_temp", "season", "EPU"))
plot(dat)

tidy(m_playtime, effects = "ran_vals") %>%
  separate(term, into = c("EPU", "year_season"), ":") %>%
  separate(EPU, into = c("junk", "EPU"), "[)]") %>%
  select(-junk) %>%
  ggplot(aes(x = as.numeric(year_season), y = estimate))+
  geom_line(aes(color = EPU))+
  facet_wrap(~EPU)


# Create EPU-level prediction grid
ndata <- expand.grid(
  year = unique(df$year),
  season = unique(df$season),
  EPU = c("GOM", "GB", "MAB"),
  bottom_temp = seq(min(df_wtemp$bottom_temp, na.rm = T), max(df_wtemp$bottom_temp, na.rm = T), length.out = 100)) %>%
  mutate(
    season_temp = ifelse(season == "SPRING", 25, 75),
    year_season = as.numeric(paste(year, season_temp, sep = "."))
  )

# Point estimates and uncertainty
est <- predict(m_playtime, newdata = ndata, re_form = NA, nsim = 100)

ndata$ci_lower <- apply(est, 1, quantile, probs = 0.025)
ndata$ci_upper <- apply(est, 1, quantile, probs = 0.975)
ndata$est_mean <- apply(est, 1, mean)


ggplot(ndata)+
  geom_line(aes(x = bottom_temp, y = exp(est_mean), color = season, group = year_season))+
  facet_wrap(~EPU)

ggplot(ndata)+
  geom_line(aes(x = year_season, y = exp(est_mean)))+
  facet_wrap(~EPU)

tv_effects <- tidy(m_playtime, effects = "ran_vals", conf.int = TRUE)
                           