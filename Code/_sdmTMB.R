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


pred_grid <- read_sf("Data/Shapefiles/prediction_grid.shp") %>%
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
mesh <- sdmTMB::make_mesh(df, xy_cols = c("X", "Y"), cutoff = 30)

plot(mesh)

m_pc1_time <- sdmTMB(
  data = df,
  formula = PC1 ~ scaled_year*EPU*season,
  mesh = mesh,
  spatial = "on",
  time = "year_season",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  silent = F
)

m_pc1_time
sanity(m_pc1_time)
write_rds(m_pc1_time, "Data/Derived/Model_output/m_pc1_time.rds")

#----------------------------------------------------
## Model PCA2
#----------------------------------------------------

# Make the mesh

m_pc2_time <- sdmTMB(
  data = df,
  formula = PC2 ~ scaled_year*EPU*season,
  mesh = mesh,
  spatial = "on",
  time = "year_season",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  silent = F
)

m_pc2_time
sanity(m_pc2_time)
write_rds(m_pc2_time, "Data/Derived/Model_output/m_pc2_time.rds")

#----------------------------------------------------
## Model PCA3
#----------------------------------------------------

# Make the mesh

m_pc3_time <- sdmTMB(
  data = df,
  formula = PC3 ~ scaled_year*EPU*season,
  mesh = mesh,
  spatial = "on",
  time = "year_season",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  silent = F
)

m_pc3_time
sanity(m_pc3_time)
write_rds(m_pc3_time, "Data/Derived/Model_output/m_pc3_time.rds")

#----------------------------------------------------
## Model length at maturity
#----------------------------------------------------

m_lm_time <- sdmTMB(
  data = df,
  formula = length_maturity ~ scaled_year*EPU*season,
  mesh = mesh,
  family = gengamma(link = "log"),
  spatial = "off",
  time = "year_season",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  # extra_time = c(2017.25, 2020.25),
  silent = F
) # Wouldn't work with the tweedie

m_lm_time
sanity(m_lm_time)
tidy(m_lm_time, effects = "ran_pars")
write_rds(m_lm_time, "Data/Derived/Model_output/m_lm_time.rds")


#----------------------------------------------------
## Model offspring size
#----------------------------------------------------

m_os_time <- sdmTMB(
  data = df,
  formula = offspring_size ~ scaled_year*EPU*season,
  mesh = mesh,
  family = gengamma(link = "log"),
  spatial = "off",
  time = "year_season",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  # extra_time = c(2017.25, 2020.25),
  silent = F
) # Tried the tweedie but had issues with convergence. gengamma worked fine.

m_os_time
sanity(m_os_time)
tidy(m_os_time, effects = "ran_pars")
write_rds(m_os_time, "Data/Derived/Model_output/m_os_time.rds")

#----------------------------------------------------
## Model trophic_level
#----------------------------------------------------

m_tl_time <- sdmTMB(
  data = df,
  formula = trophic_level ~ scaled_year*EPU*season,
  mesh = mesh,
  family = gengamma(link = "log"),
  spatial = "off",
  time = "year_season",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  # extra_time = c(2017.25, 2020.25),
  silent = F
) # Didn't work with the tweedie, reverting back to the gengamma. 

m_tl_time
sanity(m_tl_time)
tidy(m_tl_time, effects = "ran_pars")
write_rds(m_tl_time, "Data/Derived/Model_output/m_tl_time.rds")

#----------------------------------------------------
## Model fecundity
#----------------------------------------------------

# Make the mesh

m_fecundity_time <- sdmTMB(
  data = df,
  formula = fecundity ~ scaled_year*EPU*season,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "off",
  time = "year_season",
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  # extra_time = c(2017.25, 2020.25),
  silent = F
)

m_fecundity_time
sanity(m_fecundity_time)
tidy(m_fecundity_time, effects = "ran_pars")
write_rds(m_fecundity_time, "Data/Derived/Model_output/m_fecundity_time.rds")


#----------------------------------------------------
## Model age-at-maturity
#----------------------------------------------------

# Make the mesh

m_age_time <- sdmTMB(
  data = df,
  formula = age_maturity ~ scaled_year*EPU*season,
  mesh = mesh,
  spatial = "on",
  time = "year_season",
  family = tweedie(link = "log"),
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  # extra_time = c(2017.25, 2020.25),
  silent = F
)

m_age_time
sanity(m_age_time)
tidy(m_age_time, effects = "ran_pars")
write_rds(m_age_time, "Data/Derived/Model_output/m_age_time.rds")

#----------------------------------------------------
## Model k
#----------------------------------------------------

# Make the mesh

m_k_time <- sdmTMB(
  data = df,
  formula = k ~ scaled_year*EPU*season,
  mesh = mesh,
  spatial = "off",
  time = "year_season",
  family = gengamma(link = "log"),
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  # extra_time = c(2017.25, 2020.25),
  silent = F
) # Issues

m_k_time
sanity(m_k_time)
tidy(m_k_time, effects = "ran_pars")
write_rds(m_k_time, "Data/Derived/Model_output/m_k_time.rds")

#----------------------------------------------------
## Model l_inf
#----------------------------------------------------

# Make the mesh

m_linf_time <- sdmTMB(
  data = df,
  formula = l_inf ~ scaled_year*EPU*season,
  mesh = mesh,
  spatial = "off",
  time = "year_season",
  family = tweedie(link = "log"),
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  # extra_time = c(2017.25, 2020.25),
  silent = F
)

m_linf_time
sanity(m_linf_time)
tidy(m_linf_time, effects = "ran_pars")
write_rds(m_linf_time, "Data/Derived/Model_output/m_linf_time.rds")

#----------------------------------------------------
## Model max obs length
#----------------------------------------------------

# Make the mesh

hist(df$max_obs_length)
summary(df$max_obs_length)

m_maxl_time <- sdmTMB(
  data = df,
  formula = l_inf ~ scaled_year*EPU*season,
  mesh = mesh,
  spatial = "off",
  time = "year_season",
  family = tweedie(link = "log"),
  # spatial_varying = ~ scaled_year_season,
  spatiotemporal = "IID", 
  # extra_time = c(2017.25, 2020.25), 
  silent = F
)

m_maxl_time
sanity(m_maxl_time)
tidy(m_maxl_time, effects = "ran_pars")
write_rds(m_maxl_time, "Data/Derived/Model_output/m_maxl_time.rds")

  


