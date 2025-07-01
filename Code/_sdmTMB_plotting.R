library(tidyverse)
library(sdmTMB) # This is the workhorse to do the spatio-temporal modeling
library(sf)
library(cowplot)
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

# Files that are needed for plotting

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

#-----------------------------------------
## Read in models and generate plots
#-----------------------------------------

# PCA1
  m_pc1_time <- readRDS("Data/Derived/Model_output/m_pc1_time.rds")
  sanity(m_pc1_time)
  
  p1 <- plot_model_outputs(m_pc1_time, return_plot = "p1", .label = "Pace of life (PCA1)")  # for spatial map
  p2 <- plot_model_outputs(m_pc1_time, return_plot = "p2", .label = "Pace of life (PCA1)") 


#PCA2
  m_pc2_time <- readRDS("Data/Derived/Model_output/m_pc2_time.rds")
  sanity(m_pc2_time)
  
  p3 <- plot_model_outputs(m_pc2_time, return_plot = "p1", .label = "Reproduction (PCA2)")  # for spatial map
  p4 <- plot_model_outputs(m_pc2_time, return_plot = "p2", .label = "Reproduction (PCA2)") 

#PCA3
  m_pc3_time <- readRDS("Data/Derived/Model_output/m_pc3_time.rds")
  sanity(m_pc3_time)
  
  p5 <- plot_model_outputs(m_pc3_time, return_plot = "p1", .label = "Trophic level (PCA3)")  # for spatial map
  p6 <- plot_model_outputs(m_pc3_time, return_plot = "p2", .label = "Trophic level (PCA3)") 

# Plot PCA figure...
  
  library(cowplot)
  
  ggdraw(xlim = c(0, 8.5), ylim = c(0,11)) +
    draw_plot(p5, x = 0, y = 0, width = 0.7*8.5, height = 0.33*11) +
    draw_plot(p6+theme(legend.position = "none"), x = 0.7*8.5, y = 0, width = 0.3*8.5, height = 0.3*11)+
    draw_plot(p3, x = 0, y = 11*0.33, width = 0.7*8.5, height = 0.33*11)+
    draw_plot(p4+theme(legend.position = "none"), x = 0.7*8.5, y = 11*0.33, width = 0.3*8.5, height = 0.3*11)+
    draw_plot(p1, x = 0, y = 11*0.66, width = 0.7*8.5, height = 0.33*11)+
    draw_plot(p2+theme(legend.position = "none"), x = 0.7*8.5, y = 11*0.66, width = 0.3*8.5, height = 0.3*11)
  
  ggsave("Figures/fig2.svg", height = 11*1.5, width = 8.5*1.5)
  
# Length-at-maturity
  
  m_lm_time <- readRDS("Data/Derived/Model_output/m_lm_time.rds")
  m_lm_time
  sanity(m_lm_time)
  
  p7 <- plot_model_outputs(m_lm_time, return_plot = "p1", .label = "Length-at-maturity (cm)")  # for spatial map
  p8 <- plot_model_outputs(m_lm_time, return_plot = "p2", .label = "Length-at-maturity (cm)", family.log = T) 
  
  p <- ggdraw(xlim = c(0, 8.5), ylim = c(0,3.66)) +
    draw_plot(p7, x = 0, y = 0, width = 0.7*8.5, height = 3.66) +
    draw_plot(p8, x = 0.7*8.5, y = 0.1, width = 0.3*8.5, height = 3.66*0.925)
  
  ggsave("Figures/fig_lm.png", p, height = 3.66*2, width = 8.5*2, bg = "white")
  
# Offspring size
  
  m_os_time <- readRDS("Data/Derived/Model_output/m_os_time.rds")
  m_os_time
  sanity(m_os_time)
  
  p9 <- plot_model_outputs(m_os_time, return_plot = "p1", .label = "Offspring size (mm)")  # for spatial map
  p10 <- plot_model_outputs(m_os_time, return_plot = "p2", .label = "Offspring size (mm)", family.log = T) 
  
  p <- ggdraw(xlim = c(0, 8.5), ylim = c(0,3.66)) +
    draw_plot(p9, x = 0, y = 0, width = 0.7*8.5, height = 3.66) +
    draw_plot(p10, x = 0.7*8.5, y = 0.1, width = 0.3*8.5, height = 3.66*0.925)
  
  ggsave("Figures/fig_off.png", p, height = 3.66*2, width = 8.5*2, bg = "white")
  
# Trophic level
  
  m_tl_time <- readRDS("Data/Derived/Model_output/m_tl_time.rds")
  m_tl_time
  sanity(m_tl_time)
  
  p11 <- plot_model_outputs(m_tl_time, return_plot = "p1", .label = "Trophic level")  # for spatial map
  p12 <- plot_model_outputs(m_tl_time, return_plot = "p2", .label = "Trophic level", family.log = F) 
  
  p <- ggdraw(xlim = c(0, 8.5), ylim = c(0,3.66)) +
    draw_plot(p11, x = 0, y = 0, width = 0.7*8.5, height = 3.66) +
    draw_plot(p12, x = 0.7*8.5, y = 0.1, width = 0.3*8.5, height = 3.66*0.925)
  
  ggsave("Figures/fig_tl.png", p, height = 3.66*2, width = 8.5*2, bg = "white")
  
# Fecundity
  
  m_fecundity_time <- readRDS("Data/Derived/Model_output/m_fecundity_time.rds")
  m_fecundity_time
  sanity(m_fecundity_time)
  
  p13 <- plot_model_outputs(m_fecundity_time, return_plot = "p1", .label = "Fecundity (num. eggs)")  # for spatial map
  p14 <- plot_model_outputs(m_fecundity_time, return_plot = "p2", .label = "Fecundity (num. eggs)", family.log = T) 
  
  p <- ggdraw(xlim = c(0, 8.5), ylim = c(0,3.66)) +
    draw_plot(p13, x = 0, y = 0, width = 0.7*8.5, height = 3.66) +
    draw_plot(p14, x = 0.7*8.5, y = 0.1, width = 0.3*8.5, height = 3.66*0.925)
  
  ggsave("Figures/fig_fecundity.png", p, height = 3.66*2, width = 8.5*2, bg = "white")
  
# Age-at-maturity
  
  m_age_time <- readRDS("Data/Derived/Model_output/m_age_time.rds")
  m_age_time
  sanity(m_age_time)
  
  p15 <- plot_model_outputs(m_age_time, return_plot = "p1", .label = "Age-at-maturity (years)")  # for spatial map
  p16 <- plot_model_outputs(m_age_time, return_plot = "p2", .label = "Age-at-maturity (years)", family.log = T) 
  
  p <- ggdraw(xlim = c(0, 8.5), ylim = c(0,3.66)) +
    draw_plot(p15, x = 0, y = 0, width = 0.7*8.5, height = 3.66) +
    draw_plot(p16, x = 0.7*8.5, y = 0.1, width = 0.3*8.5, height = 3.66*0.925)
  
  ggsave("Figures/fig_age.png", p, height = 3.66*2, width = 8.5*2, bg = "white")
  
# k-coefficient
  
  m_k_time <- readRDS("Data/Derived/Model_output/m_k_time.rds")
  m_k_time
  sanity(m_k_time)
  
  p17 <- plot_model_outputs(m_k_time, return_plot = "p1", .label = "VB Growth (k)")  # for spatial map
  p18 <- plot_model_outputs(m_k_time, return_plot = "p2", .label = "VB Growth (k)", family.log = T) 
  
  p <- ggdraw(xlim = c(0, 8.5), ylim = c(0,3.66)) +
    draw_plot(p17, x = 0, y = 0, width = 0.7*8.5, height = 3.66) +
    draw_plot(p18, x = 0.7*8.5, y = 0.1, width = 0.3*8.5, height = 3.66*0.925)
  
  ggsave("Figures/fig_k.png", p, height = 3.66*2, width = 8.5*2, bg = "white")
  
# L_inf-coefficient
  
  m_linf_time <- readRDS("Data/Derived/Model_output/m_linf_time.rds")
  m_linf_time
  tidy(m_linf_time, effects = "ran_pars")
  sanity(m_linf_time)
  
  p19 <- plot_model_outputs(m_linf_time, return_plot = "p1", .label = "VB length-at-infinity")  # for spatial map
  p20 <- plot_model_outputs(m_linf_time, return_plot = "p2", .label = "VB length-at-infinity", family.log = T) 
  
  p <- ggdraw(xlim = c(0, 8.5), ylim = c(0,3.66)) +
    draw_plot(p19, x = 0, y = 0, width = 0.7*8.5, height = 3.66) +
    draw_plot(p20, x = 0.7*8.5, y = 0.1, width = 0.3*8.5, height = 3.66*0.925)
  
  ggsave("Figures/fig_linf.png", p, height = 3.66*2, width = 8.5*2, bg = "white")
  
  
# Max observed length
  
  m_maxl_time <- readRDS("Data/Derived/Model_output/m_maxl_time.rds")
  m_maxl_time
  tidy(m_maxl_time, effects = "ran_pars")
  sanity(m_maxl_time)
  
  p21 <- plot_model_outputs(m_maxl_time, return_plot = "p1", .label = "Maximum observed length (cm)")  # for spatial map
  p22 <- plot_model_outputs(m_maxl_time, return_plot = "p2", .label = "Maximum observed length (cm)", family.log = T) 
  
  p <- ggdraw(xlim = c(0, 8.5), ylim = c(0,3.66)) +
    draw_plot(p21, x = 0, y = 0, width = 0.7*8.5, height = 3.66) +
    draw_plot(p22, x = 0.7*8.5, y = 0.1, width = 0.3*8.5, height = 3.66*0.925)
  
  ggsave("Figures/fig_maxl.png", p, height = 3.66*2, width = 8.5*2, bg = "white")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
