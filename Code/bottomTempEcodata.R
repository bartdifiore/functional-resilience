library(sf)
library(dplyr)
library(tidyr)

# 1. Assuming your gridded data frame is loaded (e.g., from the ecodata catalog)
# Expected columns usually include: Lon, Lat, Year, Season, and Value (e.g., bottom_temp)
# if 'df_raw' contains your long-format seasonal mean bottom temperature data...

sf_use_s2(FALSE)

btg <- ecodata::bottom_temp_model_gridded

ecodata_epu <- ecodata::epu_sf %>%
  select(EPU)

# btg_sp <- btg %>% 
#   distinct(Latitude, Longitude) %>%
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(ecodata_epu)) %>%
#   st_intersection(ecodata_epu) %>%
#   select(EPU)
# 
# btg_wepu <- btg %>%
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(ecodata_epu), remove = F) %>%
#   st_join(btg_sp) %>%
#   st_drop_geometry() %>%
#   filter(Var %in% c("spring", "fall"))
  
# 2. Extract unique grid coordinate centers to define the 1/12 degree cell boxes
grid_resolution <- 1/12  # ~0.083333 degrees

st_make_grid(btg, cellsize = grid_resolution)

# Create explicit polygon bounding geometry around each point center
bt_grid_sf <- btg %>%
  filter(Var %in% c("spring", "fall"), Time >= 1970) %>%
  mutate(
    xmin = Longitude - (grid_resolution / 2),
    xmax = Longitude + (grid_resolution / 2),
    ymin = Latitude - (grid_resolution / 2),
    ymax = Latitude + (grid_resolution / 2)
  ) %>%
  rowwise() %>%
  mutate(
    geometry = st_sfc(st_polygon(list(matrix(c(
      xmin, ymin,
      xmax, ymin,
      xmax, ymax,
      xmin, ymax,
      xmin, ymin
    ), ncol = 2, byrow = TRUE))))
  ) %>%
  st_as_sf(crs = st_crs(ecodata_epu))

# 3. Join the spatial grid polygons back to your full temporal data frame
# This preserves the un-aggregated 1959-2024 time series (Season x Year)
bt_seasonal <- bt_grid_sf %>%
  select(Time, Latitude, Longitude, Var) %>%
  st_join(ecodata_epu)

meta <- readRDS("Data/NMFS_trawl/survdat_pull_20250206.rds")$survdat %>%
  janitor::clean_names() %>%
  mutate(haul_id = paste(cruise6, station, sep = "_")) %>%
  select(cruise6, station, haul_id, stratum, tow, year, season, lat, lon) %>%
  distinct()

temp <- meta %>% 
  drop_na(lat, lon) %>%
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(ecodata_epu), remove = F)

test <- st_join(temp, bt_grid_sf %>% select(-c(Longitude, Latitude)), left = T) %>%
  st_drop_geometry() %>%
  left_join(bt_seasonal)

epu_strata <- read.csv("Data/Derived/strata_to_epu.csv")

df <- read.csv("Data/Derived/CWM_trawl_dataset.csv") %>% 
  as_tibble()











