library(sf)
library(dplyr)
library(tidyr)

# 1. Assuming your gridded data frame is loaded (e.g., from the ecodata catalog)
# Expected columns usually include: Lon, Lat, Year, Season, and Value (e.g., bottom_temp)
# if 'df_raw' contains your long-format seasonal mean bottom temperature data...

sf_use_s2(FALSE)

year_filter <- c(seq(1963, 1969, by = 1), 2017, 2020)

btg <- ecodata::bottom_temp_model_gridded

ecodata_epu <- ecodata::epu_sf %>%
  select(EPU) %>%
  st_make_valid()

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
  distinct() %>%
  drop_na(lat, lon) %>%
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(ecodata_epu), remove = F) %>%
  st_make_valid() %>%
  mutate(id = paste(year, season, sep = "-")) %>% 
  filter(!year %in% year_filter)

bt_grid_sf2 <- bt_grid_sf %>%
  rename(year = Time, bottom_temp = Value) %>%
  filter(!year %in% year_filter) %>%
  mutate(season = ifelse(Var == "spring", "SPRING", "FALL")) %>%
  mutate(id = paste(year, season, sep = "-")) %>%
  select(id, season, year, bottom_temp) %>%
  st_make_valid()


index <- unique(meta$id)
out <- list()
for(i in 1:length(index)){
  out[[i]] <- meta %>% 
    filter(id == index[i]) %>%
    st_join(bt_grid_sf2 %>% filter(id == index[i]))
  
  combined_df <- data.table::rbindlist(out)
}

meta_wbt <- combined_df %>%
  select(cruise6, station, haul_id, stratum, tow, year.x, season.x, lat, lon, bottom_temp) %>%
  st_drop_geometry() %>%
  rename(year = year.x, season = season.x)

write.csv(meta_wbt, "Data/Derived/trawlmeta_wBT.csv", quote = F, row.names = F) 




polys <- read_sf("Data/Shapefiles/prediction_grid.shp") %>%
  st_transform(crs = st_crs(ecodata_epu))

btg_sp <- ecodata::bottom_temp_model_gridded %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(ecodata_epu)) %>% 
  st_join(polys) %>%
  st_drop_geometry() %>%
  filter(Var %in% c("spring", "fall")) %>% 
  mutate(season = ifelse(Var == "spring", "SPRING", "FALL")) %>% 
  rename(year = Time, fid = FID, bottom_temp = Value) %>%
  select(year, season, bottom_temp, fid) %>%
  group_by(year, season, fid) %>%
  filter(year >= 1970) %>%
  summarize(bottom_temp = mean(bottom_temp))

write_rds(btg_sp, "Data/Derived/prediction_grid_wBT.rds")






