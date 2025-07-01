library(tidyverse)
library(sf)

crs <- st_crs(4326)

epu_clipper <- ecodata::epu_sf %>% dplyr::select(EPU) %>% st_make_valid() %>% st_transform(crs = crs)

strata_sf <- sf::read_sf("Data/Shapefiles/strata.shp") %>%
  select(STRATA, geometry) %>%
  rename(stratum = STRATA) %>%
  st_make_valid() %>% 
  st_transform(crs = crs) %>%
  mutate(total_area = st_area(geometry))

intersection <- st_intersection(strata_sf, epu_clipper) %>% st_make_valid() %>%
  mutate(area = st_area(geometry)) %>%
  left_join(strata_sf %>% select(stratum, total_area) %>% st_drop_geometry()) %>%
  mutate(prop_area = area/total_area)

epu_strata <- intersection %>% 
  filter(as.numeric(prop_area) >= 0.5) %>% 
  select(stratum, EPU) %>% 
  st_drop_geometry() %>%
  distinct() %>%
  add_row(stratum = c(1630, 1640, 1170, 1180, 3500), EPU = c("MAB", "MAB", "GB", "GB", "MAB"))

write.csv(epu_strata, "Data/Derived/strata_to_epu.csv", quote = F, row.names = F)