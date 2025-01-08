library(gmRi)
library(tidyverse)
library(sf)
library(ggplot2)
library(vegan)
source("~/Github/gmRi/R/survey_prep_bd.R")


# Paths to Box folder
res_box_path<- cs_path(box_group = "Res Data")
proj_box_path<- cs_path(box_group = "Mills Lab", subfolder = "Projects/sdmTMB_sandbox")
lob_box_path <- cs_path(box_group = "Mills Lab", subfolder = "Projects/Lobster")

#-----------------------------------
## NEFSC Trawl Survey data
#-----------------------------------

df <- readRDS(paste0(res_box_path, "NMFS_trawl/SURVDAT_current/survdat_lw.rds")) # From https://github.com/NOAA-EDAB/data-requests/tree/main/EwE-menhaden-AndreBuchheister

nmfs_df_raw <- df$survdat

nmfs_df <- gmri_survdat_prep_bd(survdat = nmfs_df_raw, box_location = "cloudstorage")

traits <- read.csv("Data/trait_database.csv") %>% 
  filter(!is.na(svspp)) %>%
  dplyr::select(svspp, genus, species,
                offspring_size, 
                age_maturity,
                fecundity,
                l_inf,
                k,
                max_obs_length,
                length_maturity, 
                age_max) %>% 
  distinct(.keep_all = TRUE) %>% 
  filter_at(vars(offspring_size:age_max), any_vars(complete.cases(.)))

spp <- traits %>% 
  pull(svspp) 

df <- nmfs_df %>% 
  select(svspp:biomass_kg) %>%
  distinct() %>%
  janitor::clean_names() %>%
  as_tibble() %>%
  group_by(across(-c(catchsex, abundance, biomass_kg))) %>%
  summarize(biomass_kg = sum(biomass_kg), 
            abundance = sum(abundance)) %>%
  ungroup() %>%
  select(-c(est_towdate, avgdepth, surftemp, surfsalin, bottemp, botsalin)) %>%
  mutate(haul_id = paste(cruise6, station, sep = "_")) %>% 
  filter(svspp %in% spp,
         abundance >= 0) %>%
  rename(lat = decdeg_beglat, lon = decdeg_beglon)



#-----------------------------------------
## Build the grid
#-----------------------------------------

# Due to inconsistent sampling, we have decided to focus on 5 year periods in each hex, and discard all hex-year combinations with less than 5 tows. 

tows <- df %>% 
  dplyr::select(lat, lon, haul_id, est_year) %>% 
  distinct()

crs <- st_crs(4326)

epu_clipper <- ecodata::epu_sf %>% dplyr::select(EPU) %>% st_make_valid() %>% st_transform(crs = crs)

tows_sf <- tows %>% 
  st_as_sf(coords = c("lon", "lat"), crs = crs) %>% 
  st_transform(crs = crs) %>%
  st_join(epu_clipper) %>% 
  filter(!is.na(EPU)) %>%
  # mutate(five_year_period = est_year  %/% 5 * 5)
  mutate(five_year_period = case_when(est_year >= 1970 & est_year < 1975 ~ 1970, 
                                      est_year >= 1975 & est_year < 1980 ~ 1975, 
                                      est_year >= 1980 & est_year < 1985 ~ 1980, 
                                      est_year >= 1985 & est_year < 1990 ~ 1985, 
                                      est_year >= 1990 & est_year < 1995 ~ 1990, 
                                      est_year >= 1995 & est_year < 2000 ~ 1995, 
                                      est_year >= 2000 & est_year < 2005 ~ 2000, 
                                      est_year >= 2005 & est_year < 2010 ~ 2005, 
                                      est_year >= 2010 & est_year < 2015 ~ 2010, 
                                      est_year >= 2015 & est_year < 2020 ~ 2015))

hex_raw <- tows_sf %>% 
  st_bbox() %>% 
  st_as_sfc(crs = crs) %>% 
  st_make_grid(cellsize = 0.5, square = F) %>% 
  st_sf() %>% 
  tibble::rowid_to_column('hex_id')


hauls_in_hex <- st_join(tows_sf, hex_raw, join = st_within) %>%
  st_set_geometry(NULL) %>% 
  count(name = "hauls", hex_id, five_year_period) %>% 
  filter(!is.na(hex_id), five_year_period >= 1970)

dim(hauls_in_hex[hauls_in_hex$hauls < 5, ])/dim(hauls_in_hex)

hauls_in_hex %>% 
  ggplot(aes(x = hauls))+
  geom_histogram(color = "white", binwidth = 1, center = 0)+
  geom_vline(xintercept = 5, color = "red")+
  theme_bw()

# Based on this we are going to eliminate all hex-years with less than 5 tows.

hex <- hex_raw %>%
  left_join(hauls_in_hex, by = 'hex_id') # %>%
# st_join(ecodata::epu_sf %>% dplyr::select(EPU) %>% st_make_valid() %>% st_transform(crs = crs), join = st_nearest_feature)

write_rds(hex, "Data/Derived/hex.rds")


hex %>%
  filter(!is.na(hauls)) %>%
  filter(five_year_period >= 1970) %>%
  filter(hauls >= 5) %>%
  ggplot()+
  geom_sf(aes(fill = hauls))+
  scale_fill_binned(type = "viridis", breaks = c(5, 6, 7, 8, 9, 10, 100), alpha = 0.5)+
  geom_sf(data = ecodata::epu_sf, fill = "transparent", size = 2, color = "darkred")+
  #geom_sf(data = df_sf, size = 0.1, alpha = 0.5)+
  facet_wrap(~five_year_period)

hex_formerge <- hex %>% 
  select(hex_id, geometry) %>% 
  distinct()

formerge <- hex %>% 
  st_drop_geometry()


survey_species <- df %>%
  st_as_sf(coords = c("lon", "lat"), crs = crs) %>% 
  st_transform(crs = crs) %>%
  st_join(epu_clipper) %>% 
  filter(!is.na(EPU)) %>%
  st_join(hex_formerge, join = st_within) %>%
  st_drop_geometry() %>%
  mutate(five_year_period = case_when(est_year >= 1970 & est_year < 1975 ~ 1970, 
                                      est_year >= 1975 & est_year < 1980 ~ 1975, 
                                      est_year >= 1980 & est_year < 1985 ~ 1980, 
                                      est_year >= 1985 & est_year < 1990 ~ 1985, 
                                      est_year >= 1990 & est_year < 1995 ~ 1990, 
                                      est_year >= 1995 & est_year < 2000 ~ 1995, 
                                      est_year >= 2000 & est_year < 2005 ~ 2000, 
                                      est_year >= 2005 & est_year < 2010 ~ 2005, 
                                      est_year >= 2010 & est_year < 2015 ~ 2010, 
                                      est_year >= 2015 & est_year < 2020 ~ 2015, 
                                      est_year >= 2020 & est_year < 2025 ~ 2020)) %>%
  left_join(formerge) %>%
  drop_na(hex_id, five_year_period) %>%
  filter(hauls >= 5)

survey_mat <- survey_species %>%
  group_by(five_year_period, hex_id, EPU, haul_id, svspp) %>%
  summarize(total_abundance = sum(abundance, na.rm =T)) %>%
  pivot_wider(names_from = svspp, values_from = total_abundance, names_sort = T) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(id = paste(hex_id, five_year_period, sep = "_"))

species_in_survey <- unique(survey_species$svspp)

trait_mat <- traits %>%
  select(-genus, -species) %>% 
  filter(svspp %in% species_in_survey) %>%
  arrange(svspp) %>% 
  group_by(svspp) %>%
  summarize(across(where(is.double), \(x) mean(x, na.rm = T))) %>%
  tibble::column_to_rownames(var = "svspp") %>% 
  as.data.frame()


dim(trait_mat)
dim(survey_mat)

survey.mat <- survey_mat %>% ungroup() %>%
  select(-c(est_year, hex_id, EPU, haul_id, id)) %>% 
  as.matrix()

dim(survey.mat)


num.species <- dim(survey.mat)[2]

group_sizes <- round(num.species * c(0.1, 0.2, 0.5))
#29, 58, 146

# Run the clustering analysis to determine functional groups

neus_traits2 <-  FD::dbFD(x = trait_mat, 
                          a = survey.mat, 
                          corr = "cailliez", 
                          calc.FGR = T,
                          calc.CWM = F, 
                          print.pco = T)
str(neus_traits2)

readr::write_rds(neus_traits2, "Data/Derived/Model_output/dbFDoutput_w29groups.rds", compress = "gz")


mod29 <- readRDS("Data/Derived/Model_output/dbFDoutput_w29groups.rds")

species_in_fg <- data.frame(svspp = names(survey_mat)[-c(1:4, dim(survey_mat)[2])], spfgr = mod29$spfgr)

out <- survey_species %>%
  left_join(species_in_fg)

write_rds(out, "Data/Derived/survdat_wfG.RDS", compress = "gz")










