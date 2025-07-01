library(gmRi)
library(tidyverse)
library(sf)
library(ggplot2)
library(vegan)
library(mFD)
source("~/Github/gmRi/R/survey_prep_bd.R")


# Paths to Box folder
res_box_path<- cs_path(box_group = "Res Data")

#---------------------------------------------
## Trait database and add in trophic groups
#---------------------------------------------

trophic_groups <- ecodata::species_groupings %>% 
  janitor::clean_names() %>% 
  select(svspp, soe_20, size_cat) %>%
  # group_by(svspp) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = size_cat, values_from = soe_20) %>% 
  rename(no_size = `NA`) %>%
  mutate(temp = paste(no_size, M, S, L, XL, sep = ".")) %>% 
  mutate(trophic_group = str_remove_all(temp, ".NA|NA|[[:punct:]]"))%>%
  select(svspp, trophic_group) %>% 
  distinct()

traits <- read.csv("Data/trait_database.csv") %>% 
  left_join(trophic_groups) %>% 
  as_tibble() %>% 
  filter(!is.na(svspp)) %>%
  dplyr::select(svspp, genus, species,
                habitat, 
                # feeding_mode, 
                spawning_type,
                trophic_group,
                trophic_level,
                offspring_size, 
                age_maturity,
                fecundity,
                l_inf,
                k,
                max_obs_length,
                length_maturity, 
                # age_max # droping this one because too many gaps in the database
  ) %>% 
  # distinct(.keep_all = TRUE) %>%
  mutate(habitat = case_when(is.na(habitat) == T ~ "Other", .default = habitat), 
         spawning_type = case_when(is.na(spawning_type) == T | spawning_type == "" ~ "Other", .default = spawning_type), 
         trophic_group = case_when(is.na(trophic_group) == T ~ "Other", .default = trophic_group)) %>%
  drop_na(offspring_size:length_maturity)%>%
  arrange(svspp) %>% 
  group_by(svspp, habitat, spawning_type, trophic_group) %>%
  summarize(across(where(is.double), \(x) mean(x, na.rm = T))) %>%
  filter(!c(svspp == 573 & habitat == "reef-associated")) # Deal with one species which was listed as either reef associated or pelagic

spp <- traits %>%
  pull(svspp)



#-----------------------------------
## NEFSC Trawl Survey data
#-----------------------------------

nmfs_df_raw <- readRDS("Data/NMFS_trawl/survdat_pull_20250206.rds")$survdat

df <- nmfs_df_raw %>% 
  janitor::clean_names() %>%
  as_tibble() %>%
  group_by(across(-c(catchsex, abundance, biomass))) %>% # Species that can be sexed (like DogFish, inverts) have two or more separate entries in a tow, one for each category. Here I'm grouping all other variables in order to summerize the total (sum()) abundance or biomass. 
  summarize(abundance = sum(abundance), 
            biomass_kg = sum(biomass)) %>%
  mutate(svspp_num = as.numeric(svspp)) %>%
  drop_na(lat, lon) %>%
  mutate(haul_id = paste(cruise6, station, sep = "_")) %>%
  filter(svspp %in% spp)

crs <- st_crs(4326)

epu_clipper <- ecodata::epu_sf %>% dplyr::select(EPU) %>% st_make_valid() %>% st_transform(crs = crs)

survey_species <- df %>%
  st_as_sf(coords = c("lon", "lat")) %>% 
  st_set_crs(4269) %>%
  st_transform(crs = crs) %>%
  st_join(epu_clipper) %>% 
  filter(!is.na(EPU)) %>%
  st_drop_geometry()

survey_mat <- survey_species %>%
  select(year, season, EPU, stratum, haul_id, svspp, biomass_kg) %>%
  mutate(svspp_names = paste("X", svspp, sep = "")) %>%
  select(-svspp) %>%
  pivot_wider(names_from = svspp_names, values_from = biomass_kg, names_sort = T) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))


species_in_survey <- unique(survey_species$svspp)

dim(survey_mat)

survey.mat <- survey_mat %>% ungroup() %>%
  select(-c(year, season, EPU, stratum, haul_id)) %>% 
  as.matrix()

dim(survey.mat)

num.species <- dim(survey.mat)[2]

# So ther are 388 species in the trait database. But only 372 species after we filter species in the trait database from survdat. I suspect this is because we drop some tows outside of the EPU boundary. But may be worth confirming down the road.

#-----------------------------------
## Clean up and write out data
#-----------------------------------

trait.mat <- traits %>%
  filter(svspp %in% species_in_survey) %>% 
  mutate(svspp = as.character(svspp)) %>%
  as.data.frame()

dim(trait.mat)
dim(survey.mat)

write_rds(trait.mat, "Data/Derived/trait-df.rds", compress = "gz")
write_rds(survey.mat, "Data/Derived/survey-matrix.rds", compress = "gz")
write_rds(survey_mat, "Data/Derived/survey-df.rds", compress = "gz")



