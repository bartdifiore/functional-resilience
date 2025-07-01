#-----------------
## Libraries
#-----------------

library(tidyverse)
library(vegan)
library(ggvegan)
source("Code/1_DataSetup.R")
source("Code/theme.R")

#-----------------
## Load data
#-----------------

year_filter <- c(seq(1963, 1969, by = 1), 2017, 2020)

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
  # filter(svspp %in% spp) %>%
  filter(!year %in% year_filter)

# Get the trait data
traits <- readRDS("Data/Derived/trait-df_wPCA.rds")
spp_in_traits <- unique(traits$svspp)

trait.mat <- traits %>%
  select(-c(habitat, spawning_type, trophic_group)) %>%
  column_to_rownames("svspp") %>%
  as.matrix()

# Estimate CWM at the scale of the trawl

biomass_mat <- df %>% 
  filter(svspp %in% spp_in_traits) %>%
  ungroup() %>%
  select(year, season, stratum, haul_id, svspp, biomass_kg) %>%
  pivot_wider(names_from = svspp, values_from = biomass_kg, values_fill = 0, names_sort = T)

bio.mat <- as.matrix(biomass_mat[,-c(1:4)])

cwm_mat <- CWM(T_mat = trait.mat, S_mat = bio.mat) # Run the community weighting function. 

cwm_df <- biomass_mat[,c(1:4)] %>%
  bind_cols(as.data.frame(cwm_mat)) %>% 
  as_tibble()

# Compare with FD::functcomp() estimate
temp <- FD::functcomp(trait.mat, bio.mat)  

temp[1:6, 1] # Estimates are the same
cwm_df$trophic_level[1:6] # Estimates are the same

write.csv(cwm_df, "Data/Derived/CWM_trawl_dataset.csv", row.names = F, quote = F)



# Estimate at the whole shelf scale

shelf <- readRDS("Data/Derived/stratified_mean_biomass_whole_shelf.rds") %>% 
  filter(svspp %in% spp_in_traits) %>%
  pivot_wider(names_from = svspp, values_from = strat_mean_biomass_kg, values_fill = 0, names_sort =T)
shelf.mat <- as.matrix(shelf[,-c(1:2)])
cwm_mat <- CWM(T_mat = trait.mat, S_mat = shelf.mat)

cwm_df <- shelf[,c(1:2)] %>%
  bind_cols(as.data.frame(cwm_mat)) %>% 
  as_tibble()

write.csv(cwm_df, "Data/Derived/CWM_trawl_dataset_wholeshelf.csv", row.names = F, quote = F)








