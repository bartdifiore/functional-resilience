#--------------------------
## Libraries
#--------------------------

library(tidyverse)
library(mFD)
source("Code/1_DataSetup.R")

#------------------
## Get data
#------------------

df <- readRDS("Data/Derived/stratified_mean_biomass.rds")

survey_df <- readRDS("Data/Derived/survey-df.rds") %>% 
  mutate(five_year_period = case_when(year < 1965 ~ 1960, 
                                      year >= 1965 & year < 1970 ~ 1965,
                                      year >= 1970 & year < 1975 ~ 1970, 
                                      year >= 1975 & year < 1980 ~ 1975, 
                                      year >= 1980 & year < 1985 ~ 1980, 
                                      year >= 1985 & year < 1990 ~ 1985, 
                                      year >= 1990 & year < 1995 ~ 1990, 
                                      year >= 1995 & year < 2000 ~ 1995, 
                                      year >= 2000 & year < 2005 ~ 2000, 
                                      year >= 2005 & year < 2010 ~ 2005, 
                                      year >= 2010 & year < 2015 ~ 2010, 
                                      year >= 2015 & year < 2020 ~ 2015, 
                                      year >= 2020 ~ 2020), .before = season) %>%
  select(-year) %>%
  group_by(five_year_period, season, EPU, stratum) %>%
  mutate(num_hauls = n(), .after = haul_id) %>%
  filter(num_hauls >= 5) %>%
  mutate(id = paste(five_year_period, EPU, season, stratum, sep = "-"))

# species_not_caught <- colSums(survey_df[,-c(1:6)]) == 0 
# filt <- names(species_not_caught[species_not_caught == T])
# 
# survey_df <- survey_df %>% select(!all_of(filt))


summary(survey_df$num_hauls)
hist(survey_df$num_hauls)


trait.mat <- readRDS("Data/Derived/trait-df.rds") %>% 
  mutate(habitat = as.factor(habitat), 
         spawning_type = as.factor(spawning_type), 
         trophic_group = as.factor(trophic_group)) %>% 
  # filter(fecundity < 1.0e8) %>% # filter out the one species with the fecundity outlier
  mutate(fecundity = log(fecundity),
         offspring_size = log(offspring_size)) %>%
  select(-c(habitat, spawning_type, trophic_group)) %>%
  mutate(sp_name = paste("X", svspp, sep = "")) %>%
  arrange(sp_name) %>%
  tibble::column_to_rownames(var = "sp_name") %>%
  select(-svspp)

summary(trait.mat)

trait_cats <- data.frame(trait_name = names(trait.mat), trait_type = c("Q", "Q", "Q", "Q", "Q", "Q", "Q", "Q"))

fspace <- mFD::tr.cont.fspace(
  sp_tr        = trait.mat, 
  pca          = TRUE,
  scaling      = "scale_center",
  compute_corr = "pearson")

sp_faxes_coord <- fspace$sp_faxes_coord # These match the PCA performed with prcomp. 


#------------------------------------------
## Average and estimate metrics
#------------------------------------------
sp_faxes_coord <- fspace$sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")]

  
  mat_temp <- survey_df %>%
    group_by(five_year_period, season, EPU, stratum, id) %>%
    summarize(across(where(is.double), \(x) mean(x))) %>%
    ungroup() %>%
    select(-c(five_year_period, EPU, season, stratum)) %>% 
    column_to_rownames(var = "id") %>%
    as.matrix()
  
  alpha_fd_indices <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord,
    asb_sp_w         = mat_temp,
    # ind_vect         = c("fdis", 
    #                      "feve"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)

  fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"
  
  metrics <- fd_ind_values %>% 
    rownames_to_column(var = "id") %>%
    separate(id, into = c("year", "epu", "season", "stratum"), sep = "-") %>%
    as_tibble()

#------------------------------------------
## Plot metrics
#------------------------------------------ 
  
strata_sf <- sf::read_sf("Data/Shapefiles/strata.shp") %>%
    select(STRATA, geometry) %>%
    rename(stratum = STRATA)
  
metrics_sf <- metrics %>% 
  mutate(stratum = as.numeric(stratum)) %>%
  left_join(strata_sf) %>%
  sf::st_as_sf()

metrics_sf %>%
  # filter(year == 2015) %>%
  ggplot()+
  geom_sf(aes(fill = fric))+
  scale_fill_viridis_c()+
  facet_grid(year~season)

metrics_sf %>%
  sf::st_drop_geometry() %>%
  group_by(season, epu, year) %>%
  summarize(fric = mean(fric)) %>%
  mutate(year = as.numeric(year)) %>%
  ggplot()+
    geom_line(aes(x = year, y = fric, color = epu))+
    facet_wrap(~season)

metrics_sf %>%
  sf::st_drop_geometry() %>%
  group_by(season, epu, year) %>%
  summarize(feve = mean(feve)) %>%
  mutate(year = as.numeric(year)) %>%
  ggplot()+
  geom_line(aes(x = year, y = feve, color = epu))+
  facet_wrap(~season)

metrics_sf %>%
  sf::st_drop_geometry() %>%
  group_by(season, epu, year) %>%
  summarize(fdis = mean(fdis)) %>%
  mutate(year = as.numeric(year)) %>%
  ggplot()+
  geom_line(aes(x = year, y = fdis, color = epu))+
  facet_wrap(~season)


metrics_sf %>%
  sf::st_drop_geometry() %>%
  select(-c(fide_PC1, fide_PC2, fide_PC3)) %>%
  pivot_longer(cols = sp_richn:fspe) %>%
  group_by(season, epu, year, name) %>%
  summarize(value = mean(value)) %>%
  ggplot(aes(x = as.numeric(year), y = value))+
  geom_line(aes(color = epu))+
  facet_grid(name ~ season, scales = "free")





























