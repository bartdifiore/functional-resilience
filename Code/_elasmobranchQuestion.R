library(tidyverse)

df <- readRDS("Data/Derived/stratified_mean_biomass.rds")


library(ecodata)

?ecodata::species_groupings
View(data.frame(unique(species_groupings$EMAX)))


df %>% left_join(species_groupings, by = join_by(svspp == SVSPP)) %>% 
  mutate(elasmo_notelasmo = ifelse(EMAX %in% c("Piscivore Elasmobranchs", "Benthivore Elasmobranchs", "Omnivore Elasmobranchs", "Southern Elasmobranchs"), "elasmobranchs", "ray-finned fishes")) %>%
  filter(EPU != "SS") %>%
  group_by(year, season, EPU, elasmo_notelasmo) %>% 
  summarize(total = sum(strat_mean_biomass_kg)) %>%
  ggplot(aes(x = year, y = total))+
    geom_line(aes(color = EPU))+
    facet_wrap(~season + elasmo_notelasmo)

elasmo_plot <- df %>% left_join(species_groupings, by = join_by(svspp == SVSPP)) %>% 
  mutate(elasmo_notelasmo = ifelse(EMAX %in% c("Piscivore Elasmobranchs", "Benthivore Elasmobranchs", "Omnivore Elasmobranchs", "Southern Elasmobranchs"), "elasmobranchs", "ray-finned fishes")) %>%
  filter(EPU != "SS") %>%
  group_by(year, season, EPU, elasmo_notelasmo) %>% 
  summarize(total = sum(strat_mean_biomass_kg, na.rm = T)) %>%
  pivot_wider(names_from = elasmo_notelasmo, values_from = total) %>%
  mutate(ratio = elasmobranchs / `ray-finned fishes`, 
         season = factor(season, levels = c("SPRING", "FALL"))) %>%
  ggplot(aes(x = year, y = ratio))+
  geom_line(aes(color = EPU), size = 1)+
  scale_color_manual(values = c("GOM" = "#0072B2", "GB" = "#E69F00", "MAB" = "#009E73")) +
  facet_wrap(~season)+
  labs(x = "", y = "Elasmobranchs /\nRay-finned fishes")+
  theme_bd()

ggsave("Figures/elasmobranch_question.png", elasmo_plot)


df2 <- df %>% 
  left_join(species_groupings, by = join_by(svspp == SVSPP)) %>% 
  mutate(elasmo_notelasmo = ifelse(EMAX %in% c("Piscivore Elasmobranchs", "Benthivore Elasmobranchs", "Omnivore Elasmobranchs", "Southern Elasmobranchs"), "elasmobranchs", "ray-finned fishes")) %>%
  filter(EPU != "SS") %>%
  group_by(EPU, season, year, elasmo_notelasmo) %>%
  summarize(elasmo_biomass = sum(strat_mean_biomass_kg, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(EPU, season, year), names_from = elasmo_notelasmo, values_from = elasmo_biomass) %>%
  mutate(total_biomass = elasmobranchs + `ray-finned fishes`, 
         prop_elasmobiomass = elasmobranchs / total_biomass)

source("Code/theme.R")
ggplot(df2, aes(x = year, y = prop_elasmobiomass))+
  geom_line(aes(color = EPU), size = 1)+
  scale_color_manual(values = c("GOM" = "#0072B2", "GB" = "#E69F00", "MAB" = "#009E73")) +
  facet_wrap(~season)+
  labs(x = "", y = "Elasmobranchs /\nTotal Biomass")+
  theme_bd()







#-----------------------------------------------------
## Base it on the raw data not stratified means
#-----------------------------------------------------

year_filter <- c(seq(1963, 1969, by = 1), 2017, 2020)

nmfs_df_raw <- readRDS("Data/NMFS_trawl/survdat_pull_20250206.rds")$survdat

epu_strata <- read.csv("Data/Derived/strata_to_epu.csv")

traits <- readRDS("Data/Derived/trait-df_wPCA.rds")
spp_in_traits <- unique(traits$svspp)

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
  filter(!year %in% year_filter) %>%
  filter(svspp %in% spp_in_traits) %>%
  left_join(epu_strata)
  

forplot <- df %>% 
  left_join(species_groupings, by = join_by(svspp == SVSPP)) %>% 
  mutate(elasmo_notelasmo = ifelse(EMAX %in% c("Piscivore Elasmobranchs", "Benthivore Elasmobranchs", "Omnivore Elasmobranchs", "Southern Elasmobranchs") | svspp == 5, "elasmobranchs", "ray-finned fishes")) %>%
  filter(EPU != "SS") %>%
  group_by(EPU, season, year, elasmo_notelasmo) %>%
  summarize(elasmo_biomass = sum(biomass_kg, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(EPU, season, year), names_from = elasmo_notelasmo, values_from = elasmo_biomass) %>%
  mutate(total_biomass = elasmobranchs + `ray-finned fishes`, 
         prop_elasmobiomass = elasmobranchs / total_biomass)

ggplot(forplot, aes(x = year, y = prop_elasmobiomass))+
  geom_line(aes(color = EPU), linewidth = 0.8)+
  scale_color_manual(values = c("GOM" = "#0072B2", "GB" = "#E69F00", "MAB" = "#009E73")) +
  facet_wrap(~forcats::fct_rev(season))+
  labs(x = "", y = "Proportion elasmobranch\nbiomass")+
  theme_bd()

ggsave("Figures/elasmobranch_question_prop.png")

df %>% 
  left_join(species_groupings, by = join_by(svspp == SVSPP)) %>% 
  mutate(elasmo_notelasmo = ifelse(EMAX %in% c("Piscivore Elasmobranchs", "Benthivore Elasmobranchs", "Omnivore Elasmobranchs", "Southern Elasmobranchs"), "elasmobranchs", "ray-finned fishes")) %>%
  ungroup() %>%
  filter(stringr::str_detect(COMNAME, "UNCL")) %>%
  distinct(svspp, COMNAME, SCINAME, EMAX) %>% View()





