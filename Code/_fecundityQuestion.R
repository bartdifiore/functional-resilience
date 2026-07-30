
# Long story short, this is just to explore the extreme values of fecundity. Ocean Sunfish produce 300 million eggs. The next fish that is commonly encountered is white hake with 8 million eggs and a much smaller biomass. Many of the tows that encountered sunfish had little else in the net. So the CWM estimate is very extreme, but these are only a few outliers.  

traits <- readRDS("Data/Derived/trait-df_wPCA.rds")

traits %>% 
  select(svspp, fecundity) %>%
  mutate(svspp = as.numeric(svspp)) %>%
  left_join(ecodata::species_groupings, by = join_by(svspp == SVSPP)) %>%
  select(svspp, fecundity, COMNAME) %>%
  distinct() %>%
  arrange(desc(fecundity))
  
length(unique(species_groupings$COMNAME))
dim(species_groupings)

length(unique(traits$svspp))
dim(traits)

cwm <- read.csv("Data/Derived/CWM_trawl_dataset.csv")

hist(cwm$fecundity, breaks = 100)
hist(log(cwm$fecundity))
summary(cwm$fecundity)

cwm %>% 
  filter(fecundity == max(fecundity, na.rm = T))

cwm %>%
  arrange(desc(fecundity))


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

df %>% filter(haul_id == "200005_35") %>% View() # So this is the anomanously high fecundity tow. They caught 2 ocean sunfish and nothing else (other than a spider crab and an offshore lizardfish). 

species_groupings %>% filter(SVSPP %in% c(265, 317, 437)) %>% select(SVSPP, COMNAME)

df %>% filter(haul_id == "198306_81") %>% View() # Also 2 sunfish... almost nothing else

df %>% filter(haul_id == "202304_69") %>% View() # 1 sunfish... not much else
