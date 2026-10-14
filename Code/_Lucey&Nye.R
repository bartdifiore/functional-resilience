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

#-----------------------------------
## NEFSC Trawl Survey data
#-----------------------------------

year_filter <- c(seq(1963, 1969, by = 1), 2017, 2020)

nmfs_df_raw <- readRDS("Data/NMFS_trawl/survdat_pull_20250206.rds")$survdat

epu_strata <- read.csv("Data/Derived/strata_to_epu.csv")

df <- nmfs_df_raw %>% 
  janitor::clean_names() %>%
  as_tibble() %>%
  group_by(across(-c(catchsex, abundance, biomass))) %>% # Species that can be sexed (like DogFish, inverts) have two or more separate entries in a tow, one for each category. Here I'm grouping all other variables in order to summerize the total (sum()) abundance or biomass. 
  summarize(abundance = sum(abundance), 
            biomass_kg = sum(biomass)) %>%
  mutate(svspp_num = as.numeric(svspp)) %>%
  drop_na(lat, lon) %>%
  mutate(haul_id = paste(cruise6, station, sep = "_")) %>%
  filter(svspp %in% spp) %>%
  filter(!year %in% year_filter) %>%
  left_join(epu_strata) %>%
  drop_na(EPU)
  
formerge <- df %>%
  ungroup() %>%
  select(year, season, stratum, EPU, haul_id) %>%
  distinct() %>%
  group_by(year, season, EPU, stratum) %>%
  mutate(n_tows_per_stratum = n()) %>% # Count the number of tows conducted in each stratum in each season in each year in each EPU
  ungroup() %>%
  group_by(year, season, EPU) %>%
  mutate(n_tows_per_season = n()) # Count the number of tows conducted total in each EPU in each season in each year.
  
df_long <- df %>% 
  ungroup() %>%
  select(year, season, EPU, stratum, haul_id, svspp, biomass_kg) %>%
  left_join(formerge) %>% # Merge in the information on sampling effort
  group_by(year, season, EPU, stratum, svspp, n_tows_per_stratum, n_tows_per_season) %>%
  summarize(mean_biomass_kg = mean(biomass_kg)) %>% # Estimate the average biomass caught in each strata for each species in each season, year.
  mutate(temp_mean_biomass_kg = n_tows_per_stratum/n_tows_per_season * mean_biomass_kg) %>% # Correct each estimate by the proportional sampling effort in that strata in that season, year, EPU.
  ungroup() %>%
  group_by(year, season, EPU, svspp) %>% 
  summarize(strat_mean_biomass_kg = sum(temp_mean_biomass_kg)) # Total the contibution of each strata to the overal estimate for the EPU.

write_rds(df_long, "Data/Derived/stratified_mean_biomass.rds")

df_NEShelf <- df %>% 
  ungroup() %>%
  select(year, season, EPU, stratum, haul_id, svspp, biomass_kg) %>%
  left_join(formerge) %>% # Merge in the information on sampling effort
  group_by(year, season, EPU, stratum, svspp, n_tows_per_stratum, n_tows_per_season) %>%
  summarize(mean_biomass_kg = mean(biomass_kg)) %>% # Estimate the average biomass caught in each strata for each species in each season, year.
  mutate(temp_mean_biomass_kg = n_tows_per_stratum/n_tows_per_season * mean_biomass_kg) %>% # Correct each estimate by the proportional sampling effort in that strata in that season, year, EPU.
  ungroup() %>%
  group_by(year, season, svspp) %>% 
  summarize(strat_mean_biomass_kg = sum(temp_mean_biomass_kg, na.rm = T))

write_rds(df_NEShelf, "Data/Derived/stratified_mean_biomass_whole_shelf.rds")

df_wide <- df_long %>% 
  filter(EPU != "SS") %>%
  mutate(svspp_names = paste("X", svspp, sep = "")) %>%
  select(-svspp) %>%
  pivot_wider(names_from = "svspp_names", values_from = "strat_mean_biomass_kg", names_sort = T) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(id = paste(year, season, EPU, sep = "-"), .before = year) %>%
  column_to_rownames("id")

mat <- as.matrix(df_wide[, -c(1:3)])

# 1. Remove rows (samples) that contain only zeros
mat_clean <- mat[rowSums(mat, na.rm = TRUE) > 0, ]

# 2. Remove columns (species) that contain only zeros
mat_clean <- mat_clean[, colSums(mat_clean, na.rm = TRUE) > 0]

# 3. Handle literal NA values if any remain (Bray-Curtis cannot handle NA)
# If you have actual missing values, drop those rows or impute them first

# 4. Rerun metaMDS with the clean data
# We set previous.best to allow it to build off your last good run if desired,
# or simply increase trymax to give it more attempts to converge.
ord <- metaMDS(mat_clean, dist = "bray", trymax = 250)
ord
stressplot(ord)
plot(ord)


plot_df <- scores(ord, display = "sites") %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  separate(rowname, into = c("year", "season", "EPU"), sep = "-") %>%
  mutate(year = as.numeric(year))

p_fall <- plot_df %>% 
  filter(season == "FALL") %>% 
  arrange(EPU, year) %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = EPU, color = year), size = 3)+
  geom_path(aes(group = EPU, colour = year))+
  coord_cartesian(ylim = c(-1.1, 2), xlim = c(-1.1,2))+
  theme_bd()+
  labs(title = "Fall")

p_spring <- plot_df %>% 
  filter(season == "SPRING") %>% 
  arrange(EPU, year) %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = EPU, color = year), size = 3)+
  geom_path(aes(group = EPU, colour = year))+
  coord_cartesian(ylim = c(-1.1, 2), xlim = c(-1.1,2))+
  theme_bd()+
  labs(title = "Spring")

p_fall_smooth <- plot_df %>% 
  filter(season == "FALL") %>% 
  group_by(EPU) %>%
  mutate(NMDS1 = zoo::rollapply(NMDS1, width = 5, FUN = mean, fill = NA, align = 'center'), 
         NMDS2 = zoo::rollapply(NMDS2, width = 5, FUN = mean, fill = NA, align = "center")) %>%
  arrange(EPU, year) %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = EPU, color = year), size = 3)+
  geom_path(aes(group = EPU, colour = year), linewidth = 2)+
  coord_cartesian(ylim = c(-1.1, 2), xlim = c(-1.1,2))+
  theme_bd()+
  labs(title = "Fall")

p_spring_smooth <- plot_df %>% 
  filter(season == "SPRING") %>% 
  group_by(EPU) %>%
  mutate(NMDS1 = zoo::rollapply(NMDS1, width = 5, FUN = mean, fill = NA, align = 'center'), 
         NMDS2 = zoo::rollapply(NMDS2, width = 5, FUN = mean, fill = NA, align = "center")) %>%
  arrange(EPU, year) %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = EPU, color = year), size = 3)+
  geom_path(aes(group = EPU, colour = year), linewidth = 2)+
  coord_cartesian(ylim = c(-1.1, 2), xlim = c(-1.1,2))+
  theme_bd()+
  labs(title = "Spring")

cowplot::plot_grid(p_spring, p_spring_smooth + theme(legend.position = "none"), p_fall, p_fall_smooth+ theme(legend.position = "none"), nrow = 2, rel_widths = c(0.55, 0.45, 0.55, 0.45))

ggsave("Figures/replica_of_luceynye.png")


plot_df <- scores(ord, display = "sites") %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  separate(rowname, into = c("year", "season", "EPU"), sep = "-") %>%
  mutate(year = as.numeric(year))


p_fall <- plot_df %>% 
  filter(season == "FALL") %>% 
  arrange(EPU, year) %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = EPU, color = year), size = 3)+
  geom_path(aes(group = EPU, colour = year))+
  geom_point(data = . %>% filter(year == 2008), aes(shape = EPU), color = "red", size = 3, show.legend = F)+
  coord_cartesian(ylim = c(-0.6, 0.6), xlim = c(-1.1,1.8))+
  theme_bd()+
  labs(title = "Fall")

p_spring <- plot_df %>% 
  filter(season == "SPRING") %>% 
  arrange(EPU, year) %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = EPU, color = year), size = 3)+
  geom_path(aes(group = EPU, colour = year))+
  geom_point(data = . %>% filter(year == 2008), aes(shape = EPU), color = "red", size = 3, show.legend = F)+
  coord_cartesian(ylim = c(-0.6, 0.6), xlim = c(-1.1,1.8))+
  theme_bd()+
  labs(title = "Spring")

cowplot::plot_grid(p_spring + theme(legend.position = "none"), p_fall)

ggsave("Figures/ordi_plot.png", width = 6.5, height = 2.5, scale = 1.75)


sp_scores <- scores(ord, display = "species") %>%
  as.data.frame() %>%
  rownames_to_column("svspp_names") %>%
  mutate(SVSPP = as.numeric(str_remove(svspp_names, "^X")))

# Join in species names/common names if you have a lookup table
sp_scores <- sp_scores %>%
  left_join(ecodata::species_groupings %>% select(SVSPP, COMNAME) %>% distinct())  # adjust to your actual species metadata table

# Rank species by how extreme they are along each axis
top_nmds1_pos <- sp_scores %>% arrange(desc(NMDS1)) %>% slice_head(n = 15)
top_nmds1_neg <- sp_scores %>% arrange(NMDS1) %>% slice_head(n = 15)
top_nmds2_pos <- sp_scores %>% arrange(desc(NMDS2)) %>% slice_head(n = 15)
top_nmds2_neg <- sp_scores %>% arrange(NMDS2) %>% slice_head(n = 15)


# Keep species present in at least ~10% of samples
prevalence <- colSums(mat_clean > 0) / nrow(mat_clean)
mat_filtered <- mat_clean[, prevalence >= 0.10]

fit <- envfit(ord, mat_filtered, permutations = 999)

fit_df <- data.frame(
  svspp_names = names(fit$vectors$r),
  r2 = fit$vectors$r,
  pval = fit$vectors$pvals,
  NMDS1 = fit$vectors$arrows[,1],
  NMDS2 = fit$vectors$arrows[,2]
) %>%
  mutate(SVSPP = as.numeric(str_remove(svspp_names, "^X"))) %>%
  left_join(ecodata::species_groupings %>% select(SVSPP, COMNAME) %>% distinct()) %>%  # bring in common names
  arrange(desc(r2))

top_drivers <- fit_df %>%
  filter(pval <= 0.05) %>%
  arrange(desc(r2)) %>%
  slice_head(n = 20)

df_long %>%
  filter(svspp %in% top_drivers$SVSPP, EPU != "SS") %>%
  left_join(ecodata::species_groupings %>% select(SVSPP, COMNAME) %>% distinct() %>% rename(svspp = "SVSPP")) %>%
  ggplot(aes(x = year, y = strat_mean_biomass_kg, color = EPU)) +
  geom_line(alpha = 0.4) +
  geom_smooth(se = FALSE) +
  facet_wrap(~COMNAME, scales = "free_y") +
  theme_bd()

run_simper_comparison <- function(mat, site_ids, season_filter, 
                                  historic_epu, historic_years,
                                  recent_epu, recent_years) {
  
  historic_pattern <- paste0("^(", paste(historic_years, collapse = "|"), ")-", 
                             season_filter, "-", historic_epu, "$")
  recent_pattern <- paste0("^(", paste(recent_years, collapse = "|"), ")-", 
                           season_filter, "-", recent_epu, "$")
  
  historic_rows <- str_detect(site_ids, historic_pattern)
  recent_rows <- str_detect(site_ids, recent_pattern)
  
  target_rows <- historic_rows | recent_rows
  group_labels <- ifelse(historic_rows[target_rows], 
                         paste0("historic_", historic_epu), 
                         paste0("recent_", recent_epu))
  
  sim <- simper(mat[target_rows, ], group_labels, permutations = 999)
  
  list(sim = sim, 
       n_historic = sum(historic_rows), 
       n_recent = sum(recent_rows),
       labels = table(group_labels))
}

historic_years <- as.character(1973:1999)
recent_years <- as.character(2010:2024)
site_ids <- rownames(mat_filtered)
sum(rowSums(mat_filtered) == 0)

comp_mab_gb <- run_simper_comparison(mat_filtered, site_ids, "SPRING",
                                     "MAB", historic_years,
                                     "GB", recent_years)

comp_gb_gom <- run_simper_comparison(mat_filtered, site_ids, "SPRING",
                                     "GB", historic_years,
                                     "GOM", recent_years)

comp_mab_gom <- run_simper_comparison(mat_filtered, site_ids, "SPRING",
                                      "MAB", historic_years,
                                      "GOM", recent_years)

comp_baseline_mab <- run_simper_comparison(mat_filtered, site_ids, "SPRING",
                                           "MAB", historic_years,
                                           "MAB", recent_years)  # historic MAB vs recent MAB

comp_baseline_gb <- run_simper_comparison(mat_filtered, site_ids, "SPRING",
                                          "GB", historic_years,
                                          "GB", recent_years)  # historic GB vs recent GB

comp_baseline_gom <- run_simper_comparison(mat_filtered, site_ids, "SPRING",
                                           "GOM", historic_years,
                                           "GOM", recent_years)  # historic GOM vs recent GOM

tidy_simper <- function(sim_result, comparison_name) {
  s <- summary(sim_result)
  # summary(sim) returns a list named by the comparison, e.g. "historic_MAB_recent_GB"
  comp_df <- s[[1]] %>%
    as.data.frame() %>%
    rownames_to_column("svspp_names") %>%
    mutate(SVSPP = as.numeric(str_remove(svspp_names, "^X")),
           comparison = comparison_name) %>%
    left_join(ecodata::species_groupings %>% select(SVSPP, COMNAME) %>% distinct()) %>%
    arrange(desc(average))
  comp_df
}

mab_gb_tidy <- tidy_simper(comp_mab_gb$sim, "historic_MAB_vs_recent_GB")
gb_gom_tidy <- tidy_simper(comp_gb_gom$sim, "historic_GB_vs_recent_GOM")
mab_gom_tidy <- tidy_simper(comp_mab_gom$sim, "historic_MAB_vs_recent_GOM")

all_simper <- bind_rows(mab_gb_tidy, gb_gom_tidy, mab_gom_tidy) %>%
  mutate(similarity_type = case_when(
    abs(ava - avb) / pmax(ava, avb) < 0.25 ~ "similar_abundance_both_groups",
    TRUE ~ "differs_between_groups"
  ))

extract_overall_dissim <- function(comparison_result, comparison_name) {
  sim <- comparison_result$sim
  sim_comp <- sim[[1]]  # pull out the single pairwise comparison
  
  tibble(
    comparison = comparison_name,
    historic_epu = names(comparison_result$labels)[1],
    recent_epu = names(comparison_result$labels)[2],
    n_historic = comparison_result$n_historic,
    n_recent = comparison_result$n_recent,
    overall_dissimilarity = sim_comp$overall,
    percent_dissimilar = round(sim_comp$overall * 100, 1)
  )
}

dissim_summary <- bind_rows(
  extract_overall_dissim(comp_mab_gb, "historic_MAB_vs_recent_GB"),
  extract_overall_dissim(comp_gb_gom, "historic_GB_vs_recent_GOM"),
  extract_overall_dissim(comp_mab_gom, "historic_MAB_vs_recent_GOM"),
  extract_overall_dissim(comp_baseline_mab, "historic_MAB_vs_recent_MAB"),
  extract_overall_dissim(comp_baseline_gb, "historic_GB_vs_recent_GB"),
  extract_overall_dissim(comp_baseline_gom, "historic_GOM_vs_recent_GOM")
) %>%
  arrange(percent_dissimilar)

View(dissim_summary)

extract_top_species <- function(comparison_result, n_top = 5) {
  sim <- comparison_result$sim
  sim_comp <- sim[[1]]
  
  tibble(
    svspp_names = names(sim_comp$average),
    average = sim_comp$average,
    ava = sim_comp$ava,   # mean abundance in group A (historic)
    avb = sim_comp$avb    # mean abundance in group B (recent)
  ) %>%
    mutate(SVSPP = as.numeric(str_remove(svspp_names, "^X"))) %>%
    left_join(ecodata::species_groupings %>% select(SVSPP, COMNAME) %>% distinct(), 
              by = "SVSPP") %>%
    arrange(desc(average)) %>%
    slice_head(n = n_top)
}

top_species_string <- function(comparison_result, n_top = 5) {
  top_sp <- extract_top_species(comparison_result, n_top)
  
  top_sp <- top_sp %>%
    mutate(direction = ifelse(avb > ava, "↑", "↓"),  # ↑ = higher in recent group, ↓ = higher in historic
           label = paste0(COMNAME, " (", direction, ")"))
  
  paste(top_sp$label, collapse = ", ")
}

dissim_summary <- dissim_summary %>%
  rowwise() %>%
  mutate(top_species = case_when(
    comparison == "historic_MAB_vs_recent_GB" ~ top_species_string(comp_mab_gb),
    comparison == "historic_GB_vs_recent_GOM" ~ top_species_string(comp_gb_gom),
    comparison == "historic_MAB_vs_recent_GOM" ~ top_species_string(comp_mab_gom)
  )) %>%
  ungroup()

dissim_summary %>% View()

dissim_summary_nested <- dissim_summary %>%
  rowwise() %>%
  mutate(top_species_detail = list(case_when(
    comparison == "historic_MAB_vs_recent_GB" ~ list(extract_top_species(comp_mab_gb)),
    comparison == "historic_GB_vs_recent_GOM" ~ list(extract_top_species(comp_gb_gom)),
    comparison == "historic_MAB_vs_recent_GOM" ~ list(extract_top_species(comp_mab_gom))
  )[[1]])) %>%
  ungroup()

dissim_summary_nested %>% View()



#-----------------------------------------
## 1. Setup: distance matrix and helper indices
#-----------------------------------------

dist_full <- vegdist(mat_clean, method = "bray")
dist_mat  <- as.matrix(dist_full)

# Historic MAB reference rows (fixed comparison target)
hist_mab_rows <- str_detect(site_ids, "^(19(7[0-9]|8[0-5]))-SPRING-MAB$")

# Sanity check - should be > 0
sum(hist_mab_rows)

#-----------------------------------------
## 2. Year-by-year trend: GB distance to fixed historic MAB centroid
#-----------------------------------------

gb_to_hist_mab_dist <- function(target_year) {
  gb_row <- which(site_ids == paste0(target_year, "-SPRING-GB"))
  if (length(gb_row) == 0) return(NA_real_)
  mean(dist_mat[gb_row, hist_mab_rows])
}

gb_years <- site_ids[str_detect(site_ids, "-SPRING-GB$")] %>%
  str_extract("^[0-9]{4}") %>%
  as.numeric() %>%
  sort()

trend_df <- tibble(year = gb_years) %>%
  rowwise() %>%
  mutate(mean_dist_to_hist_MAB = gb_to_hist_mab_dist(year)) %>%
  ungroup() %>%
  drop_na(mean_dist_to_hist_MAB)

# Plot the trend
trend_df %>%
  ggplot(aes(x = year, y = mean_dist_to_hist_MAB)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(y = "Mean Bray-Curtis distance to historic MAB",
       x = "Year (GB, Spring)",
       title = "GB spring assemblage distance to historic MAB reference") +
  theme_bd()

# Formal trend test
trend_test <- cor.test(trend_df$year, trend_df$mean_dist_to_hist_MAB, method = "pearson")
trend_test

trend_lm <- lm(mean_dist_to_hist_MAB ~ year, data = trend_df)
summary(trend_lm)

#-----------------------------------------
## 3. Endpoint comparison: betadisper + adonis2
##    historic_MAB vs. historic_GB vs. recent_GB
#-----------------------------------------

historic_gb_rows <- str_detect(site_ids, "^(19(7[0-9]|8[0-5]))-SPRING-GB$")
recent_gb_rows   <- str_detect(site_ids, "^(20[1-2][0-9])-SPRING-GB$")

target_rows <- hist_mab_rows | historic_gb_rows | recent_gb_rows

group_labels <- case_when(
  hist_mab_rows[target_rows]    ~ "historic_MAB",
  historic_gb_rows[target_rows] ~ "historic_GB",
  recent_gb_rows[target_rows]   ~ "recent_GB"
)

# Sanity check group sizes
table(group_labels)

dist_sub <- vegdist(mat_filtered[target_rows, ], method = "bray")

bd <- betadisper(dist_sub, group_labels)

# Centroid-to-centroid distances (this is your key convergence number)
centroid_distances <- as.matrix(dist(bd$centroids))
centroid_distances

# Dispersion check (make sure differences aren't just a dispersion artifact)
bd
anova(bd)          # tests if groups differ in dispersion
plot(bd)           # visualize centroids + spread

# PERMANOVA - tests whether group centroids differ significantly
adonis_result <- adonis2(dist_sub ~ group_labels, permutations = 999)
adonis_result

# Pairwise PERMANOVA if adonis2's overall test is significant
# (requires pairwiseAdonis package, or manual pairwise subsetting)
# install.packages("devtools"); devtools::install_github("pmartinezarbizu/pairwiseAdonis")
# pairwiseAdonis::pairwise.adonis2(dist_sub ~ group_labels, data = data.frame(group_labels))

#-----------------------------------------
## 4. Summary table combining both pieces of evidence
#-----------------------------------------

convergence_summary <- tibble(
  comparison = c("historic_GB_vs_recent_GB", "historic_MAB_vs_recent_GB"),
  centroid_distance = c(
    centroid_distances["historic_GB", "recent_GB"],
    centroid_distances["historic_MAB", "recent_GB"]
  )
) %>%
  mutate(distance_change = centroid_distance[2] - centroid_distance[1],
         interpretation = ifelse(distance_change < 0, 
                                 "GB moved closer to historic MAB", 
                                 "GB moved farther from historic MAB"))

convergence_summary

