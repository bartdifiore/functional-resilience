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



# 
# dist_mat <- vegdist(mat, method = "bray")
# 
# # 2. Create a grouping variable (e.g., decade)
# group <- paste(df_wide$EPU, df_wide$season, sep = "-")
# 
# # 3. Run betadisper
# bd <- betadisper(dist_mat, group)
# 
# # 3. Test for differences in group centroids
# ad <- adonis2(dist_mat ~ group)
# plot(ad)
# 
# # 4. Test for homogeneity of group dispersions
# anova(bd)
# 
# # 5. Optional: Permutation test
# permutest(bd)
# 
# plot(bd)
# 
# 
# a <- anosim(ord, grouping = group)
# 
# # View results
# summary(a)
# plot(a)
# 
# 
# 
# 
# # Load required libraries
# library(vegan)
# library(tidyverse)
# 
# # Step 1: Compute distance matrix (e.g. Bray-Curtis)
# dist_mat <- vegdist(mat, method = "bray")
# 
# # Step 2: Create grouping variable (e.g. decade)
# # Replace with your own grouping variable, e.g., df_wide$decade
# group_var <- df_wide$decade
# 
# # Step 3: Run betadisper
# bd <- betadisper(dist_mat, group = group_var)
# 
# # Step 4: Extract sample scores (PCoA coordinates)
# scores_df <- as.data.frame(scores(bd, display = "sites")) %>%
#   mutate(group = group)
# 
# # Step 5: Extract group centroids
# centroids_df <- as.data.frame(scores(bd, display = "centroids")) %>%
#   rownames_to_column("group")
# 
# # Step 6: Combine sample scores with centroids
# scores_df <- left_join(scores_df, centroids_df, by = "group", suffix = c("", "_centroid"))
# 
# # Step 7: Plot with ggplot2
# ggplot(scores_df, aes(x = PCoA1, y = PCoA2, color = group)) +
#   geom_point(alpha = 0.7, size = 3) +
#   geom_segment(aes(xend = PCoA1_centroid, yend = PCoA2_centroid), 
#                arrow = arrow(length = unit(0.02, "npc")), alpha = 0.4) +
#   geom_point(data = centroids_df, aes(x = PCoA1, y = PCoA2), 
#              shape = 4, size = 5, stroke = 1.5, color = "black") +
#   labs(title = "Multivariate Dispersion (betadisper)", 
#        x = "PCoA 1", y = "PCoA 2") +
#   theme_minimal()
# 
# 
# 
# nmds_plot <- ggplot(nmds_scores, aes(x = MDS1, y = MDS2, color = period, shape = EPU)) +
#   geom_point(size = 3) +  # Plot points
#   stat_ellipse(aes(group = interaction(period, EPU)), level = 0.95, linetype = "dashed", size = 1) +  # Add ellipses
#   theme_minimal() +
#   labs(title = "NMDS: Fish Community by Period and EPU",
#        x = "NMDS1", y = "NMDS2") +
#   theme(legend.position = "bottom")
# 
# 
# 
# 
# 
# cap_model <- capscale(mat ~ year*EPU, data = df_wide, distance = "bray")
# anova(cap_model)  # Permutation test for overall model
# plot(cap_model)
# 
# 
# 
# ad <- adonis2(mat ~ year*EPU*season, data = df_wide, method = "bray")
# plot(ad)


# Simper analysis

# 1. Run SIMPER between your EPUs (assuming 'epu_vector' holds your North/South labels)
# 'species_matrix' should contain your 350+ species biomass columns
simper_res <- simper(mat, group = df_wide$EPU)

# 2. Extract the summary for your specific North vs. South comparison
simper_summary <- summary(simper_res)

# 3. Pull the top 4 species for GOM vs GB
top_gom_gb <- as.data.frame(simper_res$GB_GOM) %>% 
  arrange(desc(average)) %>% 
  head(10) %>% 
  rownames()

# 4. Pull the top 4 species for GB vs MAB
top_gb_mab <- as.data.frame(simper_res$GB_MAB) %>% 
  arrange(desc(average)) %>% 
  head(10) %>% 
  rownames()

# 5. Combine into a single list of unique driver species
all_drivers <- unique(c(top_gom_gb, top_gb_mab))
all_drivers <- str_remove(all_drivers, "X")

# 7. Plot the trajectories for all driver species
df_long %>% 
  filter(svspp %in% all_drivers, EPU != "SS", season == "FALL") %>%
  ggplot(aes(x = year, y = strat_mean_biomass_kg, color = EPU)) +
    geom_line(size = 1, alpha = 0.7) +
    geom_point(size = 2) +
    facet_wrap(~ svspp, scales = "free_y", ncol = 4) + 
    scale_color_manual(values = c("GOM" = "#1f77b4", "GB" = "#ff7f0e", "MAB" = "#2ca02c")) +
    theme_bw() +
    labs(
      title = "Key Driver Species: Transitioning Across EPUs (Fall)",
      x = "Year",
      y = "Stratified Mean Biomass",
      color = "EPU"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "bottom"
    )
