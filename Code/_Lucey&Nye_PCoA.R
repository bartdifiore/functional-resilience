# Load necessary libraries
library(vegan)
library(tidyverse)

# Assuming 'mat' is your community matrix and 'df_wide' contains metadata with 'year'

df_wide <- df_wide %>% 
  mutate(period = case_when(year <= 1985 ~ "early", 
                            year > 2009 ~ "late", 
                            .default = "middle"))

# Step 1: Perform NMDS
set.seed(123)
nmds <- metaMDS(mat, distance = "bray", trymax = 100)

# Step 1: Compute distance matrix (e.g. Bray-Curtis)
dist_mat <- vegdist(mat, method = "bray")

# Step 2: Run adonis2 to assess the effect of 'year' on community composition
adonis_results <- adonis2(mat ~ period*EPU*season, data = df_wide, permutations = 999, by = "terms")
print(adonis_results)


group_var <- paste(df_wide$period, df_wide$season, df_wide$EPU, sep = "-")
# Step 3: Run betadisper
bd <- betadisper(dist_mat, group = group_var)
permutest(bd)

plot(bd)
boxplot(bd)

mod.HSD <- TukeyHSD(bd)
mod.HSD
plot(mod.HSD)

# get betadisper dataframes ####
# have written functions to grab the necessary data from the betadisper object

# functions ####
# getting distances from betadisper() object
betadisper_distances <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(distances = unlist(model$distances))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# getting eigenvalues out of betadisper() object
betadisper_eigenvalue <- function(model){
  temp <- data.frame(eig = unlist(model$eig))
  temp$PCoA <- row.names(temp)
  row.names(temp) <- NULL
  return(temp)
}

# getting the eigenvectors out of a betadisper() object
betadisper_eigenvector <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(unlist(model$vectors))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# get centroids
betadisper_centroids <- function(model){
  temp <- data.frame(unlist(model$centroids))
  temp$group <- row.names(temp)
  temp <- dplyr::select(temp, group, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# betadisper data
get_betadisper_data <- function(model){
  temp <- list(distances = betadisper_distances(model),
               eigenvalue = betadisper_eigenvalue(model),
               eigenvector = betadisper_eigenvector(model),
               centroids = betadisper_centroids(model))
  return(temp)
}

# get betadisper data ####
betadisper_dat <- get_betadisper_data(bd)

# do some transformations on the data
betadisper_dat$eigenvalue <- mutate(betadisper_dat$eigenvalue, percent = eig/sum(eig))

# add convex hull points ####
# this could be put in a function
betadisper_dat$chull <- group_by(betadisper_dat$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(betadisper_dat$centroids, group, PCoA1, PCoA2), select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group')) %>% 
  rownames_to_column() %>% 
  separate(group, into = c("period", "season", "EPU"), sep = "-")

beta_points <- betadisper_dat$eigenvector %>% 
  select(PCoA1, PCoA2, group) %>% 
  separate(group, into = c("period", "season", "EPU"), sep = "-")

beta_centroids <- betadisper_dat$centroids %>% 
  select(PCoA1, PCoA2, group) %>% 
  separate(group, into = c("period", "season", "EPU"), sep = "-")

beta_hull <- betadisper_dat$chull %>% 
  separate(group, into = c("period", "season", "EPU"), sep = "-")

# Now the dataframes are all ready to be completely customisable in ggplot
# plot betadispersion plot

top_row <- beta_points %>% 
  bind_cols(year = df_wide$year) %>%
  ggplot()+
  geom_point(aes(PCoA1, PCoA2, col = year, shape = EPU), size = 3)+
  geom_path(aes(PCoA1, PCoA2, group = EPU, colour = year), linewidth = 1)+
  coord_cartesian(xlim = c(-0.65, 0.4), ylim = c(-0.5, 0.4))+
  scale_color_gradient()+
  facet_wrap(~season)+
  theme_bd()
  

bottom_row <- ggplot() +
  geom_path(aes(PCoA1, PCoA2, col = period, group = EPU), beta_hull, alpha = 0.5) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = period, group = EPU), betadisper_lines, alpha = 0.5) +
  geom_point(aes(PCoA1, PCoA2, col = period, shape = EPU), beta_centroids, size = 4) +
  geom_point(aes(PCoA1, PCoA2, col = period, shape = EPU), beta_points) +
  scale_color_manual(values = c(early = "#132B43", middle = "#346D9D", late = "#56B1F7"))+
  coord_cartesian(xlim = c(-0.65, 0.4), ylim = c(-0.5, 0.4))+
  facet_wrap(~season)+
  theme_bd()

cowplot::plot_grid(top_row, bottom_row, nrow =2)

# plot distances from centroid
ggplot(betadisper_dat$distances, aes(group, distances, fill = group, col = group)) +
  geom_boxplot(aes(fill = group, col = group), outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.55)) +
  stat_summary(position = position_dodge(width = 0.55), geom = 'crossbar', fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))}) +
  geom_point(aes(group, distances, col = group), shape = 21, fill ='white', position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.2)) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  # scale_color_manual('', values = c('black', 'grey'), labels = c("Grazed", 'Ungrazed')) +
  # scale_fill_manual('', values = c('black', 'grey'), labels = c("Grazed", 'Ungrazed')) +
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('') 





top_row <- beta_points %>% 
  bind_cols(year = df_wide$year) %>%
  filter(EPU != "MAB") %>%
  ggplot()+
  geom_point(aes(PCoA1, PCoA2, col = year, shape = EPU), size = 3)+
  geom_path(aes(PCoA1, PCoA2, group = EPU, colour = year), linewidth = 1)+
  # coord_cartesian(xlim = c(-0.65, 0.4), ylim = c(-0.5, 0.4))+
  # scale_color_gradient()+
  facet_wrap(~season)+
  theme_bd()


bottom_row <- ggplot() +
  geom_path(aes(PCoA1, PCoA2, col = period, group = EPU), beta_hull, alpha = 0.5) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = period, group = EPU), betadisper_lines, alpha = 0.5) +
  geom_point(aes(PCoA1, PCoA2, col = period, shape = EPU), beta_centroids, size = 4) +
  geom_point(aes(PCoA1, PCoA2, col = period, shape = EPU), beta_points) +
  scale_color_manual(values = c(early = "#132B43", middle = "#346D9D", late = "#56B1F7"))+
  coord_cartesian(xlim = c(-0.65, 0.4), ylim = c(-0.5, 0.4))+
  facet_wrap(~season)+
  theme_bd()

cowplot::plot_grid(top_row, bottom_row, nrow =2)

































# Step 3: Extract NMDS scores (coordinates)
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$period <- df_wide$period
nmds_scores$EPU <- df_wide$EPU
nmds_scores$season <- df_wide$season

# Step 4: Calculate group centroids
centroids <- nmds_scores %>%
  group_by(period, EPU, season) %>%
  summarise(across(c(NMDS1, NMDS2), mean), .groups = 'drop')

# Step 5: Plot NMDS with group centroids and dispersion ellipses
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = factor(period), shape = EPU)) +
  geom_point(alpha = 0.7, size = 3) +
  stat_ellipse(aes(group = paste(period, EPU)), level = 0.95, alpha = 0.2) +
  geom_point(data = centroids, aes(x = NMDS1, y = NMDS2), 
             shape = 4, size = 5, stroke = 1.5, color = "black") +
  labs(title = "NMDS Plot with Adonis2 Results", 
       x = "NMDS1", y = "NMDS2", color = "Year") +
  facet_wrap(~season)+
  theme_minimal()



nmds_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = period, shape = EPU)) +
  geom_point(size = 3) +  # Plot points
  stat_ellipse(aes(group = interaction(period, EPU)), level = 0.95, linetype = "dashed", size = 1) +  # Add ellipses
  theme_minimal() +
  labs(title = "NMDS: Fish Community by Period and EPU",
       x = "NMDS1", y = "NMDS2") +
  theme(legend.position = "bottom")
