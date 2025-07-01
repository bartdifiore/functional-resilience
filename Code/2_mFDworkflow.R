#--------------------------
## Libraries
#--------------------------

library(tidyverse)
library(mFD)

#--------------------------
## Get Data
#--------------------------

trait.mat <- readRDS("Data/Derived/trait-df.rds") %>% 
  mutate(habitat = as.factor(habitat), 
         spawning_type = as.factor(spawning_type), 
         trophic_group = as.factor(trophic_group)) %>% 
  filter(fecundity < 1.0e8) %>% # filter out the one species with the fecundity outlier
  mutate(fecundity = log(fecundity),
         offspring_size = log(offspring_size))

psych::pairs.panels(trait.mat)

summary(trait.mat)

survey_df <- readRDS("Data/Derived/survey-df.rds")

# So one of the things that is stumping me is that the mFD package seems to only accept a single observatio per experimental unit. I can't figure out how to deal with multiple observations of the same community. Here I will average to the EPU-year scale.
survey.mat <- survey_df %>% 
  select(-`265`) %>%
  group_by(EPU, est_year, season) %>% 
  summarize(across(where(is.double), \(x) mean(x, na.rm = T))) %>%
  mutate(id = paste(est_year, EPU, season, sep = "-")) %>% 
  ungroup() %>%
  select(-c(est_year, EPU, season)) %>% 
  column_to_rownames(var = "id") %>%
  as.matrix()

#-----------------------------------------
## Build trait designation df for mFD
#-----------------------------------------

trait_cats <- data.frame(trait_name = names(trait.mat), trait_type = c("N", "N", "N", "Q", "Q", "Q", "Q", "Q", "Q", "Q", "Q"))

#---------------------------------------------------------------------
## Computing distances between species based on functional traits
#---------------------------------------------------------------------

sp_dist <- mFD::funct.dist(
  sp_tr         = trait.mat,
  tr_cat        = trait_cats,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE) # This estimate the species x species matrix describing the distance in multivariate space between each species.

round(sp_dist, 3)   
dim(sp_dist)
str(sp_dist)

#---------------------------------------------------------------------
## Computing functional spaces & their quality
#---------------------------------------------------------------------

# Step 1. Compute multimensional functional spaces and assess their quality

fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist,
  maxdim_pcoa         = 10,
  deviation_weighting = "squared",
  fdist_scaling       = TRUE,
  fdendro             = NULL)

str(fspaces_quality)
round(fspaces_quality$"quality_fspaces", 3) 

# Step 2. Illustrating the quality of the selected functional spaces

mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

#---------------------------------------------------------------------
## Test correlation between functional axes and traits
#---------------------------------------------------------------------

sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = trait.mat,
  tr_nm          = NULL, 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3")], 
  plot           = F)

# Print traits with significant effect:
tr_faxes[which(tr_faxes$p.value < 0.05), ]
tr_faxes$"tr_faxes_stat"[which(tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]
tr_faxes$"tr_faxes_plot"
# 
# tr_faxes <- mFD::traits.faxes.cor(
#   sp_tr          = trait.mat,
#   tr_nm          = names(trait.mat)[c(1:3)], 
#   sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")], 
#   plot           = TRUE)
# tr_faxes$"tr_faxes_plot"
#---------------------------------------------------------------------
## Plot the functional space
#---------------------------------------------------------------------


sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")],
  faxes           = NULL,
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

# Plot the graph with all pairs of axes:
big_plot$patchwork



#-----------------------------------------------------------
## Compute functional diversity indices & plot them
#-----------------------------------------------------------

# Step 1. Functional alpha diversity indices in a multidimensional space

tictoc::tic()
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")],
  asb_sp_w         = survey.mat,
  ind_vect         = c("fdis", 
                       "fmpd", 
                       "fnnd", 
                       "feve", 
                       "fric",
                       "fdiv",
                       "fori", 
                       "fspe", 
                       "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)
tictoc::toc()

fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"

fd_ind_values %>% 
  rownames_to_column(var = "id") %>% 
  separate(id, into = c("year", "epu", "season"), sep = "-") %>%
  as_tibble() %>% 
  mutate(year.n = as.numeric(year)) %>% 
  select(-c(fide_PC1:fide_PC6)) %>%
  pivot_longer(cols = sp_richn:fspe) %>%
  ggplot(aes(x = year.n,y = value))+
  geom_line(aes(color = epu))+
  facet_grid(name~season, scales = "free")


plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("2021-GOM-Fall", "2021-GOM-Spring"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_vert               = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_sp                  = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_vert                = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_ch                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_ch                  = c(pool = "white", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 
plots_alpha$"fric"$"patchwork"


