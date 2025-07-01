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
         offspring_size = log(offspring_size)) %>%
  select(-c(habitat, spawning_type, trophic_group))

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

trait_cats <- data.frame(trait_name = names(trait.mat), trait_type = c("Q", "Q", "Q", "Q", "Q", "Q", "Q", "Q"))

#---------------------------------------------------------------------
## Computing distances between species based on functional traits
#---------------------------------------------------------------------

fspace <- mFD::tr.cont.fspace(
  sp_tr        = trait.mat, 
  pca          = TRUE,
  scaling      = "scale_center",
  compute_corr = "pearson")

prcomp_result <- prcomp(trait.mat, scale. = T)
print(prcomp_result)
summary(prcomp_result)
str(prcomp_result)
prcomp_result$x[1:10]

fspace$tr_correl$r

#---------------------------------------------------------------------
## Computing functional spaces & their quality
#---------------------------------------------------------------------

fspace$"quality_metrics"

fspace$"eigenvalues_percentage_var"

dist_mat <- as.matrix(fspace$sp_dist_multidim$"3D")
dist_mat[1:5, 1:5]

fspace$"tr_correl"

#---------------------------------------------------------------------
## Test correlation between functional axes and traits
#---------------------------------------------------------------------

sp_faxes_coord <- fspace$sp_faxes_coord # These match the PCA performed with prcomp. 


# This function doesn't appear to work for PCA analyses... the eigenvalues are all positive but SHOULD be negative in a PCA. Use the PCA to build out the figure.
tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = trait.mat,
  tr_nm          = names(trait.mat), 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

# Print traits with significant effect:
tr_faxes$"tr_faxes_stat"[which(tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]
tr_faxes$"tr_faxes_plot"

#---------------------------------------------------------------------
## Plot the functional space
#---------------------------------------------------------------------

big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord,
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
  sp_faxes_coord   = sp_faxes_coord[ , c("PC1", "PC2", "PC3")],
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
  pivot_longer(cols = sp_richn:fspe) %>%
  ggplot(aes(x = year.n,y = value))+
  geom_line(aes(color = epu))+
  facet_grid(name~season, scales = "free")


df_out <- fd_ind_values %>% 
  rownames_to_column(var = "id") %>% 
  separate(id, into = c("year", "epu", "season"), sep = "-") %>%
  as_tibble()

psych::pairs.panels(df_out[c(4:12)])


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

