#--------------------------
## Libraries
#--------------------------

library(tidyverse)
library(mFD)
source("Code/1_DataSetup.R")

#------------------
## Get data
#------------------

source("Code/_PCA.R")

spp_in_pca <- row.names(fspace$sp_faxes_coord)

df <- readRDS("Data/Derived/stratified_mean_biomass.rds") %>%
  mutate(svspp = paste("X", svspp, sep = "")) %>%
  filter(svspp %in% spp_in_pca) %>%
  pivot_wider(names_from = svspp, values_from = strat_mean_biomass_kg, names_sort = T) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(id = paste(year, season, EPU, sep = "-"))

mat <- df  %>%
  column_to_rownames(var = "id") %>%
  select(-c(year, season, EPU)) %>%
  as.matrix()

#------------------------------------------
## Average and estimate metrics
#------------------------------------------
sp_faxes_coord <- fspace$sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")]

colnames(mat)
row.names(sp_faxes_coord)

alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord,
  asb_sp_w         = mat,
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"

metrics <- fd_ind_values %>% 
  rownames_to_column(var = "id") %>%
  separate(id, into = c("year", "season", "epu"), sep = "-") %>%
  as_tibble()

#------------------------------------------
## Plot metrics
#------------------------------------------ 

GGally::ggpairs(metrics %>% select(fdis:fspe), diag = list(continuous = "barDiag"))

metrics %>%
  select(-c(fide_PC1, fide_PC2, fide_PC3)) %>%
  pivot_longer(cols = sp_richn:fspe) %>%
  ggplot(aes(x = as.numeric(year), y = value))+
  geom_line(aes(color = epu))+
  facet_grid(name ~ season, scales = "free")

metrics %>%
  select(c(fdis, feve, fric, fdiv, year, season, epu)) %>%
  pivot_longer(cols = fdis:fdiv) %>%
  ggplot(aes(x = as.numeric(year), y = value))+
  geom_line(aes(color = epu))+
  facet_grid(name ~ season, scales = "free")


write_rds(metrics, "Data/Derived/EPU_scale_metrics.rds")

plots <- expand.grid(year = c(1980, 1990, 2000, 2010, 2019), season = c("FALL", "SPRING"), epu = c("GOM", "GB", "MAB")) 

plots <- paste(plots$year, plots$season, plots$epu, sep = "-")


bd_multi_plot <- function(which.plot){
  plots_alpha <- mFD::alpha.multidim.plot(
    output_alpha_fd_multidim = alpha_fd_indices,
    plot_asb_nm              = which.plot,
    ind_nm                   = c("fric", "feve", "fdiv"),
    faxes                    = c("PC1", "PC2"),
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
  p1 <- plots_alpha$fric$PC1_PC2
  p2 <- plots_alpha$feve$PC1_PC2
  p3 <- plots_alpha$fdiv$PC1_PC2
  cowplot::plot_grid(p1,p2,p3, nrow = 1)
}

upper_row <- cowplot::plot_grid(p1+theme_bd(), p2+theme_bd())
lower_row <- bd_multi_plot(plots[c(1,5)])
cowplot::plot_grid(upper_row, lower_row, nrow = 2)
ggsave("Figures/multivariate_traitfigure.svg")








plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = plots[c(1,5)],
  ind_nm                   = c("fric", "feve", "fdiv"),
  faxes                    = c("PC1", "PC2"),
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
plots_alpha$fric$patchwork







