library(gmRi)
library(tidyverse)
library(sf)
library(ggplot2)
library(vegan)
library(mFD)


#---------------------
## Get data
#---------------------

trait.mat <- readRDS("Data/Derived/trait-df.rds") %>% 
  mutate(habitat = as.factor(habitat), 
         spawning_type = as.factor(spawning_type), 
         trophic_group = as.factor(trophic_group)) %>% 
  # filter(fecundity < 1.0e8) %>% # filter out the one species with the fecundity outlier
  # mutate(fecundity = log(fecundity),
  #        offspring_size = log(offspring_size)) %>%
  mutate(across(where(is.double), \(x) log(x))) %>%
  select(-c(habitat, spawning_type, trophic_group)) %>%
  as.data.frame() %>% 
  mutate(sp_name = paste("X", svspp, sep = "")) %>%
  tibble::column_to_rownames(var = "sp_name") %>%
  select(-svspp)

summary(trait.mat)
row.names(trait.mat)

#---------------------------------------------------------------------
## Computing distances between species based on functional traits
#---------------------------------------------------------------------

fspace <- mFD::tr.cont.fspace(
  sp_tr        = trait.mat, 
  pca          = TRUE,
  scaling      = "scale_center",
  compute_corr = "pearson")

#---------------------------------------------------------------------
## Computing functional spaces & their quality
#---------------------------------------------------------------------

fspace$"quality_metrics"

fspace$"eigenvalues_percentage_var"

dist_mat <- as.matrix(fspace$sp_dist_multidim$"4D")
dist_mat[1:5, 1:5]

fspace$"tr_correl"

#---------------------------------------------------------------------
## Plot functional space
#---------------------------------------------------------------------

#mFD doesn't allow easy access to the trait loadings. Switchign to prcomp (same values) to extract those values

pc <- prcomp(trait.mat, scale. = T)
print(pc)
summary(pc)

pca_df <- as.data.frame(pc$x[, 1:3]) %>%
  rownames_to_column() %>% 
  separate(rowname, into = c("junk", "svspp"), sep = "X") %>%
  select(-junk) %>%
  mutate(svspp = as.integer(svspp)) %>%
  as_tibble()

species <- ecodata::species_groupings %>% 
  select(SVSPP, COMNAME, SCINAME) %>% 
  janitor::clean_names()

examples <- pca_df %>% 
  mutate(PC1 = PC1 / sqrt(sum((PC1 - mean(PC1))^2)), 
         PC2 = PC2 / sqrt(sum((PC2 - mean(PC2))^2)), 
         PC3 = PC3 / sqrt(sum((PC3 - mean(PC3))^2))) %>%
  left_join(species) %>% # filter(PC1 > 0.1 | PC1 < -0.1) %>% V
  filter(comname %in% c("ATLANTIC HERRING", "SPINY DOGFISH", "ATLANTIC COD")) %>%
  distinct()


library(ggfortify)
p1 <- autoplot(pc,
               x = 1, 
               y = 2,
               loadings = T, 
               loadings.label = T,
               loadings.label.hjust = 1.2, 
               loadings.label.vjust = -.75,
               color = "grey")+
  annotate(geom = "point", x = examples$PC1, y = examples$PC2, size = 2, shape = 21, fill = "transparent") +
  annotate(geom = "text", x = examples$PC1 + 0.01, y = examples$PC2 + 0.01, label = examples$comname) +
  theme_classic()

autoplot(pc,
         x = 1, 
         y = 3,
         loadings = T, 
         loadings.label = T, 
         loadings.label.hjust = 1.2, 
         loadings.label.vjust = -.75,
         color = "grey")+
  theme_classic()

p2 <- autoplot(pc,
               x = 2, 
               y = 3,
               loadings = T, 
               loadings.label = T, 
               loadings.label.hjust = 1.2, 
               loadings.label.vjust = -.75,
               color = "grey")+
  annotate(geom = "point", x = examples$PC2, y = examples$PC3, size = 2, shape = 21, fill = "transparent") +
  annotate(geom = "text", x = examples$PC2 + 0.01, y = examples$PC3 + 0.01, label = examples$comname) +
  theme_classic()

cowplot::plot_grid(p1, p2, nrow = 1)
ggsave("Figures/pca_ordination.png", width = 10, height = 5)

















