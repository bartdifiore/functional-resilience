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
  # filter(fecundity < 1.0e8) %>% # filter out the one species with the fecundity outlier
  # mutate(fecundity = log(fecundity),
  #        offspring_size = log(offspring_size)) %>%
  mutate(across(where(is.double), \(x) log(x))) %>%
  select(-c(habitat, spawning_type, trophic_group)) %>%
  mutate(svspp = paste("X", svspp, sep = "")) %>%
  column_to_rownames("svspp")

GGally::ggpairs(trait.mat, diag = list(continuous = "barDiag"))


summary(trait.mat)

trait_cats <- data.frame(trait_name = names(trait.mat), trait_type = c("Q", "Q", "Q", "Q", "Q", "Q", "Q", "Q"))

species <- ecodata::species_groupings %>% select(COMNAME, SVSPP) %>% janitor::clean_names()

#---------------------------------------------------------------------
## Computing distances between species based on functional traits
#---------------------------------------------------------------------

fspace <- mFD::tr.cont.fspace(
  sp_tr        = trait.mat, 
  pca          = TRUE,
  scaling      = "scale_center",
  compute_corr = "pearson")


sp_faxes_coord <- fspace$sp_faxes_coord # These match the PCA performed with prcomp. 

tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = trait.mat,
  tr_nm          = names(trait.mat), 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

fspaces_quality<- mFD::quality.fspaces(
  sp_dist             = fspace,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")




prcomp_result <- prcomp(trait.mat, scale. = T)
print(prcomp_result)
summary(prcomp_result)
str(prcomp_result)
prcomp_result$x[1:10]
plot(prcomp_result)


pca_df <- as.data.frame(prcomp_result$x[, 1:3]) %>%
  rownames_to_column(var = "svspp") %>%
  separate(svspp, into = c("junk", "svspp"), sep = "X") %>%
  select(-junk) %>%
  mutate(svspp = as.numeric(svspp)) %>%
  as_tibble()

out <- readRDS("Data/Derived/trait-df.rds") %>% 
  mutate(habitat = as.factor(habitat), 
         spawning_type = as.factor(spawning_type), 
         trophic_group = as.factor(trophic_group)) %>% 
  left_join(pca_df %>% mutate(svspp = as.character(svspp)))

write_rds(out, "Data/Derived/trait-df_wPCA.rds")

examples <- pca_df %>% 
  mutate(PC1 = PC1 / sqrt(sum((PC1 - mean(PC1))^2)),
         PC2 = PC2 / sqrt(sum((PC2 - mean(PC2))^2)),
         PC3 = PC3 / sqrt(sum((PC3 - mean(PC3))^2))) %>%
  left_join(species) %>%
  filter(comname %in% c("ATLANTIC HERRING", "SPINY DOGFISH", "ATLANTIC COD")) %>% 
  distinct()


pca_df %>% 
  # mutate(PC1 = PC1 / sqrt(sum((PC1 - mean(PC1))^2)),
  #        PC2 = PC2 / sqrt(sum((PC2 - mean(PC2))^2)),
  #        PC3 = PC3 / sqrt(sum((PC3 - mean(PC3))^2))) %>%
  left_join(species) %>%
  arrange(PC1) %>%
  select(comname, PC1, svspp) %>%
  View()

library(ggfortify)
# autoplot(prcomp_result, 
#          x = 1, 
#          y = 2)


p1 <- autoplot(prcomp_result,
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

autoplot(prcomp_result,
         x = 1, 
         y = 3,
         loadings = T, 
         loadings.label = T, 
         loadings.label.hjust = 1.2, 
         loadings.label.vjust = -.75,
         color = "grey")+
  theme_classic()

p2 <- autoplot(prcomp_result,
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

















