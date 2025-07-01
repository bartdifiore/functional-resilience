# Load the vegan package
library(vegan)


bc_dist <- vegdist(mat, method = "bray")

# Step 3: Hierarchical clustering
hc <- hclust(bc_dist)  # UPGMA clustering

# Step 4: Plot the dendrogram
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "", sub = "")

plot(hc, labels = FALSE, main = "Dendrogram (No Site Labels)", xlab = "", sub = "")
rect.hclust(hc, k = 5, border = 2:4)

# Optional: Cut the dendrogram into groups (e.g., 3 clusters)
groups <- cutree(hc, k = 5)

# Step 5: Plot nMDS with cluster groups
ordiplot(ord, type = "n", main = "nMDS with Cluster Overlay")
points(ord, col = groups, pch = 19)
legend("topright", legend = paste("Group", 1:3), col = 1:3, pch = 19)

# Optional: Draw convex hulls around groups
test <- ordihull(ord, groups, draw = "polygon", col = 1:3, alpha = 0.3, label = TRUE)
as.data.frame(test)



plot_df <- scores(ord, display = "sites") %>% 
  as.data.frame() %>%
  bind_cols(clust = groups) %>%
  rownames_to_column() %>% 
  separate(rowname, into = c("year", "season", "EPU"), sep = "-") %>%
  mutate(year = as.numeric(year), 
         clust = as.factor(clust)) %>% 
  group_by(clust)

hulls <- plot_df %>%
  group_by(clust) %>%
  slice(chull(NMDS1, NMDS2))


# # Calculate the hulls for each group
# hull_cyl <- mtcars %>%
#   group_by(cyl) %>%
#   slice(chull(mpg, wt))
# 
# # Update the plot with a fill group, and overlay the new hulls
# p + aes(fill = factor(cyl)) + geom_polygon(data = hull_cyl, alpha = 0.5)


p_fall <- plot_df %>% 
  filter(season == "FALL") %>% 
  arrange(EPU, year) %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = EPU, fill = as.factor(clust)), size = 3)+
  scale_shape_manual(values = c(21,22,23))+
  geom_path(aes(group = EPU, colour = year))+
  geom_polygon(data = hulls, aes(fill = as.factor(clust)), alpha = 0.1)+
  coord_cartesian(ylim = c(-1.1, 1.5), xlim = c(-1.1,2))+
  theme_bd()+
  labs(title = "Fall")

p_spring <- plot_df %>% 
  filter(season == "SPRING") %>% 
  arrange(EPU, year) %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = EPU, color = year), size = 3)+
  geom_path(aes(group = EPU, colour = year))+
  coord_cartesian(ylim = c(-1.1, 1.5), xlim = c(-1.1,2))+
  theme_bd()+
  labs(title = "Spring")

cowplot::plot_grid(p_spring + theme(legend.position = "none"), p_fall)
  