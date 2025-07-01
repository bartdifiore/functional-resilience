
#----------------------------------------------------
## Estimate metrics from Flensborg et al. 2023
#----------------------------------------------------

df <- readRDS("Data/Derived/survdat_wfG.RDS")

out <- list()
for(i in 1:100){
  hauls_to_sample <- df %>% 
    distinct(haul_id, EPU, hex_id, five_year_period) %>%
    group_by(EPU, hex_id, five_year_period) %>%
    slice_sample(n=5, replace = F) %>%
    pull(haul_id)
  
  df_temp <- df %>%
    select(svspp, haul_id, EPU, hex_id, five_year_period, biomass_kg, abundance, spfgr) %>%
    filter(haul_id %in% hauls_to_sample)
  
  functional_group_richness <- df_temp %>% 
    ungroup() %>%
    group_by(EPU, hex_id, five_year_period) %>%
    summarize(functional_group_richness = n_distinct(spfgr))
  
  mean_fg_redundancy <- df_temp %>% 
    ungroup() %>%
    group_by(EPU, hex_id, five_year_period, spfgr) %>% 
    summarize(n_species_per_fg = n_distinct(svspp)) %>%
    group_by(EPU, hex_id, five_year_period) %>% 
    summarize(mean_fg_redundancy = mean(n_species_per_fg, na.rm = T))
  
  evenness.df <- df_temp %>% 
    select(-biomass_kg) %>%
    ungroup() %>%
    group_by(EPU, hex_id, five_year_period, spfgr) %>%
    pivot_wider(names_from = svspp, values_from = abundance, names_sort = T, values_fill = 0)
  
  evenness.mat <- as.matrix(evenness.df[,-c(1:5)])
  evenness.df$evennessSh <- vegan::diversity(evenness.mat, index = "shannon")
  evenness.df$richness <- vegan::specnumber(evenness.mat)
  
  evenness.df$evennessPj <- evenness.df$evennessSh/log(evenness.df$richness)
  
  evenness.df.out <- evenness.df %>%
    mutate(evennessPj = ifelse(is.nan(evennessPj) == T, NA, evennessPj)) %>%
    group_by(EPU, hex_id, five_year_period) %>%
    summarize(mean_evennessSh = mean(evennessSh, na.rm = T), 
              mean_evennessPj = mean(evennessPj, na.rm = T))
  
  out[[i]] <- functional_group_richness %>% 
    left_join(mean_fg_redundancy) %>%
    left_join(evenness.df.out) %>%
    mutate(.iteration = i) %>%
    as.data.frame()
}

fg_metrics <- do.call("rbind", out) %>% as_tibble()

fg_summary <- fg_metrics %>%
  group_by(EPU, hex_id, five_year_period) %>%
  summarize(fgr_mean = mean(functional_group_richness), 
            fgr_sd = sd(functional_group_richness), 
            fgred_mean = mean(mean_fg_redundancy), 
            fgred_sd = sd(mean_fg_redundancy), 
            fgevenSh_mean = mean(mean_evennessSh), 
            fgevenSh_sd = sd(mean_evennessSh), 
            fgevenPj_mean = mean(mean_evennessPj), 
            fgevenPj_sd = sd(mean_evennessPj)
  )


write_rds(fg_summary, "Data/Derived/fg_metrics.RDS", compress = "gz")

#------------------------------------------------
## Visualize metrics
#------------------------------------------------
hex <- readRDS("Data/Derived/hex.rds")
crs <- st_crs(4326)
fg_summary <- readRDS("Data/Derived/fg_metrics.RDS") %>%
  left_join(hex %>% select(-hauls)) %>%
  st_as_sf()

ne_countries <- rnaturalearth::ne_countries(scale = 10,
                                            continent = "North America",
                                            returnclass = "sf") %>%
  sf::st_transform(crs = crs)

ne_states <- rnaturalearth::ne_states(country = "united states of america", returnclass = "sf") %>%
  sf::st_transform(crs = crs)

xmin = -79
xmax = -63
ymin = 30
ymax = 45

xlims <- c(xmin, xmax)
ylims <- c(ymin, ymax)


fg_summary %>%
  ggplot() +
  geom_sf(aes(fill = fgr_mean)) +
  ggplot2::geom_sf(data = ne_countries, color = "grey60", size = 0.25) +
  ggplot2::geom_sf(data = ne_states, color = "grey60", size = 0.05) +
  ggplot2::coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  scale_fill_distiller(palette = "Spectral") +
  facet_wrap(~five_year_period)+
  labs(title = "Expected functional group richness per five hauls")+
  theme(legend.position = c(0.75, 0.15))

ggsave("Figures/fgr_by_5year.png")

fg_summary %>%
  ggplot() +
  geom_sf(aes(fill = fgred_mean)) +
  ggplot2::geom_sf(data = ne_countries, color = "grey60", size = 0.25) +
  ggplot2::geom_sf(data = ne_states, color = "grey60", size = 0.05) +
  ggplot2::coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  scale_fill_distiller(palette = "Spectral") +
  facet_wrap(~five_year_period)+
  labs(title = "Expected functional group redundancy per five hauls")+
  theme(legend.position = c(0.75, 0.15))

ggsave("Figures/fgred_by_5year.png")

fg_summary %>%
  ggplot() +
  geom_sf(aes(fill = fgevenSh_mean)) +
  ggplot2::geom_sf(data = ne_countries, color = "grey60", size = 0.25) +
  ggplot2::geom_sf(data = ne_states, color = "grey60", size = 0.05) +
  ggplot2::coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  scale_fill_distiller(palette = "Spectral") +
  facet_wrap(~five_year_period)+
  labs(title = "Expected Shannah functional group evenness per five hauls")+
  theme(legend.position = c(0.75, 0.15))

ggsave("Figures/fgevenSh_by_5year.png")

fg_summary %>%
  ggplot() +
  geom_sf(aes(fill = fgevenPj_mean)) +
  ggplot2::geom_sf(data = ne_countries, color = "grey60", size = 0.25) +
  ggplot2::geom_sf(data = ne_states, color = "grey60", size = 0.05) +
  ggplot2::coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  scale_fill_distiller(palette = "Spectral") +
  facet_wrap(~five_year_period)+
  labs(title = "Expected Pj functional group evenness per five hauls")+
  theme(legend.position = c(0.75, 0.15))

ggsave("Figures/fgevenPj_by_5year.png")
  
fg_summary %>% 
  group_by(EPU, five_year_period) %>% 
  summarize(fgr_mean = mean(fgr_mean)) %>%
  st_drop_geometry()%>%
  ggplot(aes(x = five_year_period, y = fgr_mean))+
  geom_line(aes(color = EPU))+
  theme_minimal()

ggsave("Figures/fgr_by_5year_timeseries.png", width = 5, height = 3)

fg_summary %>% 
  st_drop_geometry()%>%
  ggplot(aes(x = five_year_period, y = fgr_mean))+
  geom_line(aes(color = EPU, group = hex_id))+
  facet_wrap(~EPU)+
  theme_minimal()




aes(depth, exp(est),
    ymin = exp(est - 1.96 * est_se),
    ymax = exp(est + 1.96 * est_se)



