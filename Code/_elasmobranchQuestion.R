

df <- readRDS("Data/Derived/stratified_mean_biomass.rds")


library(ecodata)


?ecodata::species_groupings
View(data.frame(unique(species_groupings$EMAX)))

df %>% left_join(species_groupings, by = join_by(svspp == SVSPP)) %>% 
  mutate(elasmo_notelasmo = ifelse(EMAX %in% c("Piscivore Elasmobranchs", "Benthivore Elasmobranchs", "Omnivore Elasmobranchs", "Southern Elasmobranchs"), "elasmobranchs", "ray-finned fishes")) %>%
  filter(EPU != "SS") %>%
  group_by(year, season, EPU, elasmo_notelasmo) %>% 
  summarize(total = sum(strat_mean_biomass_kg)) %>%
  ggplot(aes(x = year, y = total))+
    geom_line(aes(color = EPU))+
    facet_wrap(~season + elasmo_notelasmo)

elasmo_plot <- df %>% left_join(species_groupings, by = join_by(svspp == SVSPP)) %>% 
  mutate(elasmo_notelasmo = ifelse(EMAX %in% c("Piscivore Elasmobranchs", "Benthivore Elasmobranchs", "Omnivore Elasmobranchs", "Southern Elasmobranchs"), "elasmobranchs", "ray-finned fishes")) %>%
  filter(EPU != "SS") %>%
  group_by(year, season, EPU, elasmo_notelasmo) %>% 
  summarize(total = sum(strat_mean_biomass_kg, na.rm = T)) %>%
  pivot_wider(names_from = elasmo_notelasmo, values_from = total) %>%
  mutate(ratio = elasmobranchs / `ray-finned fishes`, 
         season = factor(season, levels = c("SPRING", "FALL"))) %>%
  ggplot(aes(x = year, y = ratio))+
  geom_line(aes(color = EPU), size = 1)+
  scale_color_manual(values = c("GOM" = "#0072B2", "GB" = "#E69F00", "MAB" = "#009E73")) +
  facet_wrap(~season)+
  labs(x = "", y = "Elasmobranchs /\nRay-finned fishes")+
  theme_bd()

ggsave("Figures/elasmobranch_question.png", elasmo_plot)
