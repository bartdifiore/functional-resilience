library(rstan)
library(tidybayes)
library(tidyverse)


#-------------------------------
## Data and modeling
#-------------------------------

df <- readRDS("Data/Derived/fg_metrics.RDS") %>%
  ungroup() %>%
  mutate(site_id = as.numeric(as.factor(hex_id)),
         five_year.f = as.factor(five_year_period),
         epu.f = EPU,
         five_year = as.numeric(as.factor(five_year_period)), 
         epu = as.numeric(as.factor(EPU)))

five_year.f <- unique(df$five_year.f)
epu.f <- unique(df$epu.f)


ggplot(df, aes(x = fgr_mean, y = fgred_mean))+
  geom_point()

ggplot(df, aes(x = fgr_mean, y = fgevenSh_mean))+
  geom_point()


datstan <- with(df, {
  x = list(N = nrow(df),
           K = length(unique(site_id)),
           T = length(unique(five_year.f)), 
           M = length(unique(EPU)), 
           epu = epu,
           year = five_year,
           group = site_id,
           fg_richness = (fgr_mean - mean(fgr_mean))/sd(fgr_mean), 
           fg_evennessSh = (fgevenSh_mean - mean(fgevenSh_mean))/sd(fgevenSh_mean), 
           #fg_evenessPj = (mean_evennessPj - mean(mean_evennessPj))/sd(mean_evennessPj), 
           fg_redundancy = (fgred_mean - mean(fgred_mean))/sd(fgred_mean)
  )})


mod <- rstan::stan_model("Code/STAN_models/lv_mod5b.stan")

fitm1 <- rstan::sampling(mod, 
                         data = datstan,
                         iter=5000, 
                         chains=3, 
                         thin = 2, 
                         cores = 3, 
                         control=list(adapt_delta=0.99, stepsize = 0.0001)
)

shinystan::launch_shinystan(fitm1)

write_rds(fitm1, "Data/Derived/fit_heirarchical2.rds", compress = "gz")


fitm1 %>%
  tidybayes::recover_types() %>%
  tidybayes::gather_draws(beta_fg_richness, beta_fg_evennessSh, beta_fg_redundancy) %>%
  group_by(.variable) %>%
  # median_qi(.value) %>%
  ggplot(aes(y = .variable, x = .value))+
  tidybayes::geom_halfeyeh()+
  geom_vline(xintercept = 0, lty = 3)+
  theme_bw()

ggsave("figures/beta_draws.png")

samps <- rstan::extract(fitm1)

hist(samps$beta_fg_evennessSh)

# 
# nu_pred <- data.frame(year = unique(df$year),
#                       med = apply(samps$nu, 2, median),
#                       lwr = apply(samps$nu, 2, quantile, 0.05),
#                       upr = apply(samps$nu, 2, quantile, 0.95))

nu_df <- fitm1 %>%
  tidybayes::recover_types()%>%
  tidybayes::gather_draws(nu[EPU, year]) %>%
  group_by(year, EPU) %>%
  median_qi(.value) %>%
  mutate(five_year.f = (year*5)+1970, 
         epu.f = case_when(EPU == 1 ~ "GB", 
                           EPU == 2 ~ "GOM", 
                           EPU == 3 ~ "MAB", 
                           EPU == 4 ~ "SS"))

ggplot(nu_df, aes(x = five_year.f, y = .value))+
  geom_line(aes(color = epu.f))+
  geom_ribbon(aes(ymin = .lower, ymax = .upper, group = epu.f), alpha = 0.25)+
  facet_wrap(~epu.f)+
  theme_minimal()



bt <- ecodata::bottom_temp %>%
  filter(Var == "bottom temp anomaly in situ") %>%
  dplry::select(EPU, Value, Time) %>%
  rename(year.f = Time, 
         epu.f = EPU, 
         bottom_temp = Value)

nu_df %>%
  left_join(bt) %>%
  ggplot(aes(x = year.f))+
  geom_line(aes(y = .value)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),fill = "black", alpha = 0.5) +
  # geom_line(aes(y = bottom_temp))+
  labs(x = "Year", y = "Finfish community functional resilience indicator")+
  facet_wrap(~epu.f)+
  theme_bw()
ggsave("figures/resilience_indicator.png")

# Data exploration to compare our estimates of functional resilience (nu) with drivers (temperature) and responses (ecosystem services).  

# Build out a dataframe organized by year and epu with different indicators.....


bt <- ecodata::bottom_temp %>%
  filter(Var == "bottom temp anomaly in situ") %>%
  dplyr::select(EPU, Value, Time) %>%
  rename(year.f = Time, 
         epu.f = EPU, 
         bottom_temp = Value)

bt.glories <- ecodata::bottom_temp_glorys %>%
  dplyr::select(EPU, Value, Time) %>%
  rename(year.f = Time, 
         epu.f = EPU, 
         bottom_temp_glorys = Value)

fogarty <- ecodata::ppr %>%
  filter(Var == "Fogarty") %>%
  dplyr::select(EPU, Value, Time) %>%
  rename(year.f = Time, 
         epu.f = EPU, 
         ppr = Value) %>%
  mutate(ppr = scale(ppr))

# fogarty2 <- ecodata::ppr %>%
#   filter(Var == "PP") %>%
#   dplyr::select(EPU, Value, Time) %>%
#   rename(year.f = Time, 
#          epu.f = EPU, 
#          pp = Value) %>%
#   mutate(pp = scale(pp))

var_names <- unique(aggregate_biomass$Var)[1:14]

total_cpue <- ecodata::aggregate_biomass %>%
  filter(Var %in% var_names) %>%
  group_by(Time, EPU) %>%
  summarize(total_cpue = sum(Value)) %>%
  rename(year.f = Time, 
         epu.f = EPU) %>%
  filter(epu.f != "All") %>%
  ungroup() %>%
  mutate(total_cpue = scale(total_cpue))

other_indicators <- bt %>%
  left_join(bt.glories) %>%
  left_join(fogarty) %>%
  left_join(total_cpue)



# ecodata::aggregate_biomass %>%
#   filter(Var %in% var_names) %>%
#   group_by(Time, EPU) %>%
#   summarize(total_cpue = sum(Value)) %>%
#   ggplot(aes(x = Time, y = total_cpue))+
#   geom_line(aes(color = EPU))

nu_df %>%
  left_join(other_indicators) %>%
  pivot_longer(cols = bottom_temp:total_cpue) %>%
  ggplot(aes(x = year.f))+
  geom_line(aes(y = .value)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),fill = "black", alpha = 0.5) +
  geom_line(aes(y = value, color = name))+
  labs(x = "Year", y = "Fish community resilience indicator")+
  facet_wrap(~epu.f)+
  theme_bw()

cor.df <- nu_df %>%
  left_join(other_indicators)

cor.test(cor.df$.value, cor.df$bottom_temp)
cor.test(cor.df$.value, cor.df$total_cpue)
acf(cor.df$.value)
acf(cor.df$bottom_temp)
acf(cor.df$total_cpue)


cor.df %>%
  ggplot(aes(x = year.f))+
  geom_line(aes(y = .value)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),fill = "black", alpha = 0.5) +
  geom_line(aes(y = total_cpue), col = "red")+
  labs(x = "Year", y = "Fish community resilience indicator")+
  facet_wrap(~epu.f)+
  theme_bw()


df %>% 
  group_by(year.f, epu.f) %>%
  summarize(functional_group_richness = mean(functional_group_richness)) %>%
  ggplot(aes(x = year.f, y = functional_group_richness))+
  geom_line() +
  facet_wrap(~epu.f)+
  labs(title = "Mean functional group richness by EPU")+
  theme_bw()
ggsave("figures/richness_by_epu.png")

df %>% 
  group_by(year.f, epu.f) %>%
  summarize(mean_evennessSh = mean(mean_evennessSh)) %>%
  ggplot(aes(x = year.f, y = mean_evennessSh))+
  geom_line() +
  facet_wrap(~epu.f)+
  coord_cartesian(ylim = c(0,1))