library(rstan)
library(tidybayes)
library(tidyverse)


#-------------------------------
## Data and modeling
#-------------------------------

df <- readRDS("Data/Derived/EPU_scale_metrics.rds") %>%
  ungroup() %>%
  mutate(epu_season = paste(epu, season, sep = "-"), 
         epu_season.num = as.numeric(as.factor(epu_season)), 
         year.char = year, 
         year.num = as.numeric(as.factor(year))
  )
    
datstan <- with(df, {
  x = list(N = nrow(df),
           T = length(unique(year)), 
           M = length(unique(epu_season)), 
           epu = epu_season.num,
           year = year.num,
           fric = (fric - mean(fric))/sd(fric), 
           feve = (feve - mean(feve))/sd(feve), 
           fdiv = (fdiv - mean(fdiv))/sd(fdiv)
  )})


ggplot(df, aes(x = fric, y = feve))+
  geom_point()

ggplot(df, aes(x = feve, y = fdiv))+
  geom_point()

mod <- rstan::stan_model("Code/STAN_models/lv_EPUscale.stan")

fitm1 <- rstan::sampling(mod, 
                         data = datstan,
                         iter=5000, 
                         chains=3, 
                         # thin = 2, 
                         cores = 3, 
                         control=list(adapt_delta=0.99)
)

shinystan::launch_shinystan(fitm1)

write_rds(fitm1, "Data/Derived/fit_epuscale.rds", compress = "gz")

#-------------------------------
## Build figures
#-------------------------------

fitm1 %>%
  tidybayes::recover_types() %>%
  tidybayes::gather_draws(beta_fric, beta_feve, beta_fdiv) %>%
  group_by(.variable) %>%
  # median_qi(.value) %>%
  ggplot(aes(y = .variable, x = .value))+
  tidybayes::geom_halfeyeh()+
  geom_vline(xintercept = 0, lty = 3)+
  theme_bw()

ggsave("figures/beta_draws.png")


nu_df <- fitm1 %>%
  tidybayes::recover_types()%>%
  tidybayes::gather_draws(nu[epu, year]) %>%
  group_by(year, epu) %>%
  median_qi(.value) %>%
  left_join(df %>% select(year, epu, epu.f, year.char) %>% distinct())

ggplot(nu_df, aes(x = year, y = .value))+
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
