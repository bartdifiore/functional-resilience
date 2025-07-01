library(rstan)
library(tidybayes)
library(tidyverse)
source("Code/theme.R")

#-------------------------------
## Data and modeling
#-------------------------------

df <- readRDS("Data/Derived/EPU_scale_metrics.rds") %>%
  ungroup() %>%
  mutate(epu.f = as.factor(epu),
         year.char = year,
         year = as.numeric(as.factor(year)), 
         epu = as.numeric(as.factor(epu)), 
         season.f = as.factor(season), 
         season = as.numeric(as.factor(season)))

datstan <- with(df, {
  x = list(N = nrow(df),
           T = length(unique(year)), 
           M = length(unique(epu)), 
           S = length(unique(season)),
           epu = epu,
           year = year,
           season = season,
           fric = (fric - mean(fric))/sd(fric), 
           feve = (feve - mean(feve))/sd(feve), 
           fdiv = (fdiv - mean(fdiv))/sd(fdiv)
  )})


ggplot(df, aes(x = fric, y = feve))+
  geom_point()

ggplot(df, aes(x = feve, y = fdiv))+
  geom_point()

mod <- rstan::stan_model("Code/STAN_models/lv_EPUscale_wseason.stan")

fitm1 <- rstan::sampling(mod, 
                         data = datstan,
                         iter=5000, 
                         chains=4, 
                         # thin = 2, 
                         cores = 4, 
                         control=list(adapt_delta=0.99, max_treedepth = 12)
)

shinystan::launch_shinystan(fitm1)

write_rds(fitm1, "Data/Derived/fit_epuscale_wseason.rds", compress = "gz")

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
  scale_y_discrete(labels = c('Functional divergence','Functional evenness','Functional richness'))+
  labs(x = "Posterior estimate", y = "")+
  theme_bw()

ggsave("figures/beta_draws.png", width = 5, height = 4)


nu_df <- fitm1 %>%
  tidybayes::recover_types()%>%
  tidybayes::spread_draws(nu[i]) %>%
  group_by(i) %>%
  median_qi(nu) %>%
  bind_cols(df %>% ungroup() %>% select(year.char, epu.f, season.f))

ggplot(nu_df, aes(x = as.numeric(year.char), y = nu))+
  geom_line(aes(color = epu.f, linetype = season.f))+
  geom_ribbon(aes(ymin = .lower, ymax = .upper, group = season.f), alpha = 0.25)+
  facet_wrap(~epu.f)+
  labs(x = "", y = "Standardized functional resilience indicator", linetype = "Season", color = "EPU")+
  theme_bd()

ggsave("Figures/seasonal_nu.png", width = 8.5, height = 6)



bt <- ecodata::bottom_temp %>%
  filter(Var == "bottom temp anomaly in situ") %>%
  select(EPU, Value, Time) %>%
  rename(year.char = Time, 
         epu.f = EPU, 
         bottom_temp = Value) %>%
  mutate(year.char = as.character(year.char))

nu_df %>%
  left_join(bt) %>%
  ggplot(aes(x = as.numeric(year.char), y = nu))+
  geom_line(aes(color = epu.f, linetype = season.f))+
  geom_ribbon(aes(ymin = .lower, ymax = .upper, group = season.f), alpha = 0.25)+
  geom_line(aes(y = bottom_temp), color = "red")+
  facet_wrap(~epu.f)+
  labs(x = "", y = "Standardized functional resilience indicator", linetype = "Season", color = "EPU")+
  theme_bd()

ggsave("Figures/seasonal_nu_wtemp.png")

temp <- nu_df %>%
  left_join(bt)


cor.test(temp$nu,temp$bottom_temp)
