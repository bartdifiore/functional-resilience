library(tidyverse)
library(sdmTMB) # This is the workhorse to do the spatio-temporal modeling
library(sf)
library(patchwork)
source("Code/theme.R")

#-----------------------------
## Read in the data
#-----------------------------

# EPU x season x year bottom temperature from the gridded product
# (built in Code/_test.R): annual mean, 7-yr centered climate_trend, and
# annual_temp_anomaly = annual mean - climate_trend
bt_means <- read_rds("Data/Derived/bottom_temp_means.rds") %>%
  ungroup()

meta <- read.csv("Data/Derived/trawlmeta_wBT.csv")

epu_strata <- read.csv("Data/Derived/strata_to_epu.csv")

df <- read.csv("Data/Derived/CWM_trawl_dataset.csv") %>%
  as_tibble() %>%
  left_join(meta) %>%
  drop_na(lat, lon) %>%
  add_utm_columns(ll_names = c("lon", "lat"), units = "km") %>%
  left_join(epu_strata) %>%
  drop_na(EPU) %>%
  select(-bottom_temp) %>%
  filter(EPU != "SS") %>%
  left_join(bt_means, by = c("EPU", "season", "year")) %>%
  mutate(EPU = as.factor(EPU),
         season = as.factor(season),
         id = as.factor(paste(EPU, year, season, sep = "-")))

sf_use_s2(FALSE)

#-----------------------------
## Fit one model per trait
#-----------------------------

traits <- c(
  length_maturity = "Length at maturity",
  age_maturity    = "Age at maturity",
  max_obs_length  = "Max. observed length",
  k               = "Growth rate (k)",
  l_inf           = "Asymptotic length (L_inf)",
  trophic_level   = "Trophic level",
  fecundity       = "Fecundity",
  offspring_size  = "Offspring size"
)

out_dir <- "Data/Derived/Model_output"
dir.create(out_dir, showWarnings = FALSE)

fit_trait <- function(trait, refit = FALSE) {
  file <- file.path(out_dir, paste0("m_", trait, "_temp.rds"))
  if (file.exists(file) && !refit) return(read_rds(file))

  # Drop NAs before building the mesh so the response, year index, and
  # spatial projection stay row-aligned
  dat <- df %>%
    drop_na(all_of(trait), climate_trend, annual_temp_anomaly)

  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 30)

  m <- sdmTMB(
    data = dat,
    formula = as.formula(paste(trait, "~ EPU * season * (climate_trend + annual_temp_anomaly) + (1 | id)")),
    mesh = mesh,
    family = gengamma(link = "log"),
    spatial = "on",
    time = "year",
    spatiotemporal = "IID",
    extra_time = c(2017, 2020),
    silent = TRUE
  )
  write_rds(m, file)
  m
}

mods <- map(set_names(names(traits)), fit_trait)

# Convergence check for every trait
map_dfr(mods, ~ as_tibble(sanity(.x, silent = TRUE)), .id = "trait")

#-----------------------------
## Covariate slopes per EPU x season
#-----------------------------

# Slope of a covariate within each EPU x season = sum of the main effect and its
# interactions. Built as a contrast (covariate = 1 vs 0) on the model's own design
# matrix, with SEs from the fixed-effect covariance (delta method).
get_slopes <- function(m, covariates = c("climate_trend", "annual_temp_anomaly")) {
  idx <- names(m$sd_report$par.fixed) == "b_j"
  b <- m$sd_report$par.fixed[idx]
  names(b) <- colnames(m$tmb_data$X_ij[[1]])
  V <- m$sd_report$cov.fixed[idx, idx]

  f <- delete.response(terms(lme4::nobars(m$formula[[1]])))
  epu_lv <- levels(droplevels(m$data$EPU))
  season_lv <- levels(m$data$season)
  grid <- expand_grid(EPU = factor(epu_lv, epu_lv), season = factor(season_lv, season_lv))
  X <- function(g) model.matrix(f, model.frame(f, g))

  map_dfr(covariates, function(cv) {
    g0 <- grid %>% mutate(climate_trend = 0, annual_temp_anomaly = 0)
    g1 <- g0
    g1[[cv]] <- 1
    L <- X(g1) - X(g0)
    stopifnot(identical(colnames(L), names(b)))
    est <- as.vector(L %*% b)
    se <- sqrt(diag(L %*% V %*% t(L)))
    grid %>% mutate(covariate = cv, estimate = est, se = se,
                    conf.low = est - 1.96 * se, conf.high = est + 1.96 * se)
  })
}

slopes <- map_dfr(mods, get_slopes, .id = "trait") %>%
  mutate(
    # log link: convert to % change in the CWM per 1 degC
    across(c(estimate, conf.low, conf.high), ~ 100 * (exp(.x) - 1), .names = "pct_{.col}"),
    sig = conf.low > 0 | conf.high < 0,
    trait_label = factor(traits[trait], levels = rev(traits))
  )

write_csv(slopes, "Data/Derived/temperature_slopes_by_trait.csv")

slopes %>%
  filter(covariate == "climate_trend") %>%
  select(trait, EPU, season, pct_estimate, pct_conf.low, pct_conf.high, sig) %>%
  print(n = 50)

#-----------------------------
## Coefficient plot: climate_trend effects
#-----------------------------

season_cols <- c(FALL = "#eb6834", SPRING = "#2a78d6")
season_shapes <- c(FALL = 16, SPRING = 17)

epu_cols <- c(GB = "#2a78d6", GOM = "#eb6834", MAB = "#1baf7a")
epu_shapes <- c(GB = 16, GOM = 17, MAB = 15)

# Fecundity and offspring size have much larger effects, so they get their own
# panel (and x-axis) below the other traits
repro_traits <- c("fecundity", "offspring_size")

trend <- slopes %>%
  filter(covariate == "climate_trend") %>%
  mutate(trait_label = droplevels(trait_label))

coef_panel <- function(dat, color_var, facet_var, cols, shapes, dodge, show_strip = TRUE) {
  legend_name <- c(season = "Season", EPU = "EPU")[[color_var]]
  # group = reversed levels so dodged rows read top-to-bottom in legend order
  ggplot(dat, aes(x = pct_estimate, y = trait_label,
                  color = .data[[color_var]], shape = .data[[color_var]],
                  group = fct_rev(.data[[color_var]]))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_linerange(aes(xmin = pct_conf.low, xmax = pct_conf.high),
                   linewidth = 0.7, position = position_dodge(width = dodge)) +
    geom_point(size = 2.5, position = position_dodge(width = dodge)) +
    scale_color_manual(values = cols, name = legend_name) +
    scale_shape_manual(values = shapes, name = legend_name) +
    facet_wrap(as.formula(paste("~", facet_var))) +
    labs(x = NULL, y = NULL) +
    theme_bd() +
    theme(strip.background = element_blank(),
          strip.text = if (show_strip) element_text(size = 14) else element_blank(),
          panel.spacing.x = unit(1.5, "lines"),
          panel.grid.major.y = element_line(color = "gray92"),
          legend.position = "top")
}

coef_plot <- function(color_var, facet_var, cols, shapes, dodge) {
  top <- coef_panel(trend %>% filter(!trait %in% repro_traits) %>% droplevels(),
                    color_var, facet_var, cols, shapes, dodge)
  bottom <- coef_panel(trend %>% filter(trait %in% repro_traits) %>% droplevels(),
                       color_var, facet_var, cols, shapes, dodge, show_strip = FALSE) +
    labs(x = expression("Change in CWM per 1 " * degree * "C of climate trend (%)"))

  # Panel heights proportional to the number of traits in each
  (top / bottom) +
    plot_layout(heights = c(6, 2), guides = "collect") &
    theme(legend.position = "top")
}

p_epu <- coef_plot("season", "EPU", season_cols, season_shapes, dodge = 0.6)
p_epu

ggsave("Figures/fig_climate_trend_coefs_byEPU.png", p_epu, width = 10, height = 6, dpi = 300)

# Alternative layout: faceted by season, EPU as color
p_season <- coef_plot("EPU", "season", epu_cols, epu_shapes, dodge = 0.7)
p_season

ggsave("Figures/fig_climate_trend_coefs_bySeason.png", p_season, width = 8, height = 6, dpi = 300)
