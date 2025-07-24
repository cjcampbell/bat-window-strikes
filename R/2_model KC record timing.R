
source("R/0_funs.R")

df_discovery2 <- read.csv("out/data_derived/structured_surveys_bats_discovered.csv")
sd <- read.csv("out/data_derived/structured_surveys_schedule.csv")

# Find counts / survey date
dat_surv <- left_join(sd, df_discovery2) %>% 
  dplyr::filter(survey == TRUE) %>% 
  count(date, species) %>% 
  pivot_wider(names_from = species, values_from = n, names_prefix = "count_") %>% 
  dplyr::select(-`count_NA`) %>% 
  replace(is.na(.), 0) %>% 
  rowwise() %>% 
  mutate(
    count_total_bats = rowSums(across(starts_with("count_"))),
    presence_any_bats = as.numeric(count_total_bats>0),
    # count_total_bats_noEPFUVespers = rowSums(across(c("count_Evening", "count_Eastern Red", "count_Tricolored", "count_Silver-haired"))),
    # presence_noEPFUVespers = as.numeric(count_total_bats_noEPFUVespers>0),
    yday = yday(date),
    year = year(date)
  )

library(splines)
library(glmmTMB)
library(sjPlot)

# Timing models -----------------------------------------------------------
library(brms)

## Presence model (all species)
m_presence_all_bernoulli <- brm(
  bf( presence_any_bats ~ ns(yday, 5) ),
  data = dat_surv,
  family = bernoulli(),
  cores = 4,
  chains = 4,
  seed = 42,
  iter = 5000,
  warmup = 500,
  threads = threading(4),
  backend = "cmdstanr",
  file_refit = "on_change",
  file = "out/models/m_presence_all_bernoulli.rds"
)
conditional_effects(m_presence_all_bernoulli)
brms::pp_check(m_presence_all_bernoulli,ndraws = 100)
bayes_R2(m_presence_all_bernoulli)


## Abundance model (all species) -------

m_abundance_all_poi <- brm(
  bf( count_total_bats ~ ns(yday, 5)  + (1|year) ),
  data = dat_surv,
  family = poisson(),
  adapt_delta = 0.99,
  save_pars = save_pars(all = TRUE),
  cores = 4,
  chains = 4,
  seed = 42,
  iter = 5000,
  warmup = 1000,
  threads = threading(4),
  backend = "cmdstanr",
  file_refit = "on_change",
  file = "out/models/m_abundance_all_poi.rds"
)
m_abundance_all_poi <- add_criterion(
  m_abundance_all_poi, c("loo"), moment_match = TRUE, recompile = TRUE)

conditional_effects(m_abundance_all_poi)
brms::pp_check(m_abundance_all_poi,ndraws = 100)
bayes_R2(m_abundance_all_poi)

m_abundance_all_ng <- brm(
  bf( count_total_bats ~ ns(yday, 5)  + (1|year) ),
  data = dat_surv,
  family = negbinomial(),
  adapt_delta = 0.99,
  save_pars = save_pars(all = TRUE),
  cores = 4,
  chains = 4,
  seed = 42,
  iter = 5000,
  warmup = 1000,
  threads = threading(4),
  backend = "cmdstanr",
  file_refit = "on_change",
  file = "out/models/m_abundance_all_ng.rds"
)
m_abundance_all_ng <- add_criterion(
  m_abundance_all_ng, c("loo"), moment_match = TRUE, recompile = TRUE)

loo_compare(m_abundance_all_poi, m_abundance_all_ng)

plot(conditional_effects(m_abundance_all_ng, method = "posterior_predict", prob = 0.90),
     ask = FALSE,
     points = T,
     offset = T) 
plot(conditional_effects(m_abundance_all_ng,  method = "posterior_epred"),
     ask = FALSE,
     points = T,
     offset = T) 

p_epred <- conditional_effects(m_abundance_all_ng,  method = "posterior_epred")[[1]] %>% 
  ggplot() +
  aes(x = yday) +
  # geom_col(data = dat_surv, aes(y = count_total_bats)) +
  # geom_rug(data = dat_surv, aes(y = count_total_bats), sides = "b") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_path(aes(y = estimate__)) +
  geom_point(data = dat_surv, aes(y = count_total_bats), alpha = 0.25) +
  scale_y_continuous("Count of bats discovered by survey") +
  scale_x_continuous(
    "Day of year",
    breaks = monthDayYear_to_yday(monthFirsts),
    labels = format(mdy(monthFirsts), "%b")
  )

conditional_effects(m_abundance_all_ng)
brms::pp_check(m_abundance_all_ng,ndraws = 100)
bayes_R2(m_abundance_all_ng)

plot_model(m_abundance_all_poi)

pred_df <- expand.grid(
  yday = seq(min(dat_surv$yday), max(dat_surv$yday), by = 1),
  year = seq(2019, 2025, by = 1)
  )

pred_out <- predict(
  m_abundance_all_poi, 
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  newdata = pred_df
) %>% 
  cbind(pred_df, .)

ggplot(pred_out) +
  aes(x = yday) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "#ABD4E7") +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), fill = "#54A8BD") +
  geom_path(aes(y = Q50), color = "#023047", linewidth = 1) +
  geom_point(data = dat_surv, aes(y = count_total_bats)) +
  facet_wrap(~year)

plot(conditional_effects(m_abundance_all_poi,  method = "posterior_epred"),
     ask = FALSE,
     points = T,
     offset = T) 

plot(conditional_effects(m_abundance_all_poi, method = "posterior_predict", prob = 0.95),
     ask = FALSE,
     points = T,
     offset = T) 

library(tidybayes)
p_re_year <- MCMCvis::MCMCchains(m_abundance_all_ng) %>% 
  as.data.frame() %>% 
  dplyr::select(starts_with("year")) %>% 
  pivot_longer(cols = everything()) %>% 
  dplyr::mutate(name = gsub(",Intercept\\]", "", gsub("year\\[", "", name))) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.1) +
  stat_halfeye(aes(x = value, y = name), alpha = 0.5) +
  scale_y_discrete("Year of survey") +
  scale_x_continuous("Posterior estimates")


# library(patchwork)
# 
# {p_waffle | guide_area() |
#   {p_yday_species + theme(legend.position = "none")}} /
#   (p_epred | p_re_year) +
#   patchwork::plot_annotation(tag_levels = 'A') 
# 


# Add species to the model?


# dat_long <- dat_surv %>% 
#   dplyr::select(-c(count_total_bats, presence_any_bats)) %>% 
#   pivot_longer(
#     cols = c("count_Big brown", "count_Evening", "count_Vespertilionidae", 
#              "count_Eastern Red", "count_Tricolored", "count_Silver-haired"),
#     names_prefix = "count_",
#     names_to = "species",
#     values_to = "count"
#   )



## Add other predictors to model? ------

# moon_phase <- suncalc::getMoonIllumination(dat_surv$date, keep =c("fraction", "phase")) %>% 
#   rename(moon_phase = fraction)
# # moontimes0 <- suncalc::getMoonTimes(
# #   dat_surv$date, 
# #   lat = 39.10,
# #   lon = -94.58,
# #   keep = c("rise", "set"),
# #   tz = "America/Chicago"
# # )
# 
# dat_surv2 <- full_join(dat_surv, moon_phase)
# 
# m_abundance_all_poi2 <- update(
#   m_abundance_all_poi,
#   bf( count_total_bats ~ ns(yday, 5)  + (1|year) + ns(moon_phase, 2)),
#   newdata = dat_surv2
#   )
# plot(conditional_effects(m_abundance_all_poi2, method = "posterior_predict", prob = 0.95),
#      ask = FALSE,
#      points = T,
#      offset = T)

# Other Plots -------------------------------------------------------------------

# Figure -------------
# Lots of bats found close together and on the same day
# Look at the temporal clustering


# Make a map


# Lakeside nature center ------

# "updated 3-19"
