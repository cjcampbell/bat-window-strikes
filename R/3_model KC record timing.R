
source("R/0_funs.R")
library(patchwork)
library(splines)
library(glmmTMB)
library(sjPlot)

df_discovery2 <- read.csv("data/derived/structured_surveys_bats_discovered.csv") %>% 
  mutate(date = as_date(date))
sd <- read.csv("data/derived/structured_surveys_schedule.csv") %>% 
  mutate(date = as_date(date))

# Find counts / survey date
df_discovery3 <- left_join(sd, df_discovery2) %>% 
  dplyr::filter(survey == TRUE) 
dat_surv <- df_discovery3 %>% 
  count(date, species) %>% 
  pivot_wider(names_from = species, values_from = n, names_prefix = "count_") %>% 
  dplyr::select(-`count_NA`) %>% 
  replace(is.na(.), 0) %>% 
  rowwise() %>% 
  mutate(
    count_total_bats = rowSums(across(starts_with("count_"))),
    presence_any_bats = as.numeric(count_total_bats>0),
    yday = yday(date),
    year = year(date)
  ) %>% 
  ungroup

# Summarize record timing --------------

dat_surv %>% 
  dplyr::summarise(
    count_total_bats = sum(count_total_bats), .by = "yday"
  ) %>% 
  dplyr::filter(count_total_bats != 0) %>% 
  arrange(yday) %>% 
  mutate(
    month = yDay_to_Month(yday),
    season = case_when(yday < 160 ~ "spring", yday >= 160 ~ "autumn")
  ) %>% 
  group_by(season) %>% 
  uncount(count_total_bats) %>% 
  dplyr::summarise(
    q02.5 = quantile(yday, probs = c(0.025)),
    q50.0 = quantile(yday, probs = c(0.5)),
    q97.5 = quantile(yday, probs = c(0.975))
  ) %>% 
  pivot_longer(cols = -season, names_to = "quantile", values_to = "yday") %>% 
  mutate(
    monthDay = yDay_to_monthDay(yday)
  )



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
  # file_refit = "on_change",
  file = "out/models/m_presence_all_bernoulli.rds"
)
conditional_effects(m_presence_all_bernoulli)
brms::pp_check(m_presence_all_bernoulli,ndraws = 100)
bayes_R2(m_presence_all_bernoulli)


## Abundance model (all species) -------
### poisson -----
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
  # file_refit = "always",
  file = "out/models/m_abundance_all_poi.rds"
)
m_abundance_all_poi <- add_criterion(
  m_abundance_all_poi, c("loo"), moment_match = TRUE, recompile = TRUE)

conditional_effects(m_abundance_all_poi)
brms::pp_check(m_abundance_all_poi,ndraws = 100)
bayes_R2(m_abundance_all_poi)

### Neg binom ----
m_abundance_all_ng <- brm(
  bf( count_total_bats ~ ns(yday, 5)  + (1|year) ),
  data = dat_surv,
  family = negbinomial(),
  adapt_delta = 0.95,
  save_pars = save_pars(all = TRUE),
  cores = 4,
  chains = 4,
  seed = 40,
  iter = 5000,
  warmup = 1000,
  threads = threading(4),
  backend = "cmdstanr",
  # file_refit = "always",
  file = "out/models/m_abundance_all_ng.rds"
)
m_abundance_all_ng <- add_criterion(
  m_abundance_all_ng, c("loo"), moment_match = TRUE, recompile = TRUE)

### Compare ----
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
  scale_y_continuous("Bats discovered by survey") +
  scale_x_continuous(
    "Day of year",
    breaks = monthDayYear_to_yday(monthFirsts),
    labels = format(mdy(monthFirsts), "%b")
  )
p_epred
ggsave(p_epred, filename = "figs/timing_posterior.png", dpi = 600, width= 8, height = 8)

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

pp <- plot(conditional_effects(m_abundance_all_poi,  method = "posterior_epred"),
     ask = FALSE,
     points = T,
     offset = T) 

pp[[1]] + 
  ylab("Bats found")

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
p_re_year

ggsave(p_re_year, filename = "figs/timing_posterior_year.png", dpi = 600, width = 6, height = 6)



# Make combination Plot ------------

# Load plots from script 2 for inclusion in multi-panel figure w/ model results.
p_KC_records_barplot <- readRDS("tmp/p_KC_records_barplot.rds")
p_yday_species_list <- readRDS("tmp/p_yday_species_list.rds")

f3_combo22 <- 
    p_yday_species_list[[1]] +
    {p_yday_species_list[[2]]} +
    p_yday_species_list[[3]] +
    p_yday_species_list[[4]] +
    free(p_KC_records_barplot) +
    free(p_epred) + free(p_re_year) +
    plot_layout(
      design = "
    112233
    445555
    666777
    ",
      axis_titles = "collect_y",
      heights = c(1,1,1.5),
      widths = c(1,1,1)
    ) +
  # plot_annotation(
  #   tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & 
  theme(
    plot.tag.position  = c(0.1,0.95),
    plot.tag = element_text(size = 10)
    )

ggsave(f3_combo22, filename = "figs/f3_combo2.png", width = 8, height = 7, dpi = 600)
ggsave(f3_combo22, filename = "figs/f3_combo2.svg", width = 8, height = 7)
