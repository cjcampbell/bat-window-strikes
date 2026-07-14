# 8_fit iNat models.R
# Purpose: Fit the iNaturalist use-availability models of suspected bat-window
#          collisions against a target-group iNaturalist background. Reads the
#          analysis-ready table from "7_prep iNat models.R", checks predictor
#          collinearity (VIF), fits the use-availability (Bernoulli) model set,
#          compares them by LOO, and builds the manuscript figures. No raster
#          extraction or downloads happen here.
# Input:   data/derived/useavail_points.csv               (from 7_prep iNat models.R)
#          data/derived/inat_background_effort.csv         (effort reference)
# Outputs: out/models/m_useavail_*.rds                     (cached brms fits)
# Figures: figs/B_main_combo_datatop_data.{png,svg}  (main use-availability figure)
#          figs/B_contrasts.{png,svg}                (raw use-vs-available contrasts)
#          figs/SI_effort_validation.{png,svg}
#          figs/SI_effectsize_loo.{png,svg}


# Setup ----
source("R/0_funs.R")
library(sf)
library(rnaturalearth)
library(lubridate)
library(data.table)
library(brms)
library(splines)
library(patchwork)
library(ggtext)

# Fit a use-availability Bernoulli model with the shared brms settings and attach
# LOO. Cached by `file` under out/models/ (delete to refit).
brm_ua <- function(formula, data, file) {
  m <- brm(
    bf(formula), data = data, family = bernoulli(),
    cores = 4, chains = 4, seed = 42, iter = 5000, warmup = 1000,
    threads = threading(4), backend = "cmdstanr",
    file = file.path("out/models", file)
  )
  add_criterion(m, "loo")
}

# Load analysis-ready points ----
# One row per use (collision, used = 1) / available (background, used = 0) point,
# with building height, ALAN, day-of-year, and the nearest-station nightly traffic
# match (traffic = NA where no match within radar range). Built by 7_prep.
points <- fread("data/derived/useavail_points.csv")
points[, date := as.IDate(date)]

# Full-data frame (building height, ALAN, day of year at every point) and the
# radar-matched subset (points with a nightly-traffic match within 200 km). Keeping
# both from one table means the full and radar models draw from the same source.
model_data <- as.data.frame(points[, .(used, building_height, alan, yday)])
radar_data <- as.data.frame(points[!is.na(traffic) & dist_km < 200])
cat(sprintf("model_data: %d (used %d)  |  radar_data: %d (used %d)\n",
            nrow(model_data), sum(model_data$used), nrow(radar_data), sum(radar_data$used)))

# Reconstruct point geometries and the available region for the maps/figures
# (cheap; no rasters). The 100-km geodesic buffer matches 7_prep.
pts_sf        <- st_as_sf(as.data.frame(points), coords = c("longitude", "latitude"),
                          crs = proj.wgs84, remove = FALSE)
collisions_sf <- pts_sf[pts_sf$used == 1, ]
background_sf <- pts_sf[pts_sf$used == 0, ]
background_region <- collisions_sf %>%
  st_buffer(dist = 100e3) %>%
  st_union()
background_region <- st_sf(geometry = background_region)

# Collinearity check (VIF) ----
# Building height and ALAN were flagged a priori as potentially redundant (ALAN is
# an input to composite human-footprint layers). VIF is a property of the design
# matrix, not the Bayesian fit, so we read it off a plain lm() with the full model's
# right-hand side. ns(yday, 5) is a multi-column term, so car::vif() reports a
# generalized VIF (GVIF) and the scale-comparable GVIF^(1/(2*Df)); square that last
# column to compare against the usual VIF thresholds (~1-3 low, 5-10 high). Computed
# on the radar-matched subset, where all four terms coexist (as in m_useavail_radar_all).
vif_lm  <- lm(used ~ log1p(building_height) + log1p(alan) + ns(yday, 5) + log1p(traffic),
              data = radar_data)
vif_tab <- car::vif(vif_lm)
print(vif_tab)

# Pairwise Pearson correlation of the two structural terms on the modelled (log1p)
# scale, reported for both the full data and the radar-matched subset.
cat(sprintf("cor(log1p building height, log1p ALAN): full = %.2f, radar-matched = %.2f\n",
            with(model_data, cor(log1p(building_height), log1p(alan))),
            with(radar_data, cor(log1p(building_height), log1p(alan)))))

# Exploratory contrasts ----
collisionColor  <- "#0072B2"   # Okabe-Ito blue: colourblind-safe, distinct from the green buffer
backgroundColor <- "grey45"
sampleColors <- c(
  "collision site"         = collisionColor,
  "background (available)" = backgroundColor
)
labelSample <- function(used) factor(used, levels = c(0, 1),
                                     labels = c("background (available)", "collision site"))

# Building height: collisions vs the effort background.
(p_bh_density <- ggplot(model_data, aes(x = building_height, colour = labelSample(used))) +
  geom_density() +
  scale_x_continuous("Building height (m)") +
  scale_colour_manual(NULL, values = sampleColors))

# Day of year: collisions concentrate in migration windows; the background traces
# observer effort (the contrast the phenology term formalises).
(p_yday_density <- ggplot(model_data, aes(x = yday, colour = labelSample(used))) +
  geom_density() +
  scale_x_continuous("Day of year",
                     breaks = monthDayYear_to_yday(monthFirsts),
                     labels = format(mdy(monthFirsts), "%b")) +
  scale_colour_manual(NULL, values = sampleColors))

# Nightly traffic on collision vs background nights.
(p_traffic_density <- ggplot(radar_data, aes(x = traffic + 1, colour = labelSample(used))) +
  geom_density() +
  scale_x_log10("Nightly migration traffic (night prior)") +
  scale_colour_manual(NULL, values = sampleColors))

# Fit use-availability models ----
## Full data: structural term + season ----
# Use-availability logistic model: does collision favour taller buildings and
# migration-window timing, relative to where/when people observe? Building height
# is log1p-transformed (right-skewed; ~half the background is 0 m). The intercept is
# not interpretable in a use-availability design (it tracks the use:available
# ratio); coefficients express relative selection.
m_useavail_bh <- brm_ua(used ~ log1p(building_height) + ns(yday, 5),
                        model_data, "m_useavail_bh_yday")
pp_check(m_useavail_bh, ndraws = 100)
bayes_R2(m_useavail_bh)

# Shape check: is log1p adequate, or is the building-height response curved?
m_useavail_bhquad <- brm_ua(used ~ log1p(building_height) + I(log1p(building_height)^2) + ns(yday, 5),
                            model_data, "m_useavail_bhquad_yday")
loo_compare(m_useavail_bh, m_useavail_bhquad)

# ALAN (nighttime radiance; VIIRS VNL v2 2024) as a candidate structural driver in
# place of building height (parallel to m_useavail_bh, same rows for a fair LOO):
# does light alone track collisions? Building-height-plus-ALAN and the radar-matched
# comparisons follow below.
m_useavail_alan <- brm_ua(used ~ log1p(alan) + ns(yday, 5),
                          model_data, "m_useavail_alan_yday")
loo_compare(m_useavail_bh, m_useavail_alan)
# Land cover and distance to water remain candidate covariates (not yet added).

## Radar-matched subset: add night-migration traffic (and ALAN) ----
# Competing models on the radar-matched subset (same rows, for a fair LOO): season
# only, night traffic only, and both. Traffic is log1p-transformed (right-skewed).
# "Both" lets traffic express the night-level anomaly beyond the seasonal mean
# (traffic and ns(yday) share only ~27% of variance). The two ALAN models then test
# whether light, not raw height, carries the built-environment signal: m_useavail_radar_all
# adds ALAN to the top model (does the height coefficient survive?), and
# m_useavail_radar_alan swaps height for ALAN.
m_useavail_season_sub  <- brm_ua(used ~ log1p(building_height) + ns(yday, 5),
                                 radar_data, "m_useavail_season_sub")
m_useavail_radaronly   <- brm_ua(used ~ log1p(building_height) + log1p(traffic),
                                 radar_data, "m_useavail_radaronly")
m_useavail_radar_both  <- brm_ua(used ~ log1p(building_height) + ns(yday, 5) + log1p(traffic),
                                 radar_data, "m_useavail_radar_both")
m_useavail_radar_alan  <- brm_ua(used ~ log1p(alan) + ns(yday, 5) + log1p(traffic),
                                 radar_data, "m_useavail_radar_alan")
m_useavail_radar_all   <- brm_ua(used ~ log1p(building_height) + log1p(alan) + ns(yday, 5) + log1p(traffic),
                                 radar_data, "m_useavail_radar_all")

loo_compare(m_useavail_season_sub, m_useavail_radaronly, m_useavail_radar_both,
            m_useavail_radar_alan, m_useavail_radar_all)
fixef(m_useavail_radar_all)[c("log1pbuilding_height", "log1palan", "log1ptraffic"), ]
conditional_effects(m_useavail_radar_both)

# Diagnostic: nearest-station assignment (each focal point drawn to its station).
# Only runs if the station table is loaded (e.g. after sourcing 7_prep in the same
# session); it is a check, not a manuscript figure.
if (exists("stations")) {
  radar_near <- st_nearest_feature(pts_sf, stations)
  radar_crd  <- st_coordinates(pts_sf)
  radar <- data.frame(
    sample  = ifelse(pts_sf$used == 1, "collision", "background"),
    pt_lon  = radar_crd[, 1], pt_lat = radar_crd[, 2],
    st_lon  = stations$lon[radar_near], st_lat = stations$lat[radar_near],
    dist_km = as.numeric(st_distance(pts_sf, stations[radar_near, ], by_element = TRUE)) / 1000
  ) %>% filter(dist_km < 200)
  (p_station_radar <- ggplot() +
    geom_sf(data = ne_countries(scale = "medium", continent = "North America", returnclass = "sf"),
            fill = "grey97", colour = "grey80", linewidth = 0.2) +
    geom_segment(data = radar, aes(pt_lon, pt_lat, xend = st_lon, yend = st_lat, colour = sample),
                 linewidth = 0.15, alpha = 0.35) +
    geom_sf(data = stations, shape = 17, size = 0.7, colour = "black") +
    scale_colour_manual(values = c(collision = collisionColor, background = "grey45")) +
    coord_sf(xlim = c(-125, -68), ylim = c(12, 52)) +
    theme_void())
}

# Figures ----
# Report-quality drafts. Colours reuse sampleColors; "#159367" marks the buffer.

## Sampling design map ----
# North and South America land with subtle state/province boundaries and the Great
# Lakes, the 100-km "available" buffer, background points (grey), and collisions
# (blue) on top. Projected to a North America Lambert azimuthal equal-area (matches
# the geodesic buffer, so the radii read as circles). State lines and lakes are
# pulled with ne_download() and cached to tmp/. State lines are 1:10m because the
# 1:50m layer omits Mexico's states. South America is included so the northern part
# that falls in frame gets its land, borders, and admin-1 lines (coord_sf clips the rest).
na_land <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(continent %in% c("North America", "South America"))
if (!file.exists("tmp/ne_state_lines_10m.rds")) {
  saveRDS(ne_download(scale = 10, type = "admin_1_states_provinces_lines",
                      category = "cultural", returnclass = "sf"), "tmp/ne_state_lines_10m.rds")
}
na_states <- readRDS("tmp/ne_state_lines_10m.rds")
if (!file.exists("tmp/ne_lakes.rds")) {
  saveRDS(ne_download(scale = 50, type = "lakes", category = "physical", returnclass = "sf"),
          "tmp/ne_lakes.rds")
}
na_lakes <- readRDS("tmp/ne_lakes.rds")
proj_na <- "+proj=laea +lat_0=45 +lon_0=-100 +datum=WGS84 +units=m +no_defs"

(f_map <- ggplot() +
  geom_sf(data = na_land, fill = "grey95", colour = NA) +
  geom_sf(data = na_states, colour = "grey82", linewidth = 0.12) +
  geom_sf(data = na_land, fill = NA, colour = "grey65", linewidth = 0.25) +
  geom_sf(data = na_lakes, fill = "white", colour = "grey82", linewidth = 0.1) +
  geom_sf(data = background_region, fill = "#159367", colour = NA, alpha = 0.15) +
  geom_sf(data = background_sf, aes(colour = "Background observations"), size = 0.2, alpha = 0.3) +
  geom_sf(data = collisions_sf, aes(colour = "Bat collisions"), size = 0.9) +
  scale_colour_manual(NULL,
    values = c("Bat collisions" = collisionColor, "Background observations" = backgroundColor),
    breaks = c("Bat collisions", "Background observations"),
    guide = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  coord_sf(crs = proj_na, default_crs = st_crs(4326),
           xlim = st_bbox(background_region)[c("xmin", "xmax")],
           ylim = st_bbox(background_region)[c("ymin", "ymax")], expand = TRUE) +
  theme_void() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
        legend.position = "inside", legend.position.inside = c(0.01, 0.02),
        legend.justification = c(0, 0), legend.key = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(10, "pt"), legend.key.spacing.y = unit(3, "pt"),
        legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.2),
        legend.margin = margin(2, 4, 2, 3)))
# f_map is not saved on its own; it enters the combined figure below as panel (b).

## Conditional-effect panels of the top model (building height + season + night traffic) ----
ce <- conditional_effects(m_useavail_radar_both)
# Recompute the traffic effect on a log-spaced grid (other covariates held at their
# mean, as conditional_effects does) so the line is smooth on the log axis instead
# of segmented where the default linear grid is sparse.
traffic_grid <- data.frame(
  building_height = mean(radar_data$building_height),
  yday = mean(radar_data$yday),
  traffic = 10^seq(log10(min(radar_data$traffic[radar_data$traffic > 0])),
                   log10(max(radar_data$traffic)), length.out = 300)
)
traffic_epred <- posterior_epred(m_useavail_radar_both, newdata = traffic_grid)
ceData <- list(
  building_height = ce$building_height,
  yday = ce$yday,
  traffic = traffic_grid %>%
    mutate(estimate__ = colMeans(traffic_epred),
           lower__ = apply(traffic_epred, 2, quantile, 0.025),
           upper__ = apply(traffic_epred, 2, quantile, 0.975))
)
ceLabs <- c(building_height = "Building height (m)", yday = "Day of year",
            traffic = "Nightly migration traffic")

# `dat` (optional) overlays the records as rugs: collisions along the top (accent),
# background along the bottom (subtle) -- the binary-response analogue of data
# points. `y_accuracy` (optional) fixes the y-tick decimal places so several panels
# get identical axis-text width, keeping their panel regions equal when arranged.
cePanel <- function(v, logx = FALSE, dat = NULL, y_accuracy = NULL) {
  dd <- ceData[[v]]
  dd$x <- dd[[v]]
  g <- ggplot(dd, aes(x, estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__),
                fill = scales::muted(collisionColor, l = 95, c = 20)) +
    geom_line(colour = collisionColor) +
    scale_y_continuous("Relative probability of collision",
                       labels = if (is.null(y_accuracy)) waiver() else
                         scales::label_number(accuracy = y_accuracy)) +
    xlab(ceLabs[v])
  if (!is.null(dat)) {
    g <- g +
      geom_rug(data = dat[dat$used == 0, ], aes(x = .data[[v]]), inherit.aes = FALSE,
               sides = "b", colour = backgroundColor, alpha = 0.12,
               length = unit(0.03, "npc")) +
      geom_rug(data = dat[dat$used == 1, ], aes(x = .data[[v]]), inherit.aes = FALSE,
               sides = "t", colour = collisionColor, alpha = 0.5,
               length = unit(0.05, "npc"))
  }
  if (logx) g <- g + scale_x_log10()
  if (v == "yday") {
    g <- g + scale_x_continuous("Day of year",
                                breaks = monthDayYear_to_yday(monthFirsts),
                                labels = format(mdy(monthFirsts), "%b")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  g
}

## Raw use-vs-available contrasts ----
contrast_data <- radar_data %>%
  mutate(sample = factor(used, levels = c(1, 0),
                         labels = c("collision site", "background (available)")))
contrastPanel <- function(x, xlab, logx = FALSE) {
  g <- ggplot(contrast_data, aes(.data[[x]], colour = sample)) +
    geom_density() +
    scale_colour_manual(NULL, values = sampleColors) +
    scale_y_continuous("Density") +
    xlab(xlab)
  if (logx) g <- g + scale_x_log10()
  g
}
(f_contrasts <- contrastPanel("building_height", "Building height (m)") +
    contrastPanel("yday", "Day of year") +
    contrastPanel("traffic", "Nightly migration traffic (night prior)", logx = TRUE) +
    plot_layout(nrow = 1, guides = "collect") &
    theme(legend.position = "bottom"))
ggsave("figs/B_contrasts.png", f_contrasts, width = 10, height = 3.6, dpi = 600)
ggsave("figs/B_contrasts.svg", f_contrasts, width = 10, height = 3.6)

## SI: background reproduces observer effort; collisions do not ----
effort_month <- fread("data/derived/inat_background_effort.csv")[
  , .(share = sum(n_effort)), by = .(m = month(month_start))]
effort_month[, share := share / sum(share)]
month_shares <- rbindlist(list(
  cbind(effort_month, src = "effort (reference)"),
  data.table(m = month(points[used == 0, date]))[
    , .N, by = m][, .(m, share = N / sum(N), src = "background")],
  data.table(m = month(points[used == 1, date]))[
    , .N, by = m][, .(m, share = N / sum(N), src = "collisions")]
), use.names = TRUE)
(f_effort <- ggplot(month_shares, aes(m, share, colour = src)) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous("Month", breaks = 1:12, labels = month.abb) +
    scale_y_continuous("Share of records") +
    scale_colour_manual(NULL, values = c("effort (reference)" = "grey60",
                                         "background" = "#159367",
                                         "collisions" = collisionColor)))
ggsave("figs/SI_effort_validation.png", f_effort, width = 7, height = 4, dpi = 600, bg = "white")
ggsave("figs/SI_effort_validation.svg", f_effort, width = 7, height = 4)

## Species composition of collisions over the year ----
# Descriptive complement to the modeled day-of-year effect: which species collide,
# and when. Grouping and colours match the iNaturalist descriptive figure (script
# 6) -- the three commonly-identified species plus an "Other/unknown" catch-all.
# Species are assigned only from research-grade identifications; everything else
# (coarser or unverified IDs) falls into Other/unknown.
speciesFocal  <- c("Lasiurus borealis", "Lasionycteris noctivagans", "Eptesicus fuscus")
speciesLevels <- rev(c(speciesFocal, "Other/unknown"))
speciesColors <- c("Lasiurus borealis"         = "#9E2A2B",
                   "Lasionycteris noctivagans"  = "#7072A0",
                   "Eptesicus fuscus"           = "#A8541F",
                   "Other/unknown"              = "grey60")
speciesLabels <- c("Lasiurus borealis"         = "*Lasiurus borealis*",
                   "Lasionycteris noctivagans"  = "*Lasionycteris noctivagans*",
                   "Eptesicus fuscus"           = "*Eptesicus fuscus*",
                   "Other/unknown"              = "Other/unknown")
collision_species <- points[used == 1, .(yday, group = factor(
  fifelse(quality_grade == "research" & scientific_name %in% speciesFocal,
          scientific_name, "Other/unknown"),
  levels = speciesLevels))]

# Stacked day-of-year counts. `boundary = 0` anchors the bins so none straddle the
# Jan/Dec edges and get dropped. Wide-and-short by design. The legend is reversed
# (`reverse = TRUE` reorders the keys only, not the stack) and sits low and inside
# the panel, in the summer trough between the two seasonal peaks.
speciesPanel <- function() {
  ggplot(collision_species, aes(x = yday, fill = group)) +
    geom_histogram(binwidth = 14, boundary = 0, colour = "white", linewidth = 0.1) +
    scale_fill_manual(NULL, values = speciesColors, labels = speciesLabels, drop = FALSE,
                      guide = guide_legend(ncol = 1, reverse = TRUE)) +
    scale_x_continuous("Day of year",
                       breaks = monthDayYear_to_yday(monthFirsts),
                       labels = format(mdy(monthFirsts), "%b"), expand = c(0.01, 0)) +
    scale_y_continuous("Number of collision records", expand = expansion(mult = c(0, 0.08))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "inside",
          legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.title = element_blank(),
          legend.text = element_markdown(size = 8),
          legend.key.size = unit(10, "pt"),
          legend.background = element_blank(),
          legend.key.spacing.y = unit(3, "pt"))
}

## Combined main figure ----
# Top row: species composition over day of year (a) beside the sampling map (b);
# bottom row: the three modeled conditional effects (c-e). Built as two separate
# patchworks and then stacked, so the bottom three panels stay equal-width and
# aligned regardless of the map's fixed aspect ratio. Tags are set per panel because
# auto-tagging treats each nested patchwork as one unit.
p_species <- speciesPanel()                                                   + labs(tag = "(a)")
p_map     <- f_map                                                            + labs(tag = "(b)")
# y_accuracy = 0.001 gives all three the same-width y labels, so their panel
# regions come out exactly equal.
p_bh      <- cePanel("building_height", dat = radar_data, y_accuracy = 0.001) + labs(tag = "(c)")
p_yday    <- cePanel("yday", dat = radar_data, y_accuracy = 0.001)            + labs(tag = "(d)")
p_traffic <- cePanel("traffic", logx = TRUE, dat = radar_data, y_accuracy = 0.001) + labs(tag = "(e)")
datatop_top    <- p_species + p_map + plot_layout(widths = c(1, 1))
datatop_bottom <- p_bh + p_yday + p_traffic + plot_layout(nrow = 1)
(f_main_combo_datatop_data <- wrap_plots(datatop_top, datatop_bottom, ncol = 1,
                                         heights = c(1.5, 1)) &
   theme(plot.tag = element_text(size = 10), axis.text.x = element_text(size = 8)))
ggsave("figs/B_main_combo_datatop_data.png", f_main_combo_datatop_data, width = 9, height = 6, dpi = 600, bg = "white")
ggsave("figs/B_main_combo_datatop_data.svg", f_main_combo_datatop_data, width = 9, height = 6)

## SI: effect sizes and model comparison ----
# (a) posterior coefficients (relative selection strength) for the two linear
# terms; (b) LOO elpd of the three competing models (the "both" model is the ref).
coefs <- as.data.frame(fixef(m_useavail_radar_both))[c("log1pbuilding_height", "log1ptraffic"), ]
coefs$term <- c("Building height (log1p)", "Night traffic (log1p)")
p_coef <- ggplot(coefs, aes(Estimate, term)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.2) +
  geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5)) +
  scale_x_continuous("Coefficient (log-odds of use vs available)") +
  ylab(NULL)
loo_tab <- as.data.frame(loo_compare(m_useavail_season_sub, m_useavail_radaronly, m_useavail_radar_both))
loo_labs <- c(m_useavail_radar_both = "building height + season + traffic",
              m_useavail_radaronly = "building height + traffic",
              m_useavail_season_sub = "building height + season")
loo_tab$model <- factor(loo_labs[rownames(loo_tab)], levels = rev(loo_labs[rownames(loo_tab)]))
p_loo <- ggplot(loo_tab, aes(elpd_diff, model)) +
  geom_pointrange(aes(xmin = elpd_diff - se_diff, xmax = elpd_diff + se_diff)) +
  scale_x_continuous(expression(Delta * " elpd (vs best model)")) +
  ylab(NULL)
(f_effectsize_loo <- p_coef + p_loo +
   plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
   theme(plot.tag = element_text(size = 10)))
ggsave("figs/SI_effectsize_loo.png", f_effectsize_loo, width = 9, height = 3, dpi = 600, bg = "white")
ggsave("figs/SI_effectsize_loo.svg", f_effectsize_loo, width = 9, height = 3)
