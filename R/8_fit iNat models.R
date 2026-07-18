# 8_fit iNat models.R
# Purpose: Fit the iNaturalist use-availability models of suspected bat-window
#          collisions against a target-group iNaturalist background. Reads the
#          analysis-ready table from "7_prep iNat models.R", checks predictor
#          collinearity (VIF), fits the use-availability (Bernoulli) model set,
#          compares them by LOO, and builds the manuscript figures. No raster
#          extraction or downloads happen here; the Figure 5 map reads a small
#          pre-aggregated ALAN basemap baked by 7_prep.
# Input:   data/derived/useavail_points.csv               (from 7_prep iNat models.R)
#          data/derived/inat_background_effort.csv         (effort reference)
#          data/ALAN/alan_glow_log_ll.tif                  (glow basemap, from 7_prep)
# Outputs: out/models/m_useavail_*.rds                     (cached brms fits)
# Figures: figs/f5_iNaturalist results.{svg,png}   (Figure 5, both halves -- finished by hand)
#          figs/B_contrasts.{png,svg}                (raw use-vs-available contrasts)
#          figs/SI_effort_validation.{png,svg}
#          figs/SI_effectsize_loo.{png,svg}


# Setup ----
source("R/0_funs.R")
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(lubridate)
library(data.table)
library(brms)
library(splines)
library(patchwork)
library(lme4)

# Fit a use-availability Bernoulli model with the shared brms settings and attach
# LOO. The fit AND its LOO are cached by `file` under out/models/ (delete to refit);
# forwarding `file` to add_criterion() persists LOO in the .rds so it is not
# recomputed on every re-source.
brm_ua <- function(formula, data, file) {
  m <- brm(
    bf(formula), data = data, family = bernoulli(),
    cores = 4, chains = 4, seed = 42, iter = 5000, warmup = 1000,
    threads = threading(4), backend = "cmdstanr",
    file = file.path("out/models", file)
  )
  add_criterion(m, "loo", file = file.path("out/models", file))
}

# Build the design columns for a "hurdle" (two-part) version of a driver, used only as a
# robustness check on functional form (see the model comparisons below). Many points have
# the driver exactly at zero (no mapped building, or no detectable light); a hurdle splits
# the effect into two pieces: (1) a yes/no "present at all?" indicator, and (2) a smooth
# curve that shapes the response only among the above-zero values. Comparing this against a
# single smooth spline over the whole range (the form we use) tests whether treating
# "absent" as a special case buys anything. The spline basis is set to 0 where x == 0, so
# the indicator alone carries the effect there. Returned as columns to cbind onto the model
# frame; LOO only needs the in-sample fit, so precomputed basis columns are sufficient.
hurdleCols <- function(x, prefix, df = 3) {
  pos <- x > 0
  B <- matrix(0, length(x), df)
  B[pos, ] <- ns(log1p(x[pos]), df = df)
  out <- data.frame(as.numeric(pos), B)
  names(out) <- paste0(prefix, c("_pres", paste0("_sp", seq_len(df))))
  out
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
background_region <- collisions_sf |>
  st_buffer(dist = 100e3) |>
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
# Single source of truth for the two use-availability classes: each label is paired with
# its colour here, and the factor helper below derives its levels from these same names.
# So editing a label can't desync it from scale_colour_manual() -- the name mismatch that
# silently drops a class from the legend (the "factor bug"). `used == 1` is the collision
# class (first entry); `used == 0` the background class (second).
sampleColors <- c(
  "collision site"         = collisionColor,
  "background (available)" = backgroundColor
)
# Map a 0/1 `used` vector to the labelled factor. collisions_first sets legend/stack order.
sampleFactor <- function(used, collisions_first = TRUE) {
  labs <- names(sampleColors)[c(2, 1)][used + 1]   # used 0 -> background, used 1 -> collision
  factor(labs, levels = if (collisions_first) names(sampleColors) else rev(names(sampleColors)))
}
labelSample <- function(used) sampleFactor(used, collisions_first = FALSE)

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
## Full data: structural terms + season ----
# Use-availability logistic model: does collision favour taller, brighter-lit sites
# and migration-window timing, relative to where/when people observe? Building height
# and ALAN each enter as a natural spline on the log1p scale (ns(log1p(x), 3)): LOO
# supports curvature over a linear log1p term for both, and a spline is
# indistinguishable from a two-part hurdle despite the point mass at 0 (55% of points
# have no mapped building, 15% no detectable light), so the simpler spline is kept.
# Traffic (below) stays linear log1p -- no LOO support for a spline. The intercept is
# not interpretable in a use-availability design (it tracks the use:available ratio);
# effects express relative selection.
m_useavail_bh   <- brm_ua(used ~ ns(log1p(building_height), 3) + ns(yday, 5),
                          model_data, "m_useavail_bh_yday")
m_useavail_alan <- brm_ua(used ~ ns(log1p(alan), 3) + ns(yday, 5),
                          model_data, "m_useavail_alan_yday")
pp_check(m_useavail_bh, ndraws = 100)

# Functional-form check: should each structural driver enter as a straight line (on the
# log1p scale) or as a flexible curve (a spline)? We fit the straight-line version of each
# and let leave-one-out cross-validation (LOO) choose. LOO prefers the curve for both, so
# the models above use splines; these linear fits exist only for that comparison.
m_useavail_bh_lin   <- brm_ua(used ~ log1p(building_height) + ns(yday, 5),
                              model_data, "m_useavail_bh_linear")
m_useavail_alan_lin <- brm_ua(used ~ log1p(alan) + ns(yday, 5),
                              model_data, "m_useavail_alan_linear")
loo_compare(m_useavail_bh, m_useavail_bh_lin)
loo_compare(m_useavail_alan, m_useavail_alan_lin)

# Functional-form robustness for the two zero-inflated drivers: is the whole-range spline
# better than (a) a two-part hurdle (presence indicator + smooth over the positive values)
# or (b) a bare presence indicator? Built on the same "driver + season" frame as the
# spline/linear comparison, so every functional form for a driver is judged on identical
# rows. If the hurdle is indistinguishable from the spline, the simpler spline is kept; a
# presence-only indicator being clearly worse shows the graded (not just present/absent)
# response carries information.
md_bh   <- cbind(model_data, hurdleCols(model_data$building_height, "bh"))
md_alan <- cbind(model_data, hurdleCols(model_data$alan, "alan"))
m_useavail_bh_hurdle    <- brm_ua(used ~ bh_pres + bh_sp1 + bh_sp2 + bh_sp3 + ns(yday, 5),
                                  md_bh, "m_useavail_bh_hurdle")
m_useavail_bh_present   <- brm_ua(used ~ bh_pres + ns(yday, 5),
                                  md_bh, "m_useavail_bh_presence")
m_useavail_alan_hurdle  <- brm_ua(used ~ alan_pres + alan_sp1 + alan_sp2 + alan_sp3 + ns(yday, 5),
                                  md_alan, "m_useavail_alan_hurdle")
m_useavail_alan_present <- brm_ua(used ~ alan_pres + ns(yday, 5),
                                  md_alan, "m_useavail_alan_presence")
loo_compare(m_useavail_bh, m_useavail_bh_hurdle)     # hurdle vs spline (building height)
loo_compare(m_useavail_bh, m_useavail_bh_present)    # presence-only vs spline (building height)
loo_compare(m_useavail_alan, m_useavail_alan_hurdle) # hurdle vs spline (radiance)
loo_compare(m_useavail_alan, m_useavail_alan_present)# presence-only vs spline (radiance)

# Full-data co-headline model: collisions favour sites that are both tall and lit.
m_useavail_struct <- brm_ua(used ~ ns(log1p(building_height), 3) + ns(log1p(alan), 3) + ns(yday, 5),
                            model_data, "m_useavail_struct")

## Radar-matched subset: add night-migration traffic ----
# Competing models on the radar-matched subset (same rows, for a fair LOO). Traffic
# is linear log1p (right-skewed; no LOO support for a spline). Season only, traffic
# only, both, then the two light models: m_useavail_radar_alan swaps height for ALAN,
# and m_useavail_radar_all (top model) carries height, light, season and night traffic
# together -- letting us see how much of the built-environment signal is light rather
# than raw height, and whether night traffic adds beyond both.
m_useavail_season_sub <- brm_ua(used ~ ns(log1p(building_height), 3) + ns(yday, 5),
                                radar_data, "m_useavail_season_sub")
m_useavail_radaronly  <- brm_ua(used ~ ns(log1p(building_height), 3) + log1p(traffic),
                                radar_data, "m_useavail_radaronly")
m_useavail_radar_both <- brm_ua(used ~ ns(log1p(building_height), 3) + ns(yday, 5) + log1p(traffic),
                                radar_data, "m_useavail_radar_both")
m_useavail_radar_alan <- brm_ua(used ~ ns(log1p(alan), 3) + ns(yday, 5) + log1p(traffic),
                                radar_data, "m_useavail_radar_alan")
m_useavail_radar_all  <- brm_ua(used ~ ns(log1p(building_height), 3) + ns(log1p(alan), 3) + ns(yday, 5) + log1p(traffic),
                                radar_data, "m_useavail_radar_all")

# Functional-form check for traffic: does a spline beat the linear log1p term inside the
# full radar model (all other terms identical)? If not, traffic is kept linear.
m_useavail_traffic_spline <- brm_ua(
  used ~ ns(log1p(building_height), 3) + ns(log1p(alan), 3) + ns(yday, 5) + ns(log1p(traffic), 3),
  radar_data, "m_useavail_traffic_spline")
loo_compare(m_useavail_radar_all, m_useavail_traffic_spline)  # traffic spline vs linear

loo_compare(m_useavail_season_sub, m_useavail_radaronly, m_useavail_radar_both,
            m_useavail_radar_alan, m_useavail_radar_all)
conditional_effects(m_useavail_radar_all)

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
  ) |> 
    filter(dist_km < 200)
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

## Figure 5: sampling design map + conditional effects -------
### Sampling design map ----
# North and South America with state/province boundaries and the Great
# Lakes, the 100-km "available" buffer, background points (grey), and collisions
# (blue) on top. 
na_land <- ne_countries(scale = "medium", returnclass = "sf") |>
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

# Draw on the ALAN surface.
base_dk  <- "#03010a"
col_bg   <- "#E8EEF7"   # background (available) points, insets only
col_coll <- "#5AB0FF"   # collisions: luminous blue, reads against the glow
col_reg  <- "#38D6B0"   # 100-km available region outline
# Shared key for the whole figure: the records are only drawn in the insets, but the
# legend lives on the map panel, which has the room for it.
mapKeys <- c("Bat collision"           = col_coll,
             "Background (available)"  = col_bg,
             "100-km available region" = col_reg)
lim_hi   <- log1p(50)   # ~p99.9 of radiance -> palette top; cores saturate
# Figure 5 basemap (log-radiance glow), written by 7_prep iNat models.R. Not
# distributed with the repo (data/ALAN/ is untracked); re-run script 7 to regenerate.
stopifnot("data/ALAN/alan_glow_log_ll.tif not found; produced by script 7" = file.exists("data/ALAN/alan_glow_log_ll.tif"))
glow_ll  <- rast("data/ALAN/alan_glow_log_ll.tif"); names(glow_ll) <- "radiance"
glow_pr  <- project(glow_ll, proj_na, method = "bilinear")

fillGlow <- function(guide = "none") scale_fill_viridis_c(
  option = "magma", limits = c(0, lim_hi), oob = scales::squish, na.value = base_dk,
  name = "Nighttime\nradiance", breaks = log1p(c(0, 1, 10, 50)),
  labels = c("0", "1", "10", "50+"), guide = guide)
themeDark <- function() theme_void(base_size = 8) + theme(
  plot.background  = element_rect(fill = base_dk, colour = NA),
  panel.background = element_rect(fill = base_dk, colour = NA),
  legend.title = element_text(colour = "grey85", size = 7),
  legend.text  = element_text(colour = "grey75", size = 7),  # match lower-panel axis text (f5_axis)
  plot.margin  = margin(1, 1, 1, 1))

# Per-collision buffers; the union is the "available" region. Only the circles are
# drawn on the main map -- each is centred on a collision, so they mark the records
# and the available footprint at once, without a layer of overlapping points.
buffers <- st_buffer(collisions_sf, 100000)
# View = extent of all buffers plus a margin, so every circle sits fully in frame. The
# western margin is dropped: the Pacific-coast circles already reach into open ocean,
# where no background exists, so extra room there is dead space.
be   <- st_bbox(st_transform(background_region, proj_na))
mrg  <- 75000
mrg_w <- 0

# Inset windows: a dense cluster (mid-Atlantic) and an isolated collision (Denver),
# zoomed far enough that individual background points read. Denver shows the effort
# weighting directly -- available points track the Front Range glow.
winA <- list(xlim = c(-79.5, -73.5), ylim = c(38.2, 41.8), col = "#8BD3FF", tag = "(b)")
winB <- list(xlim = c(-107, -103),   ylim = c(37.8, 41.6), col = "#F5C563", tag = "(c)")
boxPoly <- function(w) st_as_sf(st_sfc(st_polygon(list(rbind(
  c(w$xlim[1], w$ylim[1]), c(w$xlim[2], w$ylim[1]), c(w$xlim[2], w$ylim[2]),
  c(w$xlim[1], w$ylim[2]), c(w$xlim[1], w$ylim[1])))), crs = 4326))

(f_map <- ggplot() +
  geom_spatraster(data = glow_pr, maxcell = 6e6) +
  geom_sf(data = na_land, fill = NA, colour = "grey30", linewidth = 0.12) +
  geom_sf(data = background_region, aes(colour = "100-km available region"),
          fill = NA, linewidth = 0.2, show.legend = "line") +

  geom_sf(data = collisions_sf[1, ], aes(colour = "Bat collision"),
          alpha = 0, size = 0.1, show.legend = "point") +
  geom_sf(data = background_sf[1, ], aes(colour = "Background (available)"),
          alpha = 0, size = 0.1, show.legend = "point") +
  geom_sf(data = boxPoly(winA), fill = NA, colour = winA$col, linewidth = 0.4) +
  geom_sf(data = boxPoly(winB), fill = NA, colour = winB$col, linewidth = 0.4) +

    
  scale_colour_manual(NULL, values = mapKeys, breaks = names(mapKeys),
                      guide = guide_legend(order = 2, override.aes = list(
                        alpha    = 1,
                        size     = c(1.6, 1.6, 0),
                        linetype = c(0, 0, 1)))) +
  fillGlow(guide = guide_colourbar(order = 1, barheight = unit(38, "pt"),
                                   barwidth = unit(5, "pt"),
                                   ticks.colour = NA, frame.colour = NA)) +
  coord_sf(crs = proj_na, xlim = c(be[["xmin"]] - mrg_w, be[["xmax"]] + mrg),
           ylim = c(be[["ymin"]] - mrg, be[["ymax"]] + mrg), expand = FALSE) +
  themeDark() +
  theme(legend.position = "inside", legend.position.inside = c(0.015, 0.30),
        legend.justification = c(0, 0.5), legend.box.just = "left",
        legend.key = element_blank(), legend.spacing.y = unit(3, "pt")))

# Inset panels: buffers merge into the union outline (as on the main map), with the
# collisions and the background points that were drawn inside them.
insetPanel <- function(w) {
  ext_w <- ext(w$xlim[1] - 0.4, w$xlim[2] + 0.4, w$ylim[1] - 0.4, w$ylim[2] + 0.4)
  pb <- st_bbox(st_transform(boxPoly(w), proj_na))
  ggplot() +
    geom_spatraster(data = project(crop(glow_ll, ext_w), proj_na), maxcell = 2e6) +
    geom_sf(data = background_region, fill = NA, colour = col_reg, linewidth = 0.25) +
    geom_sf(data = background_sf[lengths(st_intersects(background_sf, boxPoly(w))) > 0, ],
            shape = 21, fill = col_bg, colour = "grey20", stroke = 0.1, size = 0.7, alpha = 0.75) +
    geom_sf(data = collisions_sf[lengths(st_intersects(collisions_sf, boxPoly(w))) > 0, ],
            shape = 21, fill = col_coll, colour = "white", stroke = 0.3, size = 1.5) +
    fillGlow() +
    coord_sf(crs = proj_na, xlim = pb[c("xmin", "xmax")], ylim = pb[c("ymin", "ymax")],
             expand = FALSE) +
    themeDark() +
    theme(panel.border = element_rect(fill = NA, colour = w$col, linewidth = 0.9))
}
f_insetA <- insetPanel(winA); f_insetB <- insetPanel(winB)
# These are not saved alone; they enter the combined figure below as panels (a)-(c).

### Conditional-effect panels of the top model (building height + ALAN + season + night traffic) ----
ce <- conditional_effects(m_useavail_radar_all)
# Recompute the traffic effect on a log-spaced grid (other covariates held at their
# mean, as conditional_effects does) so the line is smooth on the log axis instead
# of segmented where the default linear grid is sparse.
traffic_grid <- data.frame(
  building_height = mean(radar_data$building_height),
  alan = mean(radar_data$alan),
  yday = mean(radar_data$yday),
  traffic = 10^seq(log10(min(radar_data$traffic[radar_data$traffic > 0])),
                   log10(max(radar_data$traffic)), length.out = 300)
)
traffic_epred <- posterior_epred(m_useavail_radar_all, newdata = traffic_grid)
ceData <- list(
  building_height = ce$building_height,
  alan = ce$alan,
  yday = ce$yday,
  traffic = traffic_grid |>
    mutate(estimate__ = colMeans(traffic_epred),
           lower__ = apply(traffic_epred, 2, quantile, 0.025),
           upper__ = apply(traffic_epred, 2, quantile, 0.975))
)
ceLabs <- c(building_height = "Building height (m)", alan = "Nighttime radiance (ALAN)",
            yday = "Day of year", traffic = "Nightly migration traffic")

# `dat` (optional) overlays the records as rugs: collisions along the top (accent),
# background along the bottom (subtle) -- the binary-response analogue of data
# points. `y_accuracy` (optional) fixes the y-tick decimal places so several panels
# get identical axis-text width, keeping their panel regions equal when arranged.
cePanel <- function(v, logx = FALSE, log1px = FALSE, dat = NULL, y_accuracy = NULL) {
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
               length = unit(0.05, "npc")) +
      geom_rug(data = dat[dat$used == 1, ], aes(x = .data[[v]]), inherit.aes = FALSE,
               sides = "t", colour = collisionColor, alpha = 0.5,
               length = unit(0.05, "npc"))
  }
  if (logx) g <- g + scale_x_log10()
  # ALAN is right-skewed with a point mass at 0 (log1p in the model); a pseudo-log
  # axis spreads the curve out while still showing the zero values.
  if (log1px) g <- g + scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                                          breaks = c(0, 3, 10, 30, 100, 300))
  if (v == "yday") {
    g <- g + scale_x_continuous("Day of year",
                                breaks = monthDayYear_to_yday(monthFirsts),
                                labels = format(mdy(monthFirsts), "%b")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  g
}

### Combined main figure (Figure 5) ----

p_map     <- f_map     + labs(tag = "(a)") + theme(plot.tag = element_text(colour = "grey85"))
p_insetA  <- f_insetA  + labs(tag = winA$tag) + theme(plot.tag = element_text(colour = winA$col))
p_insetB  <- f_insetB  + labs(tag = winB$tag) + theme(plot.tag = element_text(colour = winB$col))

p_bh0      <- cePanel("building_height", dat = radar_data, y_accuracy = 0.001) #+ labs(tag = "(d)")
p_alan0    <- cePanel("alan", log1px = TRUE, dat = radar_data, y_accuracy = 0.001) #+ labs(tag = "(e)")
p_yday0    <- cePanel("yday", dat = radar_data, y_accuracy = 0.001)            #+ labs(tag = "(f)")
p_traffic0 <- cePanel("traffic", logx = TRUE, dat = radar_data, y_accuracy = 1e-4) #+ labs(tag = "(g)")

f5_tag  <- theme(plot.tag = element_text(size = 9))
f5_axis <- theme(axis.text.x = element_text(size = 7))
p_map     <- p_map     + f5_tag
p_insetA  <- p_insetA  + f5_tag
p_insetB  <- p_insetB  + f5_tag
p_bh      <- p_bh0      + f5_tag + f5_axis
p_alan    <- p_alan0    + f5_tag + f5_axis
p_yday    <- p_yday0    + f5_tag + f5_axis
p_traffic <- p_traffic0 + f5_tag + f5_axis

# Whole figure is 6.5 in wide. Top-half panel widths are pinned in absolute inches so the
# map (a) exports at ~4.107 in and each inset column (b, c) at ~2.16 in, as measured in
# Inkscape; the ~0.23 in remainder falls to the inter-panel gap and outer margins.
f5_top    <- p_map + (p_insetA / p_insetB) +
  plot_layout(widths = c(4.107/2.16, 1))
# Bottom 2x2 shares one y-axis label: drop the four repeated per-panel y-titles and
# add a single centred, rotated label to their left (11 pt, matching the panel axis
# titles under theme_classic). Widths c(1, 26) give the strip ~0.24 in at 6.5 in wide.
# plot.margin (default 5.5 pt) sets the gap between patchwork panels; trim it to pull
# the 2x2 effect plots closer together.
f5_effects <- p_bh + p_alan + p_yday + p_traffic + plot_layout(nrow = 2) &
  theme(axis.title.y = element_blank(), plot.margin = margin(2, 2, 2, 2))
f5_ylab <- wrap_elements(full = grid::textGrob(
  "Relative probability of collision", rot = 90, gp = grid::gpar(fontsize = 11)))
f5_bottom <- f5_ylab + f5_effects + plot_layout(widths = c(1, 26))
ggsave("figs/f5_bottom.png", f5_bottom, width = 6.5, height = (8.06-4.226), bg = "white", dpi = 100)
ggsave("figs/f5_bottom.svg", f5_bottom, width = 6.5, height = (8.06-4.226), bg = "white")

# Stack the two halves into a single SVG. Heights track each half's native aspect.
f5_combined <- f5_top / f5_bottom + plot_layout(heights = c(1, 0.9))
ggsave("figs/f5_iNaturalist results.png", f5_combined, width = 6.5, height = 9, dpi = 200, bg = "white")
ggsave("figs/f5_iNaturalist results.svg", f5_combined, width = 6.5, height = 9, bg = "white")



## Raw use-vs-available contrasts ----
contrast_data <- radar_data |>
  mutate(sample = sampleFactor(used, collisions_first = TRUE))
contrastPanel <- function(x, xlab, logx = FALSE) {
  g <- ggplot(contrast_data, aes(.data[[x]], colour = sample)) +
    geom_density() +
    scale_colour_manual(NULL, values = sampleColors) +
    scale_y_continuous("Density") +
    xlab(xlab)
  if (logx) g <- g + scale_x_log10()
  g
}
f_contrasts <- contrastPanel("building_height", "Building height (m)") +
    contrastPanel("alan", "Nighttime radiance (ALAN)") +
    contrastPanel("yday", "Day of year") +
    contrastPanel("traffic", "Nightly migration traffic (night prior)", logx = TRUE) +
    plot_layout(nrow = 2, guides = "collect") &
    theme(legend.position = "bottom")
ggsave("figs/B_contrasts.png", f_contrasts, width = 8, height = 6, dpi = 600)

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
f_effort <- ggplot(month_shares, aes(m, share, colour = src)) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous("Month", breaks = 1:12, labels = month.abb) +
    scale_y_continuous("Share of records") +
    scale_colour_manual(NULL, values = c("effort (reference)" = "grey60",
                                         "background" = "#159367",
                                         "collisions" = collisionColor))
ggsave("figs/SI_effort_validation.png", f_effort, width = 7, height = 4, dpi = 600, bg = "white")

## SI: effect sizes and model comparison ----
# (a) comparable effect sizes: the odds ratio for an interquantile (10th -> 90th
# percentile) change in each driver from the top model, holding the others at their
# mean. Building height and ALAN are splined (no single coefficient), so each is
# summarised on a common scale as relative selection across its observed range;
# traffic is linear but shown the same way for comparability. (b) LOO elpd of the five
# radar-matched models (top model is the reference).
iqOR <- function(model, var, data, lo = 0.10, hi = 0.90) {
  base <- data.frame(building_height = mean(data$building_height), alan = mean(data$alan),
                     yday = mean(data$yday), traffic = mean(data$traffic))
  nd <- base[c(1, 1), ]; nd[[var]] <- quantile(data[[var]], c(lo, hi))
  d  <- posterior_linpred(model, newdata = nd)
  d  <- d[, 2] - d[, 1]                          # log odds ratio, 90th vs 10th pct
  data.frame(Estimate = mean(d), Q2.5 = quantile(d, .025), Q97.5 = quantile(d, .975))
}
effs <- do.call(rbind, lapply(c("building_height", "alan", "traffic"),
                              function(v) iqOR(m_useavail_radar_all, v, radar_data)))
effs$term <- factor(c("Building height", "ALAN (radiance)", "Night traffic"),
                    levels = rev(c("Building height", "ALAN (radiance)", "Night traffic")))
p_coef <- ggplot(effs, aes(exp(Estimate), term)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.2) +
  geom_pointrange(aes(xmin = exp(Q2.5), xmax = exp(Q97.5))) +
  scale_x_log10("Odds ratio (10th → 90th percentile of driver)") +
  ylab(NULL)
loo_tab <- as.data.frame(loo_compare(m_useavail_season_sub, m_useavail_radaronly,
                                     m_useavail_radar_both, m_useavail_radar_alan, m_useavail_radar_all))
loo_labs <- c(m_useavail_radar_all  = "height + light + season + traffic",
              m_useavail_radar_alan = "light + season + traffic",
              m_useavail_radar_both = "height + season + traffic",
              m_useavail_radaronly  = "height + traffic",
              m_useavail_season_sub = "height + season")
loo_tab$model <- factor(loo_labs[rownames(loo_tab)], levels = rev(loo_labs[rownames(loo_tab)]))
p_loo <- ggplot(loo_tab, aes(elpd_diff, model)) +
  geom_pointrange(aes(xmin = elpd_diff - se_diff, xmax = elpd_diff + se_diff)) +
  scale_x_continuous(expression(Delta * " elpd (vs best model)")) +
  ylab(NULL)
f_effectsize_loo <- p_coef + p_loo +
   plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
   theme(plot.tag = element_text(size = 10))
ggsave("figs/SI_effectsize_loo.png", f_effectsize_loo, width = 10, height = 3.2, dpi = 600, bg = "white")


## SI: is there enough night-to-night variation in traffic to matter? (Supplemental Figure X) ----
# This backs up a specific worry about the traffic result. The traffic effect came out
# modest and uncertain, and one innocent explanation would be that nightly migration traffic
# barely changes from night to night -- if a predictor hardly varies, no model could detect
# an effect from it. We show that is not the case: night-to-night change is in fact the
# LARGEST source of variation in migration-window traffic. We take every night's traffic at
# the weather-radar stations near our collisions and split its total variation into three
# parts:
#   - spatial:       differences between stations (some sit on busier flyways)
#   - seasonal drift: the gentle rise and fall across a migration window
#   - night-level:   what is left -- the night-to-night, weather-driven swings
# and we also report how much busier a typical "busy" night is than a "quiet" one. All
# numbers are computed here (not hard-coded), so the Results text stays in sync with the
# data. Cited as "Supplemental Figure X" pending final SI numbering.

# Collision-matched stations: those actually assigned to a collision within radar range
# (the set whose nights inform the traffic effect).
matched_stations <- unique(points[used == 1 & dist_km < 200 & !is.na(station), station])

# Full nightly-traffic series at those stations, 2019-2025 (period == night). Read
# straight from the Dark Ecology daily files (7_prep only kept the matched nights).
traffic_nightly <- rbindlist(lapply(2019:2025, function(y)
  fread(sprintf("data/darkecology/daily/%d-daily.csv", y))[
    period == "night", .(station, date = as.IDate(date), traffic)]))
traffic_nightly <- traffic_nightly[station %in% matched_stations & !is.na(traffic) & traffic > 0]
traffic_nightly[, `:=`(yday = yday(date), year = year(date), week = isoweek(date))]

# Fixed biological migration-peak windows (spring 15 Apr - 15 Jun; autumn 20 Aug - 15
# Oct), defined by migration phenology rather than by where our own records happen to
# fall, so "within-window" means within the peak, where the seasonal trend is gentle and
# night-to-night weather-driven pulses dominate.
spring_win <- yday(as.Date(c("2021-04-15", "2021-06-15")))
autumn_win <- yday(as.Date(c("2021-08-20", "2021-10-15")))
traffic_win <- traffic_nightly[
  (yday %between% spring_win) | (yday %between% autumn_win)]
traffic_win[, season := fifelse(yday %between% spring_win, "spring", "autumn")]
traffic_win[, ltraf := log1p(traffic)]

# Variance components on log1p(traffic): (1|station) spatial, (1|week) seasonal drift
# within the window, residual night-level. A separate fit adding (1|year) leaves the
# year share negligible (<1%), so it is not carried as a fourth term.
vp_mod   <- lmer(ltraf ~ 1 + (1 | station) + (1 | week), data = traffic_win, REML = TRUE)
vp       <- as.data.frame(lme4::VarCorr(vp_mod))
vp_share <- setNames(100 * vp$vcov / sum(vp$vcov), vp$grp)

# Busy:quiet ratio -- within each station x season x year with enough nights, the ratio
# of a busy night (90th pct) to a quiet one (10th pct); median across those groups.
bq <- traffic_win[, .(p10 = quantile(traffic, .10), p90 = quantile(traffic, .90), .N),
                  by = .(station, season, year)][N >= 10][, ratio := p90 / p10]
bq_median <- median(bq$ratio)

cat(sprintf(
  "Traffic variance partition (%d collision-matched stations, %d migration-window nights):\n  night-level %.0f%%, spatial %.0f%%, seasonal drift %.0f%%; busy:quiet median %.0fx (%d station-season-years)\n",
  length(matched_stations), nrow(traffic_win),
  vp_share[["Residual"]], vp_share[["station"]], vp_share[["week"]], bq_median, nrow(bq)))

# Panel (a): the three variance shares.
vp_df <- data.frame(
  component = c("Night-level\n(residual)", "Spatial\n(between-station)", "Seasonal drift\n(within-window)"),
  share     = c(vp_share[["Residual"]], vp_share[["station"]], vp_share[["week"]]))
vp_df$component <- factor(vp_df$component, levels = vp_df$component)
p_partition <- ggplot(vp_df, aes(share, component)) +
  geom_col(fill = collisionColor, width = 0.62) +
  geom_text(aes(label = sprintf("%.0f%%", share)), hjust = -0.2, size = 3.2, colour = "grey25") +
  scale_x_continuous("Share of migration-window traffic variance",
                     labels = scales::label_percent(scale = 1),
                     limits = c(0, max(vp_df$share) * 1.18), expand = c(0, 0)) +
  ylab(NULL) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, lineheight = 0.9))

# Panel (b): one representative station-season-year -- KMKX (Milwaukee) autumn 2020,
# chosen because its busy:quiet ratio sits right at the study-wide median (52x), so the
# illustration is typical rather than extreme. Each point is a night's traffic (log
# axis); the loess line is the gentle within-window seasonal drift, and the shaded band
# is this station-season-year's 10th-90th percentile (the same unit the median
# busy:quiet ratio is computed over) -- the large busy:quiet spread the night-level term
# captures, dwarfing the seasonal trend.
eg <- traffic_win[station == "KMKX" & season == "autumn" & year == 2020]
eg_band <- quantile(eg$traffic, c(.10, .90))
p_series <- ggplot(eg, aes(yday, traffic)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = eg_band[1], ymax = eg_band[2],
           fill = collisionColor, alpha = 0.10) +
  geom_point(colour = backgroundColor, alpha = 0.55, size = 1.1) +
  geom_smooth(method = "loess", formula = y ~ x, colour = collisionColor,
              fill = collisionColor, alpha = 0.2, linewidth = 0.9) +
  annotate("text", x = min(eg$yday), y = max(eg$traffic), vjust = 1, hjust = 0, size = 3,
           colour = "grey25",
           label = sprintf("busy:quiet (10th-90th pct) %.0f×\n= median across station-seasons",
                           eg_band[2] / eg_band[1])) +
  scale_y_log10("Nightly migration traffic (KMKX, autumn 2020)") +
  scale_x_continuous("Day of year",
                     breaks = monthDayYear_to_yday(monthFirsts),
                     labels = format(mdy(monthFirsts), "%b")) +
  coord_cartesian(xlim = autumn_win) +
  theme_classic()

(f_traffic_variance <- p_partition + p_series +
   plot_layout(widths = c(1, 1.15)) +
   plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
   theme(plot.tag = element_text(size = 10)))
ggsave("figs/SI_traffic_variance.png", f_traffic_variance, width = 9, height = 3.6, dpi = 600, bg = "white")
ggsave("figs/SI_traffic_variance.svg", f_traffic_variance, width = 9, height = 3.6)
