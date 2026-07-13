# 7_analyze iNat records.R
# Purpose: Compare environmental characteristics (e.g., building height)
#          at suspected bat-window collisions against a target-group
#          background of iNaturalist records.
# Inputs:  data/iNat_observations_tidy_manualChecks.csv   (retained collisions)
#          data/building height data/GBH2020_150m_GEDI.tif (Ma et al. 2023)
# Outputs: tmp/inat_background.csv          (background records, cached)
#          tmp/inat_background_effort.csv   (monthly effort / phenology reference)


# Setup ----
source("R/0_funs.R")
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(lubridate)
library(data.table)
library(rinat)
library(patchwork)
library(ggtext)  

# Extract building height at a set of points, treating unmapped (NA) cells as
# building-free ground (0 m). `points` must already be in the raster's CRS.
# `extract` is masked by tidyr, hence the terra:: qualifier.
extractBuildingHeight <- function(points, height_raster) {
  terra::extract(height_raster, points, ID = FALSE) %>%
    transmute(building_height = replace_na(building_height, 0))
}

# Bounding box of an sf/sfc geometry as the c(swlat, swlng, nelat, nelng) vector
# searchBuilder() expects (note the lat, lng, lat, lng order).
boundsVec <- function(geom) {
  bb <- st_bbox(geom)
  c(bb[["ymin"]], bb[["xmin"]], bb[["ymax"]], bb[["xmax"]])
}

# Monthly iNaturalist observation counts within a bbox (v1 API histogram; any
# quality grade, geotagged). Used as space-and-time effort weights for the
# background draw, and kept as an independent effort-phenology reference.
inatHistogram <- function(bounds, d1, d2, taxon_id = 1) {
  response <- httr::GET(
    "https://api.inaturalist.org/v1/observations/histogram",
    query = list(
      taxon_id = taxon_id, verifiable = "any", geo = "true",
      date_field = "observed", interval = "month", d1 = d1, d2 = d2,
      swlat = bounds[1], swlng = bounds[2], nelat = bounds[3], nelng = bounds[4]
    )
  )
  counts <- httr::content(response, as = "parsed", type = "application/json")$results$month
  Sys.sleep(1)
  if (length(counts) == 0) {
    return(tibble(month_start = as.Date(character()), n_effort = integer()))
  }
  tibble(month_start = as.Date(names(counts)), n_effort = as.integer(unlist(counts)))
}

# One randomly chosen page of iNaturalist observations within a bbox and month
# (any grade, geotagged). Randomising the page spreads the draw across the records
# available in that block. Returns records as a data.frame, a 0-row data.frame for a
# confirmed-empty block, or NULL on API/network failure (so the caller can retry).
fetchRandomPage <- function(bounds, d1, d2, taxon_id = 1, perpage = 200) {
  search_out <- searchBuilder(taxon_id = taxon_id, bounds = bounds, geo = TRUE, d1 = d1, d2 = d2)
  total <- howManyResults(search_out)
  if (is.null(total)) return(NULL)          # API/network failure -> caller retries
  if (total == 0) return(data.frame())      # confirmed empty block
  page <- sample.int(min(50, ceiling(total / perpage)), 1)  # API caps page*per_page at 10,000
  page_query <- paste0(search_out$search, "&per_page=", perpage, "&page=", page)
  response <- httr::GET(search_out$base_url, path = "observations.csv", query = page_query)
  text_out <- rinat:::inat_handle(response)
  if (length(text_out) == 0) return(NULL)
  read.csv(textConnection(text_out), stringsAsFactors = FALSE)
}

# Run fetch(id) for each id, appending its returned rows to out_file and logging the
# completed id to done_file, so a long download can be stopped and restarted without
# repeating finished work. On restart, ids already in done_file are skipped. fetch()
# returns a data.frame (0+ rows) on success, or NULL for a transient failure that is
# retried (up to `retries`) and, if still failing, left for a later run.
resumableDownload <- function(ids, fetch, out_file, done_file, retries = 3, sleep = 5) {
  done <- if (file.exists(done_file)) fread(done_file)$id else vector(mode = class(ids))
  todo <- setdiff(ids, done)
  cat(length(todo), "of", length(ids), "items to download\n")
  for (k in seq_along(todo)) {
    id   <- todo[k]
    rows <- NULL
    for (attempt in seq_len(retries)) {
      rows <- tryCatch(fetch(id), error = function(e) NULL)
      if (!is.null(rows)) break
      Sys.sleep(sleep)
    }
    if (is.null(rows)) next                          # still failing -> resume on a later run
    if (nrow(rows) > 0) fwrite(rows, out_file, append = file.exists(out_file))
    fwrite(data.table(id = id), done_file, append = file.exists(done_file))
    if (k %% 5 == 0) cat("\r", k, "/", length(todo), "blocks done")
  }
  cat("\n")
}

## Load collision records ----
# Manually-checked iNaturalist records (retained: CJ.manual.check == "y"),
# matched to continent via a spatial join with country polygons.
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(myproj) %>%
  select(admin, adm0_a3, continent)

collisions <- read.csv("data/iNat_observations_tidy_manualChecks.csv") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = proj.wgs84, remove = FALSE) %>%
  st_transform(myproj) %>%
  filter(CJ.manual.check == "y") %>%
  st_join(world)

## Load building height data ----
# GBH2020: global 150-m urban building height from spaceborne lidar (Ma et al.
# 2023, Scientific Data). Rename the single layer so downstream code does not
# depend on the file name.
# This is ~1.5GB and can take some time to load.
if(!exists("bh_raster")) {
  bh_raster <- rast("data/building height data/GBH2020_150m_GEDI.tif")
  names(bh_raster) <- "building_height"
}

# Precise North American collision points, restricted to 2019-2025.
# Drop six pre-2019 records are dropped from both use and background to match
# availability of iNaturalist data.
# Project to the raster's native (metric) CRS.
# Also drop iNaturalist observations w/ obscured coordinates (randomized w/in 0.2 degrees).
collisions_noam <- collisions %>%
  filter(continent == "North America", geoprivacy != "obscured") %>%
  mutate(date = as_date(as.POSIXct(datetime)),
         year = year(date),
         yday = yday(date)) %>%
  filter(year %in% 2019:2025) %>%
  st_transform(crs(bh_raster))

# Generate a 100-km buffer around collisions to define the "available" region.
# Restricting "available" to within 100 km of a collision matches the strong
# bias towards east North America records. Removes continental gradients etc.
# Buffer geodesically on lon/lat (s2): this gives a true 100-km radius. Buffering
# in the raster's Mollweide CRS instead skews the radius ~15% wider E-W away from
# the central meridian, so the "available" region would be anisotropic.
background_region <- collisions_noam %>%
  st_transform(proj.wgs84) %>%
  st_buffer(dist = 100e3) %>%
  st_union() 
background_region <- st_sf(geometry = background_region)

# Target-group background ----
# The "available" set is ~10,000 general iNaturalist animal observations (any
# quality grade, matching how collisions were retained), representing where and
# WHEN people observe. Comparing collisions against this effort surface — rather
# than against raw space and time — supports both the building-height and the
# day-of-year analyses. The whole download is cached to tmp/.
background_file   <- "tmp/inat_background.csv"
n_background      <- 10000
records_per_block <- 20    # records kept per sampled (cell x month) block
oversample        <- 3   # draw extra to survive buffer / coordinate filtering

if(!file.exists(background_file)) {

  set.seed(42)

  ## Grid the region into 2-deg cells ----
  # Querying a small cell at a time keeps per-query record counts under the API's
  # 10,000 cap and lets effort be weighted across space.
  cells <- st_make_grid(background_region, cellsize = 2) %>%
    st_sf() %>%
    st_make_valid() %>% 
    mutate(cell_id = row_number()) %>%
    st_filter(st_make_valid(background_region)) %>% 
    st_make_valid()

  ## Per-cell monthly effort ----
  # One histogram call per cell yields the space-time effort weights, saved as the
  # effort-phenology reference for validating the sample's day-of-year spread.
  # Resumable per-cell histogram download (see resumableDownload).
  resumableDownload(
    ids = cells$cell_id,
    fetch = function(id) {
      h <- inatHistogram(boundsVec(cells[cells$cell_id == id, ]), "2019-01-01", "2025-12-31")
      if (nrow(h) == 0) return(h)                    # confirmed empty cell -> logged as done
      mutate(h, cell_id = id)
    },
    out_file  = "tmp/inat_background_effort.csv",
    done_file = "tmp/inat_background_effort_done.csv"
  )
  
  # Arrange to a stable order so the block sampling below is reproducible across
  # restarts (row order feeds sample.int).
  effort <- fread("tmp/inat_background_effort.csv") %>%
    mutate(year = year(month_start)) %>%
    filter(year %in% 2019:2025, n_effort > 0) %>%
    arrange(cell_id, month_start)

  ## Year quotas matched to the collisions ----
  # Background per-year totals mirror collision per-year counts, removing
  # iNaturalist's exponential growth as a confound while leaving the within-year
  # (phenological) effort structure intact.
  year_quota <- collisions_noam %>%
    st_drop_geometry() %>%
    count(year, name = "n_collision") %>%
    mutate(quota = round(n_background * oversample * n_collision / sum(n_collision)))

  ## Sample (cell x month) blocks in proportion to effort, within each year ----
  blocks <- lapply(year_quota$year, function(y) {
    pool <- filter(effort, year == y)
    n_blocks <- ceiling(year_quota$quota[year_quota$year == y] / records_per_block)
    pool[sample.int(nrow(pool), n_blocks, replace = TRUE, prob = pool$n_effort),
         c("cell_id", "month_start")]
  }) %>%
    bind_rows()

  ## Fetch a random page per block and keep a few records each ----
  # Resumable: each block's kept records are appended to tmp/background_raw.csv as it
  # completes, and its block index is logged, so this (now much longer) download can
  # be stopped and restarted without repeating finished blocks. `blocks` is
  # deterministic given set.seed(42) + the arranged effort table, so block indices
  # are stable across restarts. Delete tmp/background_raw.csv and
  # tmp/background_blocks_done.csv to start over. `download_date` = access date.
  blocks$block_id <- seq_len(nrow(blocks))
  keep_cols <- c("id", "url", "user_login", "license", "observed_on", "datetime",
                 "latitude", "longitude", "positional_accuracy", "coordinates_obscured",
                 "quality_grade", "iconic_taxon_name")
  resumableDownload(
    ids = blocks$block_id,
    fetch = function(i) {
      cell    <- cells[cells$cell_id == blocks$cell_id[i], ]
      m_start <- blocks$month_start[i]
      m_end   <- ceiling_date(m_start, "month") - days(1)
      page <- fetchRandomPage(boundsVec(cell), as.character(m_start), as.character(m_end))
      if (is.null(page)) return(NULL)                # transient failure -> retried later
      if (nrow(page) == 0) return(page)              # confirmed empty -> logged as done
      keep <- page[sample.int(nrow(page), min(records_per_block, nrow(page))), , drop = FALSE]
      for (col in setdiff(keep_cols, names(keep))) keep[[col]] <- NA  # stable schema for append
      keep <- keep[, keep_cols]
      keep$block_id <- i
      keep
    },
    out_file  = "tmp/background_raw.csv",
    done_file = "tmp/background_blocks_done.csv"
  )
  background_raw <- fread("tmp/background_raw.csv") %>%
    mutate(download_date = as.character(Sys.Date()))

  ## Filter to precise coordinates inside the buffer, dedup, trim to target ----
  # Drop observations with obscured or missing coordinates, restrictive licenses.
  # Keep 1 row/observation, restrict to 100-km buffer.
  # Spatially restrict to the 100-km buffer, and subsample to the target region.
  # Keep  `user_login` and `license` are kept for contributor attribution 
  # and to ensure credit + reproducibility.
  
  background_clean <- background_raw %>%
    filter(
      !is.na(latitude), !is.na(longitude),
      !coordinates_obscured %in% c("true", TRUE),
      license != "") %>%
    distinct(id, .keep_all = TRUE) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = proj.wgs84, remove = FALSE) %>%
    st_filter(st_make_valid(background_region)) %>%
    st_drop_geometry() %>%
    select(any_of(c("id", "url", "user_login", "license", "observed_on", "datetime",
                    "latitude", "longitude", "positional_accuracy",
                    "coordinates_obscured", "quality_grade", "iconic_taxon_name", "download_date"))) %>%
    mutate(download_date = as.character(Sys.Date())) %>%
    slice_sample(n = min(n_background, nrow(.)))

  stopifnot("Background download returned no records" = nrow(background_clean) > 0)
  fwrite(background_clean, background_file)
}

collisions_background <- fread(background_file) %>%
  mutate(yday = yday(as_date(observed_on))) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = proj.wgs84, remove = FALSE) %>%
  st_transform(crs(bh_raster))

# Extract building height ----
# Use = collision sites; available = background records.
bh_use <- collisions_noam %>%
  extractBuildingHeight(bh_raster) %>%
  mutate(sample = "collision site")

bh_background <- collisions_background %>%
  extractBuildingHeight(bh_raster) %>%
  mutate(sample = "background (available)")

bh_compare <- bind_rows(bh_use, bh_background) %>%
  mutate(sample = factor(sample, levels = c("collision site", "background (available)")))

# Assemble modelling df ----
# Use (collision = 1) vs available (background = 0), with building height and yday.
model_data <- bind_rows(
  tibble(used = 1, building_height = bh_use$building_height,        yday = collisions_noam$yday),
  tibble(used = 0, building_height = bh_background$building_height, yday = collisions_background$yday)
)

# Exploratory comparison ----
# Okabe-Ito blue for collisions: colourblind-safe and distinct from the green
# buffer (red-on-green would fail for red-green colourblind readers).
collisionColor  <- "#0072B2"
backgroundColor <- "grey45"
sampleColors <- c(
  "collision site"         = collisionColor,
  "background (available)" = backgroundColor
)

# Building height: collisions vs the effort background.
(p_bh_density <- ggplot(bh_compare, aes(x = building_height, colour = sample)) +
  geom_density() +
  scale_x_continuous("Building height (m)") +
  scale_colour_manual(NULL, values = sampleColors))

# Day of year: collisions concentrate in migration windows; the background traces
# observer effort (the contrast the phenology term formalises).
(p_yday_density <- model_data %>%
  mutate(sample = factor(used, levels = c(0, 1),
                         labels = c("background (available)", "collision site"))) %>%
  ggplot(aes(x = yday, colour = sample)) +
  geom_density() +
  scale_x_continuous("Day of year",
                     breaks = monthDayYear_to_yday(monthFirsts),
                     labels = format(mdy(monthFirsts), "%b")) +
  scale_colour_manual(NULL, values = sampleColors))

# Map of the sampling design ----
# The 100-km buffer (the "available" region), background points, and collisions.
na_basemap <- ne_countries(scale = "medium", continent = "North America",
                           returnclass = "sf")

(p_design_map <- ggplot() +
  geom_sf(data = na_basemap, fill = "grey99", colour = "grey80", linewidth = 0.2) +
  geom_sf(data = background_region, fill = "#159367", colour = "#159367", alpha = 0.05, linewidth = 0.03) +
  geom_sf(data = st_transform(collisions_background, proj.wgs84),
          colour = "grey40", size = 0.2, alpha = 0.4) +
  geom_sf(data = st_transform(collisions_noam, proj.wgs84),
          colour = collisionColor, size = 0.45) +
  coord_sf(xlim = c(-128, -66), ylim = c(7, 53)) +
  theme_void())

# Example location ----
# A single collision point, its 1-km surroundings, and local building heights,
# for a methods/illustration figure. Swap the id as desired.
example_pt <- collisions_noam[1, ]
bh_example <- crop(bh_raster, vect(st_buffer(example_pt, dist = 3e3)))

(p_bh_example <- ggplot() +
  geom_spatraster(data = bh_example, maxcell = 5e6) +
  geom_sf(data = st_buffer(example_pt, dist = 1e3), fill = NA, colour = "red") +
  geom_sf(data = example_pt, colour = "red") +
  scale_fill_viridis_c("Building\nheight (m)", option = "turbo", na.value = NA))

# Model (Recommendation B) ----
library(brms)
library(splines)

# Use-availability logistic model: does collision favour taller buildings and
# migration-window timing, relative to where/when people observe? Building height
# is log1p-transformed (right-skewed; ~half the background is 0 m). The intercept
# is not interpretable in a use-availability design (it tracks the use:available
# ratio); coefficients express relative selection. brms conventions follow
# script 3 (ns() splines, default priors, cmdstanr, cached fits).
m_useavail <- brm(
  bf(used ~ log1p(building_height) + ns(yday, 5)),
  data = model_data,
  family = bernoulli(),
  cores = 4,
  chains = 4,
  seed = 42,
  iter = 5000,
  warmup = 1000,
  threads = threading(4),
  backend = "cmdstanr",
  file = "out/models/m_useavail_bh_yday.rds"
)
m_useavail <- add_criterion(m_useavail, "loo")

# conditional_effects(m_useavail)
pp_check(m_useavail, ndraws = 100)
bayes_R2(m_useavail)

# Shape check: is log1p adequate, or is the building-height response curved?
m_useavail_quad <- brm(
  bf(used ~ log1p(building_height) + I(log1p(building_height)^2) + ns(yday, 5)),
  data = model_data,
  family = bernoulli(),
  cores = 4,
  chains = 4,
  seed = 42,
  iter = 5000,
  warmup = 1000,
  threads = threading(4),
  backend = "cmdstanr",
  file = "out/models/m_useavail_bhquad_yday.rds"
)
m_useavail_quad <- add_criterion(m_useavail_quad, "loo")
loo_compare(m_useavail, m_useavail_quad)

conditional_effects(m_useavail_quad)

# Additional covariates (lower priority) ----
# Land cover and ALAN are expected to correlate strongly with building height, so
# they may not enter; check collinearity (e.g., VIF) before adding. Distance to
# water is a candidate third covariate.

# Radar migration term ----
# Nightly bird-migration traffic (Dark Ecology; a deliberate proxy for bat
# migration, see Discussion), matched to each record's nearest NEXRAD station on
# the NIGHT PRIOR to discovery. Bats found in the morning struck the previous
# night, so we use the local observation date minus one day; Dark Ecology "night"
# rows are dated by their evening onset, so date-1 indexes that night. We match to
# the single nearest station and keep only matches within radar range (<200 km),
# which also drops the few non-US records. Time-zone note: observed_on is the
# observer's local date; since morning discoveries are far from midnight, the
# night-prior assignment is robust to local/UTC date ambiguity.

stations <- fread("data/darkecology/nexrad-stations.csv") %>%
  st_as_sf(coords = c("lon", "lat"), crs = proj.wgs84, remove = FALSE)

radar_night <- rbindlist(lapply(2019:2025, function(y) {
  fread(sprintf("data/darkecology/daily/%d-daily.csv", y))[
    period == "night", .(station, radar_date = as.IDate(date), traffic)]
}))

# Nearest station, distance, and the night-prior date for a set of points.
radarCols <- function(points_sf, obs_date) {
  pts <- st_transform(points_sf, proj.wgs84)
  nearest <- st_nearest_feature(pts, stations)
  data.table(
    station    = stations$callsign[nearest],
    dist_km    = as.numeric(st_distance(pts, stations[nearest, ], by_element = TRUE)) / 1000,
    radar_date = as.IDate(obs_date - 1)
  )
}

radar_data <- rbind(
  cbind(data.table(used = 1, building_height = bh_use$building_height,
                   yday = collisions_noam$yday),
        radarCols(collisions_noam, collisions_noam$date)),
  cbind(data.table(used = 0, building_height = bh_background$building_height,
                   yday = collisions_background$yday),
        radarCols(collisions_background, as_date(collisions_background$observed_on)))
) %>%
  merge(radar_night, by = c("station", "radar_date"), all.x = TRUE) %>%
  filter(!is.na(traffic), dist_km < 200)

# Nightly traffic on collision vs background nights.
(p_traffic_density <- radar_data %>%
  mutate(sample = factor(used, levels = c(0, 1),
                         labels = c("background (available)", "collision site"))) %>%
  ggplot(aes(x = traffic + 1, colour = sample)) +
  geom_density() +
  scale_x_log10("Nightly migration traffic (night prior)") +
  scale_colour_manual(NULL, values = sampleColors))

# Smoke test: each focal point drawn to its assigned nearest station. Lines should
# be short and reach a nearby station (diagnostic, not a manuscript figure).
radar_pts <- rbind(
  transmute(st_transform(collisions_noam, proj.wgs84), sample = "collision"),
  transmute(st_transform(collisions_background, proj.wgs84), sample = "background")
)
radar_near <- st_nearest_feature(radar_pts, stations)
radar_crd  <- st_coordinates(radar_pts)
radar <- data.frame(
  sample  = radar_pts$sample,
  pt_lon  = radar_crd[, 1], pt_lat = radar_crd[, 2],
  st_lon  = stations$lon[radar_near], st_lat = stations$lat[radar_near],
  dist_km = as.numeric(st_distance(radar_pts, stations[radar_near, ], by_element = TRUE)) / 1000
) %>%
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
ggsave("figs/radar_station_assignment.png", p_station_radar, width = 8, height = 5, dpi = 150, bg = "white")

# Competing models on the radar-matched subset (same rows, for a fair LOO):
# seasonal spline only, night traffic only, and both. Traffic is log1p-transformed
# (highly right-skewed). "Both" lets traffic express the night-level anomaly beyond
# the seasonal mean (traffic and ns(yday) share only ~27% of variance).
m_season_sub <- brm(
  bf(used ~ log1p(building_height) + ns(yday, 5)),
  data = radar_data, family = bernoulli(),
  cores = 4, chains = 4, seed = 42, iter = 5000, warmup = 1000,
  threads = threading(4), backend = "cmdstanr",
  file = "out/models/m_useavail_season_sub.rds"
)
m_radar_only <- brm(
  bf(used ~ log1p(building_height) + log1p(traffic)),
  data = radar_data, family = bernoulli(),
  cores = 4, chains = 4, seed = 42, iter = 5000, warmup = 1000,
  threads = threading(4), backend = "cmdstanr",
  file = "out/models/m_useavail_radaronly.rds"
)
m_radar_both <- brm(
  bf(used ~ log1p(building_height) + ns(yday, 5) + log1p(traffic)),
  data = radar_data, family = bernoulli(),
  cores = 4, chains = 4, seed = 42, iter = 5000, warmup = 1000,
  threads = threading(4), backend = "cmdstanr",
  file = "out/models/m_useavail_radar_both.rds"
)
m_season_sub <- add_criterion(m_season_sub, "loo")
m_radar_only <- add_criterion(m_radar_only, "loo")
m_radar_both <- add_criterion(m_radar_both, "loo")

# "Both" is expected to win; traffic stays positive with the spline in the model,
# i.e. collisions track actual migration intensity on the night before, beyond
# seasonal timing.
loo_compare(m_season_sub, m_radar_only, m_radar_both)
fixef(m_radar_both)["log1ptraffic", ]
conditional_effects(m_radar_both)

# Figures ----
# Report-quality drafts. Colours reuse sampleColors; "#159367" marks the buffer.

## Sampling design map ----
# North and South America land with subtle state/province boundaries and the Great
# Lakes, the 100-km "available" buffer, background points (grey), and collisions
# (red) on top. Projected to a North America Lambert azimuthal equal-area (matches
# the geodesic buffer, so the radii read as circles). State lines and lakes are
# pulled with ne_download() (not ne_states(), which needs the CRAN-unavailable
# rnaturalearthhires) and cached to tmp/. State lines are 1:10m because the 1:50m
# layer omits Mexico's states. South America is included so the northern part of
# the continent that falls in frame gets its land, borders, and admin-1 lines
# (coord_sf clips the rest).
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
  geom_sf(data = st_transform(collisions_background, proj.wgs84),
          aes(colour = "Background observations"), size = 0.2, alpha = 0.3) +
  geom_sf(data = st_transform(collisions_noam, proj.wgs84),
          aes(colour = "Bat collisions"), size = 0.9) +
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
ggsave("figs/B_sampling_map.png", f_map, width = 7, height = 5, dpi = 600, bg = "white")
ggsave("figs/B_sampling_map.svg", f_map, width = 7, height = 5)

## Conditional effects of the top model (building height + season + night traffic) ----
ce <- conditional_effects(m_radar_both)
# Recompute the traffic effect on a log-spaced grid (other covariates held at their
# mean, as conditional_effects does) so the line is smooth on the log axis instead
# of segmented where the default linear grid is sparse.
traffic_grid <- data.frame(
  building_height = mean(radar_data$building_height),
  yday = mean(radar_data$yday),
  traffic = 10^seq(log10(min(radar_data$traffic[radar_data$traffic > 0])),
                   log10(max(radar_data$traffic)), length.out = 300)
)
traffic_epred <- posterior_epred(m_radar_both, newdata = traffic_grid)
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

# `dat` (optional) overlays the records as rugs: collisions along the top (accent,
# prominent), background along the bottom (subtle) — the binary-response analogue of
# the data points on the Kansas City brms panel. `y_accuracy` (optional) fixes the
# y-tick decimal places so several panels get identical axis-text width, which keeps
# their panel regions exactly equal when arranged side by side (see datatop combo).
cePanel <- function(v, logx = FALSE, dat = NULL, y_accuracy = NULL) {
  dd <- ceData[[v]]
  dd$x <- dd[[v]]
  g <- ggplot(dd, aes(x, estimate__)) +
    geom_ribbon(
      aes(ymin = lower__, ymax = upper__),
      # fill = "grey85"
      fill = scales::muted(collisionColor, l = 95, c = 20),  # "#A2CAFB"
      ) +
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
(f_condeffects <- cePanel("building_height") + cePanel("yday") +
    cePanel("traffic", logx = TRUE) +
    plot_layout(nrow = 1) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(axis.text.x = element_text(size = 7)))
ggsave("figs/B_conditional_effects.png", f_condeffects, width = 10, height = 3.2, dpi = 600)
ggsave("figs/B_conditional_effects.svg", f_condeffects, width = 10, height = 3.2)

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
effort_month <- fread("tmp/inat_background_effort.csv")[
  , .(share = sum(n_effort)), by = .(m = month(month_start))]
effort_month[, share := share / sum(share)]
month_shares <- rbindlist(list(
  cbind(effort_month, src = "effort (reference)"),
  data.table(m = month(as_date(collisions_background$observed_on)))[
    , .N, by = m][, .(m, share = N / sum(N), src = "background")],
  data.table(m = month(collisions_noam$date))[
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
# 6) — the three commonly-identified species plus an "Other/unknown" catch-all.
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
collision_species <- collisions_noam %>%
  st_drop_geometry() %>%
  transmute(yday, group = factor(
    if_else(
      quality_grade == "research" & scientific_name %in% speciesFocal,
      scientific_name,
      "Other/unknown"
    ),
    levels = speciesLevels
  ))

# Stacked day-of-year counts. `boundary = 0` anchors the bins so none straddle the
# Jan/Dec edges and get dropped. Wide-and-short by design (see layout below). The
# legend is reversed (`reverse = TRUE` reorders the keys only, not the stack, so
# Lasiurus stays at the bottom of the bars) and sits as a vertical block low and
# inside the panel, nestled in the summer trough between the two seasonal peaks.
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
          legend.text = element_markdown(size = 8),   # match the map legend size
          legend.key.size = unit(10, "pt"),           # key size + spacing match the map legend
          legend.background = element_blank(),
          legend.key.spacing.y = unit(3, "pt"))
}

## Combined main figure: map (a) over conditional effects (b-d) ----
# Layout and tagging follow the manuscript combos (e.g. F4_combo): a design-string
# with the map spanning the top row and the three effect panels beneath.
# NOTE: several layout drafts are kept below (2x2, map-on-top, map-left, and the
# data-over-inference version with the species panel) pending a final choice.

# Map on top
# (f_main_combo <- f_map + cePanel("building_height") + cePanel("yday") +
#    cePanel("traffic", logx = TRUE) +
#    plot_layout(design = "AAA\nBCD", heights = c(1.9, 1)) +
#    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
#    theme(plot.tag = element_text(size = 10), axis.text.x = element_text(size = 8)))

# 2x2
(f_main_combo <- f_map + cePanel("building_height") + cePanel("yday") +
    cePanel("traffic", logx = TRUE) +
    plot_layout(design = "AB\nCD", heights = c(1, 1)) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 10), axis.text.x = element_text(size = 8)))
ggsave("figs/B_main_combo.png", f_main_combo, width = 8, height = 8.4, dpi = 600, bg = "white")
ggsave("figs/B_main_combo.svg", f_main_combo, width = 8, height = 8.4)

# Version with the underlying records overlaid as rugs (collisions top, background bottom).
(f_main_combo_data <- f_map + cePanel("building_height", dat = radar_data) +
    cePanel("yday", dat = radar_data) + cePanel("traffic", logx = TRUE, dat = radar_data) +
    plot_layout(design = "AB\nCD", heights = c(1, 1)) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 10), axis.text.x = element_text(size = 8)))
ggsave("figs/B_main_combo_data.png", f_main_combo_data, width = 8, height = 6, dpi = 600, bg = "white")
ggsave("figs/B_main_combo_data.svg", f_main_combo_data, width = 8, height = 8)


# Another layout: top row (species bar plot + map) over a bottom row of the three
# effect panels. Built as two separate patchworks and then stacked, so the bottom
# three panels stay equal-width and aligned with each other regardless of the map's
# fixed aspect ratio (in a single shared grid the map pulls the columns it spans —
# and therefore c/d/e — out of alignment). The top-row widths are free to differ.
# Tags are set per panel (not via plot_annotation) because auto-tagging treats each
# nested patchwork as one unit and would only tag the two rows.
p_species <- speciesPanel()                                   + labs(tag = "(a)")
p_map     <- f_map                                            + labs(tag = "(b)")
# y_accuracy = 0.001 gives all three the same-width y labels (e.g. 0.600 / 0.075 /
# 0.040), so their panel regions come out exactly equal.
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
coefs <- as.data.frame(fixef(m_radar_both))[c("log1pbuilding_height", "log1ptraffic"), ]
coefs$term <- c("Building height (log1p)", "Night traffic (log1p)")
p_coef <- ggplot(coefs, aes(Estimate, term)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.2) +
  geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5)) +
  scale_x_continuous("Coefficient (log-odds of use vs available)") +
  ylab(NULL)
loo_tab <- as.data.frame(loo_compare(m_season_sub, m_radar_only, m_radar_both))
loo_labs <- c(m_radar_both = "building height + season + traffic",
              m_radar_only = "building height + traffic",
              m_season_sub = "building height + season")
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
