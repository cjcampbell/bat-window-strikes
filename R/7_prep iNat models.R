# 7_prep iNat models.R
# Purpose: Build the analysis-ready use/availability table for the iNaturalist
#          bat-window collision models. Downloads (once) and caches the target-group
#          iNaturalist background, crops the global ALAN raster to the study region
#          (once), then extracts every predictor -- building height, ALAN, and
#          nightly bird-migration traffic -- at all use (collision) and available
#          (background) points. Writes one derived table that "8_fit iNat models.R"
#          consumes; no models or figures live here. All predictor extraction (the
#          ALAN crop and radar match included) happens in this one place.
# Inputs:  data/iNat_observations_tidy_manualChecks.csv        (retained collisions)
#          data/building height data/GBH2020_150m_GEDI.tif      (Ma et al. 2023)
#          <global VIIRS VNL v2 2024 composite>                 (Elvidge et al. 2021;
#            ~11 GB, READ-ONLY, outside the repo; read only on the first run to
#            build the local crop -- see `alan_global` below)
#          data/darkecology/nexrad-stations.csv                 (NEXRAD stations)
#          data/darkecology/daily/{2019..2025}-daily.csv        (nightly traffic)
# Outputs: data/derived/inat_background.csv               (background draw, cached)
#          data/derived/inat_background_effort.csv        (effort reference)
#          data/ALAN/VNL_npp_2024_vcmslcfg_v2_NAcrop.tif  (local ALAN crop, cached)
#          data/derived/useavail_points.csv               (analysis-ready table)


# Setup ----
source("R/0_funs.R")
library(sf)
library(terra)
library(rnaturalearth)
library(lubridate)
library(data.table)
library(rinat)

# Extract building height at a set of points, treating unmapped (NA) cells as
# building-free ground (0 m). `points` must already be in the raster's CRS.
# `extract` is masked by tidyr, hence the terra:: qualifier.
extractBuildingHeight <- function(points, height_raster) {
  terra::extract(height_raster, points, ID = FALSE) %>%
    transmute(building_height = replace_na(building_height, 0))
}

# Extract VIIRS nighttime radiance (ALAN) at a set of points. Unlike building
# height, `points` are transformed to the raster's lon/lat (EPSG:4326) here, so the
# caller need not pre-project. Negative-radiance artifacts are clamped to 0, and
# any unmapped cell is treated as unlit (0), mirroring the building-height NA rule.
extractALAN <- function(points, alan_raster) {
  vals <- terra::extract(alan_raster, st_transform(points, crs(alan_raster)), ID = FALSE)$alan
  tibble(alan = pmax(replace_na(vals, 0), 0))
}

# Nearest NEXRAD station, its distance (km), and the night-prior date for a set of
# points. Bats found in the morning struck the previous night, so we index the
# local observation date minus one day; Dark Ecology "night" rows are dated by
# their evening onset, so date-1 indexes that night. Requires `stations` (below).
radarCols <- function(points_sf, obs_date) {
  pts <- st_transform(points_sf, proj.wgs84)
  nearest <- st_nearest_feature(pts, stations)
  data.table(
    station    = stations$callsign[nearest],
    dist_km    = as.numeric(st_distance(pts, stations[nearest, ], by_element = TRUE)) / 1000,
    radar_date = as.IDate(obs_date - 1)
  )
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

# Load collision records ----
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

# Load building-height raster ----
# GBH2020: global 150-m urban building height from spaceborne lidar (Ma et al.
# 2023, Scientific Data). ~1.5 GB; can take some time to load. (The ALAN raster is
# loaded further down, after the study region it is cropped to is defined.)
if (!exists("bh_raster")) {
  bh_raster <- rast("data/building height data/GBH2020_150m_GEDI.tif")
  names(bh_raster) <- "building_height"
}

# Precise North American collision points, restricted to 2019-2025. Six pre-2019
# records are dropped from both use and background to match iNaturalist data
# availability. Also drops observations with obscured (randomized within 0.2 deg)
# coordinates. Projected to the height raster's native (metric) CRS for extraction.
collisions_noam <- collisions %>%
  filter(continent == "North America", geoprivacy != "obscured") %>%
  mutate(date = as_date(as.POSIXct(datetime)),
         year = year(date),
         yday = yday(date)) %>%
  filter(year %in% 2019:2025) %>%
  st_transform(crs(bh_raster))

# 100-km geodesic buffer around collisions defines the "available" region.
# Restricting "available" to within 100 km of a collision matches the strong bias
# towards eastern North America records and removes continental gradients. Buffer
# geodesically on lon/lat (s2) for a true 100-km radius; buffering in the raster's
# metric CRS would skew the radius E-W and make the region anisotropic.
background_region <- collisions_noam %>%
  st_transform(proj.wgs84) %>%
  st_buffer(dist = 100e3) %>%
  st_union()
background_region <- st_sf(geometry = background_region)

# Load ALAN raster (crop once to the study region) ----
# VIIRS VNL v2 2024 annual composite (median-masked; Elvidge et al. 2021).
# Nighttime radiance in nW/cm2/sr; native EPSG:4326 (reprojected on extract). The
# global composite (~11 GB) is READ-ONLY and lives outside the repo. The first run
# crops it to the background region (0.5-deg margin; VIIRS is EPSG:4326, so bbox
# degrees align) and writes a small local GeoTIFF; every later run reads that crop
# directly and never touches the global file. Delete the local crop to force a
# re-crop.
alan_crop   <- "data/ALAN/VNL_npp_2024_vcmslcfg_v2_NAcrop.tif"
alan_global <- paste0("/Users/ccampbell/Library/CloudStorage/Box-Box/",
                      "- Missions & Programs/Research & Development/Data Products/ALAN/",
                      "VIIRS VNL v2 annual composites/",
                      "VNL_npp_2024_global_vcmslcfg_v2_c202502261200.median_masked.dat.tif")
if (!file.exists(alan_crop)) {
  if (!dir.exists("data/ALAN")) dir.create("data/ALAN", recursive = TRUE)
  bb       <- st_bbox(background_region)
  crop_ext <- ext(bb[["xmin"]] - 0.5, bb[["xmax"]] + 0.5, bb[["ymin"]] - 0.5, bb[["ymax"]] + 0.5)
  writeRaster(crop(rast(alan_global), crop_ext), alan_crop, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
  cat("cropped global ALAN composite to", alan_crop, "\n")
}
if (!exists("alan_raster")) {
  alan_raster <- rast(alan_crop)
  names(alan_raster) <- "alan"
}

# Target-group background ----
# The "available" set is ~10,000 general iNaturalist animal observations (any
# quality grade, matching how collisions were retained), representing where and
# WHEN people observe. Comparing collisions against this effort surface -- rather
# than against raw space and time -- supports both the structural and the
# day-of-year analyses. The download is expensive, so it runs once and is cached to
# data/derived/ (tracked, so the exact draw travels with the repo); later runs read
# that cached draw. To redraw, delete this file (and the tmp/ resumable
# intermediates) and rerun.
background_file <- "data/derived/inat_background.csv"
effort_file     <- "data/derived/inat_background_effort.csv"
n_background      <- 10000
records_per_block <- 50    # records kept per sampled (cell x month) block
oversample        <- 5     # draw extra to survive buffer / coordinate filtering

if (!file.exists(background_file)) {

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
  # Resumable: each block's kept records append to tmp/background_raw.csv, and its
  # block index is logged, so the download can stop and restart without repeating
  # finished blocks. `blocks` is deterministic given set.seed(42) + the arranged
  # effort table, so block indices are stable across restarts. `download_date` =
  # access date.
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
  # Drop observations with obscured or missing coordinates and restrictive licenses,
  # keep one row per observation, restrict to the 100-km buffer, and subsample to the
  # target region. `user_login` and `license` are kept for contributor attribution.
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
  file.copy("tmp/inat_background_effort.csv", effort_file, overwrite = TRUE)  # freeze effort ref
}

collisions_background <- fread(background_file) %>%
  mutate(date = as_date(observed_on),
         year = year(date),
         yday = yday(date)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = proj.wgs84, remove = FALSE) %>%
  st_transform(crs(bh_raster))

# Extract predictors ----
# All three predictors are extracted here, at the same use (collision) and
# available (background) points, so the modelling script never touches a raster or
# a station table. `stations` and `radar_night` are loaded first because radarCols()
# needs them.
stations <- fread("data/darkecology/nexrad-stations.csv") %>%
  st_as_sf(coords = c("lon", "lat"), crs = proj.wgs84, remove = FALSE)
radar_night <- rbindlist(lapply(2019:2025, function(y) {
  fread(sprintf("data/darkecology/daily/%d-daily.csv", y))[
    period == "night", .(station, radar_date = as.IDate(date), traffic)]
}))

## Building height and ALAN ----
bh_use          <- extractBuildingHeight(collisions_noam, bh_raster)
bh_background   <- extractBuildingHeight(collisions_background, bh_raster)
alan_use        <- extractALAN(collisions_noam, alan_raster)
alan_background <- extractALAN(collisions_background, alan_raster)

## Nightly migration traffic ----
# Nightly bird-migration traffic (Dark Ecology; a deliberate proxy for bat
# migration, see Discussion), matched to each point's nearest NEXRAD station on the
# night prior to discovery. Matches beyond radar range (>200 km) or with no traffic
# record are kept as NA here and dropped to the radar-matched subset in the
# modelling script (keeping one derived table for both the full and radar analyses).
radar_use <- cbind(collisions_noam[, c("longitude", "latitude")] %>% st_drop_geometry(),
                   radarCols(collisions_noam, collisions_noam$date))
radar_bg  <- cbind(collisions_background[, c("longitude", "latitude")] %>% st_drop_geometry(),
                   radarCols(collisions_background, collisions_background$date))

# Assemble one analysis-ready table (one row per use/available point) ----
# `used` = 1 collision, 0 background. Species fields are only meaningful for
# collisions (background lacks scientific_name); they drive the species figure.
points_use <- data.table(
  used            = 1L,
  longitude       = collisions_noam$longitude,
  latitude        = collisions_noam$latitude,
  date            = as.IDate(collisions_noam$date),
  yday            = collisions_noam$yday,
  year            = collisions_noam$year,
  building_height = bh_use$building_height,
  alan            = alan_use$alan,
  scientific_name = collisions_noam$scientific_name,
  quality_grade   = collisions_noam$quality_grade,
  iconic_taxon_name = collisions_noam$iconic_taxon_name,
  station         = radar_use$station,
  dist_km         = radar_use$dist_km,
  radar_date      = radar_use$radar_date
)
points_bg <- data.table(
  used            = 0L,
  longitude       = collisions_background$longitude,
  latitude        = collisions_background$latitude,
  date            = as.IDate(collisions_background$date),
  yday            = collisions_background$yday,
  year            = collisions_background$year,
  building_height = bh_background$building_height,
  alan            = alan_background$alan,
  scientific_name = NA_character_,
  quality_grade   = collisions_background$quality_grade,
  iconic_taxon_name = collisions_background$iconic_taxon_name,
  station         = radar_bg$station,
  dist_km         = radar_bg$dist_km,
  radar_date      = radar_bg$radar_date
)
points <- rbind(points_use, points_bg)

# Left-join nightly traffic by (station, radar_date), preserving row order and
# leaving unmatched nights as NA.
points[radar_night, on = .(station, radar_date), traffic := i.traffic]

fwrite(points, "data/derived/useavail_points.csv")
cat(sprintf("wrote data/derived/useavail_points.csv: %d rows (used %d, available %d); radar-matched %d\n",
            nrow(points), sum(points$used == 1), sum(points$used == 0),
            sum(!is.na(points$traffic) & points$dist_km < 200)))
