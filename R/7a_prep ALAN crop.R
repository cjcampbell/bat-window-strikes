# 7a_prep ALAN crop.R
# Purpose: One-time derivation of the local ALAN raster used by
#          "7_analyze iNat records.R". Crops the global VIIRS VNL v2 2024 annual
#          nighttime-light composite to the North American collision background
#          region and writes a small, portable GeoTIFF into the repo.
# Source:  VIIRS VNL v2 (median-masked) 2024 global annual composite, Earth
#          Observation Group / Payne Institute (Elvidge et al. 2021). The raw
#          global file (~11 GB, EPSG:4326, nW/cm2/sr) lives outside the repo;
#          set `src` to its location. Treated as READ-ONLY source.
# Output:  data/ALAN/VNL_npp_2024_vcmslcfg_v2_NAcrop.tif  (~60 MB)
# Run:     once, before script 7. Safe to re-run (overwrites the crop).

source("R/0_funs.R")
library(sf)
library(terra)
library(rnaturalearth)
library(lubridate)

# Raw global VIIRS composite (read-only source, outside the repo).
src <- "/Users/ccampbell/Library/CloudStorage/Box-Box/- Missions & Programs/Research & Development/Data Products/ALAN/VIIRS VNL v2 annual composites/VNL_npp_2024_global_vcmslcfg_v2_c202502261200.median_masked.dat.tif"
dst_dir <- "data/ALAN"
dst     <- file.path(dst_dir, "VNL_npp_2024_vcmslcfg_v2_NAcrop.tif")
if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE)

# Rebuild the same 100-km geodesic background region script 7 uses, so the crop
# covers every use and available point with margin. (Mirrors the collisions_noam
# and background_region construction in script 7.)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(myproj) %>%
  select(continent)

collisions_noam <- read.csv("data/iNat_observations_tidy_manualChecks.csv") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = proj.wgs84, remove = FALSE) %>%
  st_transform(myproj) %>%
  filter(CJ.manual.check == "y") %>%
  st_join(world) %>%
  filter(continent == "North America", geoprivacy != "obscured") %>%
  mutate(date = as_date(as.POSIXct(datetime))) %>%
  filter(year(date) %in% 2019:2025) %>%
  st_transform(proj.wgs84)

background_region <- collisions_noam %>%
  st_buffer(dist = 100e3) %>%
  st_union() %>%
  st_make_valid()

# Crop with a 0.5-deg margin; VIIRS is native EPSG:4326 so bbox degrees align.
bb       <- st_bbox(background_region)
crop_ext <- ext(bb["xmin"] - 0.5, bb["xmax"] + 0.5, bb["ymin"] - 0.5, bb["ymax"] + 0.5)
alan_na  <- crop(rast(src), crop_ext)

writeRaster(alan_na, dst, overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
cat("wrote", dst, sprintf("(%d x %d)\n", nrow(alan_na), ncol(alan_na)))
