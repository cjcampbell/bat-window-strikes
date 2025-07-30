
# Setup ==============
source("R/0_funs.R")
library(sf)
library(geodata)
library(rinat)
library(tidyterra)
library(data.table)


# Map locations --------------
## Load iNat records -----------
if(!file.exists("tmp/downloads.DeAnnObs")) {
  focalObs <- get_inat_obs_user("redtail5")
  
  mySearch.DeAnnObs <- searchBuilder(user_id="redtail5", taxon_id = 40268)
  howManyResults(mySearch.DeAnnObs)
  
  downloads.DeAnnObs <- downloadResults(mySearch.DeAnnObs)
  
  fwrite(downloads.DeAnnObs, file="tmp/downloads.DeAnnObs.csv", row.names = F)
}

# 39.110600348029564&nelng=-94.49173146991824&page=2&swlat=39.07396161025084&swlng=-94.65489560871218

df <- read.csv("tmp/downloads.DeAnnObs.csv") %>% 
  dplyr::filter(
    latitude <= 39.110600348029564, latitude >= 39.07396161025084,
    longitude <= -94.56, longitude >= -94.6
    )

df %>% 
  ggplot() +
  aes(x = longitude, y = latitude) +
  geom_point() +
  coord_quickmap(xlim = c(-94.55, -94.65))

## Load OSM data -------
library(osmdata)
# https://jcoliver.github.io/learn-r/017-open-street-map.html
KC_bb <- getbb("Kansas City")
available_tags(feature = "highway")
highWayTags <- available_tags(feature = "highway")
available_tags(feature = "building")


KC_water <- opq(KC_bb) %>% 
  add_osm_feature(key = "water") %>%
  osmdata_sf()

KC_highways <- KC_bb %>%
  opq() %>%
  add_osm_feature(key = "highway") %>%
  osmdata_sf()
KC_highways$osm_lines$lanes %>% unique

KC_highways.lines <- KC_highways$osm_lines %>% 
  dplyr::mutate(
    lanes = as.numeric(lanes),
    lanes = case_when(lanes == 31 ~ 1, is.na(lanes) ~ 1, TRUE ~ lanes)
  )

KC_buildings <-  KC_bb %>%
  opq() %>%
  add_osm_feature(key = "building") %>% 
  osmdata_sf()


## Generate point density estimates --------

# Pull info from polygon to convert to sf.
# https://gis.stackexchange.com/questions/444689/create-density-polygons-from-spatial-points-in-r
p <- ggplot() +
  stat_density_2d(
    data = df, 
    aes(x = longitude, y = latitude, fill = after_stat(level)), geom = "polygon")

data2d <- layer_data(p)
data2d$pol <- paste0(data2d$group, "_", data2d$subgroup)
ids <- unique(data2d$pol)
pols <- lapply(ids, function(x){
  topol <- data2d[data2d$pol == x, ]
  
  closepol <- rbind(topol, topol[1, ])
  
  pol <- st_polygon(list(as.matrix(closepol[,c("x", "y")])))
  
  # Add features
  df <- unique(topol[, grepl("level", names(topol))])
  
  tofeatures <- st_as_sf(df, geometry=st_sfc(pol))
  
  return(tofeatures)  
})

final_pols <- do.call(rbind, pols)
st_crs(final_pols) <- st_crs(proj.wgs84)

ggplot() + 
  geom_sf(final_pols, mapping = aes(fill=nlevel), color = NA)


## Plot ------------------------------

### Zoomed-in version -----------

p_dens <- ggplot() +
  geom_sf(
    data = KC_highways.lines,
    aes(linewidth = lanes),
    inherit.aes = FALSE,
    color = "grey40"
  ) +
  scale_linewidth_continuous(range = c(0.6, 4), guide = "none") +
  geom_sf(
    data = KC_buildings$osm_polygons,
    inherit.aes = FALSE,
    fill = "grey80",
    color = NA,
    size = 0.2,
  ) +
  geom_sf(final_pols, mapping = aes(fill = nlevel), color = NA) +
  scale_fill_viridis_c("Density", option = "turbo", begin = 0.6, limits = c(0,1)) +
  geom_sf(
    data = KC_buildings$osm_polygons,
    inherit.aes = FALSE,
    fill = "grey70",
    color = NA,
    size = 0.2,
    alpha = 0.5,
  ) +
  coord_sf(xlim = c(-94.57, -94.595), ylim = c(39.095, 39.108)) +
  theme(axis.title = element_blank())

ggsave(p_dens, filename = "figs/KC_density.png", dpi = 600, width = 8, height = 8)



### Don't omit lower building (another version) ------
p_dens2 <- p_dens +
  coord_sf(
    xlim = c(-94.56, -94.6),
    ylim = c(39.108, 39.08)
  )
ggsave(p_dens2, filename = "figs/KC_density2.png", dpi = 600, width = 8, height = 8)


# Make regional maps -----------
# Load urban areas
urb <- st_read("data/cb_2020_us_ua20_500k/cb_2020_us_ua20_500k.shp")

# Load states 
usa <- geodata::gadm(
  "USA", level = 1,
  path = "../../- Missions & Programs/Research & Development/Data Products"
  ) %>%
  st_as_sf()

# Load building height data (big!).
bh <- rast("data/building height data/GBH2020_150m_GEDI.tif")
usa_mol <- st_transform(usa, st_crs(bh))

bh_crop <- usa_mol %>% 
  dplyr::filter(NAME_1 %in% c("Kansas", "Missouri")) %>% 
  terra::crop(bh, .)
ggplot() +
  geom_spatraster(bh_crop, mapping = aes()) +
  scale_fill_viridis_c(option = "magma", na.value = "white")

bh_kc <- urb %>% 
  st_transform(st_crs(bh)) %>% 
  dplyr::filter(NAME20 == "Kansas City, MO--KS") %>% 
  terra::crop(bh, .)
ggplot() +
  geom_spatraster(bh_kc, mapping = aes()) +
  scale_fill_viridis_c(option = "turbo", na.value = "white")

bh_kc_sml <- final_pols %>% 
  st_transform(st_crs(bh)) %>% 
  st_buffer(dist = 5e3) %>% 
  terra::crop(bh, .) %>% 
  project(final_pols)
ggplot() +
  geom_spatraster(bh_kc_sml, mapping = aes()) +
  scale_fill_viridis_c(option = "turbo", na.value = "white") +
  geom_sf(final_pols, mapping = aes(), fill = NA, color = "white")









p_region_bh <- ggplot() +
  geom_spatraster(bh, mapping = aes()) +
  geom_sf(
    data = KC_water$osm_polygons,
    inherit.aes = FALSE,
    fill = "blue", color = NA
  ) +
  geom_sf(
    data = KC_water$osm_lines,
    inherit.aes = FALSE,
    color = "blue"
  ) +
  geom_sf(
    data = KC_water$osm_multipolygons,
    inherit.aes = FALSE,
    fill = "blue", color = NA
  ) +
  geom_sf(
    data = usa,
    mapping = aes(),
    fill = NA,
    color = "green", linewidth = 1.1
  ) +
  coord_sf(xlim = c(-94.57 - 0.1, -94.595 + 0.1),
           ylim = c(39.095 - 0.1, 39.108 + 0.1)) 

ggsave(p_region_bh, filename = "figs/p_region_bh.png", dpi = 600, width = 8, height = 8)









p_region <- ggplot() +
  geom_sf(
    data = KC_water$osm_polygons,
    inherit.aes = FALSE,
    fill = "blue", color = NA
  ) +
  geom_sf(
    data = KC_water$osm_lines,
    inherit.aes = FALSE,
    color = "blue"
  ) +
  geom_sf(
    data = KC_water$osm_multipolygons,
    inherit.aes = FALSE,
    fill = "blue", color = NA
  ) +
  geom_sf(
    urb, mapping = aes(), color = NA, fill = "grey50") +
  geom_sf(
    data = usa,
    mapping = aes(),
    fill = NA,
    color = "green", linewidth = 1.1
  ) +
  coord_sf(xlim = c(-94.57 - 0.1, -94.595 + 0.1),
           ylim = c(39.095 - 0.1, 39.108 + 0.1)) 

ggsave(p_region, filename = "figs/KC_region.png", dpi = 600, width = 8, height = 8)
