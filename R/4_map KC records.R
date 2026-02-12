
# Setup ==============
source("R/0_funs.R")
library(sf)
library(geodata)
library(rinat)
library(tidyterra)
library(data.table)


# Map locations --------------
## Load iNat records -----------
if(!file.exists("tmp/downloads.DeAnnObs.csv")) {
  focalObs <- get_inat_obs_user("redtail5")
  
  mySearch.DeAnnObs <- searchBuilder(user_id="redtail5", taxon_id = 40268)
  howManyResults(mySearch.DeAnnObs)
  
  downloads.DeAnnObs <- downloadResults(mySearch.DeAnnObs)
  
  fwrite(downloads.DeAnnObs, file="tmp/downloads.DeAnnObs.csv", row.names = F)
}

# 39.110600348029564&nelng=-94.49173146991824&page=2&swlat=39.07396161025084&swlng=-94.65489560871218

# Retain only those points observed on surveyed dates.
sd <- read.csv("data/derived/structured_surveys_schedule.csv") %>% 
  dplyr::filter(survey == TRUE) %>% 
  dplyr::select(date) %>% 
  mutate(date = date(date))


df <- read.csv("tmp/downloads.DeAnnObs.csv") %>% 
  dplyr::filter(
    latitude <= 39.110600348029564, latitude >= 39.07396161025084,
    longitude <= -94.56, longitude >= -94.6
    ) %>% 
  mutate(
    datetime = as_datetime(datetime),
    date = date(datetime)
  ) %>% 
  inner_join(sd)


df %>% 
  ggplot() +
  aes(x = longitude, y = latitude) +
  geom_point() +
  coord_quickmap(xlim = c(-94.55, -94.65))

## Load OSM data -------
# These take a minute to load so let's save and reload the outputs.
library(osmdata)
if(!file.exists("tmp/KC_osm.RData")) {
  
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
  
  save(KC_bb, KC_water, KC_highways, KC_highways.lines, KC_buildings, 
       file = "tmp/KC_osm_highRes.RData")
  
  ext_small <- st_bbox(c(xmin = -94.56, xmax = -94.6, ymin = 39.08, ymax = 39.108)) %>% 
    st_as_sfc()

  KC_highways.lines_small <- KC_highways.lines %>% 
    st_crop(ext_small) %>% 
    st_simplify(dTolerance = 0.1)
  KC_buildings_small <- KC_buildings$osm_polygons %>% 
    st_crop(ext_small) %>% 
    st_simplify(dTolerance = 0.1)
  
  save(KC_highways.lines_small, KC_buildings_small, 
       file = "tmp/KC_osm.RData")

} else if(!exists("KC_buildings")) {
  load("tmp/KC_osm.RData")
}



## Generate point density estimates of discovered bats --------

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

# Quick plot.
ggplot() + 
  geom_sf(final_pols, mapping = aes(fill=nlevel), color = NA)


## Plot ------------------------------

# Setup for plots.
arrow <- list(
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey60", "white")
  ) ,
  ggspatial::annotation_north_arrow(
    location = "br", which_north = "true",
    pad_x = unit(0, "in"), pad_y = unit(0.2, "in"),
    style = ggspatial::north_arrow_minimal(
      line_col = "grey20"
    )
  )
)

myBox_small <- st_bbox(c(xmin = -94.595-0.002, xmax = -94.57+0.0015, ymin = 39.08+0.0028, ymax = 39.108), crs =  proj.wgs84) %>% st_as_sfc() %>% st_as_sf()
myBox_medium <- st_bbox(c(xmin = -94.595-0.75, xmax = -94.57+0.75, ymin = 39.095-0.5, ymax = 39.108+0.55), crs =  proj.wgs84) %>% st_as_sfc() %>% st_as_sf()

states <- usmap::us_map(exclude = c("Alaska", "Hawaii", "Puerto Rico"))

# outlineColor1 <- "#00CC66"
outlineColor1 <- "#159367"
outlineColor2 <- "#3A0E4E"

### Make small map of city and density of occurrences ----

p_small <- ggplot() +
  geom_sf(
    data = KC_highways.lines_small,
    aes(linewidth = lanes),
    inherit.aes = FALSE,
    color = "grey80"
  ) +
  scale_linewidth_continuous(range = c(0.6, 4), guide = "none") +
  geom_sf(
    data = KC_buildings_small,
    inherit.aes = FALSE,
    fill = "grey40",
    color = NA,
    size = 0.2,
  ) +
  geom_sf(final_pols, mapping = aes(fill = nlevel), color = NA) +
  scale_fill_viridis_c("Density of\ndetections", option = "turbo", begin = 0.6, limits = c(0,1)) +
  geom_sf(
    data = KC_buildings_small,
    inherit.aes = FALSE,
    fill = "grey50",
    color = NA,
    size = 0.2,
    alpha = 0.5,
  ) +
  geom_sf(myBox_small, mapping = aes(), fill = NA, linewidth = 1, color = outlineColor2) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_sf(crs = proj.wgs84, xlim = st_bbox(myBox_small)[c(1,3)], ylim =st_bbox(myBox_small)[c(2,4)]) +
  arrow +
  theme_void() +
  theme(
    axis.title = element_blank(),
    margins = margin(l = 0, unit = "pt"),
    legend.key.width = unit(5, "mm")
    ) 

p_small_nolegend <- p_small +
  theme(legend.position = "none")

p_small_legendInset <- p_small +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.003,0.003),
    legend.justification = c(0,0),
    legend.background = element_rect(fill = "white", color = "grey40", linewidth = 0.1),
    legend.margin = margin(1,1,3,1,unit="mm"),
    )
 
# ggsave(p_small, filename = "figs/map_city.svg", width = 8,  height = 8)
ggsave(p_small_legendInset, filename = "figs/map_city.png", width = 3.08*1.5,  height = 3.25*1.5, dpi = 300)


### Make national map ----

p_big <- ggplot() +
  geom_sf(states, mapping = aes(), fill = "white", color = "grey20") +
  geom_sf(myBox_medium, mapping = aes(), fill = NA, linewidth = 1, color = outlineColor1, alpha = 0.5) +
  theme_void() +
  arrow
ggsave(p_big, filename = "figs/map_country.png", width = 1.69*2,  height = 1.11*2, dpi = 300)
ggsave(p_big, filename = "figs/map_country.svg", width = 1.69*2,  height = 1.11*2)


### Make regional map -----

# Below code uses an object, largeFileStoragePath, that must be set to the location
# of a directory containing NLCD data in a directory "nlcd".
# It will also download GADM data there in a folder `geodata` will create ("gadm").

# Load land cover dataset from NLCD.
nlcd0 <- rast(file.path(largeFileStoragePath, "nlcd/nlcd_2021_land_cover_l48_20230630/nlcd_2021_land_cover_l48_20230630.img"))
nlcd_crop_medium <- crop(nlcd0, project(terra::vect(myBox_medium),crs(nlcd0)))
myBox_medium_acea <- st_transform(myBox_medium, st_crs(nlcd_crop_medium))

geodata::geodata_path(largeFileStoragePath)
usa_48 <- geodata::gadm(country = "USA", path = "../../- Missions & Programs/Research & Development/Data Products/") %>% 
  st_as_sf() %>% 
  dplyr::filter(!NAME_1 %in% c("Alaska", "Hawaii", "Puerto Rico")) %>% 
  vect()

# Plotting info for NLCD
colorBreaks <- FedData::nlcd_colors()$Class
colorBreaks <- colorBreaks %>% 
  gsub("Pasture/Hay", "Hay/Pasture", .) %>% 
  gsub("Developed High Intensity", "Developed, High Intensity", .)%>% 
  gsub("Herbaceous", "Sedge/Herbaceous", .)
colorBreaks <- c(colorBreaks, "Herbaceous","Emergent Herbaceous Wetlands","","Barren Land", "Perennial Snow/Ice", "Unclassified")
nlcdCols <- c(
  FedData::nlcd_colors()$Color,
  "#FDE9AA",
  "#64B3D5",
  "grey50",
  "grey30",
  "grey80",
  "grey50"
)
landCovercolors <- list(
  scale_fill_manual(
    "Land cover class",
    breaks = colorBreaks,
    values = nlcdCols
  ) 
)

p_medium <- ggplot() +
  geom_spatraster(nlcd_crop_medium, mapping = aes(), maxcell = 20e5) +
  geom_sf(usa_48, mapping = aes(), fill = NA, linewidth = 0.25, color = "grey20") +
  landCovercolors +
  geom_sf(myBox_small, mapping = aes(), fill = NA, linewidth = 0.65, color = outlineColor2) +
  geom_sf(myBox_medium, mapping = aes(), fill = NA, linewidth = 1.15, color = outlineColor1) +
  theme_void() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_sf(crs = proj.wgs84, xlim = st_bbox(myBox_medium)[c(1,3)], ylim =st_bbox(myBox_medium)[c(2,4)]) +
  arrow

p_medium_noLegend <- p_medium +
  theme(legend.position = "none")
ggsave(p_medium_noLegend, filename = "figs/map_nlcd_region.png", width = 3.54*1.5,  height = 3.2*1.5, dpi = 300)
ggsave(p_medium_noLegend, filename = "figs/map_nlcd_region.svg", width = 3.54*1.5,  height = 3.2*1.5)

#### Make legend ------

myLegend <- 
  { p_medium +
    guides(fill = guide_legend(title.position="left", title.hjust = 0.5, ncol = 5)) +
    theme(
      legend.title = element_text(angle = 90, hjust = 0.5)
      )
  }%>% 
  ggpubr::get_legend()
ggsave(myLegend, filename = "figs/map_nlcd_legend.png", width = 6.5*1.5,  height = 1.5*1.5, dpi = 300)
ggsave(myLegend, filename = "figs/map_nlcd_legend.svg", width = 6.5*1.5,  height = 1.5*1.5)
