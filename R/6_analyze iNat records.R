# Load --------------
library(tidyverse)
library(tidyterra)
library(sf)
library(terra)
library(rnaturalearth)
library(lubridate)
library(patchwork)
library(data.table)

source("R/0_funs.R")

speciesColors <- c(
  "Big brown" = "#A16928",
  "Evening" = "#124559", 
  "Vespertilionidae" = "grey50", 
  "Eastern Red" = "#E46844", 
  "Tricolored" = "#f4d35e", 
  "Silver-haired" = "#484357")

myproj <- "+proj=eqearth +lon_0=0 +datum=WGS84 +units=km +no_defs"
proj.wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

taxTree <- fread("data/iNat_observations_taxTree.csv")
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  sf::st_transform(myproj) %>% 
  dplyr::select(admin, adm0_a3, continent)
mygrid <- sf::st_make_grid(world, cellsize = 500) %>%
  st_as_sf() %>% 
  dplyr::mutate(grid_id = row_number())


# Load manually-checked data
df <- read.csv("data/iNat_observations_tidy_manualChecks.csv") %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = proj.wgs84) %>% 
  st_transform(myproj) %>% 
  dplyr::filter(CJ.manual.check == "y") %>% 
  mutate(
    datetime.POSIX = as.POSIXct(datetime),
    datetime.date = as_date(datetime.POSIX),
    yday = yday(datetime.date)
  ) %>% 
  # Join with taxonomy data.
  left_join(., taxTree, by = c("scientific_name" = "search")) %>% 
  # Unite with spatial data
  st_join(., world)



## Correct  spatial data -----

df %>% dplyr::filter(is.na(continent))
# A handful don't overlap countries (probably due to obscured coordinates near the coast).
# To resolve this, find the nearest country and replace NA's.
for(i in which(is.na(df$continent)) ) {

  workingDF <- df[i,] %>% 
    st_distance(world) %>% 
    which.min() %>% 
    world[.,] %>% 
    as.data.frame() %>% 
    dplyr::select(admin, adm0_a3, continent)
  
  df$admin[i]     <- workingDF$admin
  df$adm0_a3[i]   <- workingDF$adm0_a3
  df$continent[i] <- workingDF$continent
}


# Map observations -------

df_grid <- st_join(df, mygrid) # %>% 
 #  st_join(dplyr::select(world, admin, adm0_a3, continent)) # %>% 
  # dplyr::mutate(
  #   admin = case_when(
  #     place_state_name == "New York" & place_country_name == "United States" ~ "United States of America",
  #     place_state_name == "Victoria" & place_country_name == "Australia" ~ "Australia",
  #     place_guess == "Sri Lanka" & id == "131693823" ~ "Sri Lanka",
  #     TRUE ~ admin
  #     ),
  #   adm0_a3 = case_when(
  #     place_state_name == "New York" & place_country_name == "United States" ~ "USA",
  #     place_state_name == "Victoria" & place_country_name == "Australia" ~ "AUS",
  #     place_guess == "Sri Lanka" & id == "131693823" ~ "LKA",
  #     TRUE ~ adm0_a3
  #   ),
  #   continent = case_when(
  #     place_state_name == "New York" & place_country_name == "United States" ~ "North America",
  #     place_state_name == "Victoria" & place_country_name == "Australia" ~ "Oceania",
  #     place_guess == "Sri Lanka" & id == "131693823" ~ "Asia",
  #     admin == "Russia" & place_state_name == "Krasnoyarsk" ~ "Asia",
  #     TRUE ~ continent
  #   )
  # )

df_grid_count_grid <- df_grid %>% 
  count(grid_id) %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry) %>% 
  inner_join(mygrid) %>% 
  st_as_sf()

df_grid_count <- df_grid_count_grid %>% 
  st_centroid()

ggplot() +
  geom_sf(world, mapping = aes(), fill = "white") +
  geom_sf(df, mapping = aes(color = quality_grade)) 

ggplot() +
  geom_sf(world, mapping = aes(), fill = NA) +
  geom_sf(df_grid_count, mapping = aes(size = n), alpha = 0.5, color = "magenta") +
  scale_size_continuous(
    breaks = c(1,10,40,70),
    )

ggplot() +
  geom_sf(world, mapping = aes(), fill = NA) +
  geom_sf(df_grid_count_grid, mapping = aes(fill = n), color = NA) +
  scale_fill_viridis_c(option = "turbo", begin = 0.3)


iNat_map_grid_circle <- ggplot() +
  geom_sf(world, mapping = aes(), fill = "grey95", color = "grey60", linewidth = 0.1) +
  geom_sf(df_grid_count_grid, mapping = aes(), color = "grey50", fill = NA) +
  geom_sf(df_grid_count, mapping = aes(size = n), alpha = 0.5, color = "#00B31B") +
  scale_size_continuous(
    "Number of records",
    breaks = c(1,10,40,70),
  ) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  coord_sf(ylim = c(-6792.374, 8000)) +
  theme_void() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25)
  )
ggsave(iNat_map_grid_circle, filename =  "figs/iNat_map_grid_circle.png", dpi = 600)

## By family -----
df_grid_family <- df_grid %>% 
  count(grid_id, family_name) %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry) %>% 
  inner_join(mygrid) %>% 
  st_as_sf()

df_grid_count_family <- df_grid_family %>% 
  st_centroid()


df_grid_count_family %>% 
  ggplot() +
  geom_sf(world, mapping = aes(), fill = NA) +
  geom_sf(mapping = aes(size = n, color = family_name), alpha = 0.5) +
  scale_size_continuous(
    breaks = c(1,10,40,70),
  )

ggplot() +
  geom_sf(world, mapping = aes(), fill = "grey95", color = "grey60", linewidth = 0.1) +
  geom_sf(df_grid_family, mapping = aes(), color = "grey50", fill = NA) +
  geom_sf(df_grid_count_family, mapping = aes(size = n), alpha = 0.5, color = "#00B31B") +
  scale_size_continuous(
    "Number of records",
    breaks = c(1,10,40,70),
  ) +
  facet_wrap(~family_name)

### Combine all and family-specific ----
sizeDeets <- list(
  scale_size_continuous(
    "Number of records",
    breaks = c(1,10,40,70),
    limits = c(1,70), 
    range = c(1,6)
  )
)

p_all <- ggplot() +
  geom_sf(world, mapping = aes(), fill = "grey95", color = "grey60", linewidth = 0.1) +
  geom_sf(df_grid_count_grid, mapping = aes(), color = "grey50", fill = NA) +
  geom_sf(df_grid_count, mapping = aes(size = n), alpha = 0.5, color = "green") +
  sizeDeets

p_family <- lapply(c("Vespertilionidae", "Molossidae", "Phyllostomidae"), function(x) {
  sf1 <- df_grid_family[df_grid_family$family_name == x, ]
  sf2 <- df_grid_count_family[df_grid_count_family$family_name == x, ]
  out <- ggplot() +
    geom_sf(world, mapping = aes(), fill = "grey95", color = "grey60", linewidth = 0.1) +
    geom_sf(sf1, mapping = aes(), color = "grey50", fill = NA) +
    geom_sf(sf2, mapping = aes(size = n), alpha = 0.5, color = "green") +
    sizeDeets +
    facet_wrap(~family_name) +
    theme(
      strip.background = element_blank(),
      legend.position = "none"
    )
  return(out)
})


p_all + guide_area() +
  {p_family[[1]] + p_family[[2]] + p_family[[3]] +
      plot_layout(guides = "collect") &
      theme(legend.position = "none")} +
  plot_layout(
    design = "
    12
    33
    ",
    heights = c(3,2),
    widths = c(10,1),
    ) &
  theme(legend.position = "left")


# Timing of observations by continent --------

df_grid %>% 
  ggplot() +
  geom_histogram(aes(x = yday), binwidth = 7) +
  facet_wrap(~continent) +
  theme_classic()

df_grid %>% 
  mutate(
    hemisphere = case_when(
      continent == "North America" ~ continent,
      continent %in% c("Africa", "Oceania", "South America") ~ "Southern Hemisphere",
      TRUE ~ "Northern Hemisphere"
    )
  ) %>% 
  ggplot() +
  geom_histogram(aes(x = yday), binwidth = 7) +
  facet_wrap(~hemisphere) +
  theme_classic()




## Quick checks on geography in North America ----
df_NoAm <- df_grid %>% 
  dplyr::filter(continent == "North America") %>%
  as.data.frame() %>% 
  left_join(., read.csv("data/iNat_observations_tidy_manualChecks.csv")) %>% 
  mutate(
    lat_cut = cut(latitude, breaks = seq(0,90, by = 3)),
    lon_cut = cut(longitude, breaks = seq(-180,180, by = 3))
  )

# ggplot(df_NoAm) +
#   geom_histogram(aes(x = yday)) +
#   facet_grid(lat_cut~lon_cut) +
#   theme_classic()
# 
# ggplot(df_NoAm) +
#   geom_bar(aes(x = yday)) +
#   facet_wrap(~lon_cut) +
#   theme_classic()
# 
# ggplot(df_NoAm) +
#   geom_bar(aes(x = yday)) +
#   facet_wrap(~lat_cut) +
#   theme_classic()


## Visualize species in North America -----
datebreaks <- yday(lubridate::mdy(paste(1:12, "-1-2025")))
datelabs <-  month.abb

df_NoAm %>% 
  ggplot() +
  geom_histogram(aes(x = yday, fill = scientific_name)) +
  facet_wrap(~quality_grade+scientific_name) +
  theme_classic()

df_NoAm %>% 
  dplyr::filter(
    quality_grade == "research") %>% 
  ggplot() +
  geom_histogram(aes(x = yday, fill = scientific_name)) +
  scale_fill_manual(
    "Species",
    values = 
      c(
        "Eptesicus fuscus" = "#A16928",
        "Nycticeius humeralis" = "#124559", 
        "Vespertilionidae" = "grey50", 
        "Lasiurus borealis" = "#E46844", 
        "Perimyotis subflavus" = "#f4d35e", 
        "Lasionycteris noctivagans" = "#484357")
  ) +
  scale_x_continuous(
    "Day of year",
    breaks = datebreaks,
    labels = datelabs
  ) +
  facet_wrap(~scientific_name, axes = "all_x") +
  theme_classic()

# iNaturalist_NoAm_timing <- df_NoAm %>% 
#   mutate(
#     group = case_when(
#       scientific_name %in% c("Lasiurus borealis", "Eptesicus fuscus", "Lasionycteris noctivagans") ~ scientific_name,
#       TRUE ~ "Other"
#     )
#   ) %>% 
#   dplyr::filter(
#     quality_grade == "research",
#     # scientific_name %in% c("Lasiurus borealis", "Nycticeius humeralis", "Eptesicus fuscus", "Lasionycteris noctivagans")
#   ) %>% 
#   arrange(group, scientific_name) %>% 
#   mutate(
#     # scientific_name = factor(scientific_name, levels = ),
#     group = factor(group, levels = c("Lasiurus borealis", "Lasionycteris noctivagans", "Eptesicus fuscus", "Other"))
#   ) %>% 
#   ggplot() +
#   geom_histogram(aes(x = yday, fill = scientific_name), binwidth = 14) +
#   scale_fill_manual(
#     "Species",
#     values = 
#       c(
#         "Eptesicus fuscus" = "#A16928",
#         "Nycticeius humeralis" = "#124559", 
#         "Vespertilionidae" = "grey50", 
#         "Lasiurus borealis" = "#E46844", 
#         "Perimyotis subflavus" = "#f4d35e", 
#         "Lasionycteris noctivagans" = "#484357")
#   ) +
#   scale_x_continuous(
#     "Day of year",
#     breaks = datebreaks,
#     labels = datelabs
#   ) +
#   facet_wrap(~group, axes = "all_x", scales = "free_y") +
#   theme_classic() +
#   scale_y_continuous(
#     "Count of records on iNaturalist",
#     expand = expansion(add = c(0, 2), mult = c(0, .2)),
#     breaks = seq(0,50, by = 5),
#     minor_breaks = seq(0,50, by = 1)
#   ) +
#   theme(
#     panel.grid.major.y = element_line(color = "grey80"),
#     panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.1),
#     legend.position = "none",
#     strip.background = element_blank()
#   )
# ggsave(iNaturalist_NoAm_timing, "figs/iNaturalist_NoAm_timing.png")

### Do it as a list ----


dddd <- df_NoAm %>%
  mutate(
    species_rg = case_when(
      species_name != "" & quality_grade == "research" ~ species_name,
      species_name == "Aeorestes cinereus" & quality_grade == "research" ~ "Lasiurus cinereus",
      quality_grade != "research"  ~ "Other/unknown",
      TRUE ~ "Other/unknown"
    ),
    group = case_when(
      species_rg %in% c(
        "Lasiurus borealis",
        "Eptesicus fuscus",
        "Lasionycteris noctivagans"
      ) ~ species_rg,
      TRUE ~ "Other/unknown"
    ),
    group_lab = case_when(
      group == "Other/unknown" ~ "Other/unknown",
      TRUE ~ paste0("<i>", species_rg, "</i>")
    )
  )

iNaturalist_NoAm_timing_list <- lapply(c(
  "Lasiurus borealis",
  "Lasionycteris noctivagans",
  "Eptesicus fuscus",
  "Other/unknown" ), function(x) {
  
    if(x %in% c("Eptesicus fuscus", "Other/unknown")) {
      ymax <- 5
    } else {
      ymax <- 15
    }
    
    dddd %>% 
      dplyr::filter(group == x) %>% 
      mutate(species_rg = factor(species_rg)) %>% 
    ggplot() +
      geom_histogram(aes(x = yday, fill = species_rg), binwidth = 14) +
      scale_fill_manual(
        "Species",
        values = 
          c(
            "Eptesicus fuscus" = "#A16928",
            "Nycticeius humeralis" = "#124559", 
            "Vespertilionidae" = "grey50", 
            "Lasiurus borealis" = "#E46844", 
            "Perimyotis subflavus" = "#f4d35e", 
            "Lasionycteris noctivagans" = "#484357")
      ) +
      scale_x_continuous(
        "Day of year",
        breaks = datebreaks,
        labels = datelabs
      ) +
      facet_wrap(~group_lab, axes = "all_x", scales = "free_y") +
      theme_classic() +
      scale_y_continuous(
        "Number of North American iNaturalist records",
        expand = c(0,0),
        breaks = seq(0,50, by = 5),
        minor_breaks = seq(0,50, by = 1),
        limits = c(0, ymax)
      ) +
      theme(
        strip.text = ggtext::element_markdown(),
        panel.grid.major.y = element_line(color = "grey80"), linewidth = 0.3,
        panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.15),
        legend.position = "none",
        strip.background = element_blank()
      )
})

iNaturalist_NoAm_timing2 <- patchwork::wrap_plots(iNaturalist_NoAm_timing_list, axis_titles = "collect")

iNaturalist_NoAm_timing2


# 
# # Summarize taxonomy ----------
# 
# 
# df_NoAm %>% 
#   ggplot() +
#   geom_bar(aes(y = 1, fill = species_name)) +
#   facet_wrap(~continent) +
#   scale_y_continuous(limits = c(0, NA)) +
#   coord_polar() 
# 
df_NoAm_working <- df_NoAm %>%
  mutate(
    family_name = case_when(family_name == "" ~ as.character(NA), TRUE ~ family_name),
    genus_name = case_when(genus_name == "" ~ as.character(NA), TRUE ~ genus_name),
    species_name = case_when(species_name == "" ~ as.character(NA), TRUE ~ species_name),
    fam_color = case_when(
      is.na(family_name) ~ "grey50",
      family_name == "Molossidae" ~ "orange",
      family_name == "Phyllostomidae" ~ "darkgreen",
      family_name == "Vespertilionidae" ~ "skyblue"
    )
  )

# fam.count <- df_NoAm_working %>% 
#   count(family_name) %>% 
#   mutate_each(funs(empty_as_na))  %>% 
#   mutate(
#     xmax = cumsum(n),
#     xmin = lag(xmax, 1),
#     xmin = case_when(is.na(xmin) ~ 0, TRUE ~ xmin),
#     xmid = ((xmax - xmin) / 2) + xmin,
#     ymin = 1, ymax = 2
#   ) 
# 
# gen.count <- df_NoAm_working %>% 
#   count(family_name, genus_name) %>% 
#   mutate_each(funs(empty_as_na))  %>% 
#   mutate(
#     xmax = cumsum(n),
#     xmin = lag(xmax, 1),
#     xmin = case_when(is.na(xmin) ~ 0, TRUE ~ xmin),
#     xmid = ((xmax - xmin) / 2) + xmin,
#     ymin = 2, ymax = 3
#   ) 
# 
# spp.count <- df_NoAm_working %>% 
#   count(family_name, genus_name, species_name) %>% 
#   mutate_each(funs(empty_as_na))  %>% 
#   mutate(
#     xmax = cumsum(n),
#     xmin = lag(xmax, 1),
#     xmin = case_when(is.na(xmin) ~ 0, TRUE ~ xmin),
#     xmid = ((xmax - xmin) / 2) + xmin,
#     ymin = 3, ymax = 4
#   ) 
# 
# ggplot() +
#   geom_rect(data = fam.count, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = family_name)) +
#   scale_fill_viridis_d(option = "mako", na.value = "grey") +
#   geomtextpath::geom_textpath(data = fam.count, aes(x = xmid, y = ymin + 0.5 , label = family_name) ) +
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = gen.count, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = genus_name)) +
#   scale_fill_viridis_d(option = "plasma", na.value = "grey") +
#   geomtextpath::geom_textpath(data = gen.count, aes(x = xmid, y = ymin + 0.5 , label = genus_name) ) +
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = spp.count, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = species_name)) +
#   scale_fill_viridis_d(option = "turbo", na.value = "grey") +
#   # geomtextpath::geom_textpath(data = spp.count, aes(x = xmid, y = ymin + 0.5 , label = species_name) ) +
#   # ggrepel::geom_text_repel(
#   #   data = spp.count, aes(x = xmin, y = ymin + 1 , label = species_name),
#   #   seed = 42,
#   #   nudge_x = 2,
#   #   nudge_y = 3,
#   #   min.segment.length = 0
#   #   ) +
#   geom_text(
#     data = spp.count, aes(x = xmin, y = ymin + 1 , label = species_name),
#     angle = 45, hjust = 0
#   ) +
#   coord_polar() +
#   scale_y_continuous(limits = c(0, 5)) +
#   theme_void() +
#   theme(
#     # legend.position = "none"
#   )
# ggsave("figs/iNat_count_records.png", width = 12, height = 12)
# 
# 
ddd <- df_NoAm_working %>%
  mutate(
    species_rg = case_when(
      species_name != "" & quality_grade == "research" ~ species_name,
      species_name == "Aeorestes cinereus" & quality_grade == "research" ~ "Lasiurus cinereus",
      quality_grade != "research"  ~ "Other/unknown",
      TRUE ~ "Other/unknown"
    ),
    species_rg = case_when(
      species_rg == "Aeorestes cinereus" ~ "Lasiurus cinereus",
      TRUE ~ species_rg
    ),
    species_rg_lab = 
      case_when(
        species_rg ==  "Other/unknown" ~  "Other/unknown",
        TRUE ~ paste0("<i>", species_rg, "</i>")
      )
  ) %>%
  count(species_rg, species_rg_lab) %>%
  arrange(desc(n))

iNaturalist_record_by_spp <- ddd %>%
  mutate(species_rg_lab = factor(species_rg_lab, levels = unique(.$species_rg_lab))) %>%
  ggplot() +
  geom_col(aes(x = species_rg_lab, y = n, fill = species_rg)) +
  geom_text(aes(x = species_rg_lab, y = n, label = n), vjust = -0.25, color = "grey50", size = 2) +
  scale_x_discrete("Species",) +
  scale_y_continuous(
    "Number of North American iNaturalist records",
    expand = expansion(add = c(0, 3)),
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100)
    ) +
  scale_fill_manual(
    "Species",
    values = c(
      "Eptesicus fuscus"         = speciesColors[[1]],
      "Nycticeius humeralis"     = speciesColors[[2]],
      "Lasiurus borealis"        = speciesColors[[4]],
      "Perimyotis subflavus"     = speciesColors[[5]],
      "Lasionycteris noctivagans"= speciesColors[[6]]

      # "Eptesicus fuscus"         = "#4662D7",
      # "Nycticeius humeralis"     = "#1AE4B6",
      # "Lasiurus borealis"        = "#E14209",
      # "Perimyotis subflavus"     = "#FABA39",
      # "Lasionycteris noctivagans"= "#30123B"
      )
  )  +
  theme(
    legend.position = "none",
    axis.text.x = ggtext::element_markdown(angle = 55, hjust = 1),
    axis.title.x.bottom = element_blank()
  )
ggsave(iNaturalist_record_by_spp, filename = "figs/iNaturalist_record_by_spp.png")



# Figure 4 -----------

f4_combo <- iNat_map_grid_circle  /
  iNaturalist_record_by_spp/
  iNaturalist_NoAm_timing2 +
  plot_layout(axis_titles = "collect", heights = c(3,1,2)) 
ggsave(f4_combo, filename = "figs/f4.png", width = 8, height = 10)






# Building height analyses ----------------------------------------------------
proj.aea <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs"

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

# Load NoAm
NoAm <- rnaturalearth::countries110 %>% 
  st_as_sf() %>% 
  dplyr::filter(CONTINENT == "North America") 

# Crop to NoAm
if(!file.exists("tmp/bh_NoAm.tif")) {
  bh_NoAm <- NoAm %>% 
    st_transform(st_crs(bh)) %>% 
    st_bbox() %>% 
    st_as_sfc() %>% 
    terra::crop(bh, .) %>% 
    project(proj.aea)
  writeRaster(bh_NoAm, file = "tmp/bh_NoAm.tif")
}
bh_NoAm <- rast("tmp/bh_NoAm.tif")


p_bh_usa <- ggplot() +
  geom_spatraster(bh_NoAm, mapping = aes(), maxcell = 100e5) +
  scale_fill_viridis_c(option = "turbo", na.value = "white") +
  geom_sf(NoAm, mapping = aes(), fill = NA) +
  geom_sf(usa, mapping = aes(), fill = NA) +
  coord_sf(
    xlim = c(-2236891, 2127179),
    ylim = c(-1694616, 1328895),
    crs = proj.aea
  )
ggsave("figs/p_bh_usa.png", p_bh_usa, dpi = 600, width = 12, height = 12)




df_NoAm_precise <- dplyr::filter(df, continent == "North America", geoprivacy != "obscured") %>% 
  st_transform(proj.aea)
ggplot() +
  geom_spatraster(bh_NoAm, mapping = aes()) +
  geom_sf(df_NoAm_precise, mapping = aes(), color = "red", shape = 21) +
  geom_sf(NoAm, mapping = aes(), fill = NA)

df_NoAm_precise_bh <- terra::extract(bh_NoAm, df_NoAm_precise)
ggplot() +
  geom_density(data = df_NoAm_precise_bh, aes(x = GBH2020_150m_GEDI))


sampleBoundary <- df_NoAm_precise %>% 
  st_buffer(dist = 5e3) %>% 
  st_union() %>% 
  st_transform(proj.aea) %>% 
  as_spatvector()

ggplot() +
  geom_sf(df_NoAm_precise, mapping = aes(), color = "red") +
  geom_spatvector(sampleBoundary, mapping = aes(), fill = NA)


masked_bh_NoAm <- mask(bh_NoAm, sampleBoundary)
ggplot() +
  geom_spatraster(masked_bh_NoAm, mapping = aes()) +
  geom_sf(df_NoAm_precise, mapping = aes(), color = "red", shape = 21) +
  geom_sf(NoAm, mapping = aes(), fill = NA)
# ggplot() +
#   geom_spatraster(masked_bh_NoAm, mapping = aes(), maxcell = 5e5) +
#   geom_sf(df_NoAm_precise_bh, mapping = aes(), color = "red") +
#   scale_fill_viridis_c(option = "turbo", na.value = "grey50") +
#   geom_sf(NoAm, mapping = aes(), fill = NA) +
#   geom_sf(usa, mapping = aes(), fill = NA)


ranSample <- spatSample(x = masked_bh_NoAm, size = 10000)
ranSample$GBH2020_150m_GEDI[is.na(ranSample$GBH2020_150m_GEDI)] <- 2
ggplot() +
  geom_density(data = df_NoAm_precise_bh, aes(x = GBH2020_150m_GEDI)) +
  geom_density(data = ranSample, aes(x = GBH2020_150m_GEDI))




fp <- df_NoAm_precise %>% 
  st_union %>% 
  st_transform(st_crs(bh))
fp_buffer <- st_buffer(fp, dist = 5e3)
fp_buffer_vect <- as_spatvector(fp_buffer)
fp_crop_bh <- crop(bh, ext(fp_buffer_vect))
fp_crop_bh_NA <- subst(fp_crop_bh, NA, 3)
fp_mask_bh <- mask(fp_crop_bh_NA, fp_buffer_vect)


sample_at_pts <- terra::extract(fp_mask_bh, df_NoAm_precise)
sample_random <- terra::spatSample(fp_mask_bh, size = 10000, replace = F, na.rm = T)



ggplot() +
  geom_density(sample_at_pts, mapping = aes(x = GBH2020_150m_GEDI)) +
  geom_density(sample_random, mapping = aes(x = GBH2020_150m_GEDI), color = "green")

  

