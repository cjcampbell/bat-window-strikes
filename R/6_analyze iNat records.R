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
# source("R/1_setup and load KC data.R")
# source("R/2_plot KC records.R")

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
  st_as_sf(coords = c("longitude", "latitude"), crs = proj.wgs84, remove = F) %>% 
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


# Create table with observation IDs and license info for retained obs. ------
# This is intended for sharing without violating license agreements.
df_with_info <- df  %>% 
  as.data.frame() %>% 
  dplyr::select(
    id, url, notes = CJ.notes, scientific_name, common_name, description, user_login, observed_on, license, family_name, adm0_a3, continent, 
  )

df_restrictiveLicenses <- df_with_info %>% 
  dplyr::filter(license == "") %>% 
  dplyr::select(id, url, notes, family_name, adm0_a3, continent) %>% 
  mutate(
    notes = paste0(notes, ". All rights reserved license"),
    notes = case_when(notes == ". All rights reserved license" ~ "All rights reserved license", TRUE ~ notes)
  )

df_otherLicenses <-  df_with_info %>% 
  dplyr::filter(license != "") 

df_tidy_licenses <- full_join(df_otherLicenses, df_restrictiveLicenses)
write.csv(df_tidy_licenses, "out/data_derived/iNaturalist records.csv", row.names = F)


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

## Prepare to map ----
df_grid <- st_join(df, mygrid)

df_grid_count_grid <- df_grid %>% 
  count(grid_id) %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry) %>% 
  inner_join(mygrid) %>% 
  st_as_sf()

df_grid_count <- df_grid_count_grid %>% 
  st_centroid()

## Exploratory maps -----
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

## Figure 4a: Main map of all records ----
# Plot circle in grid.
iNat_map_grid_circle <- ggplot() +
  geom_sf(world, mapping = aes(), fill = "grey95", color = "grey60", linewidth = 0.1) +
  geom_sf(df_grid_count_grid, mapping = aes(), color = "grey50", fill = NA) +
  geom_sf(df_grid_count, mapping = aes(size = n), alpha = 0.5, color = "#00B31B") +
  scale_size_continuous(
    "Number of records",
    breaks = c(1,10,40,70),
  ) +
  scale_y_continuous(expand = c(0,0)) +
 # scale_x_continuous(expand = c(0,0)) +
  coord_sf(ylim = c(-6792.374, 8000)) +
  theme_void() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25)
  )
ggsave(iNat_map_grid_circle, filename =  "figs/iNat_map_grid_circle.png", dpi = 600)

## Map by family -----
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
    heights = c(3,3),
    widths = c(10,1),
  ) 
ggsave("figs/iNaturalist_records_map_by_family.png", dpi = 600, width = 10, height = 8)

# Timing of observations by continent --------

df_grid %>% 
  ggplot() +
  geom_histogram(aes(x = yday), binwidth = 7) +
  scale_y_continuous(
    "Count of iNaturalist records",
    expand = expansion(mult = c(0, NA))
  ) +
  scale_x_continuous(
    "Day of year",
    breaks = c(seq(0,350, by = 50))
  ) +
  facet_wrap(~continent, axes = "all") +
  theme_classic() +
  theme(
    strip.background = element_blank()
  )
ggsave("figs/SI_iNatRecords_by_yday.png", dpi = 600, width = 8)



df_grid %>% 
  ggplot() +
  geom_histogram(aes(x = yday, fill = quality_grade)) +
  facet_wrap(~continent +scientific_name, axes = "all", scales = "free_y") +
  scale_y_continuous(
    "Count of iNaturalist records",
    expand = expansion(add = c(0, 4))
  ) +
  scale_fill_viridis_d(
    "iNaturalist quality grade",
    option = "mako", direction = -1, end = 0.8, begin = 0.1
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  )


df_grid %>% 
  count(continent, scientific_name) %>% 
  arrange(desc(n)) %>% 
  mutate(scientific_name=factor(scientific_name, levels = unique(.$scientific_name))) %>% 
  ggplot() +
  geom_col(aes(x = scientific_name, y = n)) +
  facet_wrap(~continent, axes = "all_x", scales = "free", space = "free_x") +
  scale_y_continuous(
    "Count",
    expand = expansion(add = c(0, 5))
  ) +
  scale_x_discrete(
    "Scientific name"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x = ggtext::element_markdown(angle = 55, hjust = 1),
  )
ggsave("figs/iNaturalist_records_by_continent_sciName.png", width = 8, dpi = 600)

# Analysis for North America -----
df_NoAm <- dplyr::filter(df_grid, continent == "North America")

datebreaks <- yday(lubridate::mdy(paste(1:12, "-1-2025")))
datelabs <-  month.abb

# df_NoAm %>% 
#   dplyr::filter(
#     quality_grade == "research") %>% 
#   ggplot() +
#   geom_histogram(aes(x = yday, fill = scientific_name)) +
#   scale_fill_manual(
#     "Species",
#     values = speciesColors2
#   ) +
#   scale_x_continuous(
#     "Day of year",
#     breaks = datebreaks,
#     labels = datelabs
#   ) +
#   facet_wrap(~scientific_name, axes = "all_x") +
#   theme_classic()

### Plot panels within a list ----

# Reclassify to match taxonomic authority.
# Classify non-research grade observations as "other/unknown".

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
    
    mySpeciesColors <- c( 
      "Eptesicus fuscus" = "#A8541F", 
      "Nycticeius humeralis" = "grey50",
      "Other/unknown" = "grey50",
      "Lasiurus borealis" = "#9E2A2B",
      "Perimyotis subflavus" = "grey50",
      "Lasionycteris noctivagans" = "#7072A0"
    )
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
        values = mySpeciesColors
      ) +
      scale_x_continuous(
        "Day of year",
        limits = c(1, 365),
        breaks = datebreaks,
        labels = datelabs,
        expand = c(0,0)
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
        axis.text.x = ggtext::element_markdown(angle = 55, hjust = 1),
        strip.text = ggtext::element_markdown(),
        panel.grid.major.y = element_line(color = "grey80"), linewidth = 0.3,
        panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.15),
        legend.position = "none",
        strip.background = element_blank()
      )
  })

iNaturalist_NoAm_timing2 <- patchwork::wrap_plots(iNaturalist_NoAm_timing_list, axis_titles = "collect")

iNaturalist_NoAm_timing2


#### Plot without colors -----
# 
# 
# iNaturalist_NoAm_timing_list_noCols <- lapply(c(
#   "Lasiurus borealis",
#   "Lasionycteris noctivagans",
#   "Eptesicus fuscus",
#   "Other/unknown" ), function(x) {
#     
#     if(x %in% c("Eptesicus fuscus", "Other/unknown")) {
#       ymax <- 5
#     } else {
#       ymax <- 15
#     }
#     
#     focalCount <- dddd %>% 
#       count(group) %>% 
#       dplyr::filter(group == x)
#     
#     dddd2 <- dddd %>%
#       dplyr::mutate(
#         group_lab = case_when(group_lab == "Other/unknown" ~ "Other", TRUE ~ group_lab),
#         group_lab = paste("<b>", group_lab, "</b><br>", "<i>n</i>=", focalCount$n)
#       )
#     
#     dddd2 %>% 
#       dplyr::filter(group == x) %>% 
#       mutate(species_rg = factor(species_rg)) %>% 
#       ggplot() +
#       geom_histogram(aes(x = yday), binwidth = 14, fill = "#337CA0") +
#       scale_x_continuous(
#         "Day of year",
#         limits = c(1, 365),
#         breaks = datebreaks,
#         labels = datelabs,
#         expand = c(0,0)
#       ) +
#       facet_wrap(~group_lab, axes = "all_x", scales = "free_y") +
#       theme_classic() +
#       scale_y_continuous(
#         "Number of North American iNaturalist records",
#         expand = c(0,0),
#         breaks = seq(0,50, by = 5),
#         minor_breaks = seq(0,50, by = 1),
#         limits = c(0, ymax)
#       ) +
#       theme(
#         strip.text = ggtext::element_markdown(),
#         panel.grid.major.y = element_line(color = "grey80"), linewidth = 0.3,
#         panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.15),
#         legend.position = "none",
#         strip.background = element_blank()
#       )
#   })
# 
# iNaturalist_NoAm_timing_list_noCols3 <- patchwork::wrap_plots(iNaturalist_NoAm_timing_list_noCols, axis_titles = "collect")
# 
# iNaturalist_NoAm_timing_list_noCols3
# # ggsave("figs/iNaturalist_timing.png",iNaturalist_NoAm_timing_list_noCols3, dpi = 600, width = 8)
# 
# 

## Summarize taxonomy ----------

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


ddd <- df_NoAm_working %>%
  mutate(
    species_rg = case_when(
      species_name != "" & quality_grade == "research" ~ species_name,
      species_name == "Aeorestes cinereus" & quality_grade == "research" ~ "Lasiurus cinereus",
      quality_grade != "research" ~ "Other/unknown",
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
  # scale_fill_manual(
  #   "Species",
  #   values = speciesColors2
  # )  +
  theme(
    legend.position = "none",
    axis.text.x = ggtext::element_markdown(angle = 55, hjust = 1),
    axis.title.x.bottom = element_blank()
  )
ggsave(iNaturalist_record_by_spp, filename = "figs/iNaturalist_record_by_spp.png")


ddd_bygroup <- ddd %>%
  as.data.frame %>% 
  mutate(
    group = case_when(
      species_rg %in% c(
        "Lasiurus borealis",
        "Lasionycteris noctivagans",
        "Eptesicus fuscus"
      ) ~ species_rg_lab,
      TRUE ~ "Other"
      ),
    species = case_when(group == "Other" ~ group, TRUE ~ species_rg)
    ) %>% 
  dplyr::summarise(n = sum(n), .by = c(`group`, `species`)) %>% 
  mutate(
    group = factor(group, levels = rev(c(
      "<i>Lasiurus borealis</i>",
      "<i>Lasionycteris noctivagans</i>",
      "<i>Eptesicus fuscus</i>",
      "Other"
    )))
  )

ddd_Others <- ddd %>%
  dplyr::filter(
    !species_rg %in% c(
      "Lasiurus borealis",
      "Lasionycteris noctivagans",
      "Eptesicus fuscus"
  )) %>% 
  mutate(
    species_rg_lab = case_when(species_rg_lab == "Other/unknown" ~ "Unknown", TRUE ~ species_rg_lab),
    species_rg_lab2 = paste0(species_rg_lab, ": ", n),
    rownum = 9 - row_number()
  )


(p_barplot <- ggplot(ddd_bygroup) +
  geom_col(aes(y = group, x = n, fill = species)) +
  geom_text(aes(y = group, x = n, label = n), vjust = 2, hjust = 1.1, color = "white", size = 3) +
  ggtext::geom_richtext(
    aes(y = group, x = 0, label = group),
    vjust = 0.5,
    hjust = 0,
    color = "white",
    size = 3,
    fill = NA, label.colour = NA
  ) +
  ggtext::geom_richtext(
    data = ddd_Others,
    aes(x = 65, y = 0.75 + (rownum / 4), label = species_rg_lab2),
    text.color = "black",
    hjust = 0,
    fill = NA, label.colour = NA
  ) +
  scale_y_discrete("Species") +
  scale_x_continuous(
    "Number of North American iNaturalist records",
    expand = expansion(add = c(0, 0)),
    breaks = seq(0, 100, by = 25),
    limits = c(0, 95)
  ) +
  scale_fill_manual(
    "Species",
    values = c(
      "Other" = "grey50", 
      "Eptesicus fuscus" = "#A8541F", 
      "Lasionycteris noctivagans" = "#7072A0", 
      "Lasiurus borealis" = "#9E2A2B"
    )
  )  +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )
)

# Time of day exploration ----------

df_time <- df_grid %>% 
  separate(datetime, into = c("date", "time", "tz"), remove = F, sep = " ") %>% 
  dplyr::mutate(
    hour = as.numeric(word(time, 1,1,sep=":")),
    min  = as.numeric(word(time, 2,2,sep=":")),
    sec  = as.numeric(word(time, 3,3,sep=":"))
  ) %>% 
  # remove 00:00:00 times.
  dplyr::filter(
    hour != 0 | min != 0 | sec != 0
  )

df_time %>% 
  count(continent, hour >= 5 & hour <= 10)

(p_timeOfDay <- df_time %>% 
  ggplot() +
  geom_bar(aes(x= hour, fill = continent)) +
  # scale_fill_viridis_d(option = "turbo", begin = 0, end = 0.5) +
    scale_fill_manual(
      "Continent",
      values = rev(c("#001219", "#0a9396", "#94d2bd", "#ABB86D", "#ee9b00",  "#9b2226"))
    ) +
  scale_y_continuous("Count of records", expand = expansion(mult=c(0,NA))) +
    scale_x_continuous("Hour of observation")
)

# Figure 4 combo ------
F4_combo <- iNat_map_grid_circle + p_timeOfDay + 
  plot_layout(
  design = "
    1
    2
    ",
  heights = c(2,1),
  axis_titles = "collect_y"
) +
  plot_annotation(
    tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & 
  theme(
    plot.tag.position  = c(0.1,0.95),
    plot.tag = element_text(size = 10)
  )

ggsave(F4_combo, filename =  "figs/F4_combo.png", dpi = 600, width = 6.5, height = 5)
ggsave(F4_combo, filename =  "figs/F4_combo.svg", width = 6.5, height = 5)


# Building height analyses ----------------------------------------------------

# Load building height data (big!).
bh <- rast("data/building height data/GBH2020_150m_GEDI.tif")

# Conduct analyses in native moll projection

# Download GADM data.
NoAm <- lapply(c("USA", "CAN", "MEX"), function(x) {
  geodata::gadm(x, level = 1,
      path = "../../- Missions & Programs/Research & Development/Data Products"
    ) 
}) %>% 
  vect %>% 
  project(terra::crs(bh))

# Interpolate 0's and then mask by expected land area.
if(!file.exists("tmp/fp_mask_bh.tif")) {
  crop(bh, NoAm, filename = "tmp/bh_NoAm0.tif")
  bh_NoAm0 <- rast("tmp/bh_NoAm0.tif")
  fp_crop_bh_NA <- subst(bh_NoAm0, NA, 0)
  fp_mask_bh <- mask(fp_crop_bh_NA, NoAm, filename = "tmp/fp_mask_bh.tif")
} else {
  fp_mask_bh <- raster("tmp/fp_mask_bh.tif")
}

if(!file.exists("tmp/built_proj.tif")) {
  geodata::landcover(var = "built") %>% 
    project(crs(bh), filename = "tmp/built_proj.tif")
} else {
  built_proj <- rast("tmp/built_proj.tif")
}



# Select example location and visualize.
collisionPoint <- df_NoAm_precise[df_NoAm_precise$id == 12539414,] 
collisionPoint_buff <- st_buffer(collisionPoint, dist = 1e3)
myOutline <- st_buffer(collisionPoint, dist = 3e3)

bh_cropped2Example <- crop(bh_NoAm0, myOutline)
bh_cropped2Example2 <- crop(fp_mask_bh, myOutline)
built_proj_resampled <- terra::resample(built_proj, bh_cropped2Example, method = "max") %>% 
  mask(., bh_cropped2Example)
# built_proj_cropped2Example <- crop(built_proj, collisionPoint_buff)

# ggplot() + 
#   geom_spatraster(built_proj_resampled, mapping = aes(), maxcell = 50e5) +
#   geom_spatvector(collisionPoint, mapping = aes(), color = "red") +
#   geom_spatvector(collisionPoint_buff, mapping = aes(), color = "red", fill = NA)

# ggplot() + 
#   geom_spatraster(bh_cropped2Example, mapping = aes(), maxcell = 50e5) +
#   geom_spatvector(collisionPoint, mapping = aes(), color = "red") +
#   geom_spatvector(collisionPoint_buff, mapping = aes(), color = "red", fill = NA)

ggplot() + 
  geom_spatraster(bh_cropped2Example2, mapping = aes(), maxcell = 50e5) +
  geom_spatvector(collisionPoint, mapping = aes(), color = "red") +
  geom_spatvector(collisionPoint_buff, mapping = aes(), color = "red", fill = NA)


# Were bats sampled at locations with greater building heights than
# the areas around them?









# # Crop to NoAm
# if(!file.exists("tmp/bh_NoAm.tif")) {
#   bh_NoAm <- NoAm %>% 
#     st_transform(st_crs(bh)) %>% 
#     st_bbox() %>% 
#     st_as_sfc() %>% 
#     terra::crop(bh, .) %>% 
#     project(proj.aea)
#   writeRaster(bh_NoAm, file = "tmp/bh_NoAm.tif")
# }
# bh_NoAm <- rast("tmp/bh_NoAm.tif")
# 

## Plot building height across USA -----
# p_bh_usa <- ggplot() +
#   geom_spatraster(bh_NoAm, mapping = aes(), maxcell = 100e5) +
#   scale_fill_viridis_c(option = "turbo", na.value = "white") +
#   geom_sf(NoAm, mapping = aes(), fill = NA) +
#   geom_sf(usa, mapping = aes(), fill = NA) +
#   coord_sf(
#     xlim = c(-2236891, 2127179),
#     ylim = c(-1694616, 1328895),
#     crs = proj.aea
#   )
# ggsave("figs/p_bh_usa.png", p_bh_usa, dpi = 600, width = 12, height = 12)


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








ranSample <- spatSample(x = masked_bh_NoAm, size = 10000)
ranSample$GBH2020_150m_GEDI[is.na(ranSample$GBH2020_150m_GEDI)] <- 2
ggplot() +
  geom_density(data = df_NoAm_precise_bh, aes(x = GBH2020_150m_GEDI)) +
  geom_density(data = ranSample, aes(x = GBH2020_150m_GEDI))


# Sample points randomly from within 5km of each collision.
fp <- df_NoAm_precise %>% 
  st_union %>% 
  st_transform(st_crs(bh))
fp_buffer <- st_buffer(fp, dist = 5e3)
fp_buffer_vect <- as_spatvector(fp_buffer)
fp_crop_bh <- crop(bh, ext(fp_buffer_vect))
fp_crop_bh_NA <- subst(fp_crop_bh, NA, 3)
fp_mask_bh <- mask(fp_crop_bh_NA, fp_buffer_vect)

# Extract building height at site of collision, and randomly from within 5km buffer.
sample_at_pts <- terra::extract(fp_mask_bh, df_NoAm_precise)
sample_random <- terra::spatSample(fp_mask_bh, size = 10000, replace = F, na.rm = T)


ggplot() +
  geom_density(sample_at_pts, mapping = aes(x = GBH2020_150m_GEDI)) +
  geom_density(sample_random, mapping = aes(x = GBH2020_150m_GEDI), color = "blue")


data.frame(
  sample = "at points", sample_at_pts
) %>% 
  full_join(
    data.frame(
      sample = "near points", sample_random
    ) 
  ) %>% 
  ggplot() +
  aes(x= sample, y = GBH2020_150m_GEDI) +
  geom_violin() +
  geom_boxplot(fill = NA) +
  scale_x_discrete(
    "Location",
    labels = c("At site of collision\n(n=207)", "Near site of collision\n(n=XX)")
  ) +
  scale_y_continuous(
    "Building height (m)"
  )




# Figure 5 ------

(f5_combo2 <- iNaturalist_NoAm_timing_list[[1]] +
   iNaturalist_NoAm_timing_list[[2]] +
   iNaturalist_NoAm_timing_list[[3]] +
   iNaturalist_NoAm_timing_list[[4]] +
   free(p_barplot) +
   plot_layout(
     design = "
    123
    455
    ",
     axis_titles = "collect_y"
   ) +
   plot_annotation(
     tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & 
   theme(
     plot.tag.position  = c(0.1,0.95),
     plot.tag = element_text(size = 10)
   )
)
ggsave(f5_combo2, filename = "figs/f5_combo2.png", width = 10, height = 6, dpi = 600)
ggsave(f5_combo2, filename = "figs/f5_combo2.svg", width = 9, height = 6)


# Tables + summaries ------
df %>% 
  as.data.frame %>% 
  count(admin) %>% 
  arrange(desc(n))
df %>% 
  as.data.frame %>% 
  count(continent) 

df %>% 
  as.data.frame() %>% 
  dplyr::filter(quality_grade == "research",species_name!="") %>% 
  count(quality_grade, family_name, species_name) %>% 
  arrange(desc(n))

df %>% 
  as.data.frame() %>% 
  dplyr::filter(continent == "North America") %>% 
  mutate(
    lon_cut = cut(longitude, breaks = c(-Inf, -100, Inf)),
    lat_cut = cut(latitude, breaks = c(-Inf, 24, Inf))
  ) %>% 
  count(lon_cut, lat_cut)