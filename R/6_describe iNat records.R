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


taxTree <- fread("data/iNat_observations_taxTree.csv")
world <- ne_countries(scale = "medium", returnclass = "sf") |> 
  sf::st_transform(myproj) |> 
  dplyr::select(admin, adm0_a3, continent)
mygrid <- sf::st_make_grid(world, cellsize = 500) |>
  st_as_sf() |> 
  dplyr::mutate(grid_id = row_number())


# Load manually-checked data
df <- read.csv("data/iNat_observations_tidy_manualChecks.csv") |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = proj.wgs84, remove = F) |> 
  st_transform(myproj) |> 
  dplyr::filter(CJ.manual.check == "y") |> 
  mutate(
    datetime.POSIX = as.POSIXct(datetime),
    datetime.date = as_date(datetime.POSIX),
    yday = yday(datetime.date)
  ) |> 
  # Join with taxonomy data.
  left_join(taxTree, by = c("scientific_name" = "search")) |>
  # Unite with spatial data
  st_join(world)


# Create table with observation IDs and license info for retained obs. ------
# This is intended for sharing without violating license agreements.
df_with_info <- df  |> 
  as.data.frame() |> 
  dplyr::select(
    id, url, notes = CJ.notes, scientific_name, common_name, description, user_login, observed_on, license, family_name, adm0_a3, continent, 
  )

df_restrictiveLicenses <- df_with_info |> 
  dplyr::filter(license == "") |> 
  dplyr::select(id, url, notes, family_name, adm0_a3, continent) |> 
  mutate(
    notes = paste0(notes, ". All rights reserved license"),
    notes = case_when(notes == ". All rights reserved license" ~ "All rights reserved license", TRUE ~ notes)
  )

df_otherLicenses <-  df_with_info |> 
  dplyr::filter(license != "") 

df_tidy_licenses <- full_join(df_otherLicenses, df_restrictiveLicenses)
write.csv(df_tidy_licenses, "data/derived/iNaturalist records.csv", row.names = F)


## Correct  spatial data -----
df |> dplyr::filter(is.na(continent))
# A handful don't overlap countries (probably due to obscured coordinates near the coast).
# To resolve this, find the nearest country and replace NA's.
for(i in which(is.na(df$continent)) ) {
  
  nearest <- df[i,] |>
    st_distance(world) |>
    which.min()
  workingDF <- world[nearest,] |>
    as.data.frame() |>
    dplyr::select(admin, adm0_a3, continent)
  
  df$admin[i]     <- workingDF$admin
  df$adm0_a3[i]   <- workingDF$adm0_a3
  df$continent[i] <- workingDF$continent
}


# Map observations -------

## Prepare to map ----
df_grid <- st_join(df, mygrid)

df_grid_count_grid <- df_grid |> 
  count(grid_id) |> 
  as.data.frame() |> 
  dplyr::select(-geometry) |> 
  inner_join(mygrid) |> 
  st_as_sf()

df_grid_count <- df_grid_count_grid |> 
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
# Box marking the North American records carried into panel (c) and Figure 5. Built
# in lon/lat and segmentized before projecting, so its edges follow the Equal Earth
# graticule rather than cutting straight across it. The southern edge sits below the
# southernmost record (Panama, 8.8 degN), making explicit that Central America falls
# within our definition of North America.
noam_box <- st_bbox(c(xmin = -126, ymin = 5, xmax = -68, ymax = 54), crs = st_crs(proj.wgs84)) |>
  st_as_sfc() |>
  st_segmentize(units::set_units(100, "km")) |>
  st_transform(myproj)

# Plot circle in grid.
iNat_map_grid_circle <- ggplot() +
  geom_sf(world, mapping = aes(), fill = "grey95", color = "grey60", linewidth = 0.1) +
  geom_sf(df_grid_count_grid, mapping = aes(), color = "grey50", fill = NA) +
  geom_sf(df_grid_count, mapping = aes(size = n), alpha = 0.5, color = "#00B31B") +
  geom_sf(noam_box, mapping = aes(), fill = NA, color = "grey25", linewidth = 0.4,
          linetype = "22") +
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
df_grid_family <- df_grid |> 
  count(grid_id, family_name) |> 
  as.data.frame() |> 
  dplyr::select(-geometry) |> 
  inner_join(mygrid) |> 
  st_as_sf()

df_grid_count_family <- df_grid_family |> 
  st_centroid()

df_grid_count_family |> 
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

df_grid |> 
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



df_grid |> 
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


df_grid |> 
  count(continent, scientific_name) |> 
  arrange(desc(n)) |> 
  mutate(scientific_name=factor(scientific_name, levels = unique(scientific_name))) |>
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


### Plot panels within a list ----

# Reclassify to match taxonomic authority.
# Classify non-research grade observations as "other/unknown".

dddd <- df_NoAm |>
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
    
    dddd |> 
      dplyr::filter(group == x) |> 
      mutate(species_rg = factor(species_rg)) |> 
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


## Summarize taxonomy ----------

df_NoAm_working <- df_NoAm |>
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


ddd <- df_NoAm_working |>
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
  ) |>
  count(species_rg, species_rg_lab) |>
  arrange(desc(n))

iNaturalist_record_by_spp <- ddd |>
  mutate(species_rg_lab = factor(species_rg_lab, levels = unique(species_rg_lab))) |>
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
  theme(
    legend.position = "none",
    axis.text.x = ggtext::element_markdown(angle = 55, hjust = 1),
    axis.title.x.bottom = element_blank()
  )
ggsave(iNaturalist_record_by_spp, filename = "figs/iNaturalist_record_by_spp.png")


# ddd_bygroup / ddd_Others / p_barplot removed: they existed only to build the
# old Figure 5 species barplot, which the use-availability driver figure replaced.


# Time of day exploration ----------

df_time <- df_grid |> 
  separate(datetime, into = c("date", "time", "tz"), remove = F, sep = " ") |> 
  dplyr::mutate(
    hour = as.numeric(word(time, 1,1,sep=":")),
    min  = as.numeric(word(time, 2,2,sep=":")),
    sec  = as.numeric(word(time, 3,3,sep=":"))
  ) |> 
  # remove 00:00:00 times.
  dplyr::filter(
    hour != 0 | min != 0 | sec != 0
  )

df_time |> 
  count(continent, hour >= 5 & hour <= 10)

(p_timeOfDay <- df_time |> 
  ggplot() +
  geom_bar(aes(x= hour, fill = continent)) +
  # scale_fill_viridis_d(option = "turbo", begin = 0, end = 0.5) +
    scale_fill_manual(
      "Continent",
      values = rev(c("#001219", "#0a9396", "#94d2bd", "#ABB86D", "#ee9b00",  "#9b2226"))
    ) +
  scale_y_continuous("Number of records", expand = expansion(mult=c(0,NA))) +
    scale_x_continuous("Hour of observation")
)

## Figure 4c: North American record timing by species ----
# Moved here from Figure 5. Restricted to exactly the records the use-availability
# models use (North America, 2019-2025, geoprivacy not obscured; n = 198), so the
# panel's N matches Figure 5 rather than the fuller descriptive set. Species are
# assigned only from research-grade identifications; grouping and colours match
# Figure 5. The box on panel (a) marks this region.
dddd_model <- dddd |>
  dplyr::filter(geoprivacy != "obscured", year(datetime.date) %in% 2019:2025)

speciesLevelsF4 <- rev(c("Lasiurus borealis", "Lasionycteris noctivagans",
                         "Eptesicus fuscus", "Other/unknown"))
speciesColorsF4 <- c("Lasiurus borealis"         = "#9E2A2B",
                     "Lasionycteris noctivagans" = "#7072A0",
                     "Eptesicus fuscus"          = "#A8541F",
                     "Other/unknown"             = "grey60")
speciesLabelsF4 <- c("Lasiurus borealis"         = "*Lasiurus borealis*",
                     "Lasionycteris noctivagans" = "*Lasionycteris noctivagans*",
                     "Eptesicus fuscus"          = "*Eptesicus fuscus*",
                     "Other/unknown"             = "Other/unknown")

# Legend sits outside on the right, matching panels (a) and (b); an inside legend
# collides with the autumn peak and with the panel tag (placed inside at top-left).
(p_NoAm_timing <- dddd_model |>
  mutate(group = factor(group, levels = speciesLevelsF4)) |>
  ggplot() +
  geom_histogram(aes(x = yday, fill = group), binwidth = 14, boundary = 0,
                 colour = "white", linewidth = 0.1) +
  scale_fill_manual(NULL, values = speciesColorsF4, labels = speciesLabelsF4, drop = FALSE,
                    guide = guide_legend(ncol = 1, reverse = TRUE)) +
  scale_x_continuous("Day of year", breaks = datebreaks, labels = datelabs,
                     expand = c(0.01, 0)) +
  scale_y_continuous("Number of records", expand = expansion(mult = c(0, 0.08))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1),
    legend.text = ggtext::element_markdown(size = 9),
    legend.key.size = unit(11, "pt"),
    legend.key.spacing.y = unit(2, "pt")
  ))


# Figure 4 combo ------
# (a) global map with the North America box, (b) hour of observation (global), then
# (c) North American record timing -- global to regional, handing off to Figure 5.
F4_combo <- iNat_map_grid_circle + p_timeOfDay + p_NoAm_timing +
  plot_layout(
  design = "
    1
    2
    3
    ",
  heights = c(2,1,1),
  axis_titles = "collect_y"
) +
  plot_annotation(
    tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") &
  theme(
    plot.tag.position  = c(0.1,0.95),
    plot.tag = element_text(size = 10)
  )

ggsave(F4_combo, filename =  "figs/F4_combo.png", dpi = 600, width = 8, height = 9.5)
ggsave(F4_combo, filename =  "figs/F4_combo.svg", width = 8, height = 9.5)




# Figure 5 is now the use-availability driver figure, built in "8_fit iNat models.R".
# The per-species timing panels that used to make up Figure 5 here were superseded by
# it; the all-species version now lives above as Figure 4c.


# Tables + summaries ------
df |> 
  as.data.frame() |>
  count(admin) |> 
  arrange(desc(n))
df |> 
  as.data.frame() |>
  count(continent) 

df |> 
  as.data.frame() |> 
  dplyr::filter(quality_grade == "research",species_name!="") |> 
  count(quality_grade, family_name, species_name) |> 
  arrange(desc(n))

df |> 
  as.data.frame() |> 
  dplyr::filter(continent == "North America") |> 
  mutate(
    lon_cut = cut(longitude, breaks = c(-Inf, -100, Inf)),
    lat_cut = cut(latitude, breaks = c(-Inf, 24, Inf))
  ) |> 
  count(lon_cut, lat_cut)