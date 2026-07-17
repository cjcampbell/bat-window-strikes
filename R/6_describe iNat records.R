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
# North America highlight for panel (c)/Figure 5. Two options are defined; the plot
# below uses `noam_shade`, with `noam_box` kept so we can swap back to the dashed box.
#
# noam_shade: the North America landmass via Natural Earth's `continent` field -- the
# same definition the analysis uses (continent == "North America"; includes Central
# America and the Caribbean, excludes South America). Already in myproj via `world`.
noam_shade <- dplyr::filter(world, continent == "North America")
# noam_box: previous rectangular region, segmentized in lon/lat so its edges follow the
# Equal Earth graticule. Overshoots into South America a little; superseded by noam_shade.
noam_box <- st_bbox(c(xmin = -126, ymin = 5, xmax = -68, ymax = 54), crs = st_crs(proj.wgs84)) |>
  st_as_sfc() |>
  st_segmentize(units::set_units(100, "km")) |>
  st_transform(myproj)

# Plot circle in grid.
iNat_map_grid_circle <- ggplot() +
  geom_sf(world, mapping = aes(), fill = "grey95", color = "grey60", linewidth = 0.1) +
  # TESTING: translucent tint over the North America landmass (ties to Figure 5's blue),
  # drawn beneath the record layers so the circles stay crisp on top. To go back to the
  # dashed box, comment this out and add, as the LAST geom:
  #   geom_sf(noam_box, mapping = aes(), fill = NA, color = "grey25", linewidth = 0.4,
  #           linetype = "22")
  geom_sf(noam_shade, mapping = aes(), fill = "#0072B2", alpha = 0.18, color = NA) +
  geom_sf(df_grid_count_grid, mapping = aes(), color = "grey50", fill = NA) +
  geom_sf(df_grid_count, mapping = aes(size = n), alpha = 0.5, color = "#00B31B") +
  scale_size_continuous(
    "Number of\nrecords",
    breaks = c(1,10,40,70),
  ) +
  scale_y_continuous(expand = c(0,0)) +
 # scale_x_continuous(expand = c(0,0)) +
  coord_sf(ylim = c(-6792.374, 8000)) +
  theme_void() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    # Legend tucked into the empty top-right corner of the map (over the NE Pacific).
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1)
  )
ggsave(iNat_map_grid_circle, filename =  "figs/iNat_map_grid_circle.png", dpi = 600)
ggsave(iNat_map_grid_circle, filename =  "figs/iNat_map_grid_circle.svg", width = 8, height = 8)

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
      values = rev(c("#001219", "#0a9396", "#94d2bd", "#ABB86D", "#ee9b00",  "#9b2226")),
      guide = guide_legend(nrow = 2, title.position = "top")
    ) +
  scale_y_continuous("Number of records", expand = expansion(mult=c(0,NA))) +
    scale_x_continuous("Hour of observation") +
  theme(legend.position = "bottom")
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

# Legend runs along the bottom as a 2x2 block (species names are long), keeping the
# plotting area clear of the autumn peak and the top-left panel tag.
(p_NoAm_timing <- dddd_model |>
  mutate(group = factor(group, levels = speciesLevelsF4)) |>
  ggplot() +
  geom_histogram(aes(x = yday, fill = group), binwidth = 14, boundary = 0,
                 colour = "white", linewidth = 0.1) +
  scale_fill_manual("Species", values = speciesColorsF4, labels = speciesLabelsF4, drop = FALSE,
                    guide = guide_legend(ncol = 2, byrow = TRUE, reverse = TRUE,
                                         title.position = "top")) +
  scale_x_continuous("Day of year", breaks = datebreaks, labels = datelabs,
                     expand = c(0.01, 0)) +
  scale_y_continuous("Number of records", expand = expansion(mult = c(0, 0.08))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1),
    legend.position = "bottom",
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
  heights = c(1.2,1,1),
  axis_titles = "collect_y"
) +
  plot_annotation(
    tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") &
  # Standardize legends across all three panels: uniform key box, uniform text size, and
  # titles trimmed to 9 pt so they no longer tower over the body / read large next to
  # Figure 5. element_markdown keeps panel (c)'s italic species names.
  theme(
    plot.tag.position  = c(0.1,0.95),
    plot.tag = element_text(size = 10),
    legend.title         = element_text(size = 9),
    legend.text          = ggtext::element_markdown(size = 9),
    legend.key.size      = unit(12, "pt"),
    legend.key.spacing.y = unit(2, "pt")
  )

ggsave(F4_combo, filename =  "figs/F4_combo.png", dpi = 400, width = 6.5, height = 5.5)
ggsave(F4_combo, filename =  "figs/F4_combo.svg", width = 6.5, height = 5.5)



# Record counts: where every number in the manuscript comes from ------
records <- df |>
  as.data.frame() |>
  mutate(
    date       = as_date(as.POSIXct(datetime)),
    year       = year(date),
    is_species = quality_grade == "research" & !is.na(species_name) & species_name != ""
  )

# Count by country.
n_country <- records |> filter(adm0_a3 %in% c("USA", "CAN", "MEX")) |> count(adm0_a3)

# Narrow to the North-American analysis set w/ multiple filters
#   step 1: keep the North America continent (this is the boundary the models use)
#   step 2: drop records whose coordinates iNaturalist has hidden ("obscured")
#   step 3: keep only 2019-2025 (the years the background covers)
na_continent  <- filter(records, continent == "North America")
na_unobscured <- filter(na_continent, geoprivacy != "obscured")
na_analysis   <- filter(na_unobscured, year %in% 2019:2025)   # == collisions_noam in 7_prep

# How many are identified to species?
# how many are the two long-distance migratory tree bats (eastern red + silver-haired)?
n_species_level <- sum(na_analysis$is_species)
n_redsilver <- sum(na_analysis$is_species &
                     na_analysis$species_name %in% c("Lasiurus borealis", "Lasionycteris noctivagans"))

# Check that no record drops out due to missing predictor.
ua_file <- "data/derived/useavail_points.csv"
if (file.exists(ua_file)) {
  ua <- fread(ua_file)
  n_used_full  <- sum(ua$used == 1)
  n_used_radar <- sum(ua$used == 1 & !is.na(ua$traffic) & ua$dist_km < 200)
} else {
  n_used_full <- nrow(na_analysis); n_used_radar <- NA_integer_
}

# The step-by-step breakdown, saved so the manuscript's counts are trackable.
count_breakdown <- tibble::tribble(
  ~record_set,                      ~n,                                    ~how_it_is_defined,
  "All retained records",           nrow(records),                        "CJ.manual.check == 'y', all continents / years",
  "United States",                  n_country$n[n_country$adm0_a3 == "USA"], "by country (admin == USA), after coastal fix",
  "Canada",                         n_country$n[n_country$adm0_a3 == "CAN"], "by country (admin == Canada)",
  "Mexico",                         n_country$n[n_country$adm0_a3 == "MEX"], "by country (admin == Mexico)",
  "North America (continent)",      nrow(na_continent),                   "step 1: continent == 'North America' (incl. Central America)",
  ".. coordinates not hidden",      nrow(na_unobscured),                  "step 2: + geoprivacy != 'obscured'",
  ".. years 2019-2025 = ANALYSED",  nrow(na_analysis),                    "step 3: + year in 2019:2025  (this is the set the models use)",
  ".. with a nearby radar station", n_used_radar,                         "+ NEXRAD station < 200 km with a night-before traffic value",
  "Identified to species",          n_species_level,                      "of the analysed set: research-grade and named to species",
  ".. red + silver-haired bats",    n_redsilver,                          "of those: Lasiurus borealis + Lasionycteris noctivagans"
)
cat("\nRecord-count breakdown (how the manuscript's numbers are derived):\n")
print(count_breakdown, n = Inf)
cat(sprintf(
  "\nNo records are lost to missing predictors: all %d analysed records have building height and\n  light (unmapped -> 0); %d are dropped only by the radar match, leaving %d with radar traffic.\n",
  n_used_full, n_used_full - n_used_radar, n_used_radar))
cat(sprintf("Species share: %d of %d analysed records identified to species (%.0f%%); red + silver = %d of %d species-level (%.0f%%), or %d of %d analysed (%.0f%%).\n",
            n_species_level, nrow(na_analysis), 100 * n_species_level / nrow(na_analysis),
            n_redsilver, n_species_level, 100 * n_redsilver / n_species_level,
            n_redsilver, nrow(na_analysis), 100 * n_redsilver / nrow(na_analysis)))

write.csv(count_breakdown, "data/derived/record_count_breakdown.csv", row.names = FALSE)