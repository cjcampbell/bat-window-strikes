# Load --------------
library(sf)
library(terra)
library(rnaturalearth)
library(lubridate)
library(patchwork)

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


ggplot() +
  geom_sf(world, mapping = aes(), fill = "grey95", color = "grey60", linewidth = 0.1) +
  geom_sf(df_grid_count_grid, mapping = aes(), color = "grey50", fill = NA) +
  geom_sf(df_grid_count, mapping = aes(size = n), alpha = 0.5, color = "#00B31B") +
  scale_size_continuous(
    "Number of records",
    breaks = c(1,10,40,70),
  )

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

ggplot(df_NoAm) +
  geom_histogram(aes(x = yday)) +
  facet_grid(lat_cut~lon_cut) +
  theme_classic()

ggplot(df_NoAm) +
  geom_bar(aes(x = yday)) +
  facet_wrap(~lon_cut) +
  theme_classic()

ggplot(df_NoAm) +
  geom_bar(aes(x = yday)) +
  facet_wrap(~lat_cut) +
  theme_classic()


## Visualize species in North America -----
df_NoAm %>% 
  ggplot() +
  geom_histogram(aes(x = yday, fill = scientific_name)) +
  facet_wrap(~quality_grade+scientific_name) +
  theme_classic()

df_NoAm %>% 
  dplyr::filter(quality_grade == "research") %>% 
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


# Summarize taxonomy ----------



