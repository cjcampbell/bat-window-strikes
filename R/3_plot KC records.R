

speciesColors <- c(
  "Big brown" = "#A16928",
  "Evening" = "#124559", 
  "Vespertilionidae" = "grey50", 
  "Eastern Red" = "#E46844", 
  "Tricolored" = "#f4d35e", 
  "Silver-haired" = "#484357")

df_discovery2 <- df_discovery2 %>% 
  mutate(
    species = factor(species, levels = c(
      "Evening", "Big brown", "Eastern Red",
      "Silver-haired", "Tricolored", "Vespertilionidae")
    ),
    plotGroup = factor(plotGroup, levels = arrange(count(df_discovery2, plotGroup), desc(n))$plotGroup),
    Building_side_general = case_when(
      Building_side %in% c("N", "S") ~ "N/S",
      Building_side %in% c("E", "W") ~ "E/W"
    ),
    season = case_when(yday<200~ "spring", TRUE ~ "fall"),
    season = factor(season, levels = c("spring", "fall")),
    Building_side = factor(Building_side, levels = c("N", "E", "S", "W")),
    locality = gsub("\\.", "", locality),
    locality = gsub("St", "", locality),
    locality = gsub("eet", "", locality),
    locality = gsub("Rd", "", locality),
    locality = gsub("Boulevard", "", locality),
    locality = stringr::str_trim(locality),
    locality = case_when(
      grepl("1100 Walnut", locality) ~ "1100 Walnut",
      locality == "Downtown KC" ~ "Downtown",
      TRUE ~ locality
    )
  ) 


# Quick exploratory plot.

rowWidth <- 10
p_waffle <- df_discovery2 %>% 
  dplyr::select(species) %>% 
  arrange(species) %>% 
  mutate(
    pos = row_number(),
    y = ceiling(pos / rowWidth),
    x = pos %% rowWidth,
    x = case_when(x == 0 ~ rowWidth, TRUE ~ x)
  ) %>% 
  ggplot() +
  geom_tile(aes(fill = species, x = x, y=y,), color = "white", size = 2) +
  scale_fill_manual(
    "Species",
    values = speciesColors,
    breaks = levels(df_discovery2$species)
    ) +
  scale_y_reverse() +
  theme_minimal() +
  coord_equal() +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1)
  ) 

ggsave(p_waffle, filename = "figs/KC_count_waffle.svg", width = 5, height = 3)

# df_discovery2 %>% 
#   count(species) %>% 
#   mutate(cumsum = cumsum(n))
#   ggplot() +
#   geom_col(aes(x = 1, y = n, fill = species)) +
#   scale_fill_manual(
#     "Species",
#     values = speciesColors,
#     breaks = levels(df_discovery2$species)
#   ) 



df_discovery2 %>% 
  count(season) %>% 
  dplyr::mutate(prop = round(n/sum(n)*100), .by = c("season"))


# Plot of counts pooled by year.
datebreaks <- yday(lubridate::mdy(paste(1:12, "-1-2025")))
datelabs <-  month.abb

p_yday_species <- df_discovery2 %>%  
  ggplot() + 
  geom_histogram(
    mapping = aes(x = yday, fill = species),
    binwidth = 14
    ) +
  scale_y_continuous(
    "Count of bats discovered by survey",
    expand = expansion(add = c(0, 2), mult = c(0, .2)),
    breaks = seq(0,50, by = 5),
    minor_breaks = seq(0,50, by = 1)
    ) +
  scale_x_continuous(
    "Day of year",
    breaks = datebreaks,
    labels = datelabs
  ) +
  facet_wrap(
    ~plotGroup,
    scales = "free_y",
    axes = "all") +
  scale_fill_manual(
    "Species",
    values = speciesColors,
    breaks = levels(df_discovery2$species)
  ) +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_line(color = "grey95", linewidth = 0.1),
    strip.background = element_blank(),
    legend.position = "none",
    axis.text.x.bottom = element_text(
      angle = 45, hjust = 1, vjust = 1
    )
  )





p_yday_species_list <- lapply(c("Evening", "Big brown", "Eastern Red", "Others"), function(x) {
  dd <- df_discovery2 %>%  
    dplyr::filter(plotGroup == x)
  
  if(x == "Evening") {
    mylimits <- c(0, NA)
  } else {
    mylimits <- c(0, 7)
  }
  
  dd %>% 
    ggplot() + 
    geom_histogram(
      mapping = aes(x = yday, fill = species),
      binwidth = 14
    ) +
    scale_y_continuous(
      "Count of bats discovered by survey",
      expand = c(0,0),
      #expand = expansion(add = c(0, 2), mult = c(0, .2)),
      limits = mylimits,
      breaks = seq(0,50, by = 5),
      minor_breaks = seq(0,50, by = 1)
    ) +
    scale_x_continuous(
      "Day of year",
      breaks = datebreaks,
      labels = datelabs, 
      limits = c(0, 365),
      expand = c(0,0)
    ) +
    facet_wrap(
      ~plotGroup,
      scales = "free_y",
      axes = "all") +
    scale_fill_manual(
      "Species",
      values = speciesColors,
      breaks = levels(df_discovery2$species)
    ) +
    theme(
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.1),
      strip.background = element_blank(),
      legend.position = "none",
      axis.text.x.bottom = element_text(
        angle = 45, hjust = 1, vjust = 1
      )
    )
})
p_yday_species2 <- patchwork::wrap_plots(p_yday_species_list, axis_titles = "collect")


ggsave(p_yday_species2 + theme(legend.position = "none"), filename = "figs/KC_yday_histogram.svg", width = 6, height = 3)


## building direction plots : -----

df_discovery2 %>% 
  dplyr::filter(!is.na(Building_side)) %>% 
  ggplot() +
  aes(x = Building_side) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=0) +
  facet_grid(species ~ season, scales = "free_y") +
  scale_y_continuous(expand = expansion(add=c(0,2))) +
  coord_polar(start = -pi/4)


df_discovery2 %>% 
  mutate(Building_side_general = case_when(
    Building_side %in% c("N", "S") ~ "N/S",
    Building_side %in% c("E", "W") ~ "E/W"
    )
  ) %>% 
  dplyr::filter(!is.na(Building_side)) %>% 
  count(species, Building_side_general) %>% 
  dplyr::mutate(prop = round(n/sum(n)*100), .by = c("species"))


df_discovery2 %>% 
  dplyr::filter(!is.na(Building_side)) %>% 
  count(Building_side_general) %>% 
  dplyr::mutate(prop = round(n/sum(n)*100))


df_discovery2 %>% 
  dplyr::filter(!is.na(Building_side)) %>% 
  count(season, Building_side_general) %>% 
  dplyr::mutate(prop = round(n/sum(n)*100), .by = c("season"))





# # Load iNat records -----------
# 
# focalObs <- get_inat_obs_user("redtail5")
# 
# mySearch.DeAnnObs <- searchBuilder(user_id="redtail5", taxon_id = 40268)
# howManyResults(mySearch.DeAnnObs)
# 
# 
# mydownloadedObs2 <- downloadResults(mysearch2)