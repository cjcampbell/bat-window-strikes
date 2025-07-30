
# Setup ==============
source("R/0_funs.R")
library(sf)
library(geodata)
library(rinat)
library(tidyterra)
library(data.table)

speciesColors <- c(
  "Big brown" = "#A16928",
  "Evening" = "#124559",
  "Vespertilionidae" = "grey50",
  "Eastern Red" = "#E46844",
  "Tricolored" = "#f4d35e",
  "Silver-haired" = "#484357"
)


# Bar plot of frequencies -------------------------
df_discovery2 <- read.csv("out/data_derived/structured_surveys_bats_discovered.csv") %>% 
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

df_discovery2 %>% 
  count(species) %>% 
  arrange(desc(n)) %>% 
  mutate(species = factor(species, levels = unique(.$species))) %>% 
  ggplot() +
  geom_col(aes(x= species, y = n, fill = species)) +
  scale_y_continuous("Number of Kansas City records",
                     expand = expansion(add = c(0, 3))) +
  scale_x_discrete("Species")+
  scale_fill_manual(
    "Species",
    values = speciesColors,
    breaks = levels(df_discovery2$species)
  ) 

## Change to scientific names --------------

p_KC_records_by_species <- df_discovery2 %>% 
  mutate(
    group = case_when(
      species == "Evening"       ~ "Nycticeius humeralis",
      species == "Big brown"     ~ "Eptesicus fuscus",
      species == "Eastern Red"   ~ "Lasiurus borealis",
      species == "Silver-haired" ~ "Nycticeius humeralis",
      species == "Tricolored"    ~ "Perimyotis subflavus",
      TRUE~ "Other/unknown"
    ),
    group_lab = case_when(
      group == "Other/unknown" ~ "Other/unknown",
      TRUE ~ paste0("<i>", group, "</i>")
    )
  ) %>% 
  count(group, group_lab) %>% 
  arrange(desc(n)) %>% 
  mutate(group_lab = factor(group_lab, levels = unique(.$group_lab))) %>% 
  ggplot() +
  geom_col(aes(x= group_lab, y = n, fill = group)) +
  geom_text(aes(x = group_lab, y = n, label = n), vjust = -0.25, color = "grey50", size = 2) +
  scale_y_continuous("Number of Kansas City records",
                     expand = expansion(add = c(0, 3))) +
  scale_x_discrete("Species") +
  scale_fill_manual(
    "",
    values = 
      c(
        "Eptesicus fuscus" = "#A16928",
        "Nycticeius humeralis" = "#124559", 
        "Other/unknown" = "grey50", 
        "Lasiurus borealis" = "#E46844", 
        "Perimyotis subflavus" = "#f4d35e", 
        "Lasionycteris noctivagans" = "#484357")
    # breaks = levels(df_discovery2$species)
  ) +
  theme(
    legend.position = "none",
    axis.title.x.bottom = element_blank(),
    axis.text.x.bottom = ggtext::element_markdown()
  )
ggsave("figs/KC_records_by_species.png", p_KC_records_by_species, dpi = 600)



# Phenology plots of yday by species/group ==============
## All together ----------
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
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.1),
    strip.background = element_blank(),
    legend.position = "none",
    axis.text.x.bottom = element_text(
      angle = 45, hjust = 1, vjust = 1
    )
  )

## Separate for consistent y axes -------------

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
      panel.grid.major.y = element_line(color = "grey80"),
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


# Counts by season ===========

df_discovery2 %>% 
  count(season) %>% 
  dplyr::mutate(prop = round(n/sum(n)*100), .by = c("season"))



# Cardinal direction of discovery plots/summaries -----

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

