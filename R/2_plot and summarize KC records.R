
# Setup ==============
source("R/0_funs.R")
library(sf)
library(geodata)
library(rinat)
library(tidyterra)
library(data.table)
library(magick)

speciesColors <- c(
  "Big brown" = "#A8541F",
  "Evening" = "#124559",
  "Vespertilionidae" = "#3B3B3B",
  "Eastern Red" = "#9E2A2B",
  "Tricolored" = "#f4d35e",
  "Silver-haired" = "#7072A0"
)
speciesColors2 <- speciesColors
names(speciesColors2) <- c( "Eptesicus fuscus", "Nycticeius humeralis", "Other/unknown", "Lasiurus borealis", "Perimyotis subflavus", "Lasionycteris noctivagans" )


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
    ),
    group = case_when(
      species == "Evening"       ~ "Nycticeius humeralis",
      species == "Big brown"     ~ "Eptesicus fuscus",
      species == "Eastern Red"   ~ "Lasiurus borealis",
      species == "Silver-haired" ~ "Lasionycteris noctivagans",
      species == "Tricolored"    ~ "Perimyotis subflavus",
      TRUE~ "Other/unknown"
    ),
    group_lab = case_when(
      group == "Other/unknown" ~ "Other/unknown",
      TRUE ~ paste0("<i>", group, "</i>")
    ),
    group_lab2 = case_when(
      group %in% c("Other/unknown", "Perimyotis subflavus", "Lasionycteris noctivagans") ~ "Other",
      TRUE ~ group_lab
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

(p_KC_records_by_species <- df_discovery2 %>% 
  count(group, group_lab) %>% 
  arrange(desc(n)) %>% 
  mutate(group_lab = factor(group_lab, levels = rev(unique(.$group_lab)))) %>% 
  ggplot() +
  geom_col(aes(y= group_lab, x = n, fill = group)) +
  geom_text(aes(y = group_lab, x = n, label = n), vjust = -0.25, color = "grey50", size = 2) +
  scale_x_continuous("Number of Kansas City records",
                     expand = expansion(add = c(0, 3))) +
  scale_y_discrete("Species") +
  scale_fill_manual(
    "",
    values = speciesColors2
  ) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = ggtext::element_markdown(),
  ))
ggsave("figs/KC_records_by_species.png", p_KC_records_by_species, dpi = 600)
saveRDS(p_KC_records_by_species, "tmp/p_KC_records_by_species.rds")


# 
ddd_bygroup <- df_discovery2 %>% 
  dplyr::mutate(species = case_when(group_lab2 != "Other" ~ species, TRUE ~ "Other")) %>% 
  count(group_lab2, species) %>% 
  arrange(n) %>% 
  mutate(group = factor(group_lab2, levels = c(.$group_lab2)))

ddd_Others <- df_discovery2 %>%
  dplyr::filter(group_lab2 == "Other") %>% 
  mutate(group_lab = case_when(group_lab == "Other/unknown"~ "Unknown", TRUE ~ group_lab)) %>% 
  count(group_lab) %>% 
  arrange(desc(n)) %>% 
  mutate(rownum = 4 - row_number())

(p_KC_records_barplot <- ggplot(ddd_bygroup) +
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
      aes(x = 50, y = 0.75 + (rownum / 3), label = paste0(group_lab, ": ", n)),
      text.color = "black",
      hjust = 0,
      fill = NA, label.colour = NA
    ) +
    scale_y_discrete("Species") +
    scale_x_continuous(
      "Number of structured survey records",
      expand = expansion(add = c(0, 0)),
      breaks = seq(0, 100, by = 25),
      limits = c(0, 95)
    ) +
    scale_fill_manual(
      "Species",
      values = speciesColors
    ) +
    theme(
      legend.position = "none",
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
)
saveRDS(p_KC_records_barplot, file = "tmp/p_KC_records_barplot.rds")



# Phenology plots of yday by species/group ==============

datebreaks <- yday(lubridate::mdy(paste(1:12, "-1-2025")))
datelabs <-  month.abb

# Separate for consistent y axes 
p_yday_species_list <- lapply(c("Evening", "Big brown", "Eastern Red", "Others"), function(x) {
  dd <- df_discovery2 %>%  
    dplyr::filter(plotGroup == x)
  
  if(x == "Evening") {
    mylimits <- c(0, NA)
  } else {
    mylimits <- c(0, 7)
  }
  
  if(x == "Others") {
    speciesColors_modified <- c(
      "Big brown" = "#A8541F",
      "Evening" = "#124559",
      "Vespertilionidae" = "grey50",
      "Eastern Red" = "#9E2A2B",
      "Tricolored" = "grey50",
      "Silver-haired" = "grey50"
    )
  } else {
    speciesColors_modified <- speciesColors
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
      ~group_lab2,
      scales = "free_y",
      axes = "all") +
    scale_fill_manual(
      "Species",
      values = speciesColors_modified,
      breaks = levels(df_discovery2$species)
    ) +
    coord_cartesian() +
    theme(
      panel.grid.major.y = element_line(color = "grey80"),
      panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.1),
      strip.background = element_blank(),
      strip.text = ggtext::element_markdown(),
      legend.position = "none",
      axis.text.x.bottom = element_text(
        angle = 45, hjust = 1, vjust = 1
      )
    )
})
p_yday_species2 <- patchwork::wrap_plots(p_yday_species_list, axis_titles = "collect")
# ggsave(p_yday_species2, filename = "figs/KC_yday_histogram.svg", width = 6, height = 3)
saveRDS(p_yday_species2, "tmp/p_yday_species2.rds")


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
  mutate(Building_side_general = case_when(
    Building_side %in% c("N", "S") ~ "N/S",
    Building_side %in% c("E", "W") ~ "E/W"
  )) %>%
  dplyr::filter(!is.na(Building_side)) %>%
  count(Building_side_general) %>%
  dplyr::mutate(prop = round(n / sum(n) * 100))

df_discovery2 %>% 
  dplyr::filter(!is.na(Building_side)) %>% 
  count(season, Building_side_general) %>% 
  dplyr::mutate(prop = round(n/sum(n)*100), .by = c("season"))


# Summaries ================================

df_discovery <- read.csv("out/data_derived/structured_surveys_bats_discovered.csv") %>% 
  mutate(date = as_date(date))
sd <- read.csv("out/data_derived/structured_surveys_schedule.csv") %>% 
  mutate(date = as_date(date))


# Bats discovered:
nrow(df_discovery)

# Surveys conducted:
sum(sd$survey)

# How many bats found alive:
df_discovery %>% count(Status)

signif(117/(117+12)*100,2)

# How many bats found in which season.
df_discovery2 %>% count(season) %>% mutate(prop = n / sum(n)*100)
