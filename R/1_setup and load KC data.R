
source("R/0_funs.R")
library(readxl)

# This script tidies records kept during structured surveys.
# Skip ahead to the next script to start with tidy, replicable csv files.

# Load discovery data -----------------------------------------------------

df_discovery <- read_excel(
    "data/Excel data from iNaturalist and cross checked with Lakeside rehab data.xlsx",
    range = "A1:P133"
  ) %>% 
  dplyr::rename(date = Date) %>% 
  dplyr::mutate(
    date = as_date(date),
    yday = yday(date),
    id = row_number()
    ) %>% 
  dplyr::filter(!is.na(date))

df_discovery2 <- df_discovery %>% 
  dplyr::select(id, starts_with("Bldg"),"Blg  C               1201 Walnut", "Other") %>% 
  pivot_longer(
    cols = c(starts_with("Bldg"), "Blg  C               1201 Walnut",  "Other"),
    names_to = "Building",
    values_to = "Building_side"
    ) %>% 
  dplyr::filter(!is.na(Building_side)) %>% 
  left_join(df_discovery, .) %>% 
  dplyr::select(-starts_with("Bldg"), -starts_with("Blg"), -"Other") %>% 
  dplyr::rename(
    locality = `Location (Corrections in green text)`,
    species = `Species (Corrections in green text)`
  ) %>% 
  # Cleaning
  dplyr::mutate(
    paired = case_when(grepl("/2", Building_side) ~ "Y", TRUE ~ "N"),
    
    species = case_when(
      species == "Vesper Bat - LNC said LBB" ~ "Vespertilionidae",
      lengths(gregexpr("\\W+", species)) + 1 == 4 ~ word(species, 3,4),
      lengths(gregexpr("\\W+", species)) + 1 == 5 ~ word(species, 4,5),
      TRUE ~ species
    ),
    species = stringr::str_trim(gsub(" Bat", "", species)),
    species = case_when(
      species == "Evening - not included in 1st manusc" ~ "Evening",
      species == "Evenings Big Brown per author" ~ "Big Brown",
      species == "Vespers" ~ "Vespertilionidae",
      startsWith("Vesper", species) ~ "Vespertilionidae",
      species == "Big Brown" ~ "Big brown",
      TRUE ~ species
    ),
    species = case_when(species == "Big Brown" ~ "Big brown", TRUE ~ species),
    
    Building = case_when(Building == "Other" ~ Building_side),
    Building = gsub("N -", "", Building),
    Building = gsub("S -", "", Building),
    Building = gsub("E -", "", Building),
    Building = gsub("E- ", "", Building),
    Building = gsub("W -", "", Building),
    Building = gsub("W-", "", Building),
    Building = gsub("1/2 8 feet ", "", Building),
    Building = gsub("2/2 8 feet ", "", Building),
    Building = str_trim(Building),
    
    Building_side = gsub("1/2 ", "", Building_side),
    Building_side = gsub("2/2 ", "", Building_side),
    Building_side = case_when(
      startsWith(Building_side, "N") ~ "N",
      startsWith(Building_side, "S") ~ "S",
      startsWith(Building_side, "E") ~ "E",
      startsWith(Building_side, "W") ~ "W",
      TRUE ~ Building_side
    ),
    Building_side = case_when(
      Building_side == "not recorded" ~ NA,
      Building_side == "unknown" ~ NA,
      Building_side == "8 feet W - 1551 McGee" ~ "W",
      TRUE ~ Building_side
    ),
    Building_side = factor(Building_side, levels = c("N", "S", "E", "W")),
    plotGroup = case_when(species %in% c("Big brown", "Eastern Red", "Evening") ~ species, TRUE ~ "Others")
  ) %>% 
  dplyr::select(
    id, 
    date, yday,
    species, plotGroup,
    locality, Status, `Description Where Found`, 
    Notes, Building_side, 
    paired
  )

write.csv(df_discovery2, "out/data_derived/structured_surveys_bats_discovered.csv", row.names = F)



# Load survey data data ---------------------------------------------------

sd <- read_excel("data/Survey dates.xlsx", sheet = "Sheet1_tidy") %>% 
  dplyr::rename(year = 1) %>% 
  pivot_longer(cols = -1, names_to = "calDate", values_to = "survey") %>% 
  dplyr::filter(!is.na(survey)) %>% 
  dplyr::mutate(
    survey = TRUE,
    date = mdy(paste(calDate, year))
    ) %>% 
  dplyr::select(-calDate, -year)

fullDates <- data.frame(date = seq.Date(ymd("2019-09-01"), ymd("2024-12-31"), by = "1 days" ))
sd <- full_join(sd, fullDates) %>% 
  arrange(date) %>% 
  replace(is.na(.), FALSE) %>% 
  dplyr::mutate(
    yday = yday(date),
    yday_bin7 = cut(yday, breaks = seq(1,365,by=7))
    ) 

write.csv(sd, "out/data_derived/structured_surveys_schedule.csv", row.names = F)

sd_yday_count <- sd %>% 
  group_by(yday_bin7) %>% 
  dplyr::summarise(n_surveys = sum(survey)) %>% 
  dplyr::mutate(
    yday_bin7_lab = as.character(yday_bin7),
    yday = as.numeric(gsub("\\(", "", stringr::word(yday_bin7_lab, 1,1, sep = ","))),
  )
