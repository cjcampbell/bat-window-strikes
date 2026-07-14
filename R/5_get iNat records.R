library(data.table)
library(rinat)
library(readxl)
# Chiroptera tax ID: 40268

all <- read_excel("data/tidy 6-2-2025 iNat bat collision Records sorted.xlsx") |> 
  dplyr::mutate(
    id = case_when(id == "112964491a" ~ "112964491", TRUE ~ id),
    id = as.integer(id),
    user_id = as.integer(user_id)
    )

# Download all info for referenced observations.
iNatObsFile <- "data/iNat downloads/iNatObsFile.csv"
iNatProjectsToCheck <- "data/iNat downloads/iNatProjectsToCheck.csv"

for(i in 1:nrow(all)) {
  
  cat(paste0("\rWorking on ", i))
  if(file.exists(iNatObsFile)) {iNatObs_downloaded <- fread(iNatObsFile)}
  if(exists("iNatObsFile")) { if(all$id[i] %in% iNatObs_downloaded$id) {next} }
  if(is.na(all$id[i])) next
  
  df <- rinat::get_inat_obs_id(all$id[i])
  Sys.sleep(1)
  # Paste in large empty dataframe to start.
  df2 <- data.frame(
    id = df$id,
    observed_on = if(is.null(df$observed_on)) NA else df$observed_on,
    description = if(is.null(df$description)) NA else df$description,
    latitude = if(is.null(df$latitude)) NA else df$latitude,
    longitude = if(is.null(df$longitude)) NA else df$longitude,
    map_scale = if(is.null(df$map_scale)) NA else df$map_scale,
    timeframe = if(is.null(df$timeframe)) NA else df$timeframe,
    species_guess = if(is.null(df$species_guess)) NA else df$species_guess,
    user_id = if(is.null(df$user_id)) NA else df$user_id,
    taxon_id = if(is.null(df$taxon_id)) NA else df$taxon_id,
    created_at = if(is.null(df$created_at)) NA else df$created_at,
    updated_at = if(is.null(df$updated_at)) NA else df$updated_at,
    place_guess = if(is.null(df$place_guess)) NA else df$place_guess,
    id_please = if(is.null(df$id_please)) NA else df$id_please,
    observed_on_string = if(is.null(df$observed_on_string)) NA else df$observed_on_string,
    iconic_taxon_id = if(is.null(df$iconic_taxon_id)) NA else df$iconic_taxon_id,
    num_identification_agreements = if(is.null(df$num_identification_agreements)) NA else df$num_identification_agreements,
    num_identification_disagreements = if(is.null(df$num_identification_disagreements)) NA else df$num_identification_disagreements,
    time_observed_at = if(is.null(df$time_observed_at)) NA else df$time_observed_at,
    time_zone = if(is.null(df$time_zone)) NA else df$time_zone,
    location_is_exact = if(is.null(df$location_is_exact)) NA else df$location_is_exact,
    delta = if(is.null(df$delta)) NA else df$delta,
    positional_accuracy = if(is.null(df$positional_accuracy)) NA else df$positional_accuracy,
    private_latitude = if(is.null(df$private_latitude)) NA else df$private_latitude,
    private_longitude = if(is.null(df$private_longitude)) NA else df$private_longitude,
    geoprivacy = if(is.null(df$geoprivacy)) NA else df$geoprivacy,
    quality_grade = if(is.null(df$quality_grade)) NA else df$quality_grade,
    positioning_method = if(is.null(df$positioning_method)) NA else df$positioning_method,
    positioning_device = if(is.null(df$positioning_device)) NA else df$positioning_device,
    out_of_range = if(is.null(df$out_of_range)) NA else df$out_of_range,
    license = if(is.null(df$license)) NA else df$license,
    uri = if(is.null(df$uri)) NA else df$uri,
    observation_photos_count = if(is.null(df$observation_photos_count)) NA else df$observation_photos_count,
    comments_count = if(is.null(df$comments_count)) NA else df$comments_count,
    zic_time_zone = if(is.null(df$zic_time_zone)) NA else df$zic_time_zone,
    oauth_application_id = if(is.null(df$oauth_application_id)) NA else df$oauth_application_id,
    observation_sounds_count = if(is.null(df$observation_sounds_count)) NA else df$observation_sounds_count,
    identifications_count = if(is.null(df$identifications_count)) NA else df$identifications_count,
    captive = if(is.null(df$captive)) NA else df$captive,
    community_taxon_id = if(is.null(df$community_taxon_id)) NA else df$community_taxon_id,
    site_id = if(is.null(df$site_id)) NA else df$site_id,
    old_uuid = if(is.null(df$old_uuid)) NA else df$old_uuid,
    public_positional_accuracy = if(is.null(df$public_positional_accuracy)) NA else df$public_positional_accuracy,
    mappable = if(is.null(df$mappable)) NA else df$mappable,
    cached_votes_total = if(is.null(df$cached_votes_total)) NA else df$cached_votes_total,
    last_indexed_at = if(is.null(df$last_indexed_at)) NA else df$last_indexed_at,
    private_place_guess = if(is.null(df$private_place_guess)) NA else df$private_place_guess,
    uuid = if(is.null(df$uuid)) NA else df$uuid,
    taxon_geoprivacy = if(is.null(df$taxon_geoprivacy)) NA else df$taxon_geoprivacy,
    tag_list = if(is.null(df$tag_list)) NA else paste(df$tag_list, collapse = "/"),
    user_login = if(is.null(df$user_login)) NA else df$user_login,
    iconic_taxon_name = if(is.null(df$iconic_taxon_name)) NA else df$iconic_taxon_name,
    captive_flag = if(is.null(df$captive_flag)) NA else df$captive_flag,
    created_at_utc = if(is.null(df$created_at_utc)) NA else df$created_at_utc,
    updated_at_utc = if(is.null(df$updated_at_utc)) NA else df$updated_at_utc,
    time_observed_at_utc = if(is.null(df$time_observed_at_utc)) NA else df$time_observed_at_utc,
    faves_count = if(is.null(df$faves_count)) NA else df$faves_count,
    owners_identification_from_vision = if(is.null(df$owners_identification_from_vision)) NA else df$owners_identification_from_vision,
    observation_field_values = if(is.null(df$observation_field_values)) NA else paste(df$observation_field_values, collapse = "/")
  )
  stopifnot(nrow(df2) == 1)
  fwrite(df2, iNatObsFile, row.names = F, append = T)
  
  if(length(df$project_observations) == 0) next
  df_proj <- data.frame(
    id = df$id, 
    projectID = unlist(df$project_observations[,10][1]),
    projectTitle = unlist(df$project_observations[,10][2])
    )
  
  fwrite(df_proj, iNatProjectsToCheck, row.names = F, append = T)
}

## Re-download the ID-based search for observations ----

for(i in 1:nrow(all)) {
  if(file.exists("data/iNat downloads/observations_from_our_searches.csv") & !exists("obs_downloaded_already")) {
    obs_downloaded_already <- fread("data/iNat downloads/observations_from_our_searches.csv")
  }
  if(all$id[i] %in% obs_downloaded_already$id) next
  cat(paste("\r Downloads complete for", round(i/nrow(all)*100), "%"))
  mysearch1 <- searchBuilder(taxon_id = 40268, id = all$id[i])
 #  howManyResults(mysearch1)
  mydownloadedObs1 <- downloadResults(mysearch1)
  fwrite(mydownloadedObs1, "data/iNat downloads/observations_from_our_searches.csv", row.names = F, append = T)
  Sys.sleep(1)
}

# Then check projects. ------

iNatProjectsToCheck_df <- fread(iNatProjectsToCheck)

projectsWithBats <- iNatProjectsToCheck_df |> 
  dplyr::select(projectID, projectTitle) |> 
  distinct() |> 
  arrange(projectTitle) |> 
  dplyr::filter(
    grepl("Collision", projectTitle) |
      grepl("collision", projectTitle) |
      grepl("Bird Safe", projectTitle) |
      grepl("Strike", projectTitle) |
      grepl("Glass City Bird", projectTitle) |
      grepl("Lights out", projectTitle) |
      grepl("Bird Strikes", projectTitle) 
  )

projNotes <- read_excel("data/Current 6-2-2025 iNat bat collision Records sorted.xlsx", sheet = 2) |> 
  dplyr::mutate(slug = gsub("https://www.inaturalist.org/projects/", "", url))

projects2Lookup <- unique(c(projectsWithBats$projectID, projNotes$slug))

## For each project, assemble a search for all bat observations. -----

# Download all specified observations from identified projects to this location:
iNatObservationsFromProjects_path <- "data/iNat downloads/iNatObservationsFromProjects.csv"

mysearch <- searchBuilder(taxon_id = 40268, projects = projects2Lookup)
howManyResults(mysearch)

mydownloadedObs <- downloadResults(mysearch)
fwrite(mydownloadedObs, "data/iNat downloads/observations_from_project_searches.csv", row.names = F)


# Search by fields ------

fields2Search <- c(
  "bird-window collision=", "Window impact=", "Window Strike=", 
  "Window Collision=", "Window kill=",
  "Building or Site Name of Collision=", 
  "Building or Site Name=",
  "Collision evidence=", "Collision Location=", "Bird Strike Location=",
  "Location of Bird Strike="
  )

for(srch in fields2Search) {
  mysearch2 <- searchBuilder(taxon_id = 40268, field = srch)
  rr <- howManyResults(mysearch2)
  if(rr == 0) next
  
  mydownloadedObs2 <- downloadResults(mysearch2)
  fwrite(
    mydownloadedObs2, 
    "data/iNat downloads/observations_from_fields_searches.csv", row.names = F,
    append = T)

}

# Search by queries and terms -----

queries2search <- c(
  "window strike", "window collision", "window impact",
  "building strike", "building collision", "building impact",
  "bird strike", "bird window", "lights out", "window search"
)

for(srch in queries2search) {
  mysearch2 <- searchBuilder(taxon_id = 40268, query = srch)
  rr <- howManyResults(mysearch2)
  if(rr == 0) next
  
  mydownloadedObs2 <- downloadResults(mysearch2)
  fwrite(
    mydownloadedObs2, 
    "data/iNat downloads/observations_from_query_searches.csv", row.names = F,
    append = T)
  Sys.sleep(1)
}



# Combine all data --------------------------------------------------------

df0 <- c("data/iNat downloads/observations_from_our_searches.csv",
  "data/iNat downloads/observations_from_project_searches.csv",
  "data/iNat downloads/observations_from_fields_searches.csv",
  "data/iNat downloads/observations_from_query_searches.csv",
  "data/iNat downloads/iNatObsFile.csv") |> 
  lapply(function(x) {
    read.csv(x) |> 
      distinct()
  }) |> 
  reduce(full_join) |> 
  full_join(all) |> 
  dplyr::select(
    id, scientific_name, datetime, place_guess, everything()
  )

# Tidy:
df <- df0 |> 
  distinct |> 
  group_by(id) |> 
  summarise(
    scientific_name = max(scientific_name,na.rm = T),
    datetime = max(datetime,na.rm = T),
    url = max(url,na.rm = T),
    place_guess = max(place_guess,na.rm = T),
    description = max(description,na.rm = T),
    latitude = max(latitude,na.rm = T),
    longitude = max(longitude,na.rm = T),
    tag_list = max(tag_list,na.rm = T),
    common_name = max(common_name,na.rm = T),
    image_url = max(image_url,na.rm = T),
    user_login = max(user_login,na.rm = T),
    species_guess = max(species_guess,na.rm = T),
    iconic_taxon_name = max(iconic_taxon_name,na.rm = T),
    taxon_id = max(taxon_id,na.rm = T),
    num_identification_agreements = max(num_identification_agreements,na.rm = T),
    num_identification_disagreements = max(num_identification_disagreements,na.rm = T),
    observed_on_string = max(observed_on_string,na.rm = T),
    observed_on = max(observed_on,na.rm = T),
    time_observed_at = max(time_observed_at,na.rm = T),
    time_zone = max(time_zone,na.rm = T),
    positional_accuracy = max(positional_accuracy,na.rm = T),
    public_positional_accuracy = max(public_positional_accuracy,na.rm = T),
    geoprivacy = max(geoprivacy,na.rm = T),
    taxon_geoprivacy = max(taxon_geoprivacy,na.rm = T),
    coordinates_obscured = max(coordinates_obscured,na.rm = T),
    positioning_method = max(positioning_method,na.rm = T),
    positioning_device = max(positioning_device,na.rm = T),
    user_id = max(user_id,na.rm = T),
    user_name = max(user_name,na.rm = T),
    created_at = max(created_at,na.rm = T),
    updated_at = max(updated_at,na.rm = T),
    quality_grade = max(quality_grade,na.rm = T),
    license = max(license,na.rm = T),
    sound_url = max(sound_url,na.rm = T),
    oauth_application_id = max(oauth_application_id,na.rm = T),
    captive_cultivated = max(captive_cultivated,na.rm = T),
    field.window.impact = max(field.window.impact,na.rm = T),
    map_scale = max(map_scale,na.rm = T),
    timeframe = max(timeframe,na.rm = T),
    id_please = max(id_please,na.rm = T),
    iconic_taxon_id = max(iconic_taxon_id,na.rm = T),
    location_is_exact = max(location_is_exact,na.rm = T),
    delta = max(delta,na.rm = T),
    private_latitude = max(private_latitude,na.rm = T),
    private_longitude = max(private_longitude,na.rm = T),
    out_of_range = max(out_of_range,na.rm = T),
    uri = max(uri,na.rm = T),
    observation_photos_count = max(observation_photos_count,na.rm = T),
    comments_count = max(comments_count,na.rm = T),
    zic_time_zone = max(zic_time_zone,na.rm = T),
    observation_sounds_count = max(observation_sounds_count,na.rm = T),
    identifications_count = max(identifications_count,na.rm = T),
    captive = max(captive,na.rm = T),
    community_taxon_id = max(community_taxon_id,na.rm = T),
    site_id = max(site_id,na.rm = T),
    old_uuid = max(old_uuid,na.rm = T),
    mappable = max(mappable,na.rm = T),
    cached_votes_total = max(cached_votes_total,na.rm = T),
    last_indexed_at = max(last_indexed_at,na.rm = T),
    private_place_guess = max(private_place_guess,na.rm = T),
    uuid = max(uuid,na.rm = T),
    captive_flag = max(captive_flag,na.rm = T),
    created_at_utc = max(created_at_utc,na.rm = T),
    updated_at_utc = max(updated_at_utc,na.rm = T),
    time_observed_at_utc = max(time_observed_at_utc,na.rm = T),
    faves_count = max(faves_count,na.rm = T),
    owners_identification_from_vision = max(owners_identification_from_vision,na.rm = T),
    observation_field_values = max(observation_field_values,na.rm = T),
    `Search word` = max(`Search word`,na.rm = T),
    place_town_name = max(place_town_name,na.rm = T),
    place_state_name = max(place_state_name,na.rm = T),
    place_country_name = max(place_country_name,na.rm = T),
    `alive?` = max(`alive?`,na.rm = T),
    `Dead?` = max(`Dead?`,na.rm = T),
    `screener_1_notes` = max(screener_1_notes,na.rm = T),
    `screener_2_notes` = max(`screener_2_notes`,na.rm = T),
  ) |> 
  dplyr::filter(
    !is.na(id),
    # Remove a user's observations:
   #  user_login != "redtail5"
    ) |> 
  arrange(id)

fwrite(df, "data/iNat_observations_tidy.csv", row.names = F)

## Check if previous checks have been completed, integrate if so. -----

previousCheckedFile <- "data/iNat_observations_tidy_manualChecks_20250711.csv"
if(file.exists(previousCheckedFile)) {
  df_checked <- fread(previousCheckedFile) |> 
    dplyr::select(id, `CJ manual check`, `CJ notes`)
  left_join(df, df_checked) |> 
    dplyr::select(id, `CJ manual check`, `CJ notes`, everything()) |> 
    fwrite("data/iNat_observations_tidy_withChecks.csv", row.names = F)
}


# Get higher taxonomic status -----------

library("taxize")

df_tax <- fread("data/iNat_observations_tidy.csv") |> 
  dplyr::select(scientific_name, taxon_id) |> 
  distinct

names2Check <- unique(df_tax$scientific_name)

# gna_data_sources() |> View # Options for data sources
names_ver <- gna_verifier(
  names2Check)
fwrite(names_ver, file = "data/iNat_observations_taxNames.csv", row.names = F)

for(searchName in names_ver$matchedCanonicalSimple) {
  
  if(exists("o")) rm(o)
  if(exists("o_df")) rm(o_df)
  
  if(is.na(searchName)) next
  # Load previously-downloaded data
  taxPath <- "data/iNat_observations_taxInfo.csv"
  if(file.exists(taxPath)) {
    taxData <- fread(taxPath, fill = T)
    # Skip if search term is already in the taxData csv
    if(searchName %in% taxData$search) next
  }
  # Download from itis.
  try({
    o <- classification(searchName, db = 'itis', rows = 1)
   # o <- classification(searchName, db = 'gbif', rows = 1)
  })
  if(!exists("o")) { # If download failed, sleep and move on.
    Sys.sleep(60)
    next
  }
  # Check download worked.
  if(names(o) == searchName & nrow(o[[1]]) != 0) {
    # Write outputs
    o_df <- as.data.frame(o[[1]]) |> 
      # pivot_wider(names_from = rank, values_from = c(name, id), names_glue = "{rank}_{.value}") |> 
      mutate(db = attr(o, "db"), search = searchName)
    fwrite(o_df, taxPath, append = T, row.names = F)
    Sys.sleep(1)
  } else {
    next
  }
}

taxInfo <- fread("data/iNat_observations_taxInfo.csv") |>
  distinct |>
  pivot_wider(
    names_from = rank,
    values_from = c(name, id),
    names_glue = "{rank}_{.value}"
  ) |>
  dplyr::select(search, db, ends_with("_name")) |>
  dplyr::filter(!is.na(kingdom_name), phylum_name == "Chordata") |>
  dplyr::select(
    c(
      "search",
      "db",
      "order_name",
      "superfamily_name",
      "family_name",
      "subfamily_name",
      "tribe_name",
      "genus_name",
      "subgenus_name",
      "species_name"
      
    )
  )

fwrite(taxInfo, "data/iNat_observations_taxTree.csv")
