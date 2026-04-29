
library(tidyverse)
theme_set(theme_classic())

# Date-related functions and objects-------------------------------------------------

yDay_to_dayMonth <- function(yday) {
  library(lubridate)
  dateOut <- as.Date(as.numeric(yday)-1, format = "%j", origin = "1.1.2018")
  format(dateOut, "%d %b")
}

yDay_to_monthDay <- function(yday) {
  library(lubridate)
  dateOut <- as.Date(as.numeric(yday)-1, format = "%j", origin = "1.1.2018")
  format(dateOut, "%b %d")
}

yDay_to_Month <- function(yday) {
  library(lubridate)
  dateOut <- as.Date(as.numeric(yday)-1, format = "%j", origin = "1.1.2018")
  format(dateOut, "%b")
}

monthDayYear_to_yday <- function(monthDayYear) {
  yday(mdy(monthDayYear))
}

monthFirsts <- paste(month.abb,1,2020, sep ="-")

# Projection strings used across spatial scripts
proj.wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"
myproj     <- "+proj=eqearth +lon_0=0 +datum=WGS84 +units=km +no_defs"
doy_labels <- format(mdy(monthFirsts), "%b")
# 
# dateGrid <- expand.grid(month.abb,c(1, 8, 15, 22, 29),2020) %>% 
#   mutate(
#     date = mdy(paste(Var1, Var2, Var3)),
#     yday = yday(date),
#     label = case_when(Var2 == 1 ~ Var1, TRUE ~ "")
#     ) %>% 
#   dplyr::filter(!is.na(date)) %>% 
#   arrange(date)


# iNaturalist download functions -----------------------------------------------



# Notes from CJ:
# These functions are directly adapted from rinat, but several contain original features.
# Those features were developed for an iNat user project (Grady et al., in review),
# New features I added in the course of this study:
# - Search arguments in searchBuilder for project, field, and observation ID.
# I would recommend citing this ms, rinat, and Grady et al. to thoroughly
# reference the source of these functions if they were used further.


emptyDF <- structure(list(scientific_name = character(0), datetime = character(0), 
                          description = logical(0), place_guess = character(0), latitude = numeric(0), 
                          longitude = numeric(0), tag_list = logical(0), common_name = character(0), 
                          url = character(0), image_url = character(0), user_login = character(0), 
                          id = integer(0), species_guess = character(0), iconic_taxon_name = character(0), 
                          taxon_id = integer(0), num_identification_agreements = integer(0), 
                          num_identification_disagreements = integer(0), observed_on_string = character(0), 
                          observed_on = character(0), time_observed_at = character(0), 
                          time_zone = character(0), positional_accuracy = integer(0), 
                          public_positional_accuracy = integer(0), geoprivacy = logical(0), 
                          taxon_geoprivacy = character(0), coordinates_obscured = character(0), 
                          positioning_method = character(0), positioning_device = character(0), 
                          user_id = integer(0), user_name = character(0), created_at = character(0), 
                          updated_at = character(0), quality_grade = character(0), 
                          license = character(0), sound_url = logical(0), oauth_application_id = integer(0), 
                          captive_cultivated = character(0)), row.names = integer(0), class = "data.frame")

# checkCons <- function() {
#   if (!curl::has_internet()) {
#     message("No Internet connection.")
#     return(invisible(NULL))
#   }
#   base_url <- "http://www.inaturalist.org/"
#   if (httr::http_error(base_url)) {
#     message("iNaturalist API is unavailable.")
#     return(invisible(NULL))
#   }
# }

searchBuilder <- function(query = NULL, taxon_name = NULL, taxon_id = NULL, place_id = NULL, 
                          quality = NULL, geo = NULL, annotation = NULL, 
                          # Observation dates.
                          d1 = NULL, d2 = NULL, 
                          year = NULL, month = NULL, day = NULL,
                          # Dates added to iNaturalist.
                          created_d1= NULL, created_d2= NULL, 
                          # Spatial bounds.
                          bounds = NULL,
                          # Projects
                          projects = NULL,
                          # Field
                          field = NULL,
                          # Observation ID
                          id = NULL,
                          # User ID
                          user_id = NULL
                          ) {
  
  arg_list <- list(query, taxon_name, taxon_id, place_id, quality, 
                   geo, 
                   year, month, day,
                   bounds, d1, d2, created_d1, created_d2, projects, field, id,
                   user_id)
  arg_vals <- lapply(arg_list, is.null)
  if (all(unlist(arg_vals))) {
    stop("All search parameters NULL. Please provide at least one.")
  }
  base_url <- "http://www.inaturalist.org/"
  search <- ""
  if (!is.null(query)) {
    search <- paste0(search, "&q=", gsub(" ", "+", query))
  }
  if (!is.null(quality)) {
    if (!sum(grepl(quality, c("casual", "research")))) {
      stop("Please enter a valid quality flag, 'casual' or 'research'.")
    }
    search <- paste0(search, "&quality_grade=", quality)
  }
  if (!is.null(taxon_name)) {
    search <- paste0(search, "&taxon_name=", gsub(" ", "+", 
                                                  taxon_name))
  }
  if (!is.null(taxon_id)) {
    search <- paste0(search, "&taxon_id=", gsub(" ", "+", 
                                                taxon_id))
  }
  if (!is.null(place_id)) {
    search <- paste0(search, "&place_id=", gsub(" ", "+", 
                                                place_id))
  }
  if (!is.null(geo) && geo) {
    search <- paste0(search, "&has[]=geo")
  }
  if (!is.null(annotation)) {
    if (length(annotation) != 2) {
      stop("annotation needs to be a vector of length 2.")
    }
    annotation <- as.character(annotation)
    if (grepl("\\D", annotation[1])) {
      stop("The annotation's term ID can only contain digits.")
    }
    if (grepl("\\D", annotation[2])) {
      stop("The annotation's value ID can only contain digits.")
    }
    search <- paste0(search, "&term_id=", annotation[1])
    search <- paste0(search, "&term_value_id=", annotation[2])
  }
  if (!is.null(year)) {
    if (length(year) > 1) {
      stop("You can only filter results by one year, please enter only one value for year.")
    }
    search <- paste0(search, "&year=", year)
  }
  if (!is.null(month)) {
    month <- as.numeric(month)
    if (is.na(month)) {
      stop("Please enter a month as a number between 1 and 12, not as a word.")
    }
    if (length(month) > 1) {
      stop("You can only filter results by one month, please enter only one value for month.")
    }
    if (month < 1 || month > 12) {
      stop("Please enter a valid month between 1 and 12")
    }
    search <- paste0(search, "&month=", month)
  }
  if (!is.null(day)) {
    day <- as.numeric(day)
    if (is.na(day)) {
      stop("Please enter a day as a number between 1 and 31, not as a word.")
    }
    if (length(day) > 1) {
      stop("You can only filter results by one day, please enter only one value for day.")
    }
    if (day < 1 || day > 31) {
      stop("Please enter a valid day between 1 and 31")
    }
    search <- paste0(search, "&day=", day)
  }
  
  if(!is.null(d1)) {
    search <- paste0(search, "&d1=", d1)
  }
  if(!is.null(d2)) {
    search <- paste0(search, "&d2=", d2)
  }
  
  
  if(!is.null(created_d1)) {
    search <- paste0(search, "&created_d1=", created_d1)
  }
  if(!is.null(created_d2)) {
    search <- paste0(search, "&created_d2=", created_d2)
  }
  if (!is.null(bounds)) {
    if (inherits(bounds, "sf")) {
      bounds_prep <- sf::st_bbox(bounds)
      bounds <- c(swlat = bounds_prep[2], swlng = bounds_prep[1], 
                  nelat = bounds_prep[4], nelng = bounds_prep[3])
    }
    if (inherits(bounds, "Spatial")) {
      bounds_prep <- sp::bbox(bounds)
      bounds <- c(swlat = bounds_prep[2, 1], swlng = bounds_prep[1, 
                                                                 1], nelat = bounds_prep[2, 2], nelng = bounds_prep[1, 
                                                                                                                    2])
    }
    if (length(bounds) != 4) {
      stop("Bounding box specifications must have 4 coordinates.")
    }
    search <- paste0(search, "&swlat=", bounds[1], "&swlng=", 
                     bounds[2], "&nelat=", bounds[3], "&nelng=", bounds[4])
  }
  
  if(!is.null(projects)) {
    
    projects2 <- paste(projects, collapse = "%2C+")
    stopifnot(length(projects2) == 1)
    
    search <- paste0(search, "&projects%5B%5D=", projects2)
  }
  
  if(!is.null(field)) {
    
    field <- gsub(" ", "%20", field)
    stopifnot(length(field) == 1)
    
    search <- paste0(search, "&field:", field)
  }
  
  if(!is.null(id)) {
    id2 <- paste0("&id=", id, collapse = "")
    search <- paste0(search, id2)
  }
  
  if(!is.null(user_id)) {
    user_id2 <- paste0("&user_id=", user_id, collapse = "")
    search <- paste0(search, user_id2)
  }

  q_path <- "observations.csv"
  ping_path <- "observations.json"
  
  return(list(base_url = base_url, search = search))
  
}


# search_out <- searchBuilder(query = "Mule Deer", bounds = c(38.44047, -125, 40.86652, -121.837),
#                             d1 = "2020-10-01", d2 = "2024-10-09")


howManyResults <- function(search_out, perpage = 1) {
  if (!curl::has_internet()) {
    message("No Internet connection.")
    return(invisible(NULL))
  }
  base_url <- "http://www.inaturalist.org/"
  if (httr::http_error(base_url)) {
    message("iNaturalist API is unavailable.")
    return(invisible(NULL))
  }
  
  # Add page = 1 for the purposes of the ping.
  page_query <- paste0(search_out$search, "&per_page=", perpage, "&page=1")
  
  ping <- httr::GET(search_out$base_url, path = "observations.json", query = page_query)
  total_res <- as.numeric(ping$headers$`x-total-entries`)
  Sys.sleep(1)
  return(total_res)
}


# howManyResults(search_out)
# howManyResults(search_out, perpage = 200)
# microbenchmark::microbenchmark(
#   times = 1,
#   howManyResults(search_out, perpage = 1),
#   howManyResults(search_out, perpage = 10), 
#   howManyResults(search_out, perpage = 200)
# )
# # Reducing the perpage value to 1 from the default of 200 decreases the ping time substantially while returning the same result.


checkResults <- function(search_out) {
  
  total_res <- howManyResults(search_out)
  
  if (total_res == 0) {
    stop("Your search returned zero results. Either your species of interest has no records or you entered an invalid search.")
  }
  else if (total_res >= 2e+05) {
    stop("Your search returned too many results, please consider breaking it up into smaller chunks by year or month.")
  }
  else if (!is.null(bounds) && total_res >= 1e+05) {
    stop("Your search returned too many results, please consider breaking it up into smaller chunks by year or month.")
  }
  
  return(total_res)
  
}

# checkResults(search_out)

downloadResults <- function(search_out, expectedResults = NULL, perpage = 200, sleepHowLong = 1) {
  
  if(is.null(expectedResults)) {
    expectedResults <- howManyResults(search_out)
  }
  data_out <- emptyDF
  if(expectedResults > 0) {
    for (i in 1:ceiling(expectedResults/200)) {
      page_query <- paste0(search_out$search, "&per_page=", perpage, "&page=", i)
      data <- httr::GET(search_out$base_url, path = "observations.csv", query = page_query)
      data2 <- c()
      counter <- 0
      while(length(data2) == 0 & counter < 5) {
        try({ data2 <- rinat:::inat_handle(data) })
        if(length(data2) == 0) { Sys.sleep(10) }
        counter <- counter + 1
      }
      data3 <- read.csv(textConnection(data2), stringsAsFactors = FALSE)
      data_out <- rbind( data_out, data3)
      Sys.sleep(sleepHowLong)
    }
  }
  if(expectedResults != nrow(data_out)) {data_out <- distinct(data_out)}
  # writeLines(paste("expecting = ", expectedResults, "Have = ", nrow(data_out)))
  stopifnot("Expected rows does not match downloaded rows" = expectedResults == nrow(data_out) )
  
  return(data_out)
}

# downloadResults(search_out)
# downloadResults(search_out, expectedResults = 478)
# microbenchmark::microbenchmark(
#   times = 1,
#   downloadResults(search_out), 
#   downloadResults(search_out, expectedResults = 478)
# )
