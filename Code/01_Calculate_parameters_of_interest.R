pacman::p_load(ggmap)
# Google API key required for creating nearby station maps
# You'll need your own and to save it as an environmental variable named `goog_key`
# or subsitute it for the `Sys.getenv("goog_key")` component of the next line
register_google(Sys.getenv("goog_key"))

source("Code/Functions/find_stations_near_deploy.R")
source("Code/Functions/create_nearby_station_maps.R")
source("Code/Functions/designate_home_stations.R")
source("Code/Functions/identify_stopovers.R")
source('Code/Functions/retrieve_detection_wind_vectors.R')
source("Code/Functions/calculate_ground_speed.R")
source("Code/Functions/retrieve_path_wind_components.R")
source("Code/Functions/create_detection_maps.R")

# Summarize detections
if (!file.exists("Data/Derived/rekn_detections.rds")) {
  # Purely academic, there were no nearby "home" stations for northbound movement
  nearby_stns <- find_stations_near_deploy(rekn, max_dist = 20, time_window = c(0, 60))
  
  rekn <- rekn %>%
    designate_home_stations(nearby_stns) %>%
    group_by(tagDeployID) %>%
    mutate(only_home = all(home) & n_distinct(recvSiteName) == 1,
           only_away = all(!home)) %>%
    ungroup()

  rekn_det <- rekn %>%
    group_by(motusTagID, spp, age, state,
             tagDeployID, tagDeployStart, tagDepLat, tagDepLon, 
             recvSiteName, recvDeployID, recvDeployLat, recvDeployLon, home, only_away) %>%
    summarize(first_det = min(ts),
              last_det = max(ts), .groups = "drop") %>%
    arrange(motusTagID, first_det) %>%
    mutate(det_wind_h = as.numeric(difftime(last_det, first_det, units = "hours"))) %>%
    identify_stopovers(solo_h = 4, adj_h = 6) %>%
    retrieve_detection_wind_vectors()
  saveRDS(rekn_det, file = "Data/Derived/rekn_detections.rds")
} else rekn_det <- readRDS("Data/Derived/rekn_detections.rds")

# Update detection maps to include *SPRING* stopovers
# Note, the code above identifies stopovers regardless of season
# That is, a bird visiting the same site in the spring and fall would be
# considered a stopover if there were no intervening detections (not correct obviously)
has_stop <- filter(rekn_det, stop_solo | stop_adj) %>% pull(motusTagID) %>% unique()
if (!file.exists("Output/REKN_detection_maps_w_stopovers.pdf"))
  create_detection_maps(det_data = filter(rekn, motusTagID %in% has_stop),
                        outfile = "Output/REKN_detection_maps_w_stopovers.pdf",
                        lat_rng = c(27, 58), lon_rng = c(-105, -50),
                        add_stops = TRUE, stopover_dat = rekn_det)

# Estimate orthodromic flight (ground) speed where appropriate
if (!file.exists("Data/Derived/rekn_speed_calculations.rds")) {
  rekn_speeds <- calculate_ground_speed(rekn_det, min_stn_sep = 150, max_time_sep = 18) %>%
    retrieve_path_wind_components()
  saveRDS(rekn_speeds, file = "Data/Derived/rekn_speed_calculations.rds")
} 
if (!file.exists("Data/Derived/rekn_speed_checks.csv")) {
  rekn_speeds %>% 
    select(-traj) %>%
    mutate(keep = TRUE) %>%
    readr::write_csv("Data/Derived/rekn_speed_checks.csv")
}

# Prepare csv for manual assignment of departure windows
if (!file.exists("Data/Derived/rekn_departures.csv")) {
  # Note that departure applies only to SC/GA coast
  # That is, REKN tagged in FL are assigned departure dates only 
  # if they were detected along the SC/GA coast
  rekn_dep <- rekn_det %>%
    group_by(motusTagID, spp, age, tagDeployID, tagDeployStart, tagDepLat, tagDepLon) %>%
    summarize(min_depart = "YYYY-MM-DD",
              max_depart = "YYYY-MM-DD",
              departure_notes = NA_character_)
  readr::write_csv(rekn_dep, "Data/Derived/rekn_departures.csv")
}
