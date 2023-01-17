create_path_map <- function(det_data, dep_data, xlim = c(-100, -60), ylim = c(25, 60), 
                            legend_pos = c(0.25, 0.2), backdrop_sf = NULL,
                            delbay_paths = NULL,
                            backdrop_label_sf = NULL,
                            path_size = 0.75, path_color = "grey60",
                            scale_loc = "bl",
                            map_crs = NULL, country_label_sf = NULL) {
  
  if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc", quiet = TRUE)
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) install.packages("rnaturalearth", quiet = TRUE)
  if (!requireNamespace("ggspatial", quietly = TRUE)) install.packages("ggspatial", quiet = TRUE)
  
  deps <- select(det_data, motusTagID, state, lat = tagDepLat, lon = tagDepLon, dt = tagDeployStart) %>%
    distinct() %>%
    mutate(recvDeployID = NA_integer_, recvSiteName = NA_character_)
  all_deps <- select(dep_data, motusTagID, state) %>%
    group_by(state) %>% tally(name = "dep_n")
  tag_loc <- group_by(deps, state) %>%
    summarize(lat = mean(lat), lon = mean(lon), det_n = n()) %>%
    left_join(all_deps, by = "state") %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  tag_labs <- tibble(lon = c(-79, -78),
                    lat = c(28, 32)) %>%
    mutate(label = sprintf("Tags (%s):\n• deployed (n = %i)\n• detected (n = %i)", tag_loc$state, tag_loc$dep_n, tag_loc$det_n)) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  tag_lines <- st_sfc(mapply(function(a,b) {st_cast(st_union(a,b),"LINESTRING")}, 
                                st_geometry(tag_labs), st_geometry(tag_loc), SIMPLIFY=FALSE)) %>%
    st_as_sf(crs = 4326)
  
  dets <- select(det_data, motusTagID, recvDeployID, recvSiteName, lat = recvDeployLat, lon = recvDeployLon, dt = first_det)

  # Restrict departures to those with detections
  deps <- filter(deps, motusTagID %in% unique(dets$motusTagID))
  
  path_dat <- bind_rows(dets, deps) %>%
    group_by(motusTagID) %>%
    arrange(motusTagID, dt) %>%
    ungroup()
  
  if (!is.null(delbay_paths)) {
    stopifnot(all(c("motusTagID", "delbay_simplified") %in% names(delbay_paths)))
    delbay_paths <- select(delbay_paths, motusTagID, delbay_simplified)
    path_dat <- left_join(path_dat, delbay_paths, by = "motusTagID")
  }
  
  # Gray palette for points
  gray_cols <- function(vals) {
    # Set palette limits
    vals <- scale_vec(vals, scale = c(0.15, 0.85))
    gray.colors(length(vals), rev = TRUE)
  }
  
  station_n <- path_dat %>%
    filter(!is.na(recvDeployID)) %>%
    group_by(recvSiteName) %>%
    summarize(n_indivs = n_distinct(motusTagID), .groups = "drop")
  
  station_locs <- select(path_dat, recvSiteName, lat, lon) %>%
    distinct() %>%
    group_by(recvSiteName) %>%
    # Keep only one location in case of slight variation in position among deployments
    slice(1) %>%
    ungroup()
  
  station_n <- left_join(station_n, station_locs, by = "recvSiteName") %>%
    mutate(pct_indivs = round(n_indivs / n_distinct(path_dat$motusTagID) * 100, 1))
  
  # Gather active receivers in network
  recvs <- get_motus_recvs(recv_path = "Data/receiver-deployments.csv", tz = attr(det_data$first_det, "tzone")) %>%
    filter(tsStart < max(det_data$last_det),
           (is.na(tsEnd) | tsEnd > min(det_data$tagDeployStart)))
  
  # Gather geographic data
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  na <- rnaturalearth::ne_states(country = c("united states of america", "canada", "mexico"), returnclass = "sf")
  if (!is.null(map_crs)) {
    na <- st_transform(na, map_crs)
    # Convert limits to new crs
    lims <- data.frame(lon = xlim, lat = ylim) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      st_transform(map_crs) %>%
      st_coordinates() %>%
      as.data.frame()
    xlim <- lims$X
    ylim <- lims$Y
  }
  
  recvs_sf <- st_as_sf(recvs, coords = c("longitude", "latitude"), crs = 4326)
  station_sf <- st_as_sf(station_n, coords = c("lon", "lat"), crs = 4326)
  det_paths_sf <- st_as_sf(path_dat, coords = c("lon", "lat"), crs = 4326)
  if (!is.null(delbay_paths))
    det_paths_sf <- group_by(det_paths_sf, motusTagID, delbay_simplified)
  else
    det_paths_sf <- group_by(det_paths_sf, motusTagID)
  det_paths_sf <- det_paths_sf %>%
    summarize(do_union=FALSE, .groups = "keep") %>%
    st_cast("LINESTRING") %>%
    ungroup()

  species <- unique(det_data$spp)
  binned_scale_lab <- paste0("% ", if (length(species) == 1) species else "individuals",
                             "\ndetected (n = ", n_distinct(path_dat$motusTagID), ")")
  
  p <- ggplot(data = na) +
    geom_sf(fill = "gray90", color = "grey50") +
    geom_sf(data = world, fill = NA, color = "grey50", lwd = 1) +
    geom_sf(data = recvs_sf, pch = 21, color = "grey5", 
               fill = "white", stroke = 0.5, size = 1) +
    geom_sf(data = tag_lines, lwd = 0.5) +
    geom_sf_label(data = tag_labs, aes(label = label), hjust = 0) 
  
  if (!is.null(delbay_paths)) {
    p <- p +
      geom_sf(data = det_paths_sf, aes(color = delbay_simplified), alpha = 0.65, size = path_size) +
      scale_color_manual("Stopover use of\nDelaware Bay", na.translate = FALSE,
                         values = c("#1b9e77", "#7570b3", "#d95f02"),
                         guide="legend") +
      guides(color = guide_legend(order = 1))
  } else {
    p <- p + geom_sf(data = det_paths_sf, color = path_color, alpha = 0.65, size = path_size)
  }
  
  p <- p + 
    geom_sf(data = station_sf, aes(fill = pct_indivs, size = pct_indivs), shape = 21, stroke = 0.5, 
            color = "black")
  
  if (!is.null(backdrop_sf)) {
    lab_lines <- tibble(lat = c(38.4, 38.75, 43, 43.684),
                        lon = c(-70.5, -74.708, -94, -89.703),
                        grp = c("DB", "DB", "GL", "GL")) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      group_by(grp) %>% summarize(grp = unique(grp), do_union=FALSE) %>% 
      st_cast("LINESTRING") %>% st_transform(st_crs(backdrop_sf))
    
    p <- p +
      geom_sf(data = lab_lines, lwd = 0.5) +
      geom_sf(data = backdrop_sf, color = "black", lwd = 1, fill = NA)
  }
  
  if (!is.null(backdrop_label_sf)) {
    p <- p +
      geom_sf_label(data = backdrop_label_sf, aes(label = label)) 
  }
    
  if (!is.null(country_label_sf)) {
    p <- p +
      geom_sf_text(data = country_label_sf, aes(label = label), size = 6,
                   fontface = "bold") 
  }
  
  p <- p +   
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_fill_binned(n.breaks = 5, palette = gray_cols) +
    scale_size_binned(n.breaks = 5, range = c(1, 6)) +
    guides(fill = guide_bins(binned_scale_lab, direction = "horizontal", order = 2),
           size = guide_bins(binned_scale_lab, direction = "horizontal", order = 2)) +
    
    labs(x = NULL, y = NULL) +
    ggspatial::annotation_scale(location = scale_loc, width_hint = 0.55, text_cex = 0.9) +
    theme_bw(base_size = 14) + 
    theme(legend.background = element_rect(color = "black"),
          legend.position = legend_pos)
  
  return(p)
}

