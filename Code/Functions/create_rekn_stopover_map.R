create_rekn_stopover_map <- function(stop_data, xlim = c(-85, -65), ylim = c(30, 51.5), 
                                     legend_pos = c(0.75, 0.15),
                                     tag_seg_end = NULL, tag_seg_lab = "Tags\ndeployed",
                                     stop_buffer = 50, map_crs = NULL,
                                     lgd_breaks = NULL, lgd_vir_end = NULL,
                                     bg_recv_start = "2017-04-29", bg_recv_end = "2019-06-15") {

  if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc", quiet = TRUE)
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) install.packages("rnaturalearth", quiet = TRUE)
  if (!requireNamespace("ggspatial", quietly = TRUE)) install.packages("ggspatial", quiet = TRUE)
  
  rekn_all_stops <- stop_data %>%
    st_as_sf(coords = c("recvDeployLon", "recvDeployLat"), crs = 4326) %>%
    buffer_ll_pts(stop_buffer)
  
  n_stop_recvs <- unlist(strsplit(paste(c(rekn_all_stops$recvDeployID, na.omit(rekn_all_stops$adj_recv_grp)), collapse = ","), ",")) %>%
    n_distinct()
  n_rekn_stopping <- n_distinct(rekn_all_stops$motusTagID)
  message("Stopovers of ", n_rekn_stopping, " individuals will be dispayed representing ", n_stop_recvs, " receivers.")
  
  rekn_all_stops <- rekn_all_stops %>%  group_by(motusTagID) %>%
    summarize() %>%
    st_buffer(0) %>% # Workaround for self-intersections
    st_intersection() %>%
    filter(st_is(., c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION")))
    
  bins <- seq(min(rekn_all_stops$n.overlaps), max(rekn_all_stops$n.overlaps))
  rekn_all_stops <- mutate(rekn_all_stops,
                           n.overlaps.f = factor(n.overlaps, 
                                               levels = if (is.null(lgd_breaks)) bins else lgd_breaks))

  # Gather active receivers in network
  recvs <- get_motus_recvs(recv_path = "Data/receiver-deployments.csv", tz = attr(stop_data$first_det, "tzone")) %>%
    filter(tsStart < lubridate::ymd(bg_recv_end),
           (is.na(tsEnd) | tsEnd > lubridate::ymd(bg_recv_start)))
  
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

  p <- 
    ggplot(data = na) +
    geom_sf(fill = "gray95", color = "black") +
    geom_sf(data = world, fill = NA, color = "black", lwd = 1) +
    geom_sf(data = rekn_all_stops, aes(fill = n.overlaps.f), shape = 21, color = NA) +
    geom_sf(data = recvs_sf, shape = 21, color = "grey5", 
            fill = "white", stroke = 0.25, size = 1.5) +
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_fill_viridis_d(paste0("# individuals\nstopping over"),
                         guide = guide_legend(label.position = "bottom",
                                              title.position = "top",
                                              direction = "horizontal", nrow = 1),
                         end = if (is.null(lgd_vir_end)) 1 else lgd_vir_end,
                         drop = FALSE) +
    ggspatial::annotation_scale(location = "br", width_hint = 0.5, text_cex = 0.9) +
    # guides(fill = guide_bins(title.position = "top", show.limits = TRUE)) +
    theme_bw(base_size = 14) + 
    theme(legend.position = legend_pos,
          # legend.direction = "horizontal",
          legend.spacing.x = unit(0, "cm"),
          legend.box.margin = margin(3, 5, 2, 2),
          legend.box.background = element_rect(color = "black"),
          legend.key = element_rect(color = "black", size = 0.5))
  return(p)
}
