get_motusTagIDs <- function(data) {
  out <- pull(data, motusTagID) %>%
    unique() %>% sort()
  return(out)
}

get_motus_recvs <- function(recv_path = NULL, tz = "GMT") {
  if (is.null(recv_path))
    recvs <- readr::read_csv("https://motus.org/data/downloads/api-proxy/receivers/deployments?fmt=csv",
                           col_types = "--c-d--c--dddd-l--")
  else 
    recvs <- readr::read_csv(recv_path, col_types = "--c-d--c--dddd-l--")
  recvs <- recvs %>%
    mutate(tsStart = as_datetime(tsStart, tz = "UTC"),
           tsEnd = as_datetime(tsEnd, tz = "UTC")) %>%
    filter(!isMobile,
           !is.na(latitude),
           !is.na(longitude))
  attr(recvs$tsStart, "tzone") <- attr(recvs$tsEnd, "tzone") <- tz
  return(recvs)
}

get_deployment_info <- function(data) {
  out <- select(data,
                motusTagID, tagDeployID, tagDepLat, tagDepLon, tagDeployStart) %>%
    distinct()
  return(out)
}

# Create geodesic buffer around points
geo_buffer <- function(data, coords = c("tagDepLon", "tagDepLat"), dist_km = 10) {
  dep_pts <- st_as_sf(data, coords = coords, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
  n_pts <- nrow(dep_pts)
  dg <- seq(from = 0, to = 360, by = 5)
  pt_buffs <- as.data.frame(geosphere::destPoint(p = as(dep_pts, "Spatial"), b = rep(dg, each = n_pts), d = dist_km * 1000))
  poly_ids <- factor(rep(dep_pts$tagDeployID, length(dg)), levels = sort(dep_pts$tagDeployID))
  lst <- split(pt_buffs, f = poly_ids)
  lst <- lapply(lst, as.matrix)
  polys <- lapply(lst, function(x) st_polygon(list(x)))
  mp <- st_sf(tagDeployID = as.integer(levels(poly_ids)),
              geometry = st_sfc(polys),
              crs = st_crs(dep_pts))
  return(mp)
}

get_active_recvs <- function(ts, recvs, time_window = c(0, 0)) {
  ts <- ts + as.difftime(time_window, units = "days")
  out <- filter(recvs,
                (is.na(tsEnd) | tsEnd > ts[1]),
                tsStart < ts[2])
  return(out)
}

find_deploy_recvs <- function(deps, recvs, dist_km = 10, coords = c("tagDepLon", "tagDepLat"), time_window = c(0, 0)) {
  recvs_sf <- st_as_sf(recvs, coords = c("longitude", "latitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
  deps_sf <- st_as_sf(deps, coords = coords, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
  out <- pbapply::pblapply(seq_len(nrow(deps)), function(i) {
    dep_i <- deps_sf[i, ]
    dep_ts <- pull(dep_i, tagDeployStart)
    
    recv_i <- get_active_recvs(dep_ts, recvs_sf, time_window = time_window)

    dists <- as.numeric(st_distance(dep_i, recv_i))
    recv_i <- mutate(recv_i,
                     tagDeployID = pull(dep_i, tagDeployID),
                     distkm = dists/1000) %>%
      filter(distkm <= dist_km) %>% as.data.frame() %>%
      select(tagDeployID, recvDeployID, distkm)
    if (nrow(recv_i) == 0) {
      message("No nearby stations for ", dep_i$motusTagID, ". Closest active station: ", round(min(dists)/1000, 1), " km")
      recv_i <- tibble::add_row(recv_i,
                                tagDeployID = pull(dep_i, tagDeployID),
                                recvDeployID = NA_integer_,
                                distkm = NA)
    }
    recv_i
  })
  out <- bind_rows(out)
  return(out)
}

get_uniq_recvs <- function(det_data) {
  out <- ungroup(det_data) %>%
    select(recvSiteName, recvDeployID, recvDeployLat, recvDeployLon) %>%
    distinct() %>%
    arrange(recvDeployID)
  return(out)
}

dates_by_month <- function(dates) {
  months <- sort(unique(lubridate::month(dates)))
  out <- sapply(months, function(m) {
    tmp_m <- sort(unique(dates[lubridate::month(dates) == m]))
    tmp_yrs <- sort(unique(lubridate::year(tmp_m)))
    tmp_out <- sapply(tmp_yrs, function(y) {
      tmp_y <- sort(unique(tmp_m[lubridate::year(tmp_m) == y]))
      tmp_y_days <- sort(unique(lubridate::day(tmp_y)))
      paste(split_long_vec(tmp_y_days), month.abb[m], y)
    })
  })
  out <- paste(unlist(out), collapse = "\n")
  out
}

split_long_vec <- function(vec, n = 8) {
  out <- split(vec, ceiling(seq_along(vec)/n))
  out <- sapply(out, function(i) paste(i, collapse = ","))
  paste(out, collapse = "\n")
}

median_wtd <- function(x, w) { 
  w <- w[order(x)]
  x <- x[order(x)]
  x[which.min(abs(stats::filter(c(0,cumsum(w)/sum(w)), c(.5,.5), sides=1)[-1] - 0.5))]
}

theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
	  plot.subtitle = element_text(size = base_size*1, color = "white"),
      plot.margin = unit(rep(1, 4), "lines")
    )
}

# Pilfered and adapted from Stack Overflow
# https://gis.stackexchange.com/a/121539/35273
buffer_ll_pts <- function(pt_ll_sf, buff_km = 15) {
  pt_crs <- st_crs(pt_ll_sf)
  n_pts <- nrow(pt_ll_sf)
  buff_list <- lapply(seq(n_pts), function(i) {
    i_coords <- st_coordinates(pt_ll_sf[i, ])
    aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                    i_coords[2], i_coords[1])
    i_proj <- st_transform(pt_ll_sf[i, ], aeqd)
    i_buff <- st_buffer(i_proj, dist = buff_km * 1000)
    i_buff <- st_transform(i_buff, pt_crs)
  })
  buffs <- do.call(rbind, buff_list)
  return(buffs)
}

st_points_to_path <- function(pt_df, coords = c("lon", "lat"), 
                              time_val = "datetime", crs = 4326) {
  out <- pt_df %>%
    st_as_sf(pt_df, coords = c("lon", "lat"), crs = 4326) %>%
    group_by(motusTagID) %>%
    arrange_at(time_val)
  summarize(do_union=FALSE, .groups = "keep") %>%
    st_cast("LINESTRING") %>%
    ungroup()
}

id_intervals <- function(intervals, current_id) {
  as.integer(cumsum(!lubridate::int_overlaps(intervals, dplyr::lag(intervals, default = dplyr::first(intervals)))) + 1 + current_id)
}


# Make date sequence
date_seq <- function(st = "2018-01-01", end = "2018-06-30", by = 1, format = "%d %b", ret_dates = FALSE) {
  out <- seq.Date(as_date(st), as_date(end), by = by)
  if (ret_dates)
    return(out)
  else
    return(format(out, format = format))
}

# Pilfered from plyr package
round_any <- function (x, accuracy, f = round) {
  f(x/accuracy) * accuracy
}

# Pilfered from stack overflow
# https://stackoverflow.com/a/45614547/2726564
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
