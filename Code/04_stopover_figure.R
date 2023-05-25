source("Code/Functions/create_rekn_stopover_map.R")

solo_stops <- filter(rekn_det, stop_solo, stop_solo_h <= 4 * 24)
adj_stops <- filter(rekn_det, stop_adj, stop_adj_h <= 6 * 24)
all_stops <- bind_rows(solo_stops, adj_stops)

# Filter SC "stopovers" of REKN tagged in SC
all_stops <- filter(all_stops, 
                    !(stop_id %in% c(1, 27, 32)),
                    age == "ASY")

# CRS for display
map_crs = "+proj=lcc +lat_1=35 +lat_2=45 +lat_0=40 +lon_0=-80 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs"


#Figure 4
# Stopover map
stops <- create_rekn_stopover_map(all_stops, map_crs = map_crs)
ggsave("Output/Fig4_rekn_all_stopovers.png", dpi = 600, height = 8, width = 5)
ggsave("Output/Fig4_rekn_all_stopovers.pdf", height = 8, width = 5)
