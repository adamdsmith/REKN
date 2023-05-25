source("Code/Functions/create_path_map.R")
source("Code/Functions/utils.R")
source("Code/Functions/plot_capture_weights.R")

rekn_det <- readRDS("Data/Derived/rekn_detections.rds") %>%
  mutate(doy = yday(first_det)) %>%
  filter(doy <= 166)

# Age distribution of all tagged birds
table(rekn_dep$age)

# Age distribution of tags detected at least once
detected_tags <- unique(rekn_det$motusTagID)
rekn_det %>% 
  select(motusTagID, age) %>% distinct() %>%
  with(., table(age, useNA = "ifany"))

# Age distribution of tags detected N of 34 degrees
rekn_det %>% filter(recvDeployLat > 34) %>% 
  select(motusTagID, age) %>% distinct() %>%
  with(., table(age, useNA = "ifany"))
rekn_w_info <- pull(filter(rekn_det, recvDeployLat > 34), motusTagID) %>% unique()

# FYI only; distribution of capture weights by detection/stopover 
plot_capture_weights(rekn_dep, detected_tags, has_stop)
ggsave("Output/rekn_capture_weights.png", dpi = 600, height = 5, width = 5)

# Figure 2
map_crs = "+proj=lcc +lat_1=30 +lat_2=55 +lat_0=42.5 +lon_0=-85 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs"
gl <- st_read("Resources/geodata/greatlakes_subbasins.shp", quiet = TRUE) %>%
  summarize()
db <- st_read("Resources/geodata/delaware-bay_HUC02040204.shp", quiet = TRUE) %>%
  st_buffer(30000) %>% summarize() %>% st_transform(st_crs(gl))
backdrop <- rbind(gl, db)
backdrop_labs <- tibble(lat = c(38.4, 43, 59, 53.75), 
                        lon = c(-70.5, -94, -85, -77), 
                        label = c("Delaware Bay", "Great Lakes\nBasin", "Hudson Bay", "James Bay")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
country_labs <- tibble(lat = c(39.45, 54.71),
                       lon = c(-98.34, -100.5),
                       label = c("United States", "Canada")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
create_path_map(filter(rekn_det, motusTagID %in% rekn_w_info, age == "ASY"), rekn_dep, map_crs = map_crs, delbay_paths = rekn_depart,
                backdrop_sf = backdrop, backdrop_label_sf = backdrop_labs,
                country_label_sf = country_labs)
ggsave("Output/Fig2_rekn_path_map.png", dpi = 600, height = 8.6, width = 6.8)
ggsave("Output/Fig2_rekn_path_map.pdf", height = 8.6, width = 6.8)

# Sample size for migration strategy categories
group_by(filter(rekn_depart, age == "ASY"), delbay_simplified) %>% tally()

# Figure 3
# What about how capture body mass related to inferred migration strategy?
rekn_depart <- left_join(rekn_depart, select(rekn_dep, motusTagID, wt_g)) %>%
  mutate(mig_strat = delbay_simplified)
levels(rekn_depart$mig_strat) <- gsub("or likely", "or\nlikely", levels(rekn_depart$mig_strat))
rekn_depart_asy <- filter(rekn_depart, age == "ASY")
rekn_depart_asy_summary <- filter(rekn_depart_asy, !is.na(wt_g)) %>%
  group_by(delbay_simplified, mig_strat) %>% tally()
levels(rekn_depart_asy_summary$mig_strat) <- gsub("or likely", "or\nlikely", levels(rekn_depart_asy_summary$mig_strat))
ggplot(rekn_depart_asy, aes(mig_strat, wt_g)) + 
  geom_boxplot(fill = "gray80") +
  geom_text(data = rekn_depart_asy_summary, aes(y = 105, label = paste0("(", n, ")"))) +
  scale_y_continuous("Capture mass in South Carolina (g)", breaks = seq(110, 160, by = 10)) +
  xlab("Northbound migration use of Delaware Bay") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
ggsave("Output/Fig3_rekn_migration_by_capture_weight.png", dpi = 600, width = 6.5, height = 4.5)

summary(wt_mod <- lm(wt_g ~ delbay_simplified, data = filter(rekn_depart, !is.na(wt_g), age == "ASY")))
par(mfrow = c(2, 2))
plot(wt_mod)
par(mfrow = c(1, 1))

# Examine durations of use of Delaware Bay and Great Lakes Basin
rekn_det_sf <- mutate(rekn_det,
                      lat = recvDeployLat,
                      lon = recvDeployLon) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

rekn_det_gl <- st_intersection(rekn_det_sf, st_transform(gl, crs = 4326)) %>%
  as.data.frame()
rekn_det_db <- st_intersection(rekn_det_sf, st_transform(db, crs = 4326)) %>%
  as.data.frame()

gl_stay_summary <- rekn_det_gl %>%
  group_by(motusTagID, age) %>%
  summarize(stay_d = as.numeric(difftime(max(last_det), min(first_det), units = "days")),
            first_day = yday(min(first_det)), .groups = "drop")
ggplot(gl_stay_summary, aes(x = first_day, y = stay_d, shape = age)) + 
  geom_point() +
  labs(x = "First day of detection in GL") +
  theme_bw()

db_stoppers <- filter(rekn_depart, delbay == "stopped") %>% pull(motusTagID)
db_stay_summary <- rekn_det_db %>%
  filter(motusTagID %in% db_stoppers) %>%
  group_by(motusTagID, age) %>%
  summarize(stay_d = as.numeric(difftime(max(last_det), min(first_det), units = "days")),
            first_day = yday(min(first_det)), .groups = "drop")

# Figure 5
# Duration of DelBay stay vs date of arrival (first detection) in DelBay?
ggplot(db_stay_summary, 
       aes(x = first_day, y = stay_d, shape = age)) + 
  geom_point(fill = "grey75", size = 4) +
  scale_shape_manual(NULL, values = c(24, 22)) +
  labs(y = "Duration of detection in Delaware Bay (days)") +
  scale_x_continuous("First day of detection in Delaware Bay",
                     breaks = seq(125, 155, by = 5), limits = c(123, 156), expand = c(0,0),
                     labels = date_seq()[seq(125, 155, by = 5)]) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.0, 1), 
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.background = element_blank())
ggsave("Output/Fig5_rekn_DB_duration_vs_arrival.png", dpi = 600, width = 6.5, height = 4.5)

with(filter(db_stay_summary, age == "ASY"), summary(stay_d))
summary(stay_mod <- lm(stay_d ~ first_day, data = filter(db_stay_summary, age == "ASY")))
par(mfrow = c(2, 2))
plot(stay_mod)
par(mfrow = c(1, 1))

# Compare first and last detections in GL and DB Basins
rekn_gl_summary <- rekn_det_gl %>%
  group_by(motusTagID, tagDeployStart) %>%
  summarize(`First detection` = min(doy),
            `Last detection` = max(doy),
            .groups = "drop") %>%
  tidyr::pivot_longer(cols = `First detection`:`Last detection`,
                      names_to = "detection",
                      values_to = "doy")

rekn_db_summary <- rekn_det_db %>%
  group_by(motusTagID, tagDeployStart) %>%
  summarize(`First detection` = min(doy),
            `Last detection` = max(doy),
            .groups = "drop") %>%
  tidyr::pivot_longer(cols = `First detection`:`Last detection`,
                      names_to = "detection",
                      values_to = "doy")

# Arctic
# Elapsed time between last GL detection and first James/Hudson Bay detection
arctic <- filter(rekn_det, recvDeployLat > 50, recvDeployLon > -100)
last_gl <- rekn_det_gl %>%
  group_by(motusTagID) %>% 
  arrange(desc(last_det)) %>%
  slice(1)
first_breed <- arctic %>%
  group_by(motusTagID) %>%
  arrange(first_det) %>%
  slice(1)

# Range of first detection dates in Arctic
label_days(range(first_breed$doy))
label_days(median(first_breed$doy))

rekn_gl_breed <- filter(last_gl, motusTagID %in% unique(first_breed$motusTagID)) %>%
  bind_rows(filter(first_breed, motusTagID %in% unique(last_gl$motusTagID))) %>%
  group_by(motusTagID) %>%
  arrange(motusTagID, first_det) %>%
  summarize(tdiff = difftime(max(first_det), min(last_det), units = "hours"))
median(rekn_gl_breed$tdiff)
range(rekn_gl_breed$tdiff)

