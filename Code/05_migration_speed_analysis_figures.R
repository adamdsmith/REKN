pacman::p_load(patchwork)
source("Code/Functions/create_path_speed_map.R")
source("Code/Functions/plot_gspd_twnd_cor.R")

rekn_speeds <- readRDS("Data/Derived/rekn_speed_calculations.rds")
rekn_speeds_redux <- readr::read_csv("Data/Derived/rekn_speed_checks.csv") %>%
  filter(keep,
         # Exclude net airspeeds < 5 m/s (Duijns et al. 2019)
         (groundspeed_net - tailwind) >= 5) 

# Number of flights
nrow(rekn_speeds_redux)

# Number of individuals
n_distinct(rekn_speeds_redux$motusTagID)

# Distance (total trajectory length) summary
median(rekn_speeds_redux$traj_len) / 1000
range(rekn_speeds_redux$traj_len) / 1000
# Distance (trajectory displacement) summary
median(rekn_speeds_redux$traj_net_disp) / 1000
range(rekn_speeds_redux$traj_net_disp) / 1000

# Figure 6
# Create map of paths highlighting speed
p <- create_path_speed_map(rekn_speeds_redux, map_crs = map_crs)
ggsave("Output/rekn_path_speed_map.pdf", height = 8.6, width = 5)
p_sm <- create_path_speed_map(rekn_speeds_redux, map_crs = map_crs, base_size = 12,
                              legend_pos = c(0.675, 0.175))
# Create additional component figures and assemble
gspd_fig <- 
  ggplot(rekn_speeds_redux, aes(x = groundspeed_net)) +
  geom_histogram(breaks = seq(0, 40, by = 2), fill = "grey75", color = "black") +
  scale_y_continuous("# Red Knot", breaks = seq(0, 10, 2)) +
  scale_x_continuous(NULL, breaks = seq(0, 40, 4), limits = c(0, 40)) +
  theme_bw(base_size = 16)

# Groundspeed and tailwind correlation, plus plot
twnd_fig <- plot_gpsd_twnd_cor(rekn_speeds_redux)

# Assemble
all <- (gspd_fig / twnd_fig) | p_sm
all + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 22))
ggsave("Output/Fig6_rekn_speed_fig.png", dpi = 600, height = 7, width = 10)

sum(rekn_speeds_redux$tailwind > 0) / nrow(rekn_speeds_redux)
