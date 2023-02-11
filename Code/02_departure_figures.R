source("Code/01_Calculate_parameters_of_interest.R")

rekn_depart <- readr::read_csv("Data/Derived/rekn_departures.csv",
                            col_types = "icciTddDDcccc") %>%
  filter(delbay != "no_detections") %>%
  mutate(tdiff = yday(max_depart) - yday(min_depart) + 1,
         delbay = factor(delbay, levels = c("stopped", "likely",
                                            "unknown", "unlikely", "skipped")),
         delbay_simplified = factor(delbay, labels = c("Stopped or likely stopped", "Stopped or likely stopped",
                                                       "Unknown", "Skipped or likely skipped",
                                                       "Skipped or likely skipped")))

# Two approaches to SC departure date figure
# 1 - use only departures with a precisely known departure date
# 2 - use all departures, but construct density that reflects our uncertainty in
#     that each individual contributes equally, but dates of departure for those
#     individuals with a possible range of dates only contribute fractionally to
#     the constructed "density" (i.e., each possible date contributes only
#     [1 / # possible days of departure])
# In both scenarios, we summarize by week of departure since data are fairly sparse

days_of_year <- seq.Date(from = as.Date("2018-01-01"), to = as.Date("2018-12-31"), by = 1) %>%
  format(format = "%d %b")
doy_labels <- format(days_of_year, format = "%d %b")
wk_starts <- seq.Date(from = as.Date("2015-01-01"), by = 7, length.out = 52)
wk_breaks <- week(wk_starts)
wk_labels <- format(wk_starts, format = "%d %b")
label_days  <- function(x, nth = 1) every_nth(doy_labels[x], nth, inverse = TRUE)
label_weeks  <- function(x, nth = 1) every_nth(wk_labels[x], nth, inverse = TRUE)
make_breaks <- function(x, by, ...) seq(floor(x[1]), ceiling(x[2]), by = by)

# Approach 1
dep_days <- filter(rekn_depart, age == "ASY", tdiff == 1) %>%
  mutate(dep_doy = yday(max_depart),
         dep_wk = week(max_depart))

# calculate median departure date for display
dep_days_sum <- dep_days %>%
  group_by(delbay_simplified) %>%
  summarize(n = n(),
            med_dep_doy = median(dep_doy),
            med_dep_wk = med_dep_doy / 7 + 1,
            med_dep_date = days_of_year[round(med_dep_doy)])

dep_fig_precise <- 
  ggplot(dep_days, aes(dep_doy)) +
  geom_vline(data = dep_days_sum, aes(xintercept = med_dep_doy), 
               color = "gray50", lwd = 4) + 
  geom_text(data = dep_days_sum, aes(x = med_dep_doy + .5, y = 2.5, 
                                      label = paste("Median:\n", med_dep_date)), 
            size = 6, hjust = 0, vjust = 1) +
  geom_bar(fill = "grey75", color = "black") +
  geom_text(data = dep_days_sum, aes(x = min(dep_days$dep_doy), y = 2.5, label = paste("n =", n)), 
            hjust = 0, size = 6, vjust = 1) +
  scale_x_continuous("Date of departure", breaks = make_breaks(range(dep_days$dep_doy), by = 3),
                     minor_breaks = make_breaks(range(dep_days$dep_doy), by = 1), 
                     labels = label_days) +
  scale_y_continuous("# Red Knots", breaks = 0:2, limits = c(0, 2.5), minor_breaks = NULL) +
  theme_bw(base_size = 16) + facet_wrap(~delbay_simplified, nrow = 1)

# Approach 2
dep_all <- rekn_depart %>%
  filter(!is.na(min_depart), age == "ASY") %>%
  mutate(date_rng = purrr::pmap(., function(min_depart, max_depart, ...) {
    seq.Date(min_depart, max_depart, by = 1)
  })) %>%
  tidyr::unnest(cols = c(date_rng)) %>%
  mutate(dt_wt = 1 / tdiff,
         doy = yday(date_rng))

# calculate median departure date by banding location for display
dep_all_sum <- dep_all %>%
  group_by(delbay_simplified) %>%
  summarize(n = n_distinct(motusTagID),
            wt_med_dep_doy = median_wtd(doy, dt_wt),
            wt_med_dep_wk = wt_med_dep_doy / 7 + 1,
            wt_med_dep_date = days_of_year[round(wt_med_dep_doy)])

dep_all_wt <- dep_all %>%
  mutate(dep_wk = week(date_rng)) %>%
  group_by(delbay_simplified, doy) %>%
  summarize(date_wt = sum(dt_wt))

# Superimpose birds with known departure dates
dep_known_wt <- dep_days %>%
  mutate(doy = dep_doy) %>%
  group_by(delbay_simplified, doy) %>%
  summarize(date_wt = sum(tdiff))

dep_fig_wt_spr <- ggplot(dep_all_wt, aes(doy)) +
  geom_vline(data = dep_all_sum, aes(xintercept = wt_med_dep_doy),
             color = "grey50", lwd = 4) +
  geom_text(data = dep_all_sum, aes(x = wt_med_dep_doy + .75, 
                                    y = 4.5, vjust = 1,
                                    label = wt_med_dep_date),
            size = 5, hjust = 0) +
  geom_bar(aes(y = date_wt), stat = "identity", fill = "grey75", color = "black") +
  geom_text(data = dep_known_wt, aes(y = 0.05, label = date_wt), vjust = 0) +
  geom_text(data = dep_all_sum, aes(x = wt_med_dep_doy + .75, y = 4, label = paste("n =", n)), 
            vjust = 1, hjust = 0, size = 5) +
  scale_x_continuous("Date of South Carolina departure", breaks = make_breaks(range(dep_all$doy) + c(-1, 1), by = 7), 
                     minor_breaks = make_breaks(range(dep_all$doy) + c(-1, 1), by = 1),
                     labels = label_days) +
  scale_y_continuous("# Red Knots", limits = c(0,4.5)) +
  facet_wrap(~ delbay_simplified, ncol = 1) +
  tag_facets(tag_levels = "A", tag_suffix = "") + 
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        tagger.panel.tag.background = element_blank(), 
        tagger.panel.tag.text = element_text(size = 16))
ggsave("Output/rekn_asy_departartures_weighted.png", dpi = 600, height = 7.5, width = 5)

