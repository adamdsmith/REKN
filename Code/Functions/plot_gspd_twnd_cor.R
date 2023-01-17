plot_gpsd_twnd_cor <- function(speed_dat, method = c("pearson", "spearman"),
                               x_breaks = seq(0, 40, by = 4)) {
  method <- match.arg(method)

  # Trajectory groundspeed association with initial tailwind to destination receiver
  gspd_twnd_cor <- with(speed_dat, cor.test(groundspeed_traj, tailwind, method = method))
  r_val <- round(gspd_twnd_cor$estimate, 2)
  p_val <- round(gspd_twnd_cor$p.value, 3)

  if (method == "pearson") 
    cor_label <- bquote("Pearson's"~italic(r) == .(r_val))
  else 
    cor_label <- bquote("Spearman's"~rho == .(r_val))

  y_rng <- c(min(floor(speed_dat$tailwind)), max(ceiling(speed_dat$tailwind)))

  p <-  
    ggplot(speed_dat, aes(x = groundspeed_net, y = tailwind)) +
    geom_point(shape = 21, fill = "grey50", color = "black", size = 5, alpha = 0.5) +
    scale_y_continuous("Initial tailwind (m/s)") +
    scale_x_continuous("Ground speed (net; m/s)",
                       breaks = x_breaks, limits = range(x_breaks)) +
    annotate(x = min(x_breaks), y = max(y_rng), label = cor_label, 
             hjust = 0, vjust = 1, geom = "text", size = 5) +
    annotate(x = min(x_breaks), y = max(y_rng) - diff(y_rng)/10, label = bquote(italic(p) == .(p_val)),
             hjust = 0, vjust = 1, geom = "text", size = 5) +
    theme_bw(base_size = 16)
  print(gspd_twnd_cor)
  return(p)
}
