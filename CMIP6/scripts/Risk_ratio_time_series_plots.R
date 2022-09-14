# annual vs 3-yr mean FAR / RR comparison

library(tidyverse)

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw())


dat <- read.csv("./CMIP6/summaries/complete_FAR_RR_time_series_with_uncertainty.csv")

dat$window_plot <- reorder(dat$window_plot, dat$window.order)
dat$region_plot <- reorder(dat$region_plot, dat$region_order)

g <- ggplot(dat) +
  geom_line(aes(x = year, y = RR, color = window_plot), size = 0.25) +
  geom_ribbon(aes(x = year, ymin = lower_RR, ymax = upper_RR, fill = window_plot), alpha = 0.15) +
  facet_wrap(~region_plot, scales = "free_y") +
  ylab("Risk Ratio") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top") +
  scale_color_manual(values = cb[c(2,6)]) +
  scale_fill_manual(values = cb[c(2,6)]) +
  coord_trans(y = "pseudo_log") +
  scale_y_continuous(breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
                     minor_breaks = NULL,
                     labels = c("10^0",
                                "10^1",
                                "10^2",
                                "10^3",
                                "10^4",
                                "10^5",
                                "10^6",
                                "10^7")) 


print(g)

ggsave("./CMIP6/figs/RR_rolling_window_time_series_annual_3yr.png", height = 4.5, width = 8, units = 'in')  
