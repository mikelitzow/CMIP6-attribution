# plots for snow crab and climate change section in council document

library(tidyverse)
theme_set(theme_bw())

# first, annual SST

dat <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv")

head(dat)

dat <- dat %>%
  filter(region == "Eastern_Bering_Sea")

ggplot(dat, aes(year, annual.unsmoothed)) +
  geom_point() +
  geom_line() +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = mean(dat$annual.unsmoothed), lty = 2) +
  labs(title = "Annual mean sea surface temperature, 1854-2021",
       subtitle = "Source: NOAA Extended Reconstructed SST v5",
       y = "°C") +
  scale_x_continuous(breaks = seq(1850, 2020, 10), minor_breaks = NULL)

ggsave("./CMIP6/figs/EBS_SST_plot.png", width = 7, height = 4, units = 'in')

dat <- dat %>%
  mutate(diff = annual.unsmoothed - mean(annual.unsmoothed))

tail(dat)

dat[dat$year %in% c(2016, 2018, 2019),]
