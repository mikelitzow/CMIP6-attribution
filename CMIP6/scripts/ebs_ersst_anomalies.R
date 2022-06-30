library(tidyverse)

dat <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv")

dat <- dat %>%
  filter(region == "Eastern_Bering_Sea") %>%
  mutate(anomaly.degrees = annual.unsmoothed - mean(annual.unsmoothed[year %in% 1854:1949]),
         anomaly.sd = (annual.unsmoothed - mean(annual.unsmoothed[year %in% 1854:1949]))/sd(annual.unsmoothed[year %in% 1854:1949]),
         anomaly.degrees.F = (9/5)*anomaly.degrees) %>%
  select(year, annual.unsmoothed, anomaly.degrees, anomaly.degrees.F, anomaly.sd)

View(dat)
