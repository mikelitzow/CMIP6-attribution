# make SI table summarizing mean and SD statistics for 1854-1949 climatology

library(tidyverse)

dat <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv")

summary <- dat %>%
  filter(year %in% 1854:1949) %>%
  group_by(region) %>%
  summarise(annual_mean = mean(annual.unsmoothed),
            annual_SD = sd(annual.unsmoothed),
            three_yr_mean = mean(annual.three.yr.running.mean, na.rm = T),
            three_year_SD = sd(annual.three.yr.running.mean, na.rm = T)) %>%
  mutate(order = c(4,2,3,1,5,6)) %>%
  arrange(order) %>%
  select(-order)

summary

write.csv(summary, "./CMIP6/summaries/regional_ERSST_mean_SD_1854-1949.csv", row.names = F)
