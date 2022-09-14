# interpret regional FAR signals

library(tidyverse)

dat <- read.csv("./CMIP6/summaries/complete_FAR_RR_time_series_with_uncertainty.csv")

# get last year FAR for each scale is indistinguishable from 0

last <- dat %>%
  filter(region != "North_Pacific",
         lower_FAR < 0) %>%
  group_by(region, window) %>%
  summarise(last_year = max(year)) %>%
  arrange(desc(window))

last

write.csv(last, "./CMIP6/summaries/last_year_FAR_zero.csv", row.names = F)

ggplot(last, aes(last_year)) +
  geom_histogram(bins = 11, fill = "grey", color = "black") +
  facet_grid(~window)


# check vs RR

check <- last %>%
  rename(year = last_year) %>%
  left_join(., dat) %>%
  arrange(desc(window))

check

# check years with three-year mean FAR >= 0.5 for EBS and GOA

FAR_0.5 <- dat %>%
  filter(region %in% c("Eastern_Bering_Sea", "Gulf_of_Alaska"),
         window == "3yr_running_mean",
         FAR >= 0.5) %>%
  arrange(region)

FAR_0.5


