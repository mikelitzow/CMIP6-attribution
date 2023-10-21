# figure of ERSST vs weighted means of CMIP6 historical, SSP245, and SSP585 runs for paper

library(tidyverse)
library(matrixStats)

theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load model weights
weights <- read.csv("./CMIP6/summaries/normalized_CMIP6_weights.csv")

# load ERSST time series
ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv") %>%
  select(region, year, annual.unsmoothed) %>%
  rename(SST = annual.unsmoothed) %>%
  mutate(group = "ERSST")

# load historical runs
historical <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv") %>%
  filter(experiment == "hist_ssp585",
         year < 2015) %>%
  select(region, model, year, annual.unsmoothed) %>%
  rename(SST = annual.unsmoothed) %>%
  left_join(., weights)

# load SSP585 temps
ssp585 <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv") %>%
  filter(experiment == "hist_ssp585",
         year > 2014) %>%
  select(region, model, year, annual.unsmoothed) %>%
  rename(SST = annual.unsmoothed) %>%
  left_join(., weights)

# load SSP245 temps
ssp245 <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.ssp245.csv") %>%
  filter(experiment == "hist_ssp245",
         year > 2014) %>%
  select(region, model, year, annual.unsmoothed) %>%
  rename(SST = annual.unsmoothed) %>%
  mutate(model = str_remove_all(model, "_245")) %>%
  left_join(., weights)


historical_mean <- historical %>%
  group_by(region, year) %>%
  summarise(weighted_mean = weighted.mean(x = SST, w = normalized_weight),
            weighted_sd = weightedSd(x = SST, w = normalized_weight)) %>%
  mutate(group = "Historical")

ssp585_mean <- ssp585 %>%
  group_by(region, year) %>%
  summarise(weighted_mean = weighted.mean(x = SST, w = normalized_weight),
            weighted_sd = weightedSd(x = SST, w = normalized_weight)) %>%
  mutate(group = "SSP585")

ssp245_mean <- ssp245 %>%
  group_by(region, year) %>%
  summarise(weighted_mean = weighted.mean(x = SST, w = normalized_weight),
            weighted_sd = weightedSd(x = SST, w = normalized_weight)) %>%
  mutate(group = "SSP245")

model_plot <- rbind(historical_mean, ssp245_mean, ssp585_mean)

g1 <- ggplot(model_plot, aes(year, weighted_mean, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = weighted_mean - 2*weighted_sd,
                    ymax = weighted_mean + 2*weighted_sd),
                alpha = 0.2, color = NA) +
  facet_wrap(~region, scales = "free_y") 

g1

g1 +
  geom_line(data = ersst, aes(year, SST), color = "black", size = 0.5) +
  theme(legend.position = "none")

## the same figure with anomalies in degrees ----------

# get pre-1950 climatological mean for each region/model
historical_climatology <- historical %>%
  filter(year < 1950) %>%
  group_by(region, model) %>%
  summarize(pre_1950_mean_temp = mean(SST))
  
historical_mean_anomaly <- left_join(historical, historical_climatology) %>%
  mutate(anomaly = SST - pre_1950_mean_temp) %>%
  group_by(region, year) %>%
  summarise(weighted_mean = weighted.mean(x = anomaly, w = normalized_weight),
            weighted_sd = weightedSd(x = anomaly, w = normalized_weight)) %>%
  mutate(group = "CMIP6 Historical")

ssp585_mean_anomaly <- left_join(ssp585, historical_climatology)  %>%
  mutate(anomaly = SST - pre_1950_mean_temp) %>%
  group_by(region, year) %>%
  summarise(weighted_mean = weighted.mean(x = anomaly, w = normalized_weight),
            weighted_sd = weightedSd(x = anomaly, w = normalized_weight)) %>%
  mutate(group = "CMIP6 SSP585")

ssp245_mean_anomaly <- left_join(ssp245, historical_climatology)  %>%
  mutate(anomaly = SST - pre_1950_mean_temp) %>%
  group_by(region, year) %>%
  summarise(weighted_mean = weighted.mean(x = anomaly, w = normalized_weight),
            weighted_sd = weightedSd(x = anomaly, w = normalized_weight)) %>%
  mutate(group = "CMIP6 SSP245")


# get pre-1950 climatological mean for each region in ersst
ersst_climatology <- ersst %>%
  filter(year < 1950) %>%
  group_by(region) %>%
  summarize(pre_1950_mean_temp = mean(SST),
            pre_1950_sd_temp = sd(SST))

ersst_mean_anomaly <- left_join(ersst, ersst_climatology) %>%
  mutate(anomaly = SST - pre_1950_mean_temp)

ersst_plot <- ersst_mean_anomaly %>%
  rename(weighted_mean = anomaly)  %>%
  mutate(group = "ERSST",
         region = str_replace_all(region, "_", " "),
         weighted_sd = NA) %>%
  select(region, year, weighted_mean, weighted_sd, group)

# all_years <- data.frame(region = rep(unique(ersst_climatology$region), each = 250),
#                         year = 1850:2099)
# 
# preindustrial_obs_envelope <- left_join(all_years, ersst_climatology) %>%
#   mutate(anomaly = 0,
#          group = "Pre-1950 envelope") %>%
#   select(region, year, anomaly, sd, group)

anomaly_plot <- rbind(ersst_plot, historical_mean_anomaly, ssp245_mean_anomaly, ssp585_mean_anomaly) 

# change <- anomaly_plot$group == "Historical"
# 
# anomaly_plot$sd[change] <- NA
# clean up and order region names
anomaly_plot <- anomaly_plot %>%
  mutate(region = str_replace_all(region, "_", " "))


region_order <- data.frame(region = unique(anomaly_plot$region),
                          plot_order = 1:6)

anomaly_plot <- left_join(anomaly_plot, region_order)

anomaly_plot$region <- reorder(anomaly_plot$region, anomaly_plot$plot_order)


g2 <- ggplot(anomaly_plot, aes(year, weighted_mean, color = group, fill = group)) +
  geom_line(size = 0.25) +
  geom_ribbon(aes(ymin = weighted_mean - 2*weighted_sd,
                  ymax = weighted_mean + 2*weighted_sd),
              alpha = 0.15, color = NA) +
  facet_wrap(~region, scales = "free_y") +
  scale_color_manual(values = c(cb[c(2,1,7)], "black")) +
  scale_fill_manual(values = c(cb[c(2,1,7)], "white")) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank()) +
  labs(y = "SST anomaly (°C)")

g2


ggsave("./CMIP6/figs/SST_time_series.png", width = 8, height = 4.5)

# add a version with uniform y-axis scale
g3 <- ggplot(anomaly_plot, aes(year, weighted_mean, color = group, fill = group)) +
  geom_line(size = 0.25) +
  geom_ribbon(aes(ymin = weighted_mean - 2*weighted_sd,
                  ymax = weighted_mean + 2*weighted_sd),
              alpha = 0.15, color = NA) +
  facet_wrap(~region) +
  scale_color_manual(values = c(cb[c(2,1,7)], "black")) +
  scale_fill_manual(values = c(cb[c(2,1,7)], "white")) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank()) +
  labs(y = "SST anomaly (°C)")


g3

ggsave("./CMIP6/figs/SST_time_series_uniform_y-axis.png", width = 8, height = 4.5)
