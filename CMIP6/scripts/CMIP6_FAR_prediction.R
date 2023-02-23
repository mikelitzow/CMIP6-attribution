# predict FAR for sst anomaly time series in each region

library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw())

ersst.dat <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv") 

ersst.dat <- ersst.dat %>%
  select(region, year, annual.anomaly.unsmoothed) %>%
  rename(annual.anomaly.1yr = annual.anomaly.unsmoothed)

# get vectors of regions and files to loop through
regions <- unique(ersst.dat$region)

# plot anomaly time series for each region
# put regions in order to plot
plot.regions <- data.frame(region = regions,
                           order = 1:6)

ersst.dat <- left_join(ersst.dat, plot.regions)

ersst.dat$region <- reorder(ersst.dat$region, ersst.dat$order)

ggplot(ersst.dat, aes(year, annual.anomaly.1yr)) +
  geom_line() +
  facet_wrap(~region) +
  ylab("Anomaly with respect to 1950-1999") +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0)

ggsave("./CMIP6/figs/observed_sst_anomalies_1950-2022_by_region.png", width = 8, height = 4)

# version with just ecosystems of focus for salmon study
just.dat <- ersst.dat %>%
  filter(region %in% c("Gulf_of_Alaska", "British_Columbia_Coast"))

ggplot(just.dat, aes(year, annual.anomaly.1yr, color = region)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)]) +
  ylab("Anomaly with respect to 1950-1999") +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0)

ggsave("./CMIP6/figs/observed_sst_anomalies_1950-2021_GOA_BC_Coast.png", width = 6, height = 3)


## loop through each region and fit predicted FAR to observed sst ------------------

file.list <- NA

for(i in 1:length(regions)){
  
  file.list[i] = paste("./CMIP6/brms_output/far_1yr_annual_base_", regions[i], ".rds", sep = "")
  
}


plot.predict <- data.frame()


for(i in 1:length(regions)){
  # i <- 1  
  model <- readRDS(file.list[i])
    
  newdata <- ersst.dat %>%
    filter(region == regions[i])

  pred.far_1yr_base <- posterior_epred(model, newdata = newdata, re_formula = NA, resp = "FAR.annual.1yr")

  # functions for different credible intervals
  f_95_l <- function(x) quantile(x, 0.025)
  f_95_u <- function(x) quantile(x, 0.975)

  f_90_l <- function(x) quantile(x, 0.05)
  f_90_u <- function(x) quantile(x, 0.95)

  f_80_l <- function(x) quantile(x, 0.1)
  f_80_u <- function(x) quantile(x, 0.9)


  # and plot
  plot.predict <- rbind(plot.predict, data.frame(
                            region = regions[i],
                            year = 1950:2021,
                            estimate__ = colMeans(pred.far_1yr_base),
                            lower_95 = apply(pred.far_1yr_base, 2, f_95_l),
                            upper_95 = apply(pred.far_1yr_base, 2, f_95_u),
                            lower_90 = apply(pred.far_1yr_base, 2, f_90_l),
                            upper_90 = apply(pred.far_1yr_base, 2, f_90_u),
                            lower_80 = apply(pred.far_1yr_base, 2, f_80_l),
                            upper_80 = apply(pred.far_1yr_base, 2, f_80_u)))
}

# put regions in order
plot.predict <- left_join(plot.predict, plot.regions)

plot.predict$region <- reorder(plot.predict$region, plot.predict$order)

ggplot(plot.predict) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  facet_wrap(~region) + 
  labs(y = "Fraction of attributable risk", x = "SST anomaly") 

ggsave("./CMIP6/figs/predicted_far_1950-2021_far_1yr_base.png", width = 6, height = 4)