# plot FAR time series for each region and each warming level
# using 15-year rolling windows of observed warming ranges to set "present"

library(tidyverse)
library(brms)

theme_set(theme_bw())

# load region names
regions <- read.csv("./CMIP6/summaries/clean_region_names.csv")
regions <- regions[1:6,1]

# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

far_pred <- data.frame()

# for(i in 1:length(regions)){ # loop through regions

for(i in 1:3){  
  # i <- 1
  
  # subset ersst.anom
  ersst.temp <- ersst.anom %>%
    filter(region == regions[i]) %>%
    select(year, annual.anomaly.unsmoothed) %>%
    rename(annual.anomaly.1yr = annual.anomaly.unsmoothed)
  
  # load regional model
  mod <- readRDS(paste("./CMIP6/brms_output/", regions[i], "_rolling_window_binomial2.rds", sep = ""))
  
  ## setup new data
  nd <- data.frame(period = c("historical", "preindustrial"),
                   ersst.year = rep(ersst.temp$year, each = 2),
                   annual.anomaly.1yr = rep(ersst.temp$annual.anomaly.1yr, each = 2),
                   N = 1000,
                   model_fac = NA)
  
  nd_pre <- nd[nd$period == "preindustrial", ]
  nd_his <- nd[nd$period == "historical", ]
  
  ## make predictions
  ## exclude random effects for model_fac
  pre_pp <- posterior_epred(mod, newdata = nd_pre, re_formula = NA)
  his_pp <- posterior_epred(mod, newdata = nd_his, re_formula = NA)
  
  ## Calc probabilities
  ## These are our posterior probabilities to use for FAR calculation
  pre_prob <- pre_pp / unique(nd$N)
  his_prob <- his_pp / unique(nd$N)
  
  
  ## Calc FAR
  far <- 1 - (pre_prob / his_prob)
  range(far, na.rm = TRUE)
  
  
  far_pred <- rbind(far_pred,
                    data.frame(region = regions[i],
                              year = nd_pre$ersst.year,
                              prob = apply(far, 2, mean),
                              lower = apply(far, 2, quantile, probs = 0.025),
                              upper = apply(far, 2, quantile, probs = 0.975)))
  
  
} # close i loop
  
# reorder areas 
region.order <- data.frame(region = regions,
                           order = 1:6)

far_pred <- left_join(far_pred, region.order)

far_pred$region <- reorder(far_pred$region, far_pred$order)

  g <- ggplot(far_pred) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_line(aes(x = year, y = prob), size = 0.25, color = "red") +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper), alpha = 0.15) +
    facet_wrap(~region) +
    ylab("Fraction of Attributable Risk") +
    theme(axis.title.x = element_blank())
  
  print(g)
  
ggsave("./CMIP6/figs/FAR_rolling_window_time_series.png", height = 2.5, width = 8, units = 'in')  
  
# print 0.5 - 1.0 degree version

plot.dat <- far_pred %>%
  filter(present == "0.5 - 1.0 degrees")

g <- ggplot(plot.dat) +
  geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
  geom_line(aes(x = year, y = prob), size = 0.3, color = "red") +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper), alpha = 0.15) +
  facet_wrap(~region) +
  theme(axis.title.x = element_blank()) +
  ylab("Fraction of Attributable Risk")

print(g)

ggsave("./CMIP6/figs/regional_FAR_0.5-1.0_degrees_warming.png", width = 6, height = 3, units = 'in')


# add risk ratios!
far_pred <- far_pred %>%
  mutate(risk_ratio = 1/(1-prob))

View(far_pred)  
