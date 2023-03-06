# plot FAR time series for each region and each warming level
# using 15-year rolling windows of observed warming ranges to set "present"
# plot annual and 3yr running mean on the same panels

library(tidyverse)
library(brms)

theme_set(theme_bw())

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# load region names
regions <- read.csv("./CMIP6/summaries/clean_region_names.csv")
regions <- regions[1:6,1]

# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# limit to 1950:2022
ersst.anom <- ersst.anom %>%
  filter(year %in% 1950:2022)

far_pred_annual <- data.frame()

for(i in 1:length(regions)){ # loop through regions

  # i <- 2
  
  # subset ersst.anom
  ersst.temp <- ersst.anom %>%
    filter(region == regions[i]) %>%
    select(year, annual.anomaly.unsmoothed) %>%
    rename(annual.anomaly.1yr = annual.anomaly.unsmoothed)
  
  # load regional model
  mod <- readRDS(paste("./CMIP6/brms_output/", regions[i], "annual_sst_rolling_window_binomial2.rds", sep = ""))
  
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
  
  
  far_pred_annual <- rbind(far_pred_annual,
                    data.frame(region = regions[i],
                              year = nd_pre$ersst.year,
                              prob = apply(far, 2, mean),
                              lower = apply(far, 2, quantile, probs = 0.025),
                              upper = apply(far, 2, quantile, probs = 0.975)))
  
  
} # close i loop


far_pred_annual$window = "annual"
  
############################
## now 3-yr running mean
  far_pred_3yr <- data.frame()
  
  for(i in 1:length(regions)){ # loop through regions
    
    # i <- 1
    
    # subset ersst.anom
    ersst.temp <- ersst.anom %>%
      filter(region == regions[i]) %>%
      select(year, annual.anomaly.three.yr.running.mean) %>%
      rename(annual.anomaly.3yr = annual.anomaly.three.yr.running.mean)
    
    # and drop NAs
    keep <- !is.na(ersst.temp$annual.anomaly.3yr)
    ersst.temp <- ersst.temp[keep,]
    
    # load regional model
    mod <- readRDS(paste("./CMIP6/brms_output/", regions[i], "_3yr_mean_annual_sst_rolling_window_binomial2.rds", sep = ""))
    
    ## setup new data
    nd <- data.frame(period = c("historical", "preindustrial"),
                     ersst.year = rep(ersst.temp$year, each = 2),
                     annual.anomaly.3yr = rep(ersst.temp$annual.anomaly.3yr, each = 2),
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
    
    
    far_pred_3yr <- rbind(far_pred_3yr,
                      data.frame(region = regions[i],
                                 year = nd_pre$ersst.year,
                                 prob = apply(far, 2, mean),
                                 lower = apply(far, 2, quantile, probs = 0.025),
                                 upper = apply(far, 2, quantile, probs = 0.975)))
    
    
  } # close i loop
  
  far_pred_3yr$window <- "3yr_running_mean"
  
  far_pred <- rbind(far_pred_annual, far_pred_3yr)
  
  # reorder areas 
  region.order <- data.frame(region = regions,
                             region_plot = c("North Pacific",
                                             "Eastern Bering Sea",
                                             "Gulf of Alaska",
                                             "British Columbia Coast",
                                             "Northern California Current",
                                             "Southern California Current"),
                             order = 1:6)
  
  far_pred <- left_join(far_pred, region.order)
  
  far_pred$region_plot <- reorder(far_pred$region_plot, far_pred$order)
  
  # reorder windows
  window.order <- data.frame(window = c("annual", "3yr_running_mean"),
                             window_plot = c("Annual SST", 
                                             "Three year mean SST"),
                             window.order = c(1, 2))
  
  far_pred <- left_join(far_pred, window.order)
  
  far_pred$window_plot <- reorder(far_pred$window_plot, far_pred$window.order)

  g <- ggplot(far_pred) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_line(aes(x = year, y = prob, color = window_plot), size = 0.25) +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = window_plot), alpha = 0.15) +
    facet_wrap(~region_plot, scales = "free_y") +
    ylab("Fraction of Attributable Risk") +
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "top") +
    scale_color_manual(values = cb[c(2,6)]) +
    scale_fill_manual(values = cb[c(2,6)])
  
  print(g)
  
ggsave("./CMIP6/figs/FAR_rolling_window_time_series_annual_3yr.png", height = 4.5, width = 8, units = 'in')  

## save version with label for paper
g <- ggplot(far_pred) +
  geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
  geom_line(aes(x = year, y = prob, color = window_plot), size = 0.25) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = window_plot), alpha = 0.15) +
  facet_wrap(~region_plot, scales = "free_y", ncol = 2) +
  ylab("Fraction of Attributable Risk") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top") +
  scale_color_manual(values = cb[c(2,6)]) +
  scale_fill_manual(values = cb[c(2,6)]) +
  labs(tag = "B") +
  theme(plot.tag.position = c(0.03, 0.95),
        plot.tag = element_text(size = 16))

print(g)

ggsave("./CMIP6/figs/FAR_rolling_window_time_series_annual_3yrlabelled.png", height = 6, width = 6, units = 'in')  

###
# calculate RR 

far_pred <- far_pred %>%
  mutate(RR = 1/(1-prob),
         lower_RR = 1/(1-lower),
         upper_RR = 1/(1-upper)) %>%
  rename(FAR = prob,
         lower_FAR = lower,
         upper_FAR = upper,
         region_order = order)


g <- ggplot(far_pred) +
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
  scale_y_continuous(breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000),
                     minor_breaks = NULL,
                     labels = c(expression(10^0),
                                expression(10^1),
                                expression(10^2),
                                expression(10^3),
                                expression(10^4),
                                expression(10^5),
                                expression(10^6),
                                expression(10^7),
                                expression(10^8))) 

print(g)

ggsave("./CMIP6/figs/Risk_Ratio_annual_3yr.png", height = 4.5, width = 8, units = 'in') 


# and save output
write.csv(far_pred, "./CMIP6/summaries/complete_FAR_RR_time_series_with_uncertainty.csv", row.names = F)

