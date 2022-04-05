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

far_pred_annual <- data.frame()

for(i in 1:length(regions)){ # loop through regions
 
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
  
  
  far_pred_annual <- rbind(far_pred_annual,
                    data.frame(region = regions[i],
                              year = nd_pre$ersst.year,
                              prob = apply(far, 2, mean),
                              lower = apply(far, 2, quantile, probs = 0.025),
                              upper = apply(far, 2, quantile, probs = 0.975)))
  
  
} # close i loop

far_pred_annual$window = "annual"
  

  
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
                             order = 1:6)
  
  far_pred <- left_join(far_pred, region.order)
  
  far_pred$region <- reorder(far_pred$region, far_pred$order)
  
  # reorder windows
  window.order <- data.frame(window = c("annual", "3yr_running_means"),
                             window.order = c(1, 2))
  
  far_pred <- left_join(far_pred, window.order)
  
  far_pred$window <- reorder(far_pred$window, far_pred$window.order)
  
  g <- ggplot(far_pred) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_line(aes(x = year, y = prob, color = window), size = 0.25) +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = window), alpha = 0.15) +
    facet_wrap(~region, scales = "free_y") +
    ylab("Fraction of Attributable Risk") +
    theme(axis.title.x = element_blank()) +
    scale_color_manual(values = cb[c(2,6)]) +
    scale_fill_manual(values = cb[c(2,6)])
  
  print(g)
  
ggsave("./CMIP6/figs/FAR_rolling_window_time_series_annual_3yr.png", height = 4.5, width = 10, units = 'in')  
  

# same plot for sst anomalies
ersst.anom

ersst.plot <- ersst.anom %>%
  select(region, year, annual.anomaly.unsmoothed, annual.anomaly.three.yr.running.mean)

names(ersst.plot)[3:4] <- c("annual", "3yr_running_mean")

ersst.plot <- ersst.plot %>%
  pivot_longer(cols = c(-region, -year), names_to = "window")

ersst.plot <- left_join(ersst.plot, region.order)

ersst.plot$region <-reorder(ersst.plot$region, ersst.plot$order)

ersst.plot <- left_join(ersst.plot, window.order)

ersst.plot$window <- reorder(ersst.plot$window, ersst.plot$window.order)


g <- ggplot(ersst.plot) +
  geom_line(aes(x = year, y = value, color = window), size = 0.25) +
  facet_wrap(~region, scales = "free_y") +
  ylab("Anomaly wrt 1950-1999") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = cb[c(2,6)]) +
  scale_fill_manual(values = cb[c(2,6)])

print(g)

ggsave("./CMIP6/figs/ersst_regional_anomalies_annual_3yr.png", height = 4.5, width = 10, units = 'in')  
