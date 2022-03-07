# evaluate CMIP6 models in terms of
# bias, correlation, and autocorrelation
# for each region

library(tidyverse)

# set theme
theme_set(theme_bw())

## load ERSST and CMIP6 sst time series ------------------------

ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv")
cmip <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv")

# confirm the regions are the same between the two
identical(unique(ersst$region), unique(cmip$region))

## model evaluation --------------------------------------------

# restricts to the reference period 
# defined by the overlap of trustworthy ERSST for GOA regions (1950-2012)
# and historical simulation period for CMIP6 (1850-2014)

# get vector of CMIP6 models to loop through
cmip.models <- unique(cmip$model)
regions <- unique(ersst$region)

# limit ersst to 1950:2014 for comparison
ersst <- ersst %>%
  filter(year %in% 1950:2014)

# limit cmip to historical runs and 1950:2014
cmip.evaluate <- cmip %>%
  filter(experiment == "hist_ssp585",
                 year %in% 1950:2014)
  
# make function for calculating first-order ar
f.ar <- function(x) stats::ar(x, order.max = 1, aic = F)$ar

# data frame for saving results
model.evaluation <- data.frame()

for(i in 1:length(cmip.models)){ # open i loop (cmip6 models)
  # i <- 1
  
  # separate model of interest, historical run, 1950:2014
  temp.cmip <- cmip.evaluate %>%
    filter(model == cmip.models[i])
  
  for(j in 1:length(regions)){
    # j <- 1
    
    # subset cmip and ersst for region of interest
    
    temp.cmip.region <- temp.cmip %>%
      filter(region == regions[j])
    
    ersst.region <- ersst %>%
      filter(region == regions[j])
    
    # calculate bias as difference in mean sst for 1950-2014: cmip - ersst 

    temp.bias.annual <- abs(mean(temp.cmip.region$annual.unsmoothed) - mean(ersst.region$annual.unsmoothed))

    temp.bias.winter <- abs(mean(temp.cmip.region$winter.unsmoothed, na.rm = T) - mean(ersst.region$winter.unsmoothed, na.rm = T))   
    
    # calculate correlation for low-frequency variability (smoothed with 10-yr running means)
    
    temp.cor.annual <- cor(zoo::rollmean(temp.cmip.region$annual.unsmoothed, 10, fill = NA),
                           zoo::rollmean(ersst.region$annual.unsmoothed, 10, fill = NA),
                           use = "p")

    temp.cor.winter <- cor(zoo::rollmean(temp.cmip.region$winter.unsmoothed, 10, fill = NA),
                           zoo::rollmean(ersst.region$winter.unsmoothed, 10, fill = NA),
                           use = "p")
    
    # calculate absolute difference in first-order autocorrelation
    
    temp.ar.annual = abs(f.ar(temp.cmip.region$annual.unsmoothed) - f.ar(ersst.region$annual.unsmoothed))
    
    temp.ar.winter = abs(f.ar(na.omit(temp.cmip.region$winter.unsmoothed)) - f.ar(na.omit(ersst.region$winter.unsmoothed)))
    
    # save result
    model.evaluation <- rbind(model.evaluation,
                              data.frame(model = cmip.models[i],
                                         region = regions[j],
                                         window = c("annual", "winter"),
                                         bias = c(temp.bias.annual, temp.bias.winter),
                                         correlation = c(temp.cor.annual, temp.cor.winter),
                                         ar = c(temp.ar.annual, temp.ar.winter)))


  } # close j loop (regions)
  
} # close i loop (cmip6 models)


# save output
write.csv(model.evaluation, "./CMIP6/summaries/model_evaluation.csv", row.names = F)


# plot
# add ordered factor
order <- data.frame(region = regions,
                    order = 1:6)

model.evaluation <- left_join(model.evaluation, order)

model.evaluation$region <- reorder(model.evaluation$region, model.evaluation$order)

plot.eval <- model.evaluation %>%
  mutate(region.window = paste(region, window, sep = "_")) %>%
  mutate(region.window = reorder(region.window, order))

ggplot(plot.eval, aes(bias, correlation)) +
  geom_point() +
  facet_wrap(~region.window, scales = "free", ncol = 2)

ggsave("./CMIP6/figs/regional_model_correlation_vs_bias.png", width = 6, height = 9, units = "in")

# really interesting - low correlation for NCC and (especially) SSC - not so predictable there
# b/c low influence of decadal Pacific variability?

# also, very high correlation/ high bias for entire North Pacific
# makes sense to see high correlation (less regional variability to capture) but higher bias at 
# larger scale seems like a surprise

ggplot(plot.eval, aes(bias, ar)) +
  geom_point() +
  facet_wrap(~region.window, scales = "free", ncol = 2)


ggsave("./CMIP6/figs/regional_model_ar1_vs_bias.png", width = 6, height = 9, units = "in")

# turn model evaluation into model weights
# first, scale

sc.bias <- model.evaluation %>%
  mutate(region.window = paste(region, window, sep = "-")) %>%
  split(.$region.window) %>%
  map_df(~ scale(.$bias)) %>%
  mutate(model = cmip.models) %>%
  pivot_longer(cols = -model, values_to = "bias.scaled", names_to = "region.window") 

sc.correlation <- model.evaluation %>%
  mutate(region.window = paste(region, window, sep = "-")) %>%
  split(.$region.window) %>%
  map_df(~ scale(.$correlation)) %>%
  mutate(model = cmip.models) %>%
  pivot_longer(cols = -model, values_to = "correlation.scaled", names_to = "region.window") 


sc.ar <- model.evaluation %>%
  mutate(region.window = paste(region, window, sep = "-")) %>%
  split(.$region.window) %>%
  map_df(~ scale(.$ar)) %>%
  mutate(model = cmip.models) %>%
  pivot_longer(cols = -model, values_to = "ar.scaled", names_to = "region.window") 
  
# recombine

scaled.model.evaluation <- model.evaluation %>%
  select(model, region, window) %>%
  mutate(region.window = paste(region, window, sep = "-")) %>%
  left_join(., sc.bias) %>%
  left_join(., sc.correlation) %>%
  left_join(., sc.ar)

# make sure this is in order for matching back to model.evaluation
 
cor(scaled.model.evaluation$bias.scaled[scaled.model.evaluation$region.window == "British_Columbia_Coast-annual"], 
    model.evaluation$bias[model.evaluation$region == "British_Columbia_Coast" & model.evaluation$window == "annual"])
# great!

# now change sign for bias and ar difference (bigger is worse)

# shorten name of df!
eval.scaled.signed <- scaled.model.evaluation %>%
  mutate(bias.scaled = -bias.scaled,
         ar.scaled = -ar.scaled)

# next step - set smallest value for each region.window - weighting variable (bias, corr, ar) combination = 1

# make function to set smallest value = 1
fmin <- function(x) (1 + (x - min(x)))

sc.bias.min <- eval.scaled.signed %>%
  split(.$region.window) %>%
  map_df(~ fmin(.$bias.scaled)) %>%
  mutate(model = cmip.models) %>%
  pivot_longer(cols = -model, values_to = "bias.scaled.min", names_to = "region.window") 

sc.correlation.min <- eval.scaled.signed %>%
  mutate(region.window = paste(region, window, sep = "-")) %>%
  split(.$region.window) %>%
  map_df(~ fmin(.$correlation.scaled)) %>%
  mutate(model = cmip.models) %>%
  pivot_longer(cols = -model, values_to = "correlation.scaled.min", names_to = "region.window") 

sc.ar.min <- eval.scaled.signed %>%
  mutate(region.window = paste(region, window, sep = "-")) %>%
  split(.$region.window) %>%
  map_df(~ fmin(.$ar.scaled)) %>%
  mutate(model = cmip.models) %>%
  pivot_longer(cols = -model, values_to = "ar.scaled.min", names_to = "region.window") 

# recombine as final df of model.weights

model.weights <- model.evaluation %>%
  mutate(region.window = paste(region, window, sep = "-")) %>% 
  left_join(., sc.bias.min) %>%
  left_join(., sc.correlation.min) %>%
  left_join(., sc.ar.min) 

# plot to check

ggplot(model.weights, aes(bias, bias.scaled.min)) +
  geom_point() +
  facet_grid(region ~ window)

ggplot(model.weights, aes(correlation, correlation.scaled.min)) +
  geom_point() +
  facet_grid(region ~ window)

ggplot(model.weights, aes(ar, ar.scaled.min)) +
  geom_point() +
  facet_grid(region ~ window)

# everything good!


# and finally, create total weight and scale for region-window combinations
model.weights <- model.weights %>%
  mutate(total.weight = bias.scaled.min^2 + correlation.scaled.min^2 + ar.scaled.min^2) 

sc.total.final <- model.weights %>%
  split(.$region.window) %>%
  map_df(~ scale(.$total.weight, center = F)) %>%
  mutate(model = cmip.models) %>%
  pivot_longer(cols = -model, values_to = "scaled.total.weight", names_to = "region.window") 

model.weights <- model.weights %>%
  left_join(., sc.total.final)

# final plot to check
ggplot(model.weights, aes(total.weight, scaled.total.weight)) +
  geom_point() +
  facet_wrap(~region.window)

# plot
order <- data.frame(region = regions,
                    order = 1:6)

plot.weights <- model.weights %>%
  left_join(., order) %>%
  mutate(region = reorder(region, order))

ggplot(plot.weights, aes(model, scaled.total.weight)) +
  geom_bar(stat = "identity", fill = "dark grey", color = "black") +
  facet_grid(region ~ window) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("./CMIP6/figs/model_weights_by_region_window.png", width = 8, height = 10, units = 'in')

ggplot(plot.weights, aes(scaled.total.weight)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(region ~ window) 


ggsave("./CMIP6/figs/model_weight_histograms_by_region_window.png", width = 6, height = 10, units = 'in')

# save

write.csv(model.weights, "./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv", row.names = F)
