## estimate probability of more extreme annual sst anomalies for different periods 

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

cb <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load extreme outcome counts
extremes <- read.csv("./CMIP6/summaries/more_extreme_annual_anomalies.csv")


# model weights for extreme events in different periods - 
# product of regional weighting (based on ar(1), correlation, bias) and 
# prediction of observed N. Pac. weighting

# load regional model weights 
regional_weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv")
  
regional_weights <- regional_weights %>%
  filter(window == "annual") %>%
  select(model, region, scaled.total.weight) %>%
  rename(regional_weight = scaled.total.weight)

# calculate region-specific model warming weights (based on prediction of experienced warming)
# save these values for future analysis
ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv")

ersst <- ersst %>%
  select(region, year, annual.unsmoothed) %>%
  mutate(model = "ersst")

models <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv")

# combine models and ersst observations into "data"
# data <- models %>% 
#   filter(experiment == "hist_ssp585",
#          year %in% 1950:2021) %>% # note that for regional warming we will calculate anomalies wrt 1950-1999 (beginning of trustworthy ERSST)
#   select(region, year, annual.unsmoothed, model)

data <- models %>% 
  filter(experiment == "hist_ssp585") %>%
  select(region, year, annual.unsmoothed, model)

data <- rbind(data, ersst) 

# calculate 1950:1999 climatology for each model and ersst
climatology <- data %>%
  filter(year %in% 1850:1949) %>%
  group_by(region, model) %>%
  summarize(climatology.mean = mean(annual.unsmoothed), climatology.sd = sd(annual.unsmoothed))

# combine climatology and data, calculate anomalies
data <- left_join(data, climatology) %>%
  mutate(anomaly = (annual.unsmoothed - climatology.mean) / climatology.sd)

# and pivot longer (ersst vs models)
ersst <- data %>%
  filter(model == "ersst") %>%
  select(region, year, anomaly) %>%
  rename(ersst.anomaly = anomaly)

data <- data %>%
  filter(model != "ersst") %>%
  left_join(., ersst)


# loop through and fit linear ersst - model regressions to get weights
regional_warming_weights <- data.frame()

regions <- unique(data$region)
models <- unique(data$model)

for(r in 1:length(regions)){ # loop through regions
  # r <- 1
  
  for(m in 1:length(models)){ # loop through models
    # m <- 1
    
    temp.dat <- data %>%
      filter(region == regions[r],
             model == models[m],
             year %in% 1972:2021)
    
    
    mod <- lm(ersst.anomaly ~ anomaly, data = temp.dat)
    
    regional_warming_weights <- rbind(regional_warming_weights,
                                      data.frame(region = regions[r],
                                                 model = models[m],
                                                 regional_warming_weight = 1 / abs(1-coefficients(mod)[2]))) # inverse of difference from 1!
    }
  
  }


weights <- left_join(regional_weights, regional_warming_weights) %>%
  mutate(total_weight = regional_weight * regional_warming_weight)


# plot to examine
ggplot(weights, aes(regional_weight, regional_warming_weight)) +
  geom_point() +
  facet_wrap(~region, scale = "free_y")

ggplot(weights, aes(total_weight)) +
  geom_histogram(fill = "grey", color = "black", bins = 20) +
  facet_wrap(~region, scale = "free_x") # I think this looks reasonable

extremes <- left_join(extremes, weights) %>%
  mutate(model_fac = as.factor(model))

# get vector of regions
regions <- unique(extremes$region)
  
## brms: setup ---------------------------------------------
  
form <-  bf(count | trials(N) + weights(total_weight, scale = TRUE) ~
                period + (1 | model_fac))

# loop through each region and fit model

# for(i in 3:length(regions)){
  for(i in 6:6){
  #   i <- 1

extremes_brms <- brm(form,
                 data = extremes[extremes$region == regions[i],],
                 family = binomial(link = "logit"),
                 seed = 1234,
                 cores = 4, chains = 4, iter = 15000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.9, max_treedepth = 15))
  
saveRDS(extremes_brms, paste("./CMIP6/brms_output/",  regions[i], "_extremes_binomial.rds", sep = ""))

}


# evaluate all six regional models

i <- 6

model <- readRDS(paste("./CMIP6/brms_output/", regions[i], "_extremes_binomial.rds", sep = ""))

check_hmc_diagnostics(model$fit)
neff_lowest(model$fit) # N. Pac, EBS and GOA are low
rhat_highest(model$fit)
summary(model)
bayes_R2(model) 
trace_plot(model$fit)

# and plot all six
new.dat <- data.frame(period = unique(extremes$period),
                      model = NA,
                      N = 1000) 

plot.dat <- data.frame()

for(i in 1:length(regions)){
# i <- 1

model <- readRDS(paste("./CMIP6/brms_output/", regions[i], "_extremes_binomial.rds", sep = ""))

probs <- posterior_epred(model, newdata = new.dat, re_formula = NA)/1000 # dive by N to get probability

plot.dat <- rbind(plot.dat,
                  data.frame(region = regions[i],
                             period = new.dat$period,
                             prob = apply(probs, 2, median),
                             lower = apply(probs, 2, quantile, probs = 0.025),
                             upper = apply(probs, 2, quantile, probs = 0.975)))
}

# calculate inverse to get expected return time
plot.dat[,c(3:5)] <- 1/plot.dat[,c(3:5)]

# and change values above 10^4 to 10^4

change <- plot.dat[,c(3:5)] > 10^4

plot.dat[,c(3:5)][change] <- 10^4

# set regions and periods in order
region.order <- data.frame(region = regions,
                         region.order = 1:6)

plot.dat <- left_join(plot.dat, region.order) %>%
  mutate(region = reorder(region, region.order))

period.order <- data.frame(period = unique(plot.dat$period),
                           period.order = 1:5)

plot.dat <- left_join(plot.dat, period.order) %>%
  mutate(period = reorder(period, period.order))


# and change labels for facets!
region_names <- c(
  North_Pacific = "North Pacific",
  Eastern_Bering_Sea = "Eastern Bering Sea",
  Gulf_of_Alaska = "Gulf of Alaska",
  British_Columbia_Coast = "British Columbia Coast",
  Northern_California_Current = "Northern California Current",
  Southern_California_Current = "Southern California Current"
)


# region_labeller <- function(variable,value){
#   return(region_names[value])
# }


extremes.plot <- ggplot(plot.dat, aes(period, prob)) +
  geom_errorbar(aes(x = period, ymin = lower, ymax = upper), width = 0.3) +
  geom_point(color = "red", size = 4) +
  facet_wrap(~region, labeller = labeller(region = region_names)) +
  scale_y_continuous(breaks=c(1,10,100,1000,10000),
                     labels = c("1", "10", "100", "1000", ">10,000"),
                     minor_breaks = c(2:9, 
                                      seq(20, 90, by = 10),
                                      seq(200, 900, by = 100),
                                      seq(2000, 9000, by = 1000))) +
  scale_x_discrete(labels = c("Preindustrial", "1950 to 0.5°", "0.5° to 1.0°", "1.0° to 1.5°", "1.5° to 2.0°")) +
  coord_trans(y = "pseudo_log") +
  ylab("Expected return time (years)") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))

extremes.plot

ggsave("./CMIP6/figs/extreme_return_time.png", width = 6, height = 8, units = 'in')        
