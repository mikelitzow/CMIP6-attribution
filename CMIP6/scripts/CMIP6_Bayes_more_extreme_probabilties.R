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

# load model warming weights (based on prediction of experienced warming)
warming_weights <- read.csv("./CMIP6/summaries/N_Pac_warming_model_weights.csv")

warming_weights <- warming_weights %>%
  rename(warming_weight = weight) %>%
  select(model, warming_weight)

weights <- left_join(regional_weights, warming_weights) %>%
  mutate(total_weight = regional_weight * warming_weight)


# plot to examine
ggplot(weights, aes(regional_weight, warming_weight)) +
  geom_point() +
  facet_wrap(~region)

ggplot(weights, aes(total_weight)) +
  geom_histogram(fill = "grey", color = "black", bins = 20) +
  facet_wrap(~region)

extremes <- left_join(extremes, weights) %>%
  mutate(model_fac = as.factor(model))

# get vector of regions
regions <- unique(extremes$region)
  
## brms: setup ---------------------------------------------
  
form <-  bf(count | trials(N) + weights(total_weight, scale = TRUE) ~
                period + (1 | model_fac))

# loop through each region and fit model

# for(i in 1:length(regions)){

for(i in 1:2){

extremes_brms <- brm(form,
                 data = extremes[extremes$region == regions[i],],
                 family = binomial(link = "logit"),
                 seed = 1234,
                 cores = 4, chains = 4, iter = 6000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.9, max_treedepth = 14))
  
saveRDS(extremes_brms, paste("./CMIP6/brms_output/",  regions[i], "_extremes_binomial.rds", sep = ""))

}


# evaluate all six regional models

i <- 6

model <- readRDS(paste("./CMIP6/brms_output/", regions[i], "_extremes_binomial.rds", sep = ""))

check_hmc_diagnostics(model$fit)
neff_lowest(model$fit) # N. Pac and EBS are low
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

# set regions and periods in order
region.order <- data.frame(region = regions,
                         region.order = 1:6)

plot.dat <- left_join(plot.dat, region.order) %>%
  mutate(region = reorder(region, region.order))

period.order <- data.frame(period = unique(plot.dat$period),
                           period.order = 1:5)

plot.dat <- left_join(plot.dat, period.order) %>%
  mutate(period = reorder(period, period.order))


ggplot(plot.dat, aes(period, prob)) +
  geom_col(fill = "grey") +
  geom_errorbar(aes(x = period, ymin = lower, ymax = upper)) +
  facet_wrap(~region, scales = "free_y")
