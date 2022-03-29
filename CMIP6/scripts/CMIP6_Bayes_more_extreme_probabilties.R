## estimate probability of more extreme annual sst anomalies for different periods 

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

# load extreme outcome counts
extremes <- read.csv("./CMIP6/summaries/more_extreme_annual_anomalies.csv")
  
# load model weights (based on ar(1), correlation, bias)
weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv")
  
weights <- weights %>%
  filter(window == "annual") %>%
  select(model, region, scaled.total.weight) %>%
  rename(model_weight = scaled.total.weight)
  
extremes <- left_join(extremes, weights) %>%
  mutate(model_fac = as.factor(model))
  
## brms: setup ---------------------------------------------
  
form <-  bf(count | trials(N) + weights(model_weight, scale = TRUE) ~
                region:period + (1 | model_fac))
  
extremes_brms <- brm(form,
                 data = extremes,
                 family = binomial(link = "logit"),
                 seed = 1234,
                 cores = 4, chains = 4, iter = 4000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 15))
  
saveRDS(extremes_brms, paste("./CMIP6/brms_output/extremes_binomial.rds", sep = ""))
