## experimental - logistic model to estimate preindustrial and historical 
## probabilities for given sst anomalies

library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
library(tidybayes)

source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## load preindustrial and historical outcome for the GOA
preindustrial <- read.csv("./CMIP6/summaries/Gulf_of_Alaska_preindustrial_outcomes.csv")
historical <- read.csv("./CMIP6/summaries/Gulf_of_Alaska_historical_outcomes.csv")

# combine
outcomes <- rbind(preindustrial, historical)

# load model weights (based on ar(1), correlation, bias)
weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv")

weights <- weights %>%
   filter(window == "annual",
          region == "Gulf_of_Alaska") %>%
   select(model, scaled.total.weight) %>%
   rename(model_weight = scaled.total.weight)

outcomes <- left_join(outcomes, weights)

## brms: setup ---------------------------------------------

# setup variables - model as factor
outcomes$model_fac <- as.factor(outcomes$model)

## fit: brms --------------------------------------

# define model formula

far_formula <-  bf(annual.1yr.events | weights(model_weight, scale = TRUE) ~
                      s(annual.anomaly.1yr, by = period, k = 6) + period + (1 | model_fac))

# run with default priors

far_brms <- brm(far_formula,
                     data = outcomes,
                     family = bernoulli(link = "logit"),
                     cores = 4, chains = 4, iter = 3000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 15))

saveRDS(far_brms, "./CMIP6/brms_output/Gulf_of_Alaska_binomial.rds")