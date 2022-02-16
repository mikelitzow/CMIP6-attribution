## estimate projected warming rate for the N. Pacific from CMIP6
## weighting models by their performance wrt 1972-2021 observed warming rate

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## load model warming wrt 1850-1949 from historical/ssp585 runs --------
model.warming.rate <- read.csv("./CMIP6/summaries/ne_pacific_annual_modeled_sst.csv", row.names = 1)

# load model weights
model.weights <- read.csv("./CMIP6/summaries/N_Pac_warming_model_weights.csv")

# wrangle data

# simplify weights
weights <- model.weights  %>%
  mutate(model_fac = as.factor(model)) %>%
  select(model_fac, weight)

dat <- model.warming.rate %>%
  mutate(model_fac = as.factor(name)) %>%
  select(-name) %>%
  rename(warming = value) 

levels(weights$model_fac); levels(dat$model_fac)

dat <- left_join(dat, weights)

# check for NA
sum(is.na(dat))

## run brms ---------------------------------

warming_formula <-  bf(warming | weights(weight) ~ s(year) + (1 | model_fac))

## Show default priors
get_prior(warming_formula, dat)

warming_brm <- brm(warming_formula,
                    data = dat,
                    cores = 4, chains = 4, iter = 4000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 12))

saveRDS(warming_brm, file = "./CMIP6/brms_output/warming_brm.rds")

warming_brm <- readRDS("./CMIP6/brms_output/warming_brm.rds")

check_hmc_diagnostics(warming_brm$fit)
neff_lowest(warming_brm$fit) # ??
rhat_highest(warming_brm$fit)
summary(warming_brm)
bayes_R2(warming_brm)
plot(warming_brm$criteria$loo, "k")
plot(conditional_effects(warming_brm), ask = FALSE)
y <- dfa$model
yrep_warming_brm  <- fitted(warming_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_warming_brm[sample(nrow(yrep_warming_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("warming_brm")
pdf("./figs/trace_warming_brm.pdf", width = 6, height = 4)
trace_plot(warming_brm$fit)
dev.off()
