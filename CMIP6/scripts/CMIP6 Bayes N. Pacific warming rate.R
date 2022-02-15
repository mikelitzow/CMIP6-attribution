## estimate projected warming rate for the N. Pacific from CMIP6
## weighting models by their performance wrt 1972-2021 observed warming rate

library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## load model warming wrt 1850-1949 from historical/ssp585 runs --------
model.warming <- read.csv("./CMIP6/summaries/N_Pac_warming_model_trends.csv")

model.weights <- read.csv("./CMIP6/summaries/N_Pac_warming_model_weights.csv")

# wrangle data

# simplify weights
weights <- model.weights %>%
  mutate(model_fac = as.factor(model)) %>%
  select(model_fac, weight)

# clean up model names in model.warming
names(model.warming) <- str_replace_all(names(model.warming), "\\.", "-")

dat <- model.warming %>%
  pivot_longer(cols = -year) %>%
  mutate(year_fac = as.factor(year),
         model_fac = as.factor(name)) %>%
  select(-year, -name) %>%
  rename(warming = value) 

levels(weights$model_fac); levels(dat$model_fac)

dat <- left_join(dat, weights)

# check for NA
sum(is.na(dat))

## run brms ---------------------------------

warming_formula <-  bf(warming | weights(weight) ~ year_fac + (1 | model_fac))

## Show default priors
get_prior(warming_formula, dat)

warming_brm <- brm(warming_formula,
                    data = dat,
                    cores = 4, chains = 4, iter = 3000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10))

saveRDS(warming_brm, file = "./CMIP6/brms_output/codR_dfa_brm.rds")

codR_dfa_brm <- readRDS("./output/codR_dfa_brm.rds")
check_hmc_diagnostics(codR_dfa_brm$fit)
neff_lowest(codR_dfa_brm$fit)
rhat_highest(codR_dfa_brm$fit)
summary(codR_dfa_brm)
bayes_R2(codR_dfa_brm)
plot(codR_dfa_brm$criteria$loo, "k")
plot(conditional_effects(codR_dfa_brm), ask = FALSE)
y <- dfa$model
yrep_codR_dfa_brm  <- fitted(codR_dfa_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR_dfa_brm[sample(nrow(yrep_codR_dfa_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("codR_dfa_brm")
pdf("./figs/trace_codR_dfa_brm.pdf", width = 6, height = 4)
trace_plot(codR_dfa_brm$fit)
dev.off()
