## estimate projected warming rate for the N. Pacific from CMIP6
## weighting models by their performance wrt 1972-2021 observed warming rate
## for ssp245

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## load model sst and create warming time series wrt 1850-1949 from historical/ssp245 runs --------
sst.245 <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.ssp245.csv")

# rearrange
sst.245 <- sst.245 %>%
  filter(region == "North_Pacific") %>%
  select(year, model, annual.unsmoothed) %>%
  pivot_wider(names_from = model, values_from = annual.unsmoothed)

model.warming.rate <- sst.245

# now get anomaly wrt 1850:1949
for(j in 2:ncol(sst.245)){
  
  model.warming.rate[,j] <- sst.245[,j] - colMeans(sst.245[1:100,])[j]
  
}

model.warming.rate <- model.warming.rate %>%
  pivot_longer(cols = -year)

# load model weights
model.weights <- read.csv("./CMIP6/summaries/N_Pac_warming_model_weights.csv")

# wrangle data

# simplify weights
weights <- model.weights  %>%
  mutate(model_fac = as.factor(prediction_model)) %>%
  rename(weight = normalized_weight) %>%
  select(model_fac, weight)

dat <- model.warming.rate %>%
  mutate(model_fac = as.factor(name)) %>%
  select(-name) %>%
  rename(warming = value) 

levels(weights$model_fac); levels(dat$model_fac)

# remove "_245" from dat$model_fac
dat <- dat %>%
  mutate (model_fac = str_remove_all(model_fac, "_245"))

dat <- left_join(dat, weights)

# check for NA
sum(is.na(dat))

## fit inverse model - year as a function of warming -----------------------------

inverse_formula <-  bf(year | weights(weight) ~ s(warming) + (1 | model_fac))

## Show default priors
get_prior(inverse_formula, dat)

inverse_warming_brm <- brm(inverse_formula,
                           data = filter(dat, year >= 1973), # limit to 1973-on to ease fitting 
                           cores = 4, chains = 4, iter = 5000,
                           save_pars = save_pars(all = TRUE),
                           control = list(adapt_delta = 0.99, max_treedepth = 16))

saveRDS(inverse_warming_brm, file = "./CMIP6/brms_output/inverse_warming_brm_ssp245.rds")

inverse_warming_brm <- readRDS("./CMIP6/brms_output/inverse_warming_brm_ssp245.rds")

check_hmc_diagnostics(inverse_warming_brm$fit)
neff_lowest(inverse_warming_brm$fit) # ??
rhat_highest(inverse_warming_brm$fit)
summary(inverse_warming_brm)
bayes_R2(inverse_warming_brm)
plot(inverse_warming_brm$criteria$loo, "k")
plot(conditional_effects(inverse_warming_brm), ask = FALSE)

trace_plot(inverse_warming_brm$fit)

# predict for 0.5, 1.0, 1.5, 2.0 degrees

new.dat <- data.frame(warming = c(0.5, 1.0, 1.5, 2.0),
                      model_fac = NA, weight = 1)


pred <- posterior_epred(inverse_warming_brm, newdata = new.dat)

## SST anomaly predictions #### 95% CI
ce1s_1 <- conditional_effects(inverse_warming_brm, effect = "warming", re_formula = NA,
                              probs = c(0.025, 0.975), resolution = 10000)

index <- ce1s_1$warming$warming

choose <- c(which.min(abs(index - 0.5)),
            which.min(abs(index - 1.0)),
            which.min(abs(index - 1.5)),
            which.min(abs(index - 2.0)))

pred.plot <- data.frame(warming = c(0.5, 1.0, 1.5, 2.0),
                        year = ce1s_1$warming$estimate__[choose],
                        UCI = ce1s_1$warming$upper__[choose],
                        LCI = ce1s_1$warming$lower__[choose])


ggplot(pred.plot, aes(warming, year)) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
  geom_point(color = "red", size = 4) +
  labs(x = "North Pacific warming (°C)",
       y = "Year reached")

ggsave("./CMIP6/figs/Bayes_estimated_warming_timing_ssp_245.png", width = 3, height = 4, units = 'in') 
