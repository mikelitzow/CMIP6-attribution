## estimate projected warming rate for the N. Pacific from observations

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

# load ersst warming wrt 1854-1949
n.pac.obs.warming <- read.csv("./CMIP6/summaries/north_pacific_annual_observed_sst.csv",
                              row.names = 1)

## fit inverse model - year as a function of warming -----------------------------

inverse_formula <-  bf(year ~ s(ersst.warming))

## Show default priors
get_prior(inverse_formula, n.pac.obs.warming)

inverse_warming_brm <- brm(inverse_formula,
                           data = filter(n.pac.obs.warming, year >= 1972), # limit to 1972-on to ease fitting 
                           cores = 4, chains = 4, iter = 4000,
                           save_pars = save_pars(all = TRUE),
                           control = list(adapt_delta = 0.99, max_treedepth = 16))

saveRDS(inverse_warming_brm, file = "./CMIP6/brms_output/inverse_warming_brm_ersst.rds")

inverse_warming_brm <- readRDS("./CMIP6/brms_output/inverse_warming_brm_ersst.rds")

check_hmc_diagnostics(inverse_warming_brm$fit)
neff_lowest(inverse_warming_brm$fit) # ??
rhat_highest(inverse_warming_brm$fit)
summary(inverse_warming_brm)
bayes_R2(inverse_warming_brm)
plot(inverse_warming_brm$criteria$loo, "k")
plot(conditional_effects(inverse_warming_brm), ask = FALSE)

trace_plot(inverse_warming_brm$fit)

# predict for 0.5, 1.0 degrees

new.dat <- data.frame(ersst.warming = c(0.5, 1.0),
                      model_fac = NA, weight = 1)


pred <- posterior_epred(inverse_warming_brm, newdata = new.dat)


## SST anomaly predictions #### 95% CI
ce1s_1 <- conditional_effects(inverse_warming_brm, effect = "ersst.warming", re_formula = NA,
                              probs = c(0.025, 0.975), resolution = 10000)

index <- ce1s_1$ersst.warming$ersst.warming

choose <- c(which.min(abs(index - 0.5)),
            which.min(abs(index - 1.0)))

pred.plot <- data.frame(ersst.warming = c(0.5, 1.0),
                        year = ce1s_1$ersst.warming$estimate__[choose],
                        UCI = ce1s_1$ersst.warming$upper__[choose],
                        LCI = ce1s_1$ersst.warming$lower__[choose])


ggplot(pred.plot, aes(ersst.warming, year)) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
  geom_point(color = "red", size = 4) +
  labs(x = "North Pacific warming (°C)",
       y = "Year reached")

# save
write.csv(pred.plot, "ERSST_Bayes_warming_timing_estimates.csv", row.names = F)
