# Bayesian models to estimate FAR across CMIP6 models

library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)

source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## set up for brms modeling ----------------------------------

# load FAR estimates
FAR.estimates <- read.csv("./CMIP6/summaries/FAR_estimates.csv")

# load model weights (based on ar(1), correlation, bias)
weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv")

# add to FAR.estimates
weights <- weights %>%
  select(model, region, window, scaled.total.weight) %>%
  rename(weight = scaled.total.weight) %>%
  pivot_wider(names_from = window, values_from = weight) %>%
  rename(annual.weight = annual,
         winter.weight = winter)

FAR <- left_join(FAR.estimates, weights)

# check
sum(is.na(FAR$annual.weight)) # ok
hist(FAR$annual.weight, breaks = 8)

sum(is.na(FAR$winter.weight)) # ok
hist(FAR$winter.weight, breaks = 8)

# look at distribution of 2016 and 2019 estimates
check <- FAR$FAR.annual.1yr[FAR$ersst.year %in% c(2016, 2019)]

hist(check, breaks = 40)

### fit Bayesian regression to estimate across CMIP models----------------------

## Check distribution --------------------------
hist(FAR$FAR.annual.1yr, breaks = 50)

## brms: setup ---------------------------------------------

# setup variables - model as factor
FAR$model_fac <- as.factor(FAR$model)


## fit: brms --------------------------------------

# first, fit separately to North Pacific, using k = 5 (k = 4 did not fit! many divergent transitions, failed convergence)
# Define model formula

# set up vector of regions
regions <- unique(FAR$region)

far_formula <-  bf(FAR.annual.1yr | weights(annual.weight, scale = TRUE) + trunc(ub = 1.01) ~ 
                     s(annual.anomaly.1yr, k = 5) + (1 | model_fac))

  temp.FAR <- FAR %>%
    filter(region == regions[1])
  
  ## base model - Gaussian distribution truncated at 1.0, each observation weighted by scaled model weight
  far_1yr_base <- brm(far_formula,
                      data = temp.FAR,
                      cores = 4, chains = 4, iter = 7000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 15))
  
  saveRDS(far_1yr_base, file = paste("./CMIP6/brms_output/far_1yr_annual_base_", regions[i], ".rds", sep = ""))
  


# fit to remaining regions

## Define model formula
far_formula <-  bf(FAR.annual.1yr | weights(annual.weight, scale = TRUE) + trunc(ub = 1.0) ~ 
                     s(annual.anomaly.1yr, k = 4) + (1 | model_fac))
for(i in 3:length(regions)){
  
  temp.FAR <- FAR %>%
    filter(region == regions[i])
  
  ## base model - Gaussian distribution truncated at 1.0, each observation weighted by scaled model weight
  far_1yr_base <- brm(far_formula,
                      data = temp.FAR,
                      cores = 4, chains = 4, iter = 7000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 15))
  
  saveRDS(far_1yr_base, file = paste("./CMIP6/brms_output/far_1yr_annual_base_", regions[i], ".rds", sep = ""))
  
}

## run model diagnostics ---------------------

# get list of model objects to load

file.list <- NA
for(i in 1:length(regions)){
  
  file.list[i] = paste("./CMIP6/brms_output/far_1yr_annual_base_", regions[i], ".rds", sep = "")
  
}

# loop through each model and run simple diagnostics

for(i in length(file.list)){
  # i <- 1
  
  model.object <- readRDS(file = file.list[i])
  
  print(regions[i])
  
  check_hmc_diagnostics(model.object$fit)
  
  neff_lowest(model.object$fit)
  
  rhat_highest(model.object$fit)
  
 }

# the following model fits have issues to be addressed:
  
  # NCC (5) and SSC (6) have divergent transitions
  
  # NPac (1), EBS(2), BC(4), NCC (5), and SCC (6) have effective sample sizes < 1000 
  
#   
# check_hmc_diagnostics(far_1yr_base$fit)
# neff_lowest(far_1yr_base$fit)
# rhat_highest(far_1yr_base$fit)
# summary(far_1yr_base)
# bayes_R2(far_1yr_base)
# 
# plot(conditional_smooths(far_1yr_base), ask = FALSE)

# y <- as.vector(na.omit(FAR$FAR.1yr)) # this does not account for weights - need to check that
# yrep_far_1yr_base  <- fitted(far_1yr_base, scale = "response", summary = FALSE)
# ppc_dens_overlay(y = y, yrep = yrep_far_1yr_base[sample(nrow(yrep_far_1yr_base), 25), ]) +
#   ggtitle("far_1yr_base.3")

## Plot predicted FAR-SST relationships ---------------------------------------
  
plot.dat <- data.frame()

for(i in 1:length(file.list)){
  # i <- 1
  
  model.object <- readRDS(file = file.list[i])
  
## SST anomaly predictions #### 95% CI
ce1s_1 <- conditional_effects(model.object, effect = "annual.anomaly.1yr", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(model.object, effect = "annual.anomaly.1yr", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(model.object, effect = "annual.anomaly.1yr", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$annual.anomaly.1yr
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$annual.anomaly.1yr[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$annual.anomaly.1yr[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$annual.anomaly.1yr[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$annual.anomaly.1yr[["lower__"]]


dat_ce$region <- regions[i]

  plot.dat <- rbind(plot.dat,
                    dat_ce)



}

# put regions in order
plot.regions <- data.frame(region = regions,
                           order = 1:6)


plot.dat <- left_join(plot.dat, plot.regions)
plot.dat$region <- reorder(plot.dat$region, plot.dat$order)

ggplot(plot.dat,
  aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(y = "Fraction of attributable risk", x = "SST anomaly") +
  facet_wrap(~region) +
  theme_bw()


ggsave("./CMIP6/figs/regional_far_annual_sst_anomaly_unsmoothed.png", width = 9, height = 6, units = 'in')




# 
# 
# ## second model - base model + ar() term -----------------------
# 
# ## Define model formula
# far_ar_formula <-  bf(FAR.1yr | weights(weight, scale = TRUE) + trunc(ub = 1.03) ~
#                         s(anomaly.1yr, k = 5) + (1 | model_fac) + ar(gr = model_fac)) 
# 
# # autocorrelation modeled within each CMIP6 model 
# 
# far_1yr_ar <- brm(far_ar_formula,
#                   data = FAR,
#                   cores = 4, chains = 4, iter = 7000,
#                   save_pars = save_pars(all = TRUE),
#                   control = list(adapt_delta = 0.999, max_treedepth = 16))
# 
# saveRDS(far_1yr_ar, file = "./CMIP6/brms_output/far_1yr_ar.rds")
# 
# # far_1yr_ar  <- add_criterion(obs_far, c("loo", "bayes_R2"), moment_match = TRUE)
# saveRDS(far_1yr_ar, file = "./CMIP6/brms_output/far_1yr_ar.rds")
# 
# far_1yr_ar <- readRDS("./CMIP6/brms_output/far_1yr_ar.rds")
# 
# check_hmc_diagnostics(far_1yr_ar$fit)
# neff_lowest(far_1yr_ar$fit)
# rhat_highest(far_1yr_ar$fit)
# summary(far_1yr_ar)
# # bayes_R2(far_1yr_ar)
# 
# # y <- as.vector(na.omit(FAR$FAR.1yr)) # this does not account for weights - need to check that
# # yrep_far_1yr_ar  <- fitted(far_1yr_ar, scale = "response", summary = FALSE)
# # ppc_dens_overlay(y = y, yrep = yrep_far_1yr_ar[sample(nrow(yrep_far_1yr_ar), 25), ]) +
# #   ggtitle("far_1yr_ar.3")
# 
# ## ar() model predicted effects ---------------------------------------
# 
# ## SST anomaly predictions #### 95% CI
# ce1s_1 <- conditional_effects(far_1yr_ar, effect = "anomaly.1yr", re_formula = NA,
#                               probs = c(0.025, 0.975))
# ## 90% CI
# ce1s_2 <- conditional_effects(far_1yr_ar, effect = "anomaly.1yr", re_formula = NA,
#                               probs = c(0.05, 0.95))
# ## 80% CI
# ce1s_3 <- conditional_effects(far_1yr_ar, effect = "anomaly.1yr", re_formula = NA,
#                               probs = c(0.1, 0.9))
# dat_ce <- ce1s_1$anomaly.1yr
# dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
# dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
# dat_ce[["upper_90"]] <- ce1s_2$anomaly.1yr[["upper__"]]
# dat_ce[["lower_90"]] <- ce1s_2$anomaly.1yr[["lower__"]]
# dat_ce[["upper_80"]] <- ce1s_3$anomaly.1yr[["upper__"]]
# dat_ce[["lower_80"]] <- ce1s_3$anomaly.1yr[["lower__"]]
# 
# ggplot(dat_ce) +
#   aes(x = effect1__, y = estimate__) +
#   geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
#   geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
#   geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
#   geom_line(size = 1, color = "red3") +
#   labs(y = "Fraction of attributable risk", x = "SST anomaly") +
#   theme_bw()
# 
# 
# ggsave("./CMIP6/figs/far_1yr_ar.png", width = 6, height = 4)
# 
# 
# ## model comparison
# loo(far_1yr_base, far_1yr_ar)
# 
# ###################
# ###################
# 
# ## 95% CI
# ce1s_1 <- conditional_effects(far_1yr_base, probs = c(0.025, 0.975))
# obs.95 <- ce1s_1$year_fac %>%
#   select(year_fac, estimate__, lower__, upper__)
# names(obs.95)[3:4] <- c("ymin.95", "ymax.95")
# 
# ## 90% CI
# ce1s_2 <- conditional_effects(obs_far, probs = c(0.05, 0.95))
# obs.90 <- ce1s_2$year_fac %>%
#   select(year_fac, estimate__, lower__, upper__)
# names(obs.90)[3:4] <- c("ymin.90", "ymax.90")
# 
# ## 80% CI
# ce1s_3 <- conditional_effects(obs_far, probs = c(0.1, 0.9))
# obs.80 <- ce1s_3$year_fac %>%
#   select(year_fac, estimate__, lower__, upper__)
# names(obs.80)[3:4] <- c("ymin.80", "ymax.80")
# 
# 
# pred.obs <- left_join(obs.95, obs.90)
# pred.obs <- left_join(pred.obs, obs.80)
# pred.obs$year <- as.numeric(as.character(pred.obs$year_fac))
# 
# theme_set(theme_bw())
# 
# g1 <- ggplot(pred.obs) +
#   aes(x = year, y = estimate__) +
#   geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90") +
#   geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
#   geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
#   geom_line(size = 0.5, color = "red3") +
#   theme(axis.title.x = element_blank()) +
#   ylab("Fraction of Attributable Risk") +
#   ggtitle("FAR for observed 3-year running mean SST") +
#   scale_x_continuous(breaks=seq(1970, 2020, 10)) 
# 
# print(g1)
# 
# ggsave("./figs/year_predicted_effect_obs_far_3.yr_mean_sst.png", width = 6, height = 4)
# 
# 
# ## save output
# 
# # first, retrieve se__
# pred.obs$se__ <- ce1s_1$year_fac$se__
# 
# write.csv(pred.obs, "./summaries/observed_3yr_running_mean_sst_FAR_bayes_estimates.csv")
# 
# ## now fit Bayes model to projected FAR as a function of NE Pacific warming
# FAR.warming <- read.csv("./CMIP6/summaries/projected_FAR.csv", row.names = 1)
# 
# FAR.warming$level <- as.factor(FAR.warming$level)
# 
# ggplot(FAR.warming, aes(y = FAR.1, group = level)) +
#   geom_boxplot()
# 
# 
# 
# # change FAR = 1 to FAR = 0.9999 for beta distribution
# 
# for(j in 3:5){
#   change <- FAR.warming[,j] == 1
#   FAR.warming[change, j] <- 0.9999
# }
# 
# for(j in 3:5){
#   change <- FAR.warming[,j] <= 0
#   FAR.warming[change, j] <- 0.0001
# }
# 
# # and set up explanatory variables as factors
# FAR.warming$model_fac <- as.factor(FAR.warming$model)
# FAR.warming$level_fac <- as.factor(FAR.warming$level)
# 
# ## Check distribution --------------------------
# hist(FAR.warming$FAR.1, breaks = 50)
# 
# ## brms: setup ---------------------------------------------
# 
# ## Define model formulas
# far_formula <-  bf(FAR.1 ~ level_fac + (1 | model_fac))
# 
# 
# ## fit: brms --------------------------------------
# 
# ## observed time series
# proj_far <- brm(far_formula,
#                 data = FAR.warming,
#                 family = Beta(),
#                 cores = 4, chains = 4, iter = 2500,
#                 save_pars = save_pars(all = TRUE),
#                 control = list(adapt_delta = 0.99, max_treedepth = 12))
# 
# saveRDS(proj_far, file = "brms_output/proj_far.1.rds")
# 
# 
# check_hmc_diagnostics(proj_far$fit)
# neff_lowest(proj_far$fit)
# rhat_highest(proj_far$fit)
# summary(proj_far)
# bayes_R2(proj_far)
# y <- FAR.warming$FAR.1
# yrep_proj_far  <- fitted(proj_far, scale = "response", summary = FALSE)
# ppc_dens_overlay(y = y, yrep = yrep_proj_far[sample(nrow(yrep_proj_far), 25), ]) +
#   ggtitle("proj_far.3")
# 
# 
# 
# 
# ########################-----------------
# 
# ## THE FOLLOWING SHOULD BE RE-FIT TO HISTORICAL AND SSP_585!
# 
# # starting with annual SST FAR - for each model, calculate FAR from hist_ssp585 vs preindustrial
# 
# # need to select FAR values based on warming rather than year
# 
# # warming <- read.csv("./CMIP6/summaries/model.ne.pacific.warming.timing.csv", row.names = 1)
# #  
# # FAR.1 <- data.frame()
# # 
# # models <- unique(warming$model)
# # levels <- unique(warming$level)
# # 
# # for(i in 1:length(models)){
# #   # i <- 1
# #   for(j in 1:length(levels)){
# #     # j <- 1
# #     # select the model and warming level of interest
# #     temp <- warming %>%
# #       dplyr::filter(model == models[i], level == levels[j])
# #     # identify the 15 years after a warming level was reached for a particular model
# #     time.frame <- temp$year:(temp$year+14) 
# #     
# #     temp.FAR.3 <- observed.FAR %>%
# #       dplyr::filter(model == models[i],
# #                     year %in% time.frame)
# #     
# #     FAR.3 <- rbind(FAR.3, 
# #                    data.frame(model = models[i],
# #                               level = levels[j],
# #                               FAR.3 = temp.FAR.3$FAR.3))
# #   }
# #   
# # }
# 
# 
# 
# 
# temp.df <- data.frame(year = 1900:2020, 
#                       model = models[i],
#                       FAR = temp.out)
# 
# observed.FAR <- rbind(observed.FAR, temp.df)
# 
# }
# 
# 
# ## plot to check
# 
# ggplot(observed.FAR, aes(year, FAR, color = model)) +
#   geom_line()
# 
# 
# obs$year_fac <- as.factor(obs$year)
# obs$model_fac <- as.factor(obs$model)
# 
# mod$year_fac <- as.factor(mod$year)
# mod$model_fac <- as.factor(mod$source)
# 
# # change FAR = 1 to Far = 0.9999 to allow beta distribution
# change <- obs$FAR == 1 
# obs$FAR[change] <- 0.9999
# 
# change <- mod$FAR == 1 
# mod$FAR[change] <- 0.9999
# 
# ## Check distribution --------------------------
# hist(obs$FAR, breaks = 50)
# 
# ## brms: setup ---------------------------------------------
# ## This is the observed FAR time series for Fig. 1a
# ## Define model formulas
# far_formula_fixef <-  bf(FAR ~ year_fac + (1 | model_fac))
# 
# ## fit: brms --------------------------------------
# 
# ## observed time series
# obs_far_fixef <- brm(far_formula_fixef,
#                      data = obs,
#                      family = Beta(),
#                      cores = 4, chains = 4, iter = 15000,
#                      save_pars = save_pars(all = TRUE),
#                      control = list(adapt_delta = 0.999, max_treedepth = 15))
# obs_far_fixef  <- add_criterion(obs_far_fixef, c("loo", "bayes_R2"), moment_match = TRUE)
# saveRDS(obs_far_fixef, file = "output/obs_far_fixef.rds")
# 
# obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")
# check_hmc_diagnostics(obs_far_fixef$fit)
# neff_lowest(obs_far_fixef$fit)
# rhat_highest(obs_far_fixef$fit)
# summary(obs_far_fixef)
# bayes_R2(obs_far_fixef)
# y <- obs$FAR
# yrep_obs_far_fixef  <- fitted(obs_far_fixef, scale = "response", summary = FALSE)
# ppc_dens_overlay(y = y, yrep = yrep_obs_far_fixef[sample(nrow(yrep_obs_far_fixef), 25), ]) +
#   xlim(0, 500) +
#   ggtitle("obs_far_fixef")
# 
# ## Predicted effects ---------------------------------------
# 
# ## Year predictions ##
# 
# ## 95% CI
# ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
# obs.95 <- ce1s_1$year_fac %>%
#   select(year_fac, estimate__, lower__, upper__)
# names(obs.95)[3:4] <- c("ymin.95", "ymax.95")
# 
# ## 90% CI
# ce1s_2 <- conditional_effects(obs_far_fixef, probs = c(0.05, 0.95))
# obs.90 <- ce1s_2$year_fac %>%
#   select(year_fac, estimate__, lower__, upper__)
# names(obs.90)[3:4] <- c("ymin.90", "ymax.90")
# 
# ## 80% CI
# ce1s_3 <- conditional_effects(obs_far_fixef, probs = c(0.1, 0.9))
# obs.80 <- ce1s_3$year_fac %>%
#   select(year_fac, estimate__, lower__, upper__)
# names(obs.80)[3:4] <- c("ymin.80", "ymax.80")
# 
# 
# pred.obs <- left_join(obs.95, obs.90)
# pred.obs <- left_join(pred.obs, obs.80)
# pred.obs$year <- as.numeric(as.character(pred.obs$year_fac))
# 
# theme_set(theme_bw())
# 
# g1 <- ggplot(pred.obs) +
#   aes(x = year, y = estimate__) +
#   geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90") +
#   geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
#   geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
#   geom_line(size = 0.5, color = "red3") +
#   geom_hline(yintercept = 0.95, lty=2) +
#   theme(axis.title.x = element_blank()) +
#   ylab("FAR") +
#   scale_x_continuous(breaks=seq(1960, 2020, 10)) 
# 
# print(g1)
# 
# ggsave("./figs/year_predicted_effect_obs_far.png", width = 4.5, height = 2)
