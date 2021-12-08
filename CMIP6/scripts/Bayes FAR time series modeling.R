## Cod FAR
## model the FAR time series to account for model error

library(ggplot2)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
source("./scripts/stan_utils.R")

theme_set(theme_bw())

## Read in data --------------------------------------------
obs <- read.csv("./CMIP6/summaries/ersst_scaled_GOA_output.csv", row.names=1)

# add 2-yr and 3-yr running means
obs$two.yr.mean <- zoo::rollmean(obs$anomaly, 2, fill = NA, align = "right")
obs$three.yr.mean <- zoo::rollmean(obs$anomaly, 3, fill = NA, align = "right")

mod <- read.csv("./CMIP6/summaries/selected_model_scaled_GOA_output.csv", row.names=1)

## calculate FAR for observed SST from each model -----------

# separate preindustrial and hist_ssp585 runs for convenience

preindustrial <- mod %>%
  dplyr::filter(name == "piControl") 

hist_ssp585 <- mod %>%
  dplyr::filter(name == "hist_ssp585")

# get vector of model names
models <- unique(preindustrial$model)

# object to catch FAR values for each model
observed.FAR <- data.frame()


# loop through each model
for(i in 1:length(models)){
  # i <- 1
  
  # separate model of interest
  pre.temp <- preindustrial %>% 
    dplyr::filter(model == models[i])
  
  # add 2-yr and 3-yr running means
  pre.temp$two.yr.mean <- zoo::rollmean(pre.temp$anomaly, 2, fill = NA, align = "right")
  pre.temp$three.yr.mean <- zoo::rollmean(pre.temp$anomaly, 3, fill = NA, align = "right")
  
 
  for(j in 1:nrow(obs)){
    # j <- 120
    
    # calculate FAR for annual, 2-yr, and 3-yr
    preind.prob <- sum(pre.temp$anomaly >= obs$anomaly[j])/length(pre.temp$anomaly)
    obs.prob <- sum(obs$anomaly >=  obs$anomaly[j])/length(obs$anomaly)
    FAR.1 <- 1 - preind.prob/obs.prob

    
    preind.prob <- sum(na.omit(pre.temp$two.yr.mean) >= obs$two.yr.mean[j])/length(na.omit(pre.temp$two.yr.mean))
    obs.prob <- sum(na.omit(obs$two.yr.mean) >=  obs$two.yr.mean[j])/length(na.omit(obs$two.yr.mean))
    FAR.2 <- 1 - preind.prob/obs.prob
    
    preind.prob <- sum(na.omit(pre.temp$three.yr.mean) >= obs$three.yr.mean[j])/length(na.omit(pre.temp$three.yr.mean))
    obs.prob <- sum(na.omit(obs$three.yr.mean) >=  obs$three.yr.mean[j])/length(na.omit(obs$three.yr.mean))
    FAR.3 <- 1 - preind.prob/obs.prob
        
    observed.FAR <- rbind(observed.FAR, 
                      data.frame(model = models[i],
                                 year = obs$year[j],
                                 FAR.1 = FAR.1,
                                 FAR.2 = FAR.2,
                                 FAR.3 = FAR.3))
    
  }
}
 

ggplot(observed.FAR, aes(year, FAR.3, color = model)) +
  geom_line()

ggsave("./CMIP6/figs/obs.FAR.3.by.model.png", width = 8, height = 4)

# same plot for annual SST FAR values
ggplot(observed.FAR, aes(year, FAR.1, color = model)) +
  geom_line()

ggsave("./CMIP6/figs/obs.FAR.1.by.model.png", width = 8, height = 4)

# and fit Bayesian regression to estimate across CMIP models

# remove na
observed.FAR <- na.omit(observed.FAR)

# change FAR = 1 to FAR = 0.9999  and FAR <= 0 to FAR = 0.0001 for beta distribution

change <- observed.FAR$FAR.3 == 1
observed.FAR$FAR.3[change] <- 0.9999

change <- observed.FAR$FAR.3 <= 0
observed.FAR$FAR.3[change] <- 0.0001

# and set up explanatory variables as factors
observed.FAR$year_fac <- as.factor(observed.FAR$year)
observed.FAR$model_fac <- as.factor(observed.FAR$model)

## Check distribution --------------------------
hist(observed.FAR$FAR.3, breaks = 50)

## brms: setup ---------------------------------------------

## Define model formulas
far_formula <-  bf(FAR.3 ~ year_fac + (1 | model_fac))

## limit to 1970-2020
## (historical period for pollock assessment)
observed.FAR.1970.2020 <- observed.FAR %>%
  dplyr::filter(year >= 1970)

## fit: brms --------------------------------------

## observed time series
obs_far <- brm(far_formula,
                     data = observed.FAR.1970.2020,
                     family = Beta(),
                     cores = 4, chains = 4, iter = 15000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 15))

saveRDS(obs_far, file = "brms_output/obs_far.rds")

obs_far  <- add_criterion(obs_far, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(obs_far, file = "output/obs_far.rds")

obs_far <- readRDS("./brms_output/obs_far.rds")
check_hmc_diagnostics(obs_far$fit)
neff_lowest(obs_far$fit)
rhat_highest(obs_far$fit)
summary(obs_far)
bayes_R2(obs_far)
y <- observed.FAR.1970.2020$FAR.3
yrep_obs_far  <- fitted(obs_far, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_obs_far[sample(nrow(yrep_obs_far), 25), ]) +
  ggtitle("obs_far.3")

## Predicted effects ---------------------------------------

## Year predictions ##

## 95% CI
ce1s_1 <- conditional_effects(obs_far, probs = c(0.025, 0.975))
obs.95 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.95)[3:4] <- c("ymin.95", "ymax.95")

## 90% CI
ce1s_2 <- conditional_effects(obs_far, probs = c(0.05, 0.95))
obs.90 <- ce1s_2$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.90)[3:4] <- c("ymin.90", "ymax.90")

## 80% CI
ce1s_3 <- conditional_effects(obs_far, probs = c(0.1, 0.9))
obs.80 <- ce1s_3$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.80)[3:4] <- c("ymin.80", "ymax.80")


pred.obs <- left_join(obs.95, obs.90)
pred.obs <- left_join(pred.obs, obs.80)
pred.obs$year <- as.numeric(as.character(pred.obs$year_fac))

theme_set(theme_bw())

g1 <- ggplot(pred.obs) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90") +
  geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
  geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  theme(axis.title.x = element_blank()) +
  ylab("Fraction of Attributable Risk") +
  ggtitle("FAR for observed 3-year running mean SST") +
  scale_x_continuous(breaks=seq(1970, 2020, 10)) 

print(g1)

ggsave("./figs/year_predicted_effect_obs_far_3.yr_mean_sst.png", width = 6, height = 4)


## save output

# first, retrieve se__
pred.obs$se__ <- ce1s_1$year_fac$se__

write.csv(pred.obs, "./summaries/observed_3yr_running_mean_sst_FAR_bayes_estimates.csv")


########################-----------------

## THE FOLLOWING SHOULD BE RE-FIT TO HISTORICAL AND SSP_585!

# need to select FAR values based on warming rather than year

warming <- read.csv("./summaries/model.ne.pacific.warming.timing.csv", row.names = 1)
 
FAR.3 <- data.frame()

models <- unique(warming$model)
levels <- unique(warming$level)

for(i in 1:length(models)){
  # i <- 1
  for(j in 1:length(levels)){
    # j <- 1
    # select the model and warming level of interest
    temp <- warming %>%
      dplyr::filter(model == models[i], level == levels[j])
    # identify the 15 years after a warming level was reached for a particular model
    time.frame <- temp$year:(temp$year+14) 
    
    temp.FAR.3 <- observed.FAR %>%
      dplyr::filter(model == models[i],
                    year %in% time.frame)
    
    FAR.3 <- rbind(FAR.3, 
                   data.frame(model = models[i],
                              level = levels[j],
                              FAR.3 = temp.FAR.3$FAR.3))
  }
  
}




temp.df <- data.frame(year = 1900:2020, 
                      model = models[i],
                      FAR = temp.out)
  
observed.FAR <- rbind(observed.FAR, temp.df)

}


## plot to check

ggplot(observed.FAR, aes(year, FAR, color = model)) +
  geom_line()


obs$year_fac <- as.factor(obs$year)
obs$model_fac <- as.factor(obs$model)

mod$year_fac <- as.factor(mod$year)
mod$model_fac <- as.factor(mod$source)

# change FAR = 1 to Far = 0.9999 to allow beta distribution
change <- obs$FAR == 1 
obs$FAR[change] <- 0.9999

change <- mod$FAR == 1 
mod$FAR[change] <- 0.9999

## Check distribution --------------------------
hist(obs$FAR, breaks = 50)

## brms: setup ---------------------------------------------
## This is the observed FAR time series for Fig. 1a
## Define model formulas
far_formula_fixef <-  bf(FAR ~ year_fac + (1 | model_fac))

## fit: brms --------------------------------------

## observed time series
obs_far_fixef <- brm(far_formula_fixef,
                     data = obs,
                     family = Beta(),
                     cores = 4, chains = 4, iter = 15000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 15))
obs_far_fixef  <- add_criterion(obs_far_fixef, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(obs_far_fixef, file = "output/obs_far_fixef.rds")

obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")
check_hmc_diagnostics(obs_far_fixef$fit)
neff_lowest(obs_far_fixef$fit)
rhat_highest(obs_far_fixef$fit)
summary(obs_far_fixef)
bayes_R2(obs_far_fixef)
y <- obs$FAR
yrep_obs_far_fixef  <- fitted(obs_far_fixef, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_obs_far_fixef[sample(nrow(yrep_obs_far_fixef), 25), ]) +
  xlim(0, 500) +
  ggtitle("obs_far_fixef")

## Predicted effects ---------------------------------------

## Year predictions ##

## 95% CI
ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
obs.95 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.95)[3:4] <- c("ymin.95", "ymax.95")

## 90% CI
ce1s_2 <- conditional_effects(obs_far_fixef, probs = c(0.05, 0.95))
obs.90 <- ce1s_2$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.90)[3:4] <- c("ymin.90", "ymax.90")

## 80% CI
ce1s_3 <- conditional_effects(obs_far_fixef, probs = c(0.1, 0.9))
obs.80 <- ce1s_3$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.80)[3:4] <- c("ymin.80", "ymax.80")


pred.obs <- left_join(obs.95, obs.90)
pred.obs <- left_join(pred.obs, obs.80)
pred.obs$year <- as.numeric(as.character(pred.obs$year_fac))

theme_set(theme_bw())

g1 <- ggplot(pred.obs) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90") +
  geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
  geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  geom_hline(yintercept = 0.95, lty=2) +
  theme(axis.title.x = element_blank()) +
  ylab("FAR") +
  scale_x_continuous(breaks=seq(1960, 2020, 10)) 

print(g1)

ggsave("./figs/year_predicted_effect_obs_far.png", width = 4.5, height = 2)