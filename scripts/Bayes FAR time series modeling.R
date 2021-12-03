## Cod FAR
## model the FAR time series to account for model error

library(ggplot2)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
source("./scripts/stan_utils.R")

## Read in data --------------------------------------------
obs <- read.csv("./summaries/ersst_scaled_GOA_output.csv", row.names=1)
mod <- read.csv("./summaries/selected_model_scaled_GOA_output.csv", row.names=1)

## calculate FAR for observed SST from each model -----------

# separate preindustrial and hist_ssp585 runs for convenience

preindustrial <- mod %>%
  dplyr::filter(name == "piControl") 

hist_ssp585 <- preindustrial <- mod %>%
  dplyr::filter(name == "hist_ssp585")

# get vector of model names
models <- unique(preindustrial$model)

# object to catch FAR values for each model
observed.FAR <- data.frame()

# loop through each model
for(i in 1:length(models)){
  # i <- 1
  
  temp.out <- vector() # place to cache output
  
  # separate model of interest
  pre.temp <- preindustrial %>% 
    dplyr::filter(model == models[i])
  
  hist_585.temp <- hist_ssp585 %>%
    dplyr::filter(model == models[i], year %in% 1981:2020) # so here we're using 1981-2020 for the present day period for FAR
  
  for(j in 1:nrow(obs)){
    # j <- 120
    FAR <- 1 - (sum(pre.temp$anomaly >= obs$anomaly[j])/length(pre.temp$anomaly))/
      (sum(hist_585.temp$anomaly >=  obs$anomaly[j])/length(hist_585.temp$anomaly))
    
    temp.out <- c(temp.out, FAR)
    
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