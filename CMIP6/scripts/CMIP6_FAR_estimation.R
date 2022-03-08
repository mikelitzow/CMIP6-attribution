## Bayesian models to estimate FAR across different models
## for each region

library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)

source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## STEP 1 ----------------------------------------------------
# calculate ERSST anomalies wrt 1950-1999 for each region
# this was done ERSST summaries.R"

## STEP 2 -----------------------
# calculate anomaly wrt 1950-1999 for each model preindustrial and hist_ssp585
# N.B. - historical 1950-1999 period is climatology for both preindustrial and hist_ssp585

# this was done in "CMIP6 processing.R"

## STEP 3 ---------------------------------------------
# Calculate FAR for each anomaly in the observed time series, 
# for each model using years between 1950 and warming = 1.0 as "present"

# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv")

# load CMIP6 model weights
model.weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv") 

# load estimated warming level timing for each model
warming.timing <- read.csv("./CMIP6/summaries/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# get vector of regions
regions <- unique(cmip.anom$region)

# create df to catch probabilites for preindustrial runs
FAR.estimates <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1

  
  # loop through each region
  for(j in 1:length(regions)) { # start j loop (regions)
    # j <- 1
    
    # separate model and region of interest
    pre.temp <- cmip.anom %>% 
      filter(experiment == "piControl",
             model == models[i],
             region == regions[j])
    
    ersst.temp <- ersst.anom %>%
      filter(region == regions[j])

      
    # loop through each year of observation
    for(k in 1:nrow(ersst.temp)){ # start k loop (years)
    # k <- 2
    
    # calculate prob for annual unsmoothed, annual 2-yr running mean, and annual 3-yr running mean
    
    annual.prob <- sum(pre.temp$annual.unsmoothed >= ersst.temp$annual.anomaly.unsmoothed[k])/length(pre.temp$annual.unsmoothed)
    
    two.yr.prob <- NA
    
    ifelse(is.na(ersst.temp$annual.anomaly.two.yr.running.mean[k]), two.yr.prob <- NA, 
           two.yr.prob <- sum(pre.temp$annual.two.yr.running.mean >= ersst.temp$annual.anomaly.two.yr.running.mean[k], na.rm = T)/length(na.omit(pre.temp$annual.two.yr.running.mean)))
    
    three.yr.prob <- NA
    
    ifelse(is.na(ersst.temp$annual.anomaly.three.yr.running.mean[k]), three.yr.prob <- NA, 
           three.yr.prob <- sum(pre.temp$annual.three.yr.running.mean >= ersst.temp$annual.anomaly.three.yr.running.mean[k], na.rm = T)/length(na.omit(pre.temp$annual.three.yr.running.mean)))

    
    # calculate prob for winter unsmoothed, winter 2-yr running mean, and winter 3-yr running mean
    
    winter.prob <- NA
      
    ifelse(is.na(ersst.temp$winter.anomaly.unsmoothed[k]), winter.prob <- NA,
           winter.prob <- sum(pre.temp$winter.unsmoothed >= ersst.temp$winter.anomaly.unsmoothed[k], na.rm = T)/length(na.omit(pre.temp$winter.unsmoothed)))
    
    two.yr.winter.prob <- NA
    
    ifelse(is.na(ersst.temp$winter.anomaly.two.yr.running.mean[k]), two.yr.winter.prob <- NA, 
           two.yr.winter.prob <- sum(pre.temp$winter.two.yr.running.mean >= ersst.temp$winter.anomaly.two.yr.running.mean[k], na.rm = T)/length(na.omit(pre.temp$winter.two.yr.running.mean)))
    
    three.yr.winter.prob <- NA
    
    ifelse(is.na(ersst.temp$winter.anomaly.three.yr.running.mean[k]), three.yr.winter.prob <- NA, 
           three.yr.winter.prob <- sum(pre.temp$winter.three.yr.running.mean >= ersst.temp$winter.anomaly.three.yr.running.mean[k], na.rm = T)/length(na.omit(pre.temp$winter.three.yr.running.mean)))
    
  
  # add to df
  FAR.estimates <- rbind(FAR.estimates,
                              data.frame(model = models[i],
                                         region = regions[j],
                                         ersst.year = ersst.temp$year[k],
                                         
                                         annual.anomaly.1yr = ersst.temp$annual.anomaly.unsmoothed[k],
                                         preind.prob.annual.1yr = annual.prob,
                                         
                                         annual.anomaly.2yr = ersst.temp$annual.anomaly.two.yr.running.mean[k],
                                         preind.prob.annual.2yr = two.yr.prob,
                                         
                                         annual.anomaly.3yr = ersst.temp$annual.anomaly.three.yr.running.mean[k],
                                         preind.prob.annual.3yr = three.yr.prob,
                                         
                                         winter.anomaly.1yr = ersst.temp$winter.anomaly.unsmoothed[k],
                                         preind.prob.winter.1yr = winter.prob,
                                         
                                         winter.anomaly.2yr = ersst.temp$winter.anomaly.two.yr.running.mean[k],
                                         preind.prob.winter.2yr = two.yr.winter.prob,
                                         
                                         winter.anomaly.3yr = ersst.temp$winter.anomaly.three.yr.running.mean[k],
                                         preind.prob.winter.3yr = three.yr.winter.prob))

      } # close k loop (ersst years)
    
    } # close j loop (regions)
  
  } # close i loop (models)


# Calculate present probability using 1950 through 1.0 degree warming from hist.585 as "present"

# create addition of historical probabilities to add to FAR.estimates

hist.addition <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # loop through each region
  # j <- 1
  for(j in 1:length(regions)){ # start j loop (regions)
    
  # separate model and region of interest
  hist.temp <- cmip.anom %>% 
      filter(experiment == "hist_ssp585",
             model == models[i],
             region == regions[j])
    
  # pull "present" years (1950 to 1.0 degrees warming)
  use = 1950:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp <- hist.temp %>%
    filter(year %in% use)
  
  # break out the relevant chunk of Far.estimates (with the model and region of interest)
  FAR.temp <- FAR.estimates %>%
    filter(model == models[i],
           region == regions[j])
  
  # now loop through each observation and calculate present probability
  for(k in 1:nrow(FAR.temp)){ # start k loop (years)
    # k <- 1
    
    # calculate prob for annual unsmoothed, annual 2-yr running mean, and annual 3-yr running mean
    # names are inconsistent across different dataframes!
    
    annual.prob <- sum(hist.temp$annual.unsmoothed >= FAR.temp$annual.anomaly.1yr[k])/length(hist.temp$annual.unsmoothed)
    
    two.yr.prob <- NA
    
    ifelse(is.na(FAR.temp$annual.anomaly.2yr[k]), two.yr.prob <- NA, 
           two.yr.prob <- sum(hist.temp$annual.two.yr.running.mean >= FAR.temp$annual.anomaly.2yr[k], na.rm = T)/length(na.omit(hist.temp$annual.two.yr.running.mean)))
    
    three.yr.prob <- NA
    
    ifelse(is.na(FAR.temp$annual.anomaly.3yr[k]), three.yr.prob <- NA, 
           three.yr.prob <- sum(hist.temp$annual.three.yr.running.mean >= FAR.temp$annual.anomaly.3yr[k], na.rm = T)/length(na.omit(hist.temp$annual.three.yr.running.mean)))
    
    
    # calculate prob for winter unsmoothed, winter 2-yr running mean, and winter 3-yr running mean
    
    winter.prob <- NA
    
    ifelse(is.na(FAR.temp$winter.anomaly.1yr[k]), winter.prob <- NA,
           winter.prob <- sum(hist.temp$winter.unsmoothed >= FAR.temp$winter.anomaly.1yr[k])/length(hist.temp$winter.unsmoothed))
    
    two.yr.winter.prob <- NA
    
    ifelse(is.na(FAR.temp$winter.anomaly.2yr[k]), two.yr.winter.prob <- NA, 
           two.yr.winter.prob <- sum(hist.temp$winter.two.yr.running.mean >= FAR.temp$winter.anomaly.2yr[k], na.rm = T)/length(na.omit(hist.temp$winter.two.yr.running.mean)))
    
    three.yr.winter.prob <- NA
    
    ifelse(is.na(FAR.temp$winter.anomaly.3yr[k]), three.yr.winter.prob <- NA, 
           three.yr.winter.prob <- sum(hist.temp$winter.three.yr.running.mean >= FAR.temp$winter.anomaly.3yr[k], na.rm = T)/length(na.omit(hist.temp$winter.three.yr.running.mean)))
    
    
    hist.addition <- rbind(hist.addition,
                           data.frame(model = models[i],
                                      region = regions[j],
                                      ersst.year = FAR.temp$ersst.year[k],
                                      
                                      hist.prob.annual.1yr = annual.prob,
                                      
                                      hist.prob.annual.2yr = two.yr.prob,
                                      
                                      hist.prob.annual.3yr = three.yr.prob,
                                      
                                      hist.prob.winter.1yr = winter.prob,
                                      
                                      hist.prob.winter.2yr = two.yr.winter.prob,
                                      
                                      hist.prob.winter.3yr = three.yr.winter.prob))
   } # close k loop (years)
  
  } # close j loop (regions)
  
} # close i loop (models)

# combine preindustrial and historical anomalies in one df

FAR.estimates <- left_join(FAR.estimates, hist.addition)

# now add FAR values

# annual
FAR.estimates$FAR.annual.1yr <- 1 - FAR.estimates$preind.prob.annual.1yr / FAR.estimates$hist.prob.annual.1yr

FAR.estimates$FAR.annual.2yr <- 1 - FAR.estimates$preind.prob.annual.2yr / FAR.estimates$hist.prob.annual.2yr
  
FAR.estimates$FAR.annual.3yr <- 1 - FAR.estimates$preind.prob.annual.3yr / FAR.estimates$hist.prob.annual.3yr

# winter
FAR.estimates$FAR.winter.1yr <- 1 - FAR.estimates$preind.prob.winter.1yr / FAR.estimates$hist.prob.winter.1yr

FAR.estimates$FAR.winter.2yr <- 1 - FAR.estimates$preind.prob.winter.2yr / FAR.estimates$hist.prob.winter.2yr

FAR.estimates$FAR.winter.3yr <- 1 - FAR.estimates$preind.prob.winter.3yr / FAR.estimates$hist.prob.winter.3yr


## NEED TO REMOVE FAR = -Inf!! (i.e., undefined because present probability = 0)
change <- FAR.estimates == -Inf
sum(change, na.rm = T) # 529 instances

# percent of -Inf instances
529 / (nrow(FAR.estimates)*6) # 0.9% 

FAR.estimates[change] <- NA # replace with NA

# plot to check
plot.dat <- FAR.estimates %>%
  select(model, region, ersst.year, FAR.annual.1yr, FAR.annual.2yr, FAR.annual.3yr,
         FAR.winter.1yr, FAR.winter.2yr, FAR.winter.3yr) %>%
  pivot_longer(cols = c(-ersst.year, -model, -region))

# looks good - some negatives!

ggplot(plot.dat, aes(ersst.year, value, color = model)) +
  geom_line() +
  facet_wrap(name~region, scales = "free_y", ncol = 6)


ggsave("./CMIP6/figs/FAR_0.5-1.0_warming_by.model.png", width = 10, height = 8)

# and save
write.csv(FAR.estimates, "./CMIP6/summaries/FAR_estimates.csv", row.names = F)

## set up for brms modeling ----------------------------------

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

## Define model formula
far_formula <-  bf(FAR.annual.1yr | weights(annual.weight, scale = TRUE) + trunc(ub = 1.04) ~ 
                     s(annual.anomaly.1yr, k = 5) + (1 | model_fac))


## fit: brms --------------------------------------


# loop through regions
regions <- unique(FAR$region)

for(i in 1:length(regions)){
  
  temp.FAR <- FAR %>%
    filter(region == regions[i])

## base model - Gaussian distribution truncated at 1.04, each observation weighted by scaled model weight
far_1yr_base <- brm(far_formula,
                     data = temp.FAR,
                     cores = 4, chains = 4, iter = 6000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 16))

saveRDS(far_1yr_base, file = paste("./CMIP6/brms_output/far_1yr_annual_base_", regions[i], ".rds", sep = ""))

}

# far_1yr_base  <- add_criterion(obs_far, c("loo", "bayes_R2"), moment_match = TRUE)
# saveRDS(far_1yr_base, file = "./CMIP6/brms_output/far_1yr_base.rds")

far_1yr_base <- readRDS("./CMIP6/brms_output/far_1yr_base.rds")

check_hmc_diagnostics(far_1yr_base$fit)
neff_lowest(far_1yr_base$fit)
rhat_highest(far_1yr_base$fit)
summary(far_1yr_base)
bayes_R2(far_1yr_base)

plot(conditional_smooths(far_1yr_base), ask = FALSE)

# y <- as.vector(na.omit(FAR$FAR.1yr)) # this does not account for weights - need to check that
# yrep_far_1yr_base  <- fitted(far_1yr_base, scale = "response", summary = FALSE)
# ppc_dens_overlay(y = y, yrep = yrep_far_1yr_base[sample(nrow(yrep_far_1yr_base), 25), ]) +
#   ggtitle("far_1yr_base.3")

## Base model predicted effects ---------------------------------------

## SST anomaly predictions #### 95% CI
ce1s_1 <- conditional_effects(far_1yr_base, effect = "anomaly.1yr", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(far_1yr_base, effect = "anomaly.1yr", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(far_1yr_base, effect = "anomaly.1yr", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$anomaly.1yr
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$anomaly.1yr[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$anomaly.1yr[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$anomaly.1yr[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$anomaly.1yr[["lower__"]]

ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(y = "Fraction of attributable risk", x = "SST anomaly") +
  theme_bw()


ggsave("./CMIP6/figs/far_1yr_base.png", width = 6, height = 4)

# predict

newdata <- obs.sst %>%
  select(year, sc.sst) %>%
  rename(anomaly.1yr = sc.sst)

pred.far_1yr_base <- posterior_epred(far_1yr_base, newdata = newdata, re_formula = NA, resp = "FAR.1yr")

# functions for different credible intervals
f_95_l <- function(x) quantile(x, 0.025)
f_95_u <- function(x) quantile(x, 0.975)

f_90_l <- function(x) quantile(x, 0.05)
f_90_u <- function(x) quantile(x, 0.95)

f_80_l <- function(x) quantile(x, 0.1)
f_80_u <- function(x) quantile(x, 0.9)


# and plot
plot.predict <- data.frame(year = 1950:2021,
                           estimate__ = colMeans(pred.far_1yr_base),
                           lower_95 = apply(pred.far_1yr_base, 2, f_95_l),
                           upper_95 = apply(pred.far_1yr_base, 2, f_95_u),
                           lower_90 = apply(pred.far_1yr_base, 2, f_90_l),
                           upper_90 = apply(pred.far_1yr_base, 2, f_90_u),
                           lower_80 = apply(pred.far_1yr_base, 2, f_80_l),
                           upper_80 = apply(pred.far_1yr_base, 2, f_80_u))


ggplot(plot.predict) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(y = "Fraction of attributable risk", x = "SST anomaly") +
  theme_bw()

ggsave("./CMIP6/figs/predicted_far_1950-2021_far_1yr_base.png", width = 6, height = 4)


## second model - base model + ar() term -----------------------

## Define model formula
far_ar_formula <-  bf(FAR.1yr | weights(weight, scale = TRUE) + trunc(ub = 1.03) ~
                        s(anomaly.1yr, k = 5) + (1 | model_fac) + ar(gr = model_fac)) 

# autocorrelation modeled within each CMIP6 model 

far_1yr_ar <- brm(far_ar_formula,
                    data = FAR,
                    cores = 4, chains = 4, iter = 7000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.999, max_treedepth = 16))

saveRDS(far_1yr_ar, file = "./CMIP6/brms_output/far_1yr_ar.rds")

# far_1yr_ar  <- add_criterion(obs_far, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(far_1yr_ar, file = "./CMIP6/brms_output/far_1yr_ar.rds")

far_1yr_ar <- readRDS("./CMIP6/brms_output/far_1yr_ar.rds")

check_hmc_diagnostics(far_1yr_ar$fit)
neff_lowest(far_1yr_ar$fit)
rhat_highest(far_1yr_ar$fit)
summary(far_1yr_ar)
# bayes_R2(far_1yr_ar)

# y <- as.vector(na.omit(FAR$FAR.1yr)) # this does not account for weights - need to check that
# yrep_far_1yr_ar  <- fitted(far_1yr_ar, scale = "response", summary = FALSE)
# ppc_dens_overlay(y = y, yrep = yrep_far_1yr_ar[sample(nrow(yrep_far_1yr_ar), 25), ]) +
#   ggtitle("far_1yr_ar.3")

## ar() model predicted effects ---------------------------------------

## SST anomaly predictions #### 95% CI
ce1s_1 <- conditional_effects(far_1yr_ar, effect = "anomaly.1yr", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(far_1yr_ar, effect = "anomaly.1yr", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(far_1yr_ar, effect = "anomaly.1yr", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$anomaly.1yr
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$anomaly.1yr[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$anomaly.1yr[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$anomaly.1yr[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$anomaly.1yr[["lower__"]]

ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(y = "Fraction of attributable risk", x = "SST anomaly") +
  theme_bw()


ggsave("./CMIP6/figs/far_1yr_ar.png", width = 6, height = 4)


## model comparison
loo(far_1yr_base, far_1yr_ar)

###################
###################

## 95% CI
ce1s_1 <- conditional_effects(far_1yr_base, probs = c(0.025, 0.975))
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

## now fit Bayes model to projected FAR as a function of NE Pacific warming
FAR.warming <- read.csv("./CMIP6/summaries/projected_FAR.csv", row.names = 1)

FAR.warming$level <- as.factor(FAR.warming$level)

ggplot(FAR.warming, aes(y = FAR.1, group = level)) +
  geom_boxplot()



# change FAR = 1 to FAR = 0.9999 for beta distribution

for(j in 3:5){
change <- FAR.warming[,j] == 1
FAR.warming[change, j] <- 0.9999
}

for(j in 3:5){
  change <- FAR.warming[,j] <= 0
  FAR.warming[change, j] <- 0.0001
}

# and set up explanatory variables as factors
FAR.warming$model_fac <- as.factor(FAR.warming$model)
FAR.warming$level_fac <- as.factor(FAR.warming$level)

## Check distribution --------------------------
hist(FAR.warming$FAR.1, breaks = 50)

## brms: setup ---------------------------------------------

## Define model formulas
far_formula <-  bf(FAR.1 ~ level_fac + (1 | model_fac))


## fit: brms --------------------------------------

## observed time series
proj_far <- brm(far_formula,
               data = FAR.warming,
               family = Beta(),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.99, max_treedepth = 12))

saveRDS(proj_far, file = "brms_output/proj_far.1.rds")


check_hmc_diagnostics(proj_far$fit)
neff_lowest(proj_far$fit)
rhat_highest(proj_far$fit)
summary(proj_far)
bayes_R2(proj_far)
y <- FAR.warming$FAR.1
yrep_proj_far  <- fitted(proj_far, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_proj_far[sample(nrow(yrep_proj_far), 25), ]) +
  ggtitle("proj_far.3")




########################-----------------

## THE FOLLOWING SHOULD BE RE-FIT TO HISTORICAL AND SSP_585!

# starting with annual SST FAR - for each model, calculate FAR from hist_ssp585 vs preindustrial

# need to select FAR values based on warming rather than year

# warming <- read.csv("./CMIP6/summaries/model.ne.pacific.warming.timing.csv", row.names = 1)
#  
# FAR.1 <- data.frame()
# 
# models <- unique(warming$model)
# levels <- unique(warming$level)
# 
# for(i in 1:length(models)){
#   # i <- 1
#   for(j in 1:length(levels)){
#     # j <- 1
#     # select the model and warming level of interest
#     temp <- warming %>%
#       dplyr::filter(model == models[i], level == levels[j])
#     # identify the 15 years after a warming level was reached for a particular model
#     time.frame <- temp$year:(temp$year+14) 
#     
#     temp.FAR.3 <- observed.FAR %>%
#       dplyr::filter(model == models[i],
#                     year %in% time.frame)
#     
#     FAR.3 <- rbind(FAR.3, 
#                    data.frame(model = models[i],
#                               level = levels[j],
#                               FAR.3 = temp.FAR.3$FAR.3))
#   }
#   
# }




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