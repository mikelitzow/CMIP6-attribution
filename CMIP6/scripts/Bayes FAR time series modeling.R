## Bayesian models to estimate GOA FAR

library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)
source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## annual SST means --------------------

# load
nc <- nc_open("./CMIP6/data/nceiErsstv5_5d38_fa81_4a70.nc")

# process
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- years(d)

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

SST <- aperm(SST, 3:1)  

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check

# need to drop Bristol Bay cells
BB <- c("N58E200", "N58E202", "N56E200")
SST[,BB] <- NA

# trim spatial area

drop <-  lat < 56 
SST[,drop] <- NA

drop <-  lat > 55 & lon == 198
SST[,drop] <- NA

# and check
temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(195,230), ylim=c(52,62))
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")

# calculate monthly mean
obs.sst <- rowMeans(SST, na.rm = T)

# and annual observed means
ann.sst <- tapply(obs.sst, as.numeric(as.character(yr)), mean)

## STEP 1 ----------------------------------------------------

# calculate GOA anomalies wrt 1950-1999

# limit to 1950-2021 (years we're confident in)  and scale wrt 1950-1999

climatology.mean <- mean(ann.sst[names(ann.sst) %in% 1950:1999])
climatology.sd <- sd(ann.sst[names(ann.sst) %in% 1950:1999])

obs.sst <- data.frame(year = 1950:2021,
                      sc.sst = (ann.sst[names(ann.sst) %in% 1950:2021] - climatology.mean) / climatology.sd)

# combine with 2-year and 3-yr rolling means

obs.sst$sc.sst2 <- rollmean(obs.sst$sc.sst, 2, fill = NA, align = "right") # corresponds to year before and year of ocean entry!

obs.sst$sc.sst3 <- rollmean(obs.sst$sc.sst, 3, fill = NA) # year before, year of, year after ocean entry

plot.sst <- obs.sst %>%
  pivot_longer(cols = -year)

ggplot(plot.sst, aes(year, value, color = name)) +
  geom_line() +
  geom_hline(yintercept = 0)


## STEP 2 -----------------------
# calculate probability wrt 1950-1999 for each model preindustrial and hist_ssp585
# N.B. - historical 1950-1999 period is climatology for both preindustrial and hist_ssp585

# load monthly sst model runs

mod.dat <- read.csv("./CMIP6/summaries/GOA_monthly_sst_piControl_hist585.csv")

# get annual means
yr <- rep(1850:2099, each = 12)

ff <- function(x) tapply(x, yr, mean)

mod.dat <- mod.dat %>%
  select(-date)

annual.sst <- as.data.frame(apply(mod.dat, 2, ff))

# change to C if needed
change <- colMeans(annual.sst) > 200

annual.sst[,change] <- annual.sst[,change] - 273.15

# separate preindustrial and hist_ssp585 runs 
keep <- grep("piControl", names(annual.sst))

preindustrial <- (annual.sst[,keep])

# clean up names
names(preindustrial) <- str_replace(names(preindustrial), ".nc_piControl", "")

keep <- grep("hist", names(annual.sst))

hist_ssp585 <- annual.sst[,keep]

# clean up names
names(hist_ssp585) <- str_replace(names(hist_ssp585), ".nc_hist_ssp585", "")


# change to anomalies to wrt historical 1950-1999
yr <- 1850:2099

# calculate climatology (mean and SD) for 1950-1999 historical runs
ff <- function(x) mean(x[yr %in% 1950:1999])

hist.clim.mean <- apply(hist_ssp585, 2, ff)

ff <- function(x) sd(x[yr %in% 1950:1999])

hist.clim.sd <- apply(hist_ssp585, 2, ff)

# confirm names are lined up
identical(names(hist.clim.mean), names(hist_ssp585))
identical(names(hist.clim.sd), names(hist_ssp585))

# now calculate anomalies

for(i in 1:length(hist.clim.mean)){
  
  hist_ssp585[,i] <- (hist_ssp585[,i] - hist.clim.mean[i]) / hist.clim.sd[i]
  
}

# check 
colMeans(hist_ssp585[yr %in% 1950:1999,]) # perfecto

# now the same with preindustrial runs

# confirm names are lined up
identical(names(hist.clim.mean), names(preindustrial))
identical(names(hist.clim.sd), names(preindustrial))

# now calculate anomalies

for(i in 1:length(hist.clim.mean)){
  
  preindustrial[,i] <- (preindustrial[,i] - hist.clim.mean[i]) / hist.clim.sd[i]
  
}

hist(colMeans(preindustrial[yr %in% 1950:1999,]))

# pivot longer
preindustrial <- preindustrial %>%
  mutate(year = 1850:2099) %>%
  pivot_longer(cols = -year) %>%
  rename(model = name,
         goa.sst = value)

hist_ssp585 <- hist_ssp585 %>%
  mutate(year = 1850:2099) %>%
  pivot_longer(cols = -year) %>%
  rename(model = name,
         goa.sst = value)

# plot to check
plot.1 <- preindustrial %>%
  mutate(run = "preindustrial")

plot.2 <- hist_ssp585 %>%
  mutate(run = "historical/ssp585")

plot.dat <- rbind(plot.1, plot.2)

ggplot(plot.dat, aes(year, goa.sst, color = run)) +
  geom_line() +
  facet_wrap(~model) +
  labs(y = "Anomaly") +
  theme(axis.title.x = element_blank())

ggsave("./CMIP6/figs/hist585_vs_preindustrial_anomalies_wrt_1950-1999.png", width = 9, height = 7, units = 'in')

## STEP 3 ---------------------------------------------
# Calculate FAR for each anomaly in the observed time series, 
# for each model using years between 1950 and warming = 1.0 as present

# get vector of model names
models <- unique(preindustrial$model)
preindustrial.prob <- data.frame()

# loop through each model
for(i in 1:length(models)){
  # i <- 1
  
  # separate model of interest
  pre.temp <- preindustrial %>% 
    dplyr::filter(model == models[i])
  
  # calculate 2-yr and 3-yr running means
  pre.temp$two.yr.mean <- zoo::rollmean(pre.temp$goa.sst, 2, fill = NA, align = "right")
  pre.temp$three.yr.mean <- zoo::rollmean(pre.temp$goa.sst, 3, fill = NA, align = "center") # again - year before, year of, year after
  
  # set up temporary objects to catch probability and save relevant anomaly
  preind.1 <- preind.2 <- preind.3 <- preind.anom.1 <- preind.anom.2 <- preind.anom.3 <- NA
  
  for(j in 1:nrow(obs.sst)){
    # j <- 1
    
    # calculate prob for annual, 2-yr, and 3-yr
    
    preind.anom.1[j] <- obs.sst$sc.sst[j]
    preind.1[j] <- sum(pre.temp$goa.sst >= obs.sst$sc.sst[j])/length(pre.temp$goa.sst)
    
    preind.anom.2[j] <- obs.sst$sc.sst2[j] 
    ifelse(is.na(preind.anom.2[j]), preind.2[j] <- NA, preind.2[j] <- sum(pre.temp$two.yr.mean >= obs.sst$sc.sst2[j], na.rm = T)/length(na.omit(pre.temp$two.yr.mean)))
    
    preind.anom.3[j] <- obs.sst$sc.sst3[j]
    ifelse(is.na(preind.anom.3[j]), preind.3[j] <- NA, preind.3[j] <- sum(pre.temp$three.yr.mean >= obs.sst$sc.sst3[j], na.rm = T)/length(na.omit(pre.temp$three.yr.mean)))
    
  }
  
  # add to df
  preindustrial.prob <- rbind(preindustrial.prob,
                              data.frame(model = models[i],
                                         anomaly.1yr = preind.anom.1,
                                         prob.1yr = preind.1,
                                         anomaly.2yr = preind.anom.2,
                                         prob.2yr = preind.2,
                                         anomaly.3yr = preind.anom.3,
                                         prob.3yr = preind.3))
}


# Calculate probability using 1950 through 1.0 degree warming from hist.585 as "present"
# then calculate FAR

# load warming timing for each model
timing <- read.csv("./CMIP6/summaries/model.north.pacific.warming.timing.csv", row.names = 1)

# fix model names to match timing
hist_ssp585$model <- str_replace_all(hist_ssp585$model, "\\.", "-")

# calculate 2-yr and 3-yr rolling means
hist_ssp585$two.yr.mean <- zoo::rollmean(hist_ssp585$goa.sst, 2, fill = NA, align = "right")
hist_ssp585$three.yr.mean <- zoo::rollmean(hist_ssp585$goa.sst, 3, fill = NA, align = "center")

# get vector of model names
models <- unique(hist_ssp585$model)

# check that these match with model names from timing
identical(models, unique(timing$model)) # hot dog


# rename preind.prob
FAR.estimates <- preindustrial.prob

# and fix model names to match others
FAR.estimates$model <- str_replace_all(FAR.estimates$model, "\\.", "-")

# add columns to save present probabilities
FAR.estimates$present.prob.1yr <- FAR.estimates$present.prob.2yr <- FAR.estimates$present.prob.3yr <- NA

# and create df for final FAR results
FAR.final <- data.frame()

# loop through each model
for(i in 1:length(models)){
  # i <- 1
  
  # separate model of interest
  hist.temp <- hist_ssp585 %>% 
    dplyr::filter(model == models[i])
  
  # pull "present" years (0.5 - 1.0 degrees warming)
  use = 1950:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp <- hist.temp %>%
    filter(year %in% use)
  
  # break out the relevant chunk of Far.estimates (with the model of interest)
  FAR.temp <- FAR.estimates %>%
    filter(model == models[i])
  
  # now loop through each observation and calculate present probability
  for(j in 1:nrow(FAR.temp)){
    # j <- 1
    
    # calculate prob for annual, 2-yr, and 3-yr
    
    FAR.temp$present.prob.1yr[j] <- sum(hist.temp$goa.sst >= FAR.temp$anomaly.1yr[j])/length(hist.temp$goa.sst)
    
    FAR.temp$present.prob.2yr[j] <- ifelse(is.na(FAR.temp$anomaly.2yr[j]), NA,
                                           sum(hist.temp$two.yr.mean >= FAR.temp$anomaly.2yr[j], na.rm = T)/length(na.omit(hist.temp$two.yr.mean)))
    
    FAR.temp$present.prob.3yr[j] <- ifelse(is.na(FAR.temp$anomaly.3yr[j]), NA,
                                           sum(hist.temp$three.yr.mean >= FAR.temp$anomaly.3yr[j], na.rm = T)/length(na.omit(hist.temp$three.yr.mean)))
    
  }
  
  # add year!
  FAR.temp$year <- 1950:2021
  
  # and add to final
  FAR.final <- rbind(FAR.final, FAR.temp)
  
} 

# now add FAR values
FAR.final$FAR.1yr <- 1 - FAR.final$prob.1yr / FAR.final$present.prob.1yr

FAR.final$FAR.2yr <- 1 - FAR.final$prob.2yr / FAR.final$present.prob.2yr

FAR.final$FAR.3yr <- 1 - FAR.final$prob.3yr / FAR.final$present.prob.3yr

## NEED TO REMOVE FAR = -Inf!! (i.e., undefined because present probability = 0)
change <- FAR.final == -Inf
sum(change, na.rm = T) # 18 instances
FAR.final[change] <- NA # replace with NA


# plot to check
plot.dat <- FAR.final %>%
  select(model, year, FAR.1yr, FAR.2yr, FAR.3yr) %>%
  pivot_longer(cols = c(-year, -model))

# looks good - some negatives!


ggplot(plot.dat, aes(year, value, color = name)) +
  geom_line() +
  facet_wrap(~model, scales = "free_y")


ggsave("./CMIP6/figs/FAR_0.5-1.0_warming_by.model.png", width = 8, height = 4)


## set up for brms modeling ----------------------------------

# load model weights (based on ar(1), correlation, bias)
weight <- read.csv("./CMIP6/summaries/GOA_model_weights_bias_corr_ar1.csv", row.names = 1)

# add to FAR.final

weight <- weight %>%
  select(name, total) %>%
  rename(model = name,
         weight = total) %>%
  mutate(weight = weight / sd(weight)) # scaling by SD of weights

FAR <- left_join(FAR.final, weight)

# check
sum(is.na(FAR$weight)) # ok
hist(FAR$weight, breaks = 8)

# look at distribution of 2016 and 2019 estimates
check <- FAR.final$FAR.1yr[FAR.final$year %in% c(2016, 2019)]

hist(check, breaks = 40)

### fit Bayesian regression to estimate across CMIP models----------------------

## Check distribution --------------------------
hist(FAR$FAR.1yr, breaks = 50)

## brms: setup ---------------------------------------------

# setup variables - model as factor
FAR$model_fac <- as.factor(FAR$model)

## Define model formula
far_formula <-  bf(FAR.1yr | weights(weight, scale = TRUE) + trunc(ub = 1.04) ~ s(anomaly.1yr, k = 5) + (1 | model_fac))


## fit: brms --------------------------------------

## base model - Gaussian distribution truncated at 1.04, each observation weighted by scaled model weight
far_1yr_base <- brm(far_formula,
                     data = FAR,
                     cores = 4, chains = 4, iter = 6000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 16))

saveRDS(far_1yr_base, file = "./CMIP6/brms_output/far_1yr_base.rds")

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