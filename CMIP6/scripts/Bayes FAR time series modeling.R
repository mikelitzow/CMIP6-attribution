## Cod FAR
## model the FAR time series to account for model error

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

# load monthly sst model runs

mod.dat <- read.csv("./CMIP6/summaries/GOA_monthly_sst_piControl_hist585.csv")

# get annual means
yr <- rep(1850:2099, each = 12)

ff <- function(x) tapply(x, yr, mean)

mod.dat <- mod.dat %>%
  select(-date)

annual.sst <- apply(mod.dat, 2, ff)

# change to C if needed
change <- colMeans(annual.sst) > 200

annual.sst[,change] <- annual.sst[,change] - 273.15

# change to anomalies to wrt 1950-1999
yr <- 1850:2099

ff <- function(x) (x - mean(x[yr %in% 1950:1999])) / sd(x[yr %in% 1950:1999])

annual.anom <- apply(annual.sst, 2, ff)
  
# separate preindustrial and hist_ssp585 runs for convenience
keep <- grep("piControl", names(mod.dat))

preindustrial <- as.data.frame(annual.anom[,keep])

# clean up names
names(preindustrial) <- str_replace(names(preindustrial), ".nc_piControl", "")

keep <- grep("hist", names(mod.dat))

hist_ssp585 <- as.data.frame(annual.anom[,keep])

# clean up names
names(hist_ssp585) <- str_replace(names(hist_ssp585), ".nc_hist_ssp585", "")


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


## STEP 3 ---------------------------------------------
# calculate preindustrial probability

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
  
  # set up temporary objects to catch probability
  preind.1 <- preind.2 <- preind.3 <- NA
 
  for(j in 1:nrow(obs.sst)){
    # j <- 1
    
    # calculate prob for annual, 2-yr, and 3-yr
    preind.1[j] <- sum(pre.temp$goa.sst >= obs.sst$sc.sst[j])/length(pre.temp$goa.sst)

    preind.2[j] <- sum(pre.temp$two.yr.mean >= obs.sst$sc.sst2[j], na.rm = T)/length(pre.temp$two.yr.mean)
    
    preind.3[j] <- sum(pre.temp$three.yr.mean >= obs.sst$sc.sst3[j], na.rm = T)/length(pre.temp$three.yr.mean)
    
  }
  
  # add to df
  preindustrial.prob <- rbind(preindustrial.prob,
                              data.frame(model = models[i],
                                         prob.1yr = preind.1,
                                         prob.2yr = preind.2,
                                         prob.3yr = preind.3))
}
 
## STEP 4 -------------------------------------
# Calculate probability at 0.5, 1.0, 1.5, 2.0 warming from hist.585

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

hist_ssp585_prob <- data.frame()

# warming levels
levels <- seq(0.5, 2, by = 0.5)

# make an object to save probabilites
warming.level.prob.temp.1yr <- warming.level.prob.temp.2yr <- warming.level.prob.temp.3yr <- NA

# loop through each model
for(i in 1:length(models)){
  # i <- 1
  
  # separate model of interest
  hist.temp <- hist_ssp585 %>% 
    dplyr::filter(model == models[i])
  
  # pull warming years
  temp <- NA
  
  for(j in 1:length(levels)){
  j <- 1
  temp <- timing %>%
    filter(model == models[i],
           level == levels[j])
  
  # pull out year warming level is reached, and following 9 years
  hist.temp <- hist.temp %>%
    filter(year %in% temp$year:(temp$year + 9))
    
  for(k in 1:nrow(obs.sst)){
    # k <- 1
    
    # calculate prob for annual, 2-yr, and 3-yr
    warming.level.prob.temp.1yr[j] <- sum(hist.temp$goa.sst >= obs.sst$sc.sst[k])/length(hist.temp$goa.sst)
    
    warming.level.prob.temp.2yr[j] <- sum(hist.temp$two.yr.mean >= obs.sst$sc.sst[k])/length(hist.temp$goa.sst)
    
    preind.1[j] <- sum(pre.temp$goa.sst >= obs.sst$sc.sst[j])/length(pre.temp$goa.sst)
    
    preind.2[j] <- sum(pre.temp$two.yr.mean >= obs.sst$sc.sst2[j], na.rm = T)/length(pre.temp$two.yr.mean)
    
    preind.3[j] <- sum(pre.temp$three.yr.mean >= obs.sst$sc.sst3[j], na.rm = T)/length(pre.temp$three.yr.mean)
    
  }
  
  # add to df
  preindustrial.prob <- rbind(preindustrial.prob,
                              data.frame(model = models[i],
                                         prob.1yr = preind.1,
                                         prob.2yr = preind.2,
                                         prob.3yr = preind.3))
}
















ggplot(observed.FAR, aes(year, FAR.3, color = model)) +
  geom_line()

ggsave("./CMIP6/figs/obs.FAR.3.by.model.png", width = 8, height = 4)

# same plot for annual SST FAR values
ggplot(observed.FAR, aes(year, FAR.1, color = model)) +
  geom_line()

ggsave("./CMIP6/figs/obs.FAR.1.by.model.png", width = 8, height = 4)

## now calculate historical / projected FAR for each model

# object to catch FAR values for each model
projected.FAR <- data.frame()


# loop through each model
for(i in 1:length(models)){
  # i <- 1
  
  # separate model of interest
  pre.temp <- preindustrial %>% 
    dplyr::filter(model == models[i])
  
  proj.temp <- hist_ssp585 %>%
    dplyr::filter(model == models[i])
  
  # add 2-yr and 3-yr running means
  pre.temp$two.yr.mean <- zoo::rollmean(pre.temp$anomaly, 2, fill = NA, align = "right")
  pre.temp$three.yr.mean <- zoo::rollmean(pre.temp$anomaly, 3, fill = NA, align = "right")
  
  proj.temp$two.yr.mean <- zoo::rollmean(proj.temp$anomaly, 2, fill = NA, align = "right")
  proj.temp$three.yr.mean <- zoo::rollmean(proj.temp$anomaly, 3, fill = NA, align = "right")
  
  
  for(j in 1:nrow(proj.temp)){
    # j <- 200
    
    # calculate FAR for annual, 2-yr, and 3-yr
    preind.prob <- sum(pre.temp$anomaly >= proj.temp$anomaly[j])/length(pre.temp$anomaly)
    proj.prob <- sum(proj.temp$anomaly >=  proj.temp$anomaly[j])/length(proj.temp$anomaly)
    FAR.1 <- 1 - preind.prob/proj.prob
    
    preind.prob <- sum(na.omit(pre.temp$two.yr.mean) >= na.omit(proj.temp$two.yr.mean[j]))/length(na.omit(pre.temp$two.yr.mean))
    proj.prob <- sum(na.omit(proj.temp$two.yr.mean) >=  na.omit(proj.temp$two.yr.mean[j]))/length(na.omit(proj.temp$two.yr.mean))
    FAR.2 <- 1 - preind.prob/proj.prob
    
    preind.prob <- sum(na.omit(pre.temp$three.yr.mean) >= na.omit(proj.temp$three.yr.mean[j]))/length(na.omit(pre.temp$three.yr.mean))
    proj.prob <- sum(na.omit(proj.temp$three.yr.mean) >=  na.omit(proj.temp$three.yr.mean[j]))/length(na.omit(proj.temp$three.yr.mean))
    FAR.3 <- 1 - preind.prob/proj.prob
    
    projected.FAR <- rbind(projected.FAR, 
                          data.frame(model = models[i],
                                     year = proj.temp$year[j],
                                     FAR.1 = FAR.1,
                                     FAR.2 = FAR.2,
                                     FAR.3 = FAR.3))
    
  }
}

# now sort based on warming timing

warming <- read.csv("./CMIP6/summaries/model.ne.pacific.warming.timing.csv", row.names = 1)

FAR.warming <- data.frame()

models <- unique(warming$model)
levels <- unique(warming$level)

for(i in 1:length(models)){
  # i <- 1
  for(j in 1:length(levels)){
    # j <- 1
    
    # select the model and warming level of interest
    temp <- warming %>%
      dplyr::filter(model == models[i], level == levels[j])
    
    # identify the 10 years after a warming level was reached for a particular model
    time.frame <- temp$year:(temp$year + 9) 
    
    temp.FAR <- projected.FAR %>%
      dplyr::filter(model == models[i],
                    year %in% time.frame)
    
    FAR.warming <- rbind(FAR.warming, 
                   data.frame(model = models[i],
                              level = levels[j],
                              FAR.1 = temp.FAR$FAR.1,
                              FAR.2 = temp.FAR$FAR.2,
                              FAR.3 = temp.FAR$FAR.3))
  }
  
}

# save
write.csv(FAR.warming, "./CMIP6/summaries/projected_FAR.csv")

### fit Bayesian regression to estimate across CMIP models----------------------
## first, observations
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