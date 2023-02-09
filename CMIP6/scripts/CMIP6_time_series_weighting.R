# weight CMIP6 by climatology, change in temperature, annual SD, and annual autocorrelation
# using the approach of Knutti et al. 2017 and time series of annual mean SST, rather than
# gridded / interpolated data
#
# step one - calculate distance from ERSST for each variable

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)

# set palette
new.col <- oceColorsPalette(64)

# set theme
theme_set(theme_bw())

## load model and observation time series ---------------------------
cmip <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv") %>%
  filter(experiment == "hist_ssp585") %>%
  select(region, model, year, annual.unsmoothed) %>%
  rename(annual.sst = annual.unsmoothed)

ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv") %>%
  select(region, year, annual.unsmoothed) %>%
  rename(annual.sst = annual.unsmoothed)


# plot to check
plot <- cmip %>%
  rename(model_sst = annual.sst) %>%
  filter(year %in% 1950:2022) %>%
  left_join(.,ersst) %>%
  rename(observed_sst = annual.sst) 

ggplot(plot, aes(year, model_sst, color = model)) +
  geom_line() +
  facet_wrap(~region, scales = "free_y") +
  geom_line(aes(year, observed_sst), color = "black", size = 1)


## model evaluation -------------------------------------------------

# get vector or models to roll through
models <- unique(cmip$model)

# get regions to roll through
regions <- unique(ersst$region)


## exploratory plots---------------------------------------
plot.out <- data.frame()

for(r in 1:length(regions)){
  # r <- 2
  
  ersst.temp <- ersst %>%
    filter(region == regions[r],
           year %in% 1950:2022)
  
  for(m in 1:length(models)){
   # m <- 1 
  cmip.temp <- cmip %>%
    filter(region == regions[r],
           model == models[m],
           year %in% 1950:2022)  
    
  plot.out <- rbind(plot.out,
                    data.frame(region = regions[r],
                               year = cmip.temp$year,
                               model = models[m],
                               model_sst = cmip.temp$annual.sst,
                               observed_sst = ersst.temp$annual.sst)) 
    
  }
  
}

temp <- plot.out %>%
  filter(region == "Eastern_Bering_Sea") %>%
  pivot_longer(cols = c(model_sst, observed_sst))

ggplot(temp, aes(year, value, color = name)) +
  geom_line() +
  facet_wrap(~model)

# calculate model weighting differences ----------------
# create object to catch results
output <- data.frame()

for(r in 1:length(regions)){
  # r <- 2 
  
  temp.ersst <- ersst %>%
    filter(region == regions[r])
  
  # calculate regional ERSST climatology, SD, AR(1), and trend
  obs.climatology <- mean(temp.ersst$annual.sst[temp.ersst$year %in% 1950:2014])
  obs.sd <- sd(temp.ersst$annual.sst[temp.ersst$year %in% 1950:2014])
  obs.ar <- ar(temp.ersst$annual.sst[temp.ersst$year %in% 1950:2014], order.max = 1, aic = F)$ar
  obs.trend <- summary(lm(annual.sst ~ year, 
                          data = temp.ersst[temp.ersst$year %in% 1973:2022,]))$coefficients[2,1]
  
  for(m in 1:length(models)){
    # m <- 1
    
    temp.cmip <- cmip %>%
      filter(region == regions[r],
             model == models[m]) 
    
    # calculate regional CMIP climatology, SD, AR(1), and trend
    cmip.climatology <- mean(temp.cmip$annual.sst[temp.cmip$year %in% 1950:2014])
    cmip.sd <- sd(temp.cmip$annual.sst[temp.cmip$year %in% 1950:2014])
    cmip.ar <- ar(temp.cmip$annual.sst[temp.cmip$year %in% 1950:2014], order.max = 1, aic = F)$ar
    cmip.trend <- summary(lm(annual.sst ~ year, 
                            data = temp.cmip[temp.cmip$year %in% 1973:2022,]))$coefficients[2,1]
    
    

    # add to output
    output <- rbind(output,
                    data.frame(region = regions[r],
                               model = models[m],
                               climatology_diff = abs(obs.climatology - cmip.climatology),
                               sd_diff = abs(obs.sd - cmip.sd),
                               ar_diff = abs(obs.ar - cmip.ar),
                               trend_diff = abs(obs.trend - cmip.trend)))
  
 
  } # close m loop (models)
} # close r loop (regions) 


# plot to check
plot <- output %>%
  pivot_longer(cols = c(-region, -model))

ggplot(plot, aes(value)) +
  geom_histogram(bins = 12, fill = "grey", color = "black") +
  facet_grid(region~name, scales = "free")

# now save output
write.csv(output, "./CMIP6/summaries/CMIP6_time_series_weights.csv", row.names = F)

