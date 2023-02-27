## summarize events across each model-region combination
## to estimate preindustrial probability

library(tidyverse)

## STEP 1 ----------------------------------------------------
# calculate ERSST anomalies wrt 1854-1949 for each region
# this was done ERSST summaries.R"

## STEP 2 -----------------------
# calculate anomaly wrt 1850-1949 for each model preindustrial and hist_ssp585
# this was done in "CMIP6 processing.R"

## STEP 3 ---------------------------------------------
# Record the positive and negative outcomes for each preindustrial run
# for each model-anomaly combination (>= this anomaly or not)
# using years between 1950 and warming = 1.0 as "present"

# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# limit ERSST anomalies to 1950-2021 (period for attribution runs)
ersst.anom <- ersst.anom %>%
  filter(year %in% 1950:2021)

# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv")

# # load CMIP6 model weights
# model.weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv") 

# # load estimated warming level timing for each model
# timing <- read.csv("./CMIP6/summaries/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# get vector of regions
regions <- unique(cmip.anom$region)

# create df to catch outcomes for preindustrial runs
preindustrial.outcomes <- data.frame()

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
    
    # separate this region from ersst.anom
    
    ersst.temp <- ersst.anom %>%
      filter(region == regions[j])

      
    # loop through each year of observation
    for(k in 1:nrow(ersst.temp)){ # start k loop (years)
    # k <- 5
    
    # record outcome for annual unsmoothed, annual 2-yr running mean, and annual 3-yr running mean
    
    annual.1yr <- ifelse(pre.temp$annual.unsmoothed >= ersst.temp$annual.anomaly.unsmoothed[k], 1, 0)
    
    annual.2yr <- ifelse(pre.temp$annual.two.yr.running.mean >= ersst.temp$annual.anomaly.two.yr.running.mean[k], 1, 0)
    
    annual.3yr <- ifelse(pre.temp$annual.three.yr.running.mean >= ersst.temp$annual.anomaly.three.yr.running.mean[k], 1, 0)
    
    # calculate prob for winter unsmoothed, winter 2-yr running mean, and winter 3-yr running mean
    
    winter.1yr <- ifelse(pre.temp$winter.unsmoothed >= ersst.temp$winter.anomaly.unsmoothed[k], 1, 0)
    
    winter.2yr <- ifelse(pre.temp$winter.two.yr.running.mean >= ersst.temp$winter.anomaly.two.yr.running.mean[k], 1, 0)
    
    winter.3yr <- ifelse(pre.temp$winter.three.yr.running.mean >= ersst.temp$winter.anomaly.three.yr.running.mean[k], 1, 0) 
  
    # add to df
    preindustrial.outcomes <- rbind(preindustrial.outcomes,
                              data.frame(model = models[i],
                                         period = "preindustrial",
                                         region = regions[j],
                                         ersst.year = ersst.temp$year[k],
                                         
                                         annual.anomaly.1yr = ersst.temp$annual.anomaly.unsmoothed[k],
                                         annual.1yr.events = annual.1yr,
                                         
                                         annual.anomaly.2yr = ersst.temp$annual.anomaly.two.yr.running.mean[k],
                                         annual.2yr.events = annual.2yr,
                                         
                                         annual.anomaly.3yr = ersst.temp$annual.anomaly.three.yr.running.mean[k],
                                         annual.3yr.events = annual.3yr,
                                         
                                         winter.anomaly.1yr = ersst.temp$winter.anomaly.unsmoothed[k],
                                         winter.1yr.events = winter.1yr,
                                         
                                         winter.anomaly.2yr = ersst.temp$winter.anomaly.two.yr.running.mean[k],
                                         winter.2yr.events = winter.2yr,
                                         
                                         winter.anomaly.3yr = ersst.temp$winter.anomaly.three.yr.running.mean[k],
                                         winter.3yr.events = winter.3yr))

      } # close k loop (ersst years)
    
    } # close j loop (regions)
  
  } # close i loop (models)

nrow(preindustrial.outcomes) # 2.48 M rows!


head(preindustrial.outcomes) # this is right - each model/period/region/year associated with all 200 outcomes 

# break into separate objects for each region and save

for(i in 1:length(regions)){
  
  temp <- preindustrial.outcomes %>%
    filter(region == regions[i]) 
  
  write.csv(temp, file = paste("./CMIP6/summaries/", regions[i], "_preindustrial_outcomes.csv", sep = ""), row.names = F)
  
  
}

## add 2021 and 2022 data (2021 was previously incomplete for running means) --------------------

## reload data
# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# get vector of regions
regions <- unique(cmip.anom$region)

# reload for each region and combine
preindustrial.outcomes <- data.frame()

for(i in 1:length(regions)){
  
temp <-  read.csv(file = paste("./CMIP6/summaries/", regions[i], "_preindustrial_outcomes.csv", sep = ""))

preindustrial.outcomes <- rbind(preindustrial.outcomes, temp)  
  
}

## loop through and calculate 2021 / 2022 values

# create df to catch outcomes for preindustrial runs
preindustrial.outcomes.21.22 <- data.frame()

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
    
    # separate this region from ersst.anom
    
    ersst.temp <- ersst.anom %>%
      filter(region == regions[j])
    
    
    # loop through each year of observation
    for(k in (nrow(ersst.temp)-1):nrow(ersst.temp)){ # start k loop (final two years)
      # k <- 5
      
      # record outcome for annual unsmoothed, annual 2-yr running mean, and annual 3-yr running mean
      
      annual.1yr <- ifelse(pre.temp$annual.unsmoothed >= ersst.temp$annual.anomaly.unsmoothed[k], 1, 0)
      
      annual.2yr <- ifelse(pre.temp$annual.two.yr.running.mean >= ersst.temp$annual.anomaly.two.yr.running.mean[k], 1, 0)
      
      annual.3yr <- ifelse(pre.temp$annual.three.yr.running.mean >= ersst.temp$annual.anomaly.three.yr.running.mean[k], 1, 0)
      
      # calculate prob for winter unsmoothed, winter 2-yr running mean, and winter 3-yr running mean
      
      winter.1yr <- ifelse(pre.temp$winter.unsmoothed >= ersst.temp$winter.anomaly.unsmoothed[k], 1, 0)
      
      winter.2yr <- ifelse(pre.temp$winter.two.yr.running.mean >= ersst.temp$winter.anomaly.two.yr.running.mean[k], 1, 0)
      
      winter.3yr <- ifelse(pre.temp$winter.three.yr.running.mean >= ersst.temp$winter.anomaly.three.yr.running.mean[k], 1, 0) 
      
      # add to df
      preindustrial.outcomes.21.22 <- rbind(preindustrial.outcomes.21.22,
                                      data.frame(model = models[i],
                                                 period = "preindustrial",
                                                 region = regions[j],
                                                 ersst.year = ersst.temp$year[k],
                                                 
                                                 annual.anomaly.1yr = ersst.temp$annual.anomaly.unsmoothed[k],
                                                 annual.1yr.events = annual.1yr,
                                                 
                                                 annual.anomaly.2yr = ersst.temp$annual.anomaly.two.yr.running.mean[k],
                                                 annual.2yr.events = annual.2yr,
                                                 
                                                 annual.anomaly.3yr = ersst.temp$annual.anomaly.three.yr.running.mean[k],
                                                 annual.3yr.events = annual.3yr,
                                                 
                                                 winter.anomaly.1yr = ersst.temp$winter.anomaly.unsmoothed[k],
                                                 winter.1yr.events = winter.1yr,
                                                 
                                                 winter.anomaly.2yr = ersst.temp$winter.anomaly.two.yr.running.mean[k],
                                                 winter.2yr.events = winter.2yr,
                                                 
                                                 winter.anomaly.3yr = ersst.temp$winter.anomaly.three.yr.running.mean[k],
                                                 winter.3yr.events = winter.3yr))
      
    } # close k loop (ersst years)
    
  } # close j loop (regions)
  
} # close i loop (models)

# now remove 2021 from original and join with the two new years
preindustrial_original <- preindustrial.outcomes %>%
  filter(ersst.year != 2021)

preindustrial_updated <- rbind(preindustrial_original, preindustrial.outcomes.21.22)

# break into separate objects for each region and save

for(i in 1:length(regions)){
  
  temp <- preindustrial_updated %>%
    filter(region == regions[i]) 
  
  write.csv(temp, file = paste("./CMIP6/summaries/", regions[i], "_preindustrial_outcomes_updated.csv", sep = ""), row.names = F)
  
  
}

## recalculate North Pacific values to ensure that SST for latitude south of 20N was excluded -------------------------

## reload data
# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# limit ERSST anomalies to 1950-2022 (period for attribution runs)
ersst.anom <- ersst.anom %>%
  filter(year %in% 1950:2022)

# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# get vector of regions
regions <- unique(cmip.anom$region)

# create df to catch outcomes for preindustrial runs
preindustrial.NP <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # loop through each region
  for(j in 6) { # start j loop (regions) # locking at region == "North Pacific
    # j <- 1
    
    # separate model and region of interest
    pre.temp <- cmip.anom %>% 
      filter(experiment == "piControl",
             model == models[i],
             region == regions[j])
    
    # separate this region from ersst.anom
    
    ersst.temp <- ersst.anom %>%
      filter(region == regions[j])
    
    
    # loop through each year of observation
    for(k in 1:nrow(ersst.temp)){ # start k loop (years)
      # k <- 5
      
      # record outcome for annual unsmoothed, annual 2-yr running mean, and annual 3-yr running mean
      
      annual.1yr <- ifelse(pre.temp$annual.unsmoothed >= ersst.temp$annual.anomaly.unsmoothed[k], 1, 0)
      
      annual.2yr <- ifelse(pre.temp$annual.two.yr.running.mean >= ersst.temp$annual.anomaly.two.yr.running.mean[k], 1, 0)
      
      annual.3yr <- ifelse(pre.temp$annual.three.yr.running.mean >= ersst.temp$annual.anomaly.three.yr.running.mean[k], 1, 0)
      
      # calculate prob for winter unsmoothed, winter 2-yr running mean, and winter 3-yr running mean
      
      winter.1yr <- ifelse(pre.temp$winter.unsmoothed >= ersst.temp$winter.anomaly.unsmoothed[k], 1, 0)
      
      winter.2yr <- ifelse(pre.temp$winter.two.yr.running.mean >= ersst.temp$winter.anomaly.two.yr.running.mean[k], 1, 0)
      
      winter.3yr <- ifelse(pre.temp$winter.three.yr.running.mean >= ersst.temp$winter.anomaly.three.yr.running.mean[k], 1, 0) 
      
      # add to df
      preindustrial.NP <- rbind(preindustrial.NP,
                                      data.frame(model = models[i],
                                                 period = "preindustrial",
                                                 region = regions[j],
                                                 ersst.year = ersst.temp$year[k],
                                                 
                                                 annual.anomaly.1yr = ersst.temp$annual.anomaly.unsmoothed[k],
                                                 annual.1yr.events = annual.1yr,
                                                 
                                                 annual.anomaly.2yr = ersst.temp$annual.anomaly.two.yr.running.mean[k],
                                                 annual.2yr.events = annual.2yr,
                                                 
                                                 annual.anomaly.3yr = ersst.temp$annual.anomaly.three.yr.running.mean[k],
                                                 annual.3yr.events = annual.3yr,
                                                 
                                                 winter.anomaly.1yr = ersst.temp$winter.anomaly.unsmoothed[k],
                                                 winter.1yr.events = winter.1yr,
                                                 
                                                 winter.anomaly.2yr = ersst.temp$winter.anomaly.two.yr.running.mean[k],
                                                 winter.2yr.events = winter.2yr,
                                                 
                                                 winter.anomaly.3yr = ersst.temp$winter.anomaly.three.yr.running.mean[k],
                                                 winter.3yr.events = winter.3yr))
      
    } # close k loop (ersst years)
    
  } # close j loop (regions)
  
} # close i loop (models)


# check


check <- preindustrial.NP %>%
  group_by(period, region, ersst.year) %>%
  summarise(count = n())
  

# and save 
write.csv(preindustrial.NP, "./CMIP6/summaries/North_Pacific_preindustrial_outcomes_updated.csv", row.names = F)

check.ebs <- read.csv("./CMIP6/summaries/Eastern_Bering_Sea_preindustrial_outcomes_updated.csv")

check2 <- check.ebs %>%
  group_by(period, region, ersst.year) %>%
  summarise(count = n())
