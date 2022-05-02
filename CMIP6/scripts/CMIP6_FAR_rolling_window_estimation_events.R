## summarize events across each model-region combination
## to estimate historical probability

library(tidyverse)

theme_set(theme_bw())

## STEP 1 ----------------------------------------------------
# calculate ERSST anomalies wrt 1854-1949 for each region
# this was done ERSST summaries.R"

## STEP 2 -----------------------
# calculate anomaly wrt 1850-1949 for each model preindustrial and hist_ssp585
# this was done in "CMIP6 processing.R"

## STEP 3 ---------------------------------------------
# Record the positive and negative outcomes for each preindustrial run
# for each model-anomaly combination (>= this anomaly or not)
# this was done in "CMIP6_FAR_estimation_events.R"

## STEP 4 ---------------------------------------------
# Record the positive and negative outcomes for each hist_ssp585 run
# for each model-anomaly combination (>= this anomaly or not)
# using 15-year rolling window of estimated warming as "present"
# THIS SCRIPT


# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# limit to 1950-2021 (period of interest for attribution)
ersst.anom <- ersst.anom %>%
  filter(year %in% 1950:2021)

# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# get vector of regions
regions <- unique(cmip.anom$region)

# load CMIP6 model weights
model.weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv") 

# load smoothed warming trend for each model
model.warming.trends <- read.csv("./CMIP6/summaries/loess_smoothed_model_warming_rates.csv")

# check model names
check <- data.frame(models = models,
                    trend.names = colnames(model.warming.trends)[2:ncol(model.warming.trends)])
check # line up, just have the -/. difference

# fix
colnames(model.warming.trends)[2:ncol(model.warming.trends)] <- models

# and pivot longer 
model.warming.trends <- model.warming.trends %>%
  pivot_longer(cols = -year)


# and load predicted warming rate across all models (weighted by 1972:2021 predictions)
predicted.warming <- read.csv("./CMIP6/summaries/brms_predicted_North_Pac_warming.csv")


## record historical outcomes ----------------
# approach for each observation:
# -select window of year-7 to year+7 for each observation
# -find range of predicted warming during that window in predicted.warming
# -limit each model to the relevant warming range in model.warming.trends
# calculate proportion as big as or larger than ersst anomaly

# create df of historical outcomes
historical.rolling.window.outcomes <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # loop through each region

  for(j in 1:length(regions)){ # start j loop (regions)
  # j <- 3 # (GOA only for now!)
  
  # separate model and region of interest
  hist.temp <- cmip.anom %>% 
      filter(experiment == "hist_ssp585",
             model == models[i],
             region == regions[j])

  # break out the relevant region from ersst
  ersst.temp <- ersst.anom %>%
    filter(region == regions[j])
  
  # loop through each year of observation
  for(k in 1:nrow(ersst.temp)){ # start k loop (years)
    # k <- 5
    
    # define 15-year window
    window <- (ersst.temp$year[k] - 7) : (ersst.temp$year[k] + 7)
    
    # define range of predicted warming values for this window
    warming.range <- range(predicted.warming$pred_mean[predicted.warming$year %in% window])
    
    # find years for the model of interest that fall into this warming range
    use <- model.warming.trends %>%
      filter(name == models[i],
             value >= warming.range[1] & value <= warming.range[2])
    
    
    hist.temp.use <- hist.temp %>%
      filter(year %in% use$year)
    
    # record outcome for annual unsmoothed, annual 2-yr running mean, and annual 3-yr running mean
    
    annual.1yr <- ifelse(hist.temp.use$annual.unsmoothed >= ersst.temp$annual.anomaly.unsmoothed[k], 1, 0)
    
    annual.2yr <- ifelse(hist.temp.use$annual.two.yr.running.mean >= ersst.temp$annual.anomaly.two.yr.running.mean[k], 1, 0)
    
    annual.3yr <- ifelse(hist.temp.use$annual.three.yr.running.mean >= ersst.temp$annual.anomaly.three.yr.running.mean[k], 1, 0)
    
    # calculate prob for winter unsmoothed, winter 2-yr running mean, and winter 3-yr running mean
    
    winter.1yr <- ifelse(hist.temp.use$winter.unsmoothed >= ersst.temp$winter.anomaly.unsmoothed[k], 1, 0)
    
    winter.2yr <- ifelse(hist.temp.use$winter.two.yr.running.mean >= ersst.temp$winter.anomaly.two.yr.running.mean[k], 1, 0)
    
    winter.3yr <- ifelse(hist.temp.use$winter.three.yr.running.mean >= ersst.temp$winter.anomaly.three.yr.running.mean[k], 1, 0) 
    
    # add to df
    historical.rolling.window.outcomes <- rbind(historical.rolling.window.outcomes,
                                    data.frame(model = models[i],
                                               period = "historical",
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

# break into separate objects for each region and save

for(i in 1:length(regions)){
  
  temp <- historical.rolling.window.outcomes %>%
    filter(region == regions[i]) 
  
  write.csv(temp, file = paste("./CMIP6/summaries/", regions[i], "_historical_outcomes_rolling_window.csv", sep = ""), row.names = F)
  
}

