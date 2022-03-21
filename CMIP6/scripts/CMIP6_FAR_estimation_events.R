## summarize events across each model-region combination
## to estimate preindustrial and historical probability

library(tidyverse)

theme_set(theme_bw())

## STEP 1 ----------------------------------------------------
# calculate ERSST anomalies wrt 1950-1999 for each region
# this was done ERSST summaries.R"

## STEP 2 -----------------------
# calculate anomaly wrt 1950-1999 for each model preindustrial and hist_ssp585
# this was done in "CMIP6 processing.R"

## STEP 3 ---------------------------------------------
# Record the positive and negative outcomes for each preindustrial and historical run
# for each model-anomaly combination (>= this anomaly or not)
# using years between 1950 and warming = 1.0 as "present"

# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv")

# load CMIP6 model weights
model.weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv") 

# load estimated warming level timing for each model
timing <- read.csv("./CMIP6/summaries/model.north.pacific.warming.timing.csv")

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

## record outcomes using 1950 through 0.5 degree warming from hist.585 as "present" ----------------

# create df of historical outcomes
historical.outcomes <- data.frame()

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
    
  # pull "present" years (0.5 to 1.0 degrees warming OR 1950 - 0.5 degrees warming)
  
  # use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
  
  
  # and limit hist.temp to these years
  hist.temp <- hist.temp %>%
    filter(year %in% use)
  
  # break out the relevant region from ersst
  ersst.temp <- ersst.anom %>%
    filter(region == regions[j])
  
  
  # loop through each year of observation
  for(k in 1:nrow(ersst.temp)){ # start k loop (years)
    # k <- 5
    
    # record outcome for annual unsmoothed, annual 2-yr running mean, and annual 3-yr running mean
    
    annual.1yr <- ifelse(hist.temp$annual.unsmoothed >= ersst.temp$annual.anomaly.unsmoothed[k], 1, 0)
    
    annual.2yr <- ifelse(hist.temp$annual.two.yr.running.mean >= ersst.temp$annual.anomaly.two.yr.running.mean[k], 1, 0)
    
    annual.3yr <- ifelse(hist.temp$annual.three.yr.running.mean >= ersst.temp$annual.anomaly.three.yr.running.mean[k], 1, 0)
    
    # calculate prob for winter unsmoothed, winter 2-yr running mean, and winter 3-yr running mean
    
    winter.1yr <- ifelse(hist.temp$winter.unsmoothed >= ersst.temp$winter.anomaly.unsmoothed[k], 1, 0)
    
    winter.2yr <- ifelse(hist.temp$winter.two.yr.running.mean >= ersst.temp$winter.anomaly.two.yr.running.mean[k], 1, 0)
    
    winter.3yr <- ifelse(hist.temp$winter.three.yr.running.mean >= ersst.temp$winter.anomaly.three.yr.running.mean[k], 1, 0) 
    
    # add to df
    historical.outcomes <- rbind(historical.outcomes,
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
  
  temp <- historical.outcomes %>%
    filter(region == regions[i]) 
  
  # write.csv(temp, file = paste("./CMIP6/summaries/", regions[i], "_historical_outcomes.csv", sep = ""), row.names = F)
  write.csv(temp, file = paste("./CMIP6/summaries/", regions[i], "_historical_outcomes_1950-0.5_degrees_warming.csv", sep = ""), row.names = F)
  
}




#########################  
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
