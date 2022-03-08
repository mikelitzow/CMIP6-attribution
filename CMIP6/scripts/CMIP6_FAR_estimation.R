## Estimate FAR across each model-region combination

library(tidyverse)

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