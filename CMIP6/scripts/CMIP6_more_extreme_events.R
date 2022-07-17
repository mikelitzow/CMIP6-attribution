## summarize proportion of annual sst anomalies greater 
## than yet observed at different warming levels

library(tidyverse)

theme_set(theme_bw())

# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# plot distributions to check
ggplot(filter(ersst.anom, year %in% 1950:2021), aes(annual.anomaly.unsmoothed)) +
  geom_density(fill = "grey") +
  facet_wrap(~region) +
  xlim(-3, 6.5)

# find max annual sst anomaly for each time series
ersst.max <- ersst.anom %>%
  group_by(region) %>%
  summarise(max.anomaly = max(annual.anomaly.unsmoothed)) 


# and define the mean maximum observed across all six time series
mean.max <- round(mean(ersst.max$max.anomaly), 1) # 4.02, rounding to 4 SD


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

# create df to catch outcomes for extreme runs
extreme.outcomes <- data.frame()

# loop through two extreme levels - SD = 4 and SD = 5

extreme.sd <- c(4, 5)

for(x in 1:length(extreme.sd)) { # start x loop (extremes)

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
    
      
    # record how many model years are more extreme
    extreme.outcomes <- rbind(extreme.outcomes,
                                    data.frame(extreme.level = extreme.sd[x],
                                    model = models[i],
                                    region = regions[j],
                                    period = "preindustrial",
                                    count = sum(pre.temp$annual.unsmoothed >= extreme.sd[x], na.rm = T),
                                    N = length(!is.na(pre.temp$annual.unsmoothed))))
 
    
    } # close j loop (regions)
  
  } # close i loop (models)
  
} # close x loop (extreme level)


head(extreme.outcomes) 

colSums(extreme.outcomes[, c(5:6)]) # only two 4SD or 5SD occurrences in preindustrial

# look at GOA
checkGOA <- extreme.outcomes %>%
  filter(region == "Gulf_of_Alaska") %>%
  group_by(extreme.level) %>%
  summarise(count = sum(count),
            N = sum(N))

checkGOA

# find positive occurences
filter(extreme.outcomes, count > 0)

# two CCE occurrences of 4SD 

## record outcomes using different warming levels from hist.585 ----------------

# loop through two extreme levels - SD = 4 and SD = 5
for(x in 1:length(extreme.sd)) { # start x loop (extremes)

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # loop through each region

  for(j in 1:length(regions)){ # start j loop (regions)
  # j <- 1
  
  # separate model and region of interest
  hist.temp <- cmip.anom %>% 
      filter(experiment == "hist_ssp585",
             model == models[i],
             region == regions[j])
  
  # # separate this region from ersst.max
  # ersst.temp <- ersst.max %>%
  #   filter(region == regions[j])
  
  ## pull 1950 - 0.5 degrees warming
  
  use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)

  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(extreme.level = extreme.sd[x],
                                       model = models[i],
                                       region = regions[j],
                                       period = "1950_to_0.5",
                                       count = sum(hist.temp.use$annual.unsmoothed >= extreme.sd[x], na.rm = T),
                                       N = length(!is.na(hist.temp.use$annual.unsmoothed))))
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]

  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(extreme.level = extreme.sd[x],
                                       model = models[i],
                                       region = regions[j],
                                       period = "0.5_to_1.0",
                                       count = sum(hist.temp.use$annual.unsmoothed >= extreme.sd[x], na.rm = T),
                                       N = length(!is.na(hist.temp.use$annual.unsmoothed))))
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(extreme.level = extreme.sd[x],
                                       model = models[i],
                                       region = regions[j],
                                       period = "1.0_to_1.5",
                                       count = sum(hist.temp.use$annual.unsmoothed >= extreme.sd[x], na.rm = T),
                                       N = length(!is.na(hist.temp.use$annual.unsmoothed))))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(extreme.level = extreme.sd[x],
                                       model = models[i],
                                       region = regions[j],
                                       period = "1.5_to_2.0",
                                       count = sum(hist.temp.use$annual.unsmoothed >= extreme.sd[x], na.rm = T),
                                       N = length(!is.na(hist.temp.use$annual.unsmoothed))))
  

  } # close j loop (regions)

} # close i loop (models)

} # close x loop (extreme levels)

# check
check <- extreme.outcomes %>%
  group_by(extreme.level, region, period) %>%
  summarise(count = sum(count),
            N = sum(N)) %>%
  mutate(prop = count/N)

View(check)

# and save
write.csv(extreme.outcomes, "./CMIP6/summaries/more_extreme_annual_anomalies.csv", row.names = F)
