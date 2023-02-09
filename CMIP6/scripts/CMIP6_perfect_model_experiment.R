# compare each model with every other model to tune sigma-d for model weighting

library(tidyverse)

# set theme
theme_set(theme_bw())

## load model time series ---------------------------
cmip <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv") %>%
  filter(experiment == "hist_ssp585") %>%
  select(region, model, year, annual.unsmoothed) %>%
  rename(annual.sst = annual.unsmoothed)

## model evaluation -------------------------------------------------

# get vector of models to roll through
models <- unique(cmip$model)

# get regions to roll through
regions <- unique(ersst$region)


# loop through each model as "truth" and calculate differences for all other models
# create object to catch results
differences <- data.frame()


for(m in 1:length(models)){
  # m <- 2 
  # create vector of models to compare with
  compare.cmip <- models[-m]
  
  for(r in 1:length(regions)){
    # r <- 1
    true.cmip <- cmip %>%
      filter(region == regions[r],
             model == models[m]) 
  
    # calculate "true" regional climatology, SD, AR(1), and trend
    true.climatology <- mean(true.cmip$annual.sst[true.cmip$year %in% 1950:2014])
    true.sd <- sd(true.cmip$annual.sst[true.cmip$year %in% 1950:2014])
    true.ar <- ar(true.cmip$annual.sst[true.cmip$year %in% 1950:2014], order.max = 1, aic = F)$ar
    true.trend <- summary(lm(annual.sst ~ year, 
                          data = true.cmip[true.cmip$year %in% 1973:2022,]))$coefficients[2,1]
    
    for(c in 1: length(compare.cmip)){
      # c <- 1
      comp.cmip <- cmip %>%
        filter(region == regions[r],
               model == compare.cmip[c]) 
      
      # calculate regional CMIP climatology, SD, AR(1), and trend
      comp.climatology <- mean(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014])
      comp.sd <- sd(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014])
      comp.ar <- ar(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014], order.max = 1, aic = F)$ar
      comp.trend <- summary(lm(annual.sst ~ year, 
                               data = comp.cmip[comp.cmip$year %in% 1973:2022,]))$coefficients[2,1]

      # add to output
      differences <- rbind(differences,
                      data.frame(target_model = models[m],
                                 region = regions[r],
                                 prediction_model = compare.cmip[c],
                                 climatology_diff = abs(true.climatology - comp.climatology),
                                 sd_diff = abs(true.sd - comp.sd),
                                 ar_diff = abs(true.ar - comp.ar),
                                 trend_diff = abs(true.trend - comp.trend)))
      
    } # close c loop (comparison models)
  } # close m loop (models)
} # close r loop (regions) 

# now save output
write.csv(differences, "./CMIP6/summaries/CMIP6_prefect_model_differences.csv", row.names = F)

# load model similarities
# these were calculated in CMIP6_time_series_weighting.R
similarities <- read.csv("./CMIP6/summaries/CMIP6_time_series_model_similarities.csv")


## step 2: calculate model similarities ----------------------------------------

# calculate model weighting differences ----------------
# create object to catch results
similarity <- data.frame()

for(r in 1:length(regions)){
  # r <- 2 
  
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
    
    ## loop through every other model and compare
    
    # create a vector of models to compare with
    compare.cmip <- models[-m]
    
    # and loop through 
    for(c in 1:length(compare.cmip)){
      # c <- 1
      comp.cmip <- cmip %>%
        filter(region == regions[r],
               model == compare.cmip[c]) 
      
      # calculate regional CMIP climatology, SD, AR(1), and trend
      comp.climatology <- mean(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014])
      comp.sd <- sd(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014])
      comp.ar <- ar(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014], order.max = 1, aic = F)$ar
      comp.trend <- summary(lm(annual.sst ~ year, 
                               data = comp.cmip[comp.cmip$year %in% 1973:2022,]))$coefficients[2,1]
    
    
    # add to output
    similarity <- rbind(similarity,
                    data.frame(region = regions[r],
                               model = models[m],
                               comparison = compare.cmip[c],
                               climatology_diff = abs(comp.climatology - cmip.climatology),
                               sd_diff = abs(comp.sd - cmip.sd),
                               ar_diff = abs(comp.ar - cmip.ar),
                               trend_diff = abs(comp.trend - cmip.trend)))
    
    } # close c loop (compare)
  } # close m loop (models)
} # close r loop (regions) 

# now save output
write.csv(similarity, "./CMIP6/summaries/CMIP6_time_series_model_similarities.csv", row.names = F)
