# trying to predict from gam while ignoring model_fac

library(tidyverse)
library(mgcv)

# load binomial climate model outputs
historical <- read.csv("./CMIP6/summaries/Gulf_of_Alaska_historical_outcomes.csv")

# add climate model weights
weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv")

weights <- weights %>%
  filter(region == "Gulf_of_Alaska",
         window == "annual") %>%
  select(model, scaled.total.weight) %>%
  rename(weight = scaled.total.weight)

historical <- left_join(historical, weights) %>%
  mutate(model_fac = as.factor(model))

# model
hist.mod <- bam(annual.1yr.events ~ 
                  model_fac + s(annual.anomaly.1yr) + s(annual.anomaly.1yr, by = model_fac, bs = "re"), weights = weight, data = historical)

summary(hist.mod1)

# and predict probabilities from this model for ersst anomaly time series
# we want to predict from the model while remaining agnostic as to model differences

ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv") 

new.dat1 <- ersst %>%
  filter(region == "Gulf_of_Alaska") %>%
  select(year, annual.anomaly.unsmoothed) %>%
  rename(annual.anomaly.1yr = annual.anomaly.unsmoothed) %>%
  mutate(model_fac = unique(historical$model_fac)[1]) # we have to specify model_fac even if excluding from the predict statement!

# try with a different model_fac level to confirm that exclude is working in predict
new.dat2 <- ersst %>%
  filter(region == "Gulf_of_Alaska") %>%
  select(year, annual.anomaly.unsmoothed) %>%
  rename(annual.anomaly.1yr = annual.anomaly.unsmoothed) %>%
  mutate(model_fac = unique(historical$model_fac)[2])

test1 <- predict(hist.mod1, 
                newdata = new.dat1, 
                type = "response", 
                exclude = c("model_fac", rownames(summary(hist.mod1)$s.table)[2:24]), 
                discrete = FALSE)

test2 <- predict(hist.mod1, 
                 newdata = new.dat2,
                 type = "response", 
                 exclude = c("model_fac", rownames(summary(hist.mod1)$s.table)[2:24]), 
                 discrete = FALSE)


plot(test1, test2)
abline(a = 0, b = 1, col = "red")

# model_fac effects haven't been excluded
