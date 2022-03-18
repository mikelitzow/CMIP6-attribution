# fit binomial model to preindustrial/historical outcomes in mgcv;
# estimate FAR from the resulting probabilities;
# bootstrap the climate model outputs to estimate uncertainty around the resulting FAR estimate

library(tidyverse)
library(mgcv)

historical <- read.csv("./CMIP6/summaries/Gulf_of_Alaska_historical_outcomes.csv")
preindustrial <- read.csv("./CMIP6/summaries/Gulf_of_Alaska_preindustrial_outcomes.csv")

head(preindustrial)

# add model weights
weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv")

weights <- weights %>%
  filter(region == "Gulf_of_Alaska",
         window == "annual") %>%
  select(model, scaled.total.weight) %>%
  rename(weight = scaled.total.weight)

preindustrial <- left_join(preindustrial, weights) %>%
  select(annual.anomaly.1yr, annual.1yr.events, weight) %>%
  rename(sst.anomaly = annual.anomaly.1yr,
         event = annual.1yr.events)

historical <- left_join(historical, weights) %>%
  select(annual.anomaly.1yr, annual.1yr.events, weight) %>%
  rename(sst.anomaly = annual.anomaly.1yr,
         event = annual.1yr.events)

# observed ersst anomalies to predict under preindustrial/historical conditions
ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv") 

new.dat <- ersst %>%
  filter(region == "Gulf_of_Alaska") %>%
  select(year, annual.anomaly.unsmoothed) %>%
  rename(sst.anomaly = annual.anomaly.unsmoothed)

# object to catch FAR estimates
FAR.est <- data.frame()


for(i in 1:1000){
  # i <- 1
  # resample CMIP6 output
  draws <- sample(1:nrow(preindustrial), 100000, replace = T)
  preind.temp <- preindustrial[draws,]
  
  draws <- sample(1:nrow(historical), 100000, replace = T)
  hist.temp <- historical[draws,]

  # model
  preind.mod <- bam(event ~ 
                  s(sst.anomaly), weights = weight, data = preind.temp)
  
  hist.mod <- bam(event ~ 
                    s(sst.anomaly), weights = weight, data = hist.temp)
  
  FAR.est <- rbind(FAR.est,
                   data.frame(iteration = i,
                              year = new.dat$year,
                              anomaly = new.dat$sst.anomaly,
                              preindustrial.probability = predict(preind.mod, newdata = new.dat, type = "response"),
                              historical.probability = predict(hist.mod, newdata = new.dat, type = "response")))
}

# model
hist.mod2 <- bam(annual.1yr.events ~ 
                   model_fac + s(annual.anomaly.1yr), weights = weight, data = historical)

MuMIn::AICc(hist.mod1, hist.mod2) # mod 1 better

summary(hist.mod1)

# and predict

ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv") 

new.dat <- ersst %>%
  filter(region == "Gulf_of_Alaska") %>%
  select(year, annual.anomaly.unsmoothed) %>%
  rename(annual.anomaly.1yr = annual.anomaly.unsmoothed) %>%
  mutate(model_fac = unique(historical$model_fac)[2])

str(summary(hist.mod1))

rownames(summary(hist.mod1)$s.table)


hist.mod1$s.table

test1 <- predict(hist.mod1, 
                newdata = new.dat, 
                type = "response", 
                exclude = c("model_fac", rownames(summary(hist.mod1)$s.table)[2:24]), 
                discrete = FALSE)

test2 <- predict(hist.mod1, 
                newdata = new.dat, 
                type = "response", 
                newdata.guaranteed = TRUE, 
                discrete = FALSE)

preind.mod <- bam(annual.1yr.events ~ 
                    s(annual.anomaly.1yr) + s(model_fac, bs = "re"), weights = weight, data = preindustrial, family = binomial())
 
preind.mod2 <- bam(annual.1yr.events ~ 
                    s(annual.anomaly.1yr), weights = weight, data = preindustrial, family = binomial())


ranef(preind.mod)

summary(preind.mod)
plot(preind.mod)


hist.mod <- bam(annual.1yr.events ~ 
                  s(annual.anomaly.1yr) + s(model_fac, bs = "re"), weights = weight, data = historical)


summary(hist.mod)
plot(hist.mod)





