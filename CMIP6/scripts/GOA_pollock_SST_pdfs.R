# plot sst pdfs relative to critical temps for GOA pollock

library(tidyverse)

# load ersst anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

unique(ersst.anom$region)

ersst.anom <- ersst.anom %>%
  dplyr::filter(region == "Gulf_of_Alaska")

# plot distributions to check
ggplot(filter(ersst.anom, year %in% 1950:1999), aes(annual.anomaly.three.yr.running.mean)) +
  geom_density(fill = "grey")

# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv")

# load estimated warming level timing for each model
timing <- read.csv("./CMIP6/summaries/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# create df to catch outcomes for extreme runs
extreme.outcomes <- data.frame()

# plot hindcast and projected pdfs for sst anomalies --------------------------
# create df to catch outcomes for extreme runs
anomaly.pdfs <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  pre.temp <- cmip.anom %>% 
    filter(experiment == "piControl",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "preindustrial",
                                   anomaly = na.omit(pre.temp$annual.three.yr.running.mean)))
  
  
  
} # close i loop (models)


## record outcomes using different warming levels from hist.585 ----------------


# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  hist.temp <- cmip.anom %>% 
    filter(experiment == "hist_ssp585",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  # # separate this region from ersst.max
  # ersst.temp <- ersst.max %>%
  #   filter(region == regions[j])
  
  ## pull 1950 - 0.5 degrees warming
  
  use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1950_to_0.5",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "0.5_to_1.0",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.0_to_1.5",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.5_to_2.0",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  
} # close i loop (models)


# model weights for anomalies in different periods - 
# product of regional weighting (based on ar(1), correlation, bias) and 
# prediction of observed N. Pac. weighting

# load CMIP6 model weights
model.weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv") 

# clean up model weights 
model.weights <- model.weights %>%
  filter(window == "annual",
         region == "Gulf_of_Alaska") %>%
  select(model, scaled.total.weight) 

# calculate GOA-specific model warming weights (based on prediction of experienced warming)

ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv")

ersst <- ersst %>%
  select(year, annual.unsmoothed) %>%
  mutate(model = "ersst")

models <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv")

# combine models and ersst observations into "data"
data <- models %>% 
  filter(experiment == "hist_ssp585",
         region == "Gulf_of_Alaska",
         year %in% 1850:2021) %>% # note that for regional warming we will calculate anomalies wrt 1950-1999 (beginning of trustworthy ERSST)
  select(year, annual.unsmoothed, model)

data <- rbind(data, ersst) 

# calculate 1850:1949 climatology for each model and ersst
climatology <- data %>%
  filter(year %in% 1850:1949) %>%
  group_by(model) %>%
  summarize(climatology.mean = mean(annual.unsmoothed), climatology.sd = sd(annual.unsmoothed))

# combine climatology and data, calculate anomalies
data <- left_join(data, climatology) %>%
  mutate(anomaly = (annual.unsmoothed - climatology.mean) / climatology.sd)

# and pivot longer (ersst vs models)
ersst <- data %>%
  filter(model == "ersst") %>%
  select(year, anomaly) %>%
  rename(ersst.anomaly = anomaly)

data <- data %>%
  filter(model != "ersst") %>%
  left_join(., ersst)

# loop through and fit linear ersst - model regressions to get weights
regional_warming_weights <- data.frame()

models <- unique(data$model)


for(m in 1:length(models)){ # loop through models
  # m <- 1
  
  temp.dat <- data %>%
    filter(model == models[m],
           year %in% 1972:2021)
  
  
  mod <- lm(ersst.anomaly ~ anomaly, data = temp.dat)
  
  regional_warming_weights <- rbind(regional_warming_weights,
                                    data.frame(model = models[m],
                                               regional_warming_weight = 1 / abs(1-coefficients(mod)[2]))) # inverse of difference from 1!
}




weights <- left_join(model.weights, regional_warming_weights) %>%
  mutate(total_weight = scaled.total.weight * regional_warming_weight)


# plot to examine
ggplot(weights, aes(scaled.total.weight, regional_warming_weight)) +
  geom_point() 

ggplot(weights, aes(total_weight)) +
  geom_histogram(fill = "grey", color = "black", bins = 20) 

anomaly.pdfs <- left_join(anomaly.pdfs, weights) 


# resample to weight models
resample.pdf <- data.frame()

periods <- unique(anomaly.pdfs$period)

for(i in 1:length(periods)){
  # i <- 1
  
  temp <- anomaly.pdfs[anomaly.pdfs$period == periods[i],]
  
  resample.pdf <- rbind(resample.pdf,
                        data.frame(period = periods[i],
                                   anomaly = sample(temp$anomaly, 1000, replace = T, prob = temp$total_weight)))
  
  
  
}

# reorder
plot.order <- data.frame(period = unique(resample.pdf$period),
                         order = 1:5)



resample.pdf <- left_join(resample.pdf, plot.order) %>%
  mutate(period =  reorder(period, order))

# and plot
theme_set(theme_bw())
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

shann_pdf <- ggplot(resample.pdf, aes(period, anomaly)) +
  geom_violin(fill = cb[6], lty = 0, alpha = 0.5) +
  coord_flip() +
  xlab("North Pacific warming") +
  ylab("Anomaly (Std. Dev.)") +
  geom_hline(yintercept = 1.915, lty = 5) +
  ggtitle("Shannon age diversity - 3yr mean annual SST")


###########
### annual sst for weights -------------------------


anomaly.pdfs <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  pre.temp <- cmip.anom %>% 
    filter(experiment == "piControl",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "preindustrial",
                                   anomaly = na.omit(pre.temp$annual.unsmoothed)))
  
  
  
} # close i loop (models)


## record outcomes using different warming levels from hist.585 ----------------


# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  hist.temp <- cmip.anom %>% 
    filter(experiment == "hist_ssp585",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  # # separate this region from ersst.max
  # ersst.temp <- ersst.max %>%
  #   filter(region == regions[j])
  
  ## pull 1950 - 0.5 degrees warming
  
  use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1950_to_0.5",
                                   anomaly = na.omit(hist.temp.use$annual.unsmoothed)))
  
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "0.5_to_1.0",
                                   anomaly = na.omit(hist.temp.use$annual.unsmoothed)))
  
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.0_to_1.5",
                                   anomaly = na.omit(hist.temp.use$annual.unsmoothed)))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.5_to_2.0",
                                   anomaly = na.omit(hist.temp.use$annual.unsmoothed)))
  
  
  
} # close i loop (models)


# model weights for anomalies in different periods - 
# product of regional weighting (based on ar(1), correlation, bias) and 
# prediction of observed N. Pac. weighting

# load CMIP6 model weights
model.weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv") 

# clean up model weights 
model.weights <- model.weights %>%
  filter(window == "annual",
         region == "Gulf_of_Alaska") %>%
  select(model, scaled.total.weight) 

# calculate GOA-specific model warming weights (based on prediction of experienced warming)

ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv")

ersst <- ersst %>%
  select(year, annual.unsmoothed) %>%
  mutate(model = "ersst")

models <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv")

# combine models and ersst observations into "data"
data <- models %>% 
  filter(experiment == "hist_ssp585",
         region == "Gulf_of_Alaska",
         year %in% 1850:2021) %>% # note that for regional warming we will calculate anomalies wrt 1950-1999 (beginning of trustworthy ERSST)
  select(year, annual.unsmoothed, model)

data <- rbind(data, ersst) 

# calculate 1850:1949 climatology for each model and ersst
climatology <- data %>%
  filter(year %in% 1850:1949) %>%
  group_by(model) %>%
  summarize(climatology.mean = mean(annual.unsmoothed), climatology.sd = sd(annual.unsmoothed))

# combine climatology and data, calculate anomalies
data <- left_join(data, climatology) %>%
  mutate(anomaly = (annual.unsmoothed - climatology.mean) / climatology.sd)

# and pivot longer (ersst vs models)
ersst <- data %>%
  filter(model == "ersst") %>%
  select(year, anomaly) %>%
  rename(ersst.anomaly = anomaly)

data <- data %>%
  filter(model != "ersst") %>%
  left_join(., ersst)

# loop through and fit linear ersst - model regressions to get weights
regional_warming_weights <- data.frame()

models <- unique(data$model)


for(m in 1:length(models)){ # loop through models
  # m <- 1
  
  temp.dat <- data %>%
    filter(model == models[m],
           year %in% 1972:2021)
  
  
  mod <- lm(ersst.anomaly ~ anomaly, data = temp.dat)
  
  regional_warming_weights <- rbind(regional_warming_weights,
                                    data.frame(model = models[m],
                                               regional_warming_weight = 1 / abs(1-coefficients(mod)[2]))) # inverse of difference from 1!
}




weights <- left_join(model.weights, regional_warming_weights) %>%
  mutate(total_weight = scaled.total.weight * regional_warming_weight)


# plot to examine
ggplot(weights, aes(scaled.total.weight, regional_warming_weight)) +
  geom_point() 

ggplot(weights, aes(total_weight)) +
  geom_histogram(fill = "grey", color = "black", bins = 20) 

anomaly.pdfs <- left_join(anomaly.pdfs, weights) 


# resample to weight models
resample.pdf <- data.frame()

periods <- unique(anomaly.pdfs$period)

for(i in 1:length(periods)){
  # i <- 1
  
  temp <- anomaly.pdfs[anomaly.pdfs$period == periods[i],]
  
  resample.pdf <- rbind(resample.pdf,
                        data.frame(period = periods[i],
                                   anomaly = sample(temp$anomaly, 1000, replace = T, prob = temp$total_weight)))
  
  
  
}

# reorder
plot.order <- data.frame(period = unique(resample.pdf$period),
                         order = 1:5)



resample.pdf <- left_join(resample.pdf, plot.order) %>%
  mutate(period =  reorder(period, order))

# and plot
theme_set(theme_bw())
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

weight_pdf <- ggplot(resample.pdf, aes(period, anomaly)) +
  geom_violin(fill = cb[6], lty = 0, alpha = 0.5) +
  coord_flip() +
  xlab("North Pacific warming") +
  ylab("Anomaly (Std. Dev.)") +
  geom_hline(yintercept = 2.602, lty = 5) +
  ggtitle("Spawning weight - annual SST")

##########
### 2yr winter sst for recruitment -------------------------


anomaly.pdfs <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  pre.temp <- cmip.anom %>% 
    filter(experiment == "piControl",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "preindustrial",
                                   anomaly = na.omit(pre.temp$winter.two.yr.running.mean)))
  
  
  
} # close i loop (models)


## record outcomes using different warming levels from hist.585 ----------------


# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  hist.temp <- cmip.anom %>% 
    filter(experiment == "hist_ssp585",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  # # separate this region from ersst.max
  # ersst.temp <- ersst.max %>%
  #   filter(region == regions[j])
  
  ## pull 1950 - 0.5 degrees warming
  
  use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1950_to_0.5",
                                   anomaly = na.omit(hist.temp.use$winter.two.yr.running.mean)))
  
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "0.5_to_1.0",
                                   anomaly = na.omit(hist.temp.use$winter.two.yr.running.mean)))
  
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.0_to_1.5",
                                   anomaly = na.omit(hist.temp.use$winter.two.yr.running.mean)))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.5_to_2.0",
                                   anomaly = na.omit(hist.temp.use$winter.two.yr.running.mean)))
  
  
  
} # close i loop (models)


# model weights for anomalies in different periods - 
# product of regional weighting (based on ar(1), correlation, bias) and 
# prediction of observed N. Pac. weighting

# load CMIP6 model weights
model.weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv") 

# clean up model weights 
model.weights <- model.weights %>%
  filter(window == "annual",
         region == "Gulf_of_Alaska") %>%
  select(model, scaled.total.weight) 

# calculate GOA-specific model warming weights (based on prediction of experienced warming)

ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv")

ersst <- ersst %>%
  select(year, annual.unsmoothed) %>%
  mutate(model = "ersst")

models <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv")

# combine models and ersst observations into "data"
data <- models %>% 
  filter(experiment == "hist_ssp585",
         region == "Gulf_of_Alaska",
         year %in% 1850:2021) %>% # note that for regional warming we will calculate anomalies wrt 1950-1999 (beginning of trustworthy ERSST)
  select(year, annual.unsmoothed, model)

data <- rbind(data, ersst) 

# calculate 1850:1949 climatology for each model and ersst
climatology <- data %>%
  filter(year %in% 1850:1949) %>%
  group_by(model) %>%
  summarize(climatology.mean = mean(annual.unsmoothed), climatology.sd = sd(annual.unsmoothed))

# combine climatology and data, calculate anomalies
data <- left_join(data, climatology) %>%
  mutate(anomaly = (annual.unsmoothed - climatology.mean) / climatology.sd)

# and pivot longer (ersst vs models)
ersst <- data %>%
  filter(model == "ersst") %>%
  select(year, anomaly) %>%
  rename(ersst.anomaly = anomaly)

data <- data %>%
  filter(model != "ersst") %>%
  left_join(., ersst)

# loop through and fit linear ersst - model regressions to get weights
regional_warming_weights <- data.frame()

models <- unique(data$model)


for(m in 1:length(models)){ # loop through models
  # m <- 1
  
  temp.dat <- data %>%
    filter(model == models[m],
           year %in% 1972:2021)
  
  
  mod <- lm(ersst.anomaly ~ anomaly, data = temp.dat)
  
  regional_warming_weights <- rbind(regional_warming_weights,
                                    data.frame(model = models[m],
                                               regional_warming_weight = 1 / abs(1-coefficients(mod)[2]))) # inverse of difference from 1!
}




weights <- left_join(model.weights, regional_warming_weights) %>%
  mutate(total_weight = scaled.total.weight * regional_warming_weight)


# plot to examine
ggplot(weights, aes(scaled.total.weight, regional_warming_weight)) +
  geom_point() 

ggplot(weights, aes(total_weight)) +
  geom_histogram(fill = "grey", color = "black", bins = 20) 

anomaly.pdfs <- left_join(anomaly.pdfs, weights) 


# resample to weight models
resample.pdf <- data.frame()

periods <- unique(anomaly.pdfs$period)

for(i in 1:length(periods)){
  # i <- 1
  
  temp <- anomaly.pdfs[anomaly.pdfs$period == periods[i],]
  
  resample.pdf <- rbind(resample.pdf,
                        data.frame(period = periods[i],
                                   anomaly = sample(temp$anomaly, 1000, replace = T, prob = temp$total_weight)))
  
  
  
}

# reorder
plot.order <- data.frame(period = unique(resample.pdf$period),
                         order = 1:5)



resample.pdf <- left_join(resample.pdf, plot.order) %>%
  mutate(period =  reorder(period, order))

# and plot
theme_set(theme_bw())
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

recr_pdf <- ggplot(resample.pdf, aes(period, anomaly)) +
  geom_violin(fill = cb[6], lty = 0, alpha = 0.5) +
  coord_flip() +
  xlab("North Pacific warming") +
  ylab("Anomaly (Std. Dev.)") +
  geom_hline(yintercept = 1.655, lty = 5) +
  ggtitle("Recruitment - 2yr mean winter SST")

png("./CMIP6/figs/pollock_sst_pdfs.png", width = 5, height = 10, units = 'in', res = 300)

ggpubr::ggarrange(recr_pdf, shann_pdf, weight_pdf, ncol = 1, labels = "auto")

dev.off()
