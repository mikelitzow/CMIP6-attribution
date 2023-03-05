# EBS SST pdf plots for SWAMC outreach talk
library(tidyverse)

# set theme
theme_set(theme_bw())

# custom colorblind colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## different approach - plot hindcast and projected pdfs for sst anomalies --------------------------

# load ersst anomalies
ersst.anom <- read.csv("./data/regional_north_pacific_ersst_anomaly_time_series.csv")

unique(ersst.anom$region)

ersst.anom <- ersst.anom %>%
  dplyr::filter(region == "Eastern_Bering_Sea")

# plot distributions to check
ggplot(filter(ersst.anom, year %in% 1950:1999), aes(annual.anomaly.three.yr.running.mean)) +
  geom_density(fill = "grey")

# sst anomaly threshold predicted to correspond with highest 5 FAR values (2017 value is threshold) = 2.075
ersst.max <- 4

# load CMIP6 anomalies
cmip.anom <- read.csv("./data/CMIP6.anomaly.time.series.csv")

# load estimated warming level timing for each model
timing <- read.csv("./data/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# create df to catch outcomes for extreme runs
anomaly.pdfs <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  pre.temp <- cmip.anom %>% 
    filter(experiment == "piControl",
           model == models[i],
           region == "Eastern_Bering_Sea")
  
  
  
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
           region == "Eastern_Bering_Sea")
  
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
model.weights <- read.csv("./Data/CMIP6_model_weights_by_region_window.csv") 

# clean up model weights 
model.weights <- model.weights %>%
  filter(window == "annual",
         region == "Eastern_Bering_Sea") %>%
  select(model, scaled.total.weight) 

# calculate EBS-specific model warming weights (based on prediction of experienced warming)

ersst <- read.csv("./Data/regional_north_pacific_ersst_time_series.csv")

ersst <- ersst %>%
  select(year, annual.unsmoothed) %>%
  mutate(model = "ersst")

models <- read.csv("./Data/CMIP6.sst.time.series.csv")

# combine models and ersst observations into "data"
data <- models %>% 
  filter(experiment == "hist_ssp585",
         region == "Eastern_Bering_Sea",
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

# get good labels for each warming period
labs <- data.frame(period = unique(resample.pdf$period),
                   plot_period = c("Preindustrial",
                                   "1950 to 0.5°",
                                   "0.5° to 1.0°",
                                   "1.0° to 1.5°",
                                   "1.5° to 2.0°"),
                   plot_order = 1:5)

resample.pdf <- left_join(resample.pdf, labs) %>%
  mutate(plot_period = reorder(plot_period, plot_order))

# get percentage of each period exceeding 4SD 
percentage_4sd <- resample.pdf %>%
  group_by(plot_period) %>%
  summarise(percentage_4sd = sum(anomaly > 4)/n())

# get additional summary statistics
summary__stats <- resample.pdf %>%
  group_by(plot_period) %>%
  summarise(percentage_4sd = sum(anomaly > 4)/n(),
            mean_anomaly = mean(anomaly),
            sd_anomaly = sd(anomaly),
            percentage_4sd = sum(anomaly > 4)/n())

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf_plot <- ggplot(resample.pdf, aes(plot_period, anomaly)) +
  geom_violin(fill = cb[6], lty = 0, alpha = 0.5) +
  coord_flip() +
  xlab("North Pacific warming") +
  ylab("SST anomaly (Std. Dev.)") +
  geom_hline(yintercept = ersst.max, lty = 2)

pdf_plot

ggsave("./Figs/SST_pdfs_by_warming_period.png", width = 3.5, height = 3.5, units = 'in')
