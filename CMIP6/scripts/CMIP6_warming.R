## evaluate North Pacific-wide warming rate for CMIP6 models

## compare historical / ssp585 with pre-industrial 

library(tidyverse)

# set theme
theme_set(theme_bw())

# load CMIP6 and ERSST sst

cmip <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv")

# select historical_ssp585 runs for North Pacific, remove un-needed columns

cmip <- cmip %>%
  filter(region == "North_Pacific",
         experiment == "hist_ssp585") %>%
  select(model, year, annual.unsmoothed)

ersst <- read.csv("./CMIP6/summaries/North_Pacific_ersst_1854-2021.csv")

# steps:

# 1) fit loess to full time series for each model and get estimated time
# for 0.5, 1, 1.5, 2 degrees of warming wrt 1850-1949

# 2) fit loess to ersst to calculate warming wrt 1854-1949 through 2021

# 3) regress annual observed on annual modeled to weight each model

# 4) fit Bayesian regression model to estimate overall warming event


# step 1
# fit loess to full time series for each model and get estimated time
# for 0.5, 1, 1.5, 2 degrees of warming wrt 1850-1949

cmip.models <- unique(cmip$model)

cmip$warming.wrt.1850.1949 <- NA

for(i in 1:length(cmip.models)){

  temp <- cmip %>%
    filter(model == cmip.models[i])
  
  climatology.mean <- mean(temp$annual.unsmoothed[temp$year %in% 1850:1949])
  
  temp.warming <- temp$annual.unsmoothed - climatology.mean
  
  cmip$warming.wrt.1850.1949[cmip$model == cmip.models[i]] <- temp.warming
  
}


ersst$warming.wrt.1854.1949 <- ersst$annual.unsmoothed - mean(ersst$annual.unsmoothed[ersst$year %in% 1854:1949])

plot.ersst <- ersst %>%
  mutate(model = "ersst")

ggplot(cmip, aes(year, warming.wrt.1850.1949, color = model)) +
  geom_line() +
  geom_line(data = plot.ersst, aes(year, warming.wrt.1854.1949, color = model), color = "black") +
  ggtitle("North Pacific warming (ERSST in black)")

ggsave("./CMIP6/figs/n_pac_model_estimated_warming_rate_by_model.png", width = 7, height = 4, units = 'in')

warming.rate <- warming.rate %>%
  pivot_wider(names_from = name, values_from = value)
warming.rate <- as.matrix(warming.rate)

# create an object to catch timing of warming
smoothed.warming <- data.frame()

models <- unique(cmip$model)

for(i in 1:length(models)){

  temp <- cmip %>%
    filter(model == models[i])
  
  mod <- loess(temp$warming.wrt.1850.1949 ~ temp$year)
  
  smoothed.warming <- rbind(smoothed.warming,
                                data.frame(model = models[i],
                                           year = 1850:2099,
                                           smoothed.warming = predict(mod)))
  
} 

# plot to check
ggplot(smoothed.warming, aes(year, smoothed.warming)) +
  geom_line() + 
  facet_wrap(~model) +
  geom_hline(yintercept = c(0.5, 1, 1.5, 2), lty = 2, color = "red") +
  coord_cartesian(xlim = c(2000,2050), ylim = c(0,4))

# looks good!

# get the year that each warming level is reached for each model
levels <- c(0.5, 1, 1.5, 2)

timing <- data.frame()


for(i in 1:length(models)){

  temp <- smoothed.warming %>%
    filter(model == models[i])
  
  temp.timing <- NA
  
  for(j in 1:length(levels)){
 
    
    temp.timing[j] <- min(temp$year[temp$smoothed.warming >= levels[j]])
    
  }
  
  timing <- rbind(timing,
                  data.frame(
                    model = models[i],
                    level = levels,
                    year = temp.timing))
}

# and plot
ggplot(timing, aes(year)) +
  geom_histogram(bins = 6) +
  facet_wrap(~level, scales = "free_x")

ggplot(timing, aes(level, year, color = model)) +
  geom_point() +
  labs(x = "N. Pacific warming wrt 1850-1949 (°C)",
       y = "Year first reached") +
  scale_y_continuous(breaks = seq(1960, 2090, by = 10))

ggsave("./CMIP6/figs/N_Pac_warming_rate_by_model.png", width = 6, height = 4, units = 'in')

timing$level <- as.factor(timing$level)

obs.timing <- data.frame()


# fit loess to ersst to get observed warming rate
mod <- loess(ersst$warming.wrt.1854.1949 ~ ersst$year)
ersst$trend <- predict(mod)

# plot to check
ggplot(ersst, aes(year, trend)) +
  geom_line()

# looks right


# now get observed timing of different warming levels - 
# we've only reached 1, which is 0.5 

i <- 1
temp <- min(ersst$year[ersst$trend >= levels[i]]) 

obs.timing <- data.frame(level = levels,
                         timing = c(2003, NA, NA, NA))

# 2003


obs.timing$level <- as.factor(obs.timing$level)

ggplot(timing, aes(y = year, level)) +
  geom_boxplot() +
  labs(x = "N Pacific warming wrt 1850-1949 (°C)",
       y = "Year first reached",
       title = "Red dot = ERSSTv5 warming") +
  geom_point(data = obs.timing, aes(y = timing, as.factor(level)), color = "red") 

ggsave("./CMIP6/figs/ne_pacific_warming_rate_models_obs.png", width = 5, height = 4, units = 'in')



# save timing
write.csv(timing, "./CMIP6/summaries/model.north.pacific.warming.timing.csv")

###########################
# STOPPED HERE            #
###########################
# evaluate ability of different models to predict warming for 1950-2021

# reload unsmoothed warming data
model.warming.rate <- read.csv("./CMIP6/summaries/ne_pacific_annual_modeled_sst.csv", row.names = 1)

obs.warming.rate <- read.csv("./CMIP6/summaries/ne_pacific_annual_observed_sst.csv", row.names = 1)

predict.timing <- model.warming.rate %>%
  filter(year %in% 1972:2021) 

response.timing <- obs.warming.rate %>%
  filter(year %in% 1972:2021) %>%
  rename(ersst = ersst.warming)

predict.timing <- left_join(predict.timing, response.timing)

# plot to check
ggplot(predict.timing, aes(value, ersst)) +
  geom_point() +
  facet_wrap(~name) # avoids a 'fishook' around declining rate of 1950s


model.warming.evaluation <- data.frame()

models <- unique(predict.timing$name)

for(i in 1:length(models)){
  
  # i <- 1
  
  temp <- predict.timing %>%
    filter(name == models[i])
  
  linear.fit <- lm(ersst ~ value, data = temp)
  
  RSS <- c(crossprod(linear.fit$residuals))
  Pearson.resid <- RSS / linear.fit$df.residual
  
  
  model.warming.evaluation <- rbind(model.warming.evaluation,
                                    data.frame(model = models[i],
                                               coeff = coefficients(linear.fit)[2],
                                               Pearson.resid = Pearson.resid))
  
}


model.warming.evaluation # should use coefficients! (inverse of difference from 1)

model.warming.evaluation$coeff.from.one <- 1-model.warming.evaluation$coeff
model.warming.evaluation$weight <- abs(1/model.warming.evaluation$coeff.from.one)

ggplot(model.warming.evaluation, aes(weight)) +
  geom_histogram(bins = 8, fill = "grey", color = "black") # seems right!


ggplot(model.warming.evaluation, aes(abs(coeff.from.one), weight)) +
  geom_point()

# save 
write.csv(model.warming.evaluation, "./CMIP6/summaries/N_Pac_warming_model_weights.csv", row.names = F)


# next, brms estimates of warming timing
