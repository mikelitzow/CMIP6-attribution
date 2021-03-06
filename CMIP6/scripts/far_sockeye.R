## FAR-catch relationships for GOA sockeye

library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## Read in data --------------------------------------------
data <- read.csv("./CMIP6/data/GOA_sockeye_catch_far.csv")

# plot time series and experienced FAR
plot.dat <- data %>%
  select(year, annual_far_3, log_catch) %>%
  rename('Ocean entry FAR' = annual_far_3,
         'Log(catch)' = log_catch) %>%
  pivot_longer(cols = -year) %>%
  mutate(order = if_else(name == "Log(catch)", 2, 1),
         name = reorder(name, order))

ggplot(plot.dat, aes(year, value)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~name, scale = "free_y", ncol = 1)

## brms: setup ---------------------------------------------

sockeye_form <- bf(log_catch_stnd ~ 1 + s(annual_far_3) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term --------------------------------------
sockeye_model1 <- brm(sockeye_form,
                   data = data,
                   prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(sockeye_model1, file = "./CMIP6/brms_output/sockeye_far.rds")

sockeye_model1 <- readRDS("./CMIP6/brms_output/sockeye_far.rds")

check_hmc_diagnostics(sockeye_model1$fit)

neff_lowest(sockeye_model1$fit)

rhat_highest(sockeye_model1$fit)

summary(sockeye_model1)

bayes_R2(sockeye_model1)

plot(conditional_smooths(sockeye_model1), ask = FALSE)

plot(sockeye_model1, ask = FALSE)

y <- data$log_catch_stnd

yrep_sockeye_model1  <- fitted(sockeye_model1, scale = "response", summary = FALSE)

ppc_dens_overlay(y = y, yrep = yrep_sockeye_model1[sample(nrow(yrep_sockeye_model1), 25), ]) +
  ggtitle("sockeye_model1")

trace_plot(sockeye_model$fit)

## refit with cor, not cov in ar() ---------------------

sockeye_form2 <- bf(log_catch_stnd ~ 1 + s(annual_far_3) + ar(time = year, p = 1))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit --------------------------------------
sockeye_model2 <- brm(sockeye_form2,
                     data = data,
                     prior = priors,
                     cores = 4, chains = 4, iter = 3000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(sockeye_model2, file = "./CMIP6/brms_output/sockeye_far2.rds")

sockeye_model2 <- readRDS("./CMIP6/brms_output/sockeye_far2.rds")

check_hmc_diagnostics(sockeye_model2$fit)

neff_lowest(sockeye_model2$fit)

rhat_highest(sockeye_model2$fit)

summary(sockeye_model2)

bayes_R2(sockeye_model2)

plot(conditional_smooths(sockeye_model2), ask = FALSE)

plot(sockeye_model2, ask = FALSE)

y <- data$log_catch_stnd

yrep_sockeye_model2  <- fitted(sockeye_model2, scale = "response", summary = FALSE)

ppc_dens_overlay(y = y, yrep = yrep_sockeye_model2[sample(nrow(yrep_sockeye_model2), 25), ]) +
  ggtitle("sockeye_model2") # this looks better than model1

trace_plot(sockeye_model2$fit)

## plot model2 ------------------------

## 95% CI
ce1s_1 <- conditional_effects(sockeye_model2, effect = "annual_far_3", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(sockeye_model2, effect = "annual_far_3", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(sockeye_model2, effect = "annual_far_3", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$annual_far_3

#########################

dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$annual_far_3[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$annual_far_3[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$annual_far_3[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$annual_far_3[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(data$annual_far_3), amount = 0.0051),
                          rep(NA, 100-length(unique(data$annual_far_3))))


g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Fraction of Attributable Risk", y = "log(catch)") +
  theme_bw() +
  geom_point(data = data, aes(x = annual_far_3, y = log_catch_stnd), color = "grey40") +
  geom_hline(yintercept = 0, lty = 2) +
  geom_rug(aes(x=rug.anom, y=NULL)) 

print(g2)

## fit third model with categorical FAR -----------------------------


# add categorical variable - is FAR above or below 0.95?
data <- data %>%
  mutate(far_fac = as.factor(if_else(annual_far_3 > 0.95, "above", "below")))

sockeye_form3 <- bf(log_catch_stnd ~ 1 + far_fac + ar(time = year, p = 1))


get_prior(sockeye_form3, data = data)

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit
sockeye_model3 <- brm(sockeye_form3,
                      data = data,
                      prior = priors,
                      cores = 4, chains = 4, iter = 3000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(sockeye_model3, file = "./CMIP6/brms_output/sockeye_far3.rds")

sockeye_model3 <- readRDS("./CMIP6/brms_output/sockeye_far3.rds")

check_hmc_diagnostics(sockeye_model3$fit)

neff_lowest(sockeye_model3$fit)

rhat_highest(sockeye_model3$fit)

summary(sockeye_model3)

bayes_R2(sockeye_model3)

plot(sockeye_model3, ask = FALSE)

y <- data$log_catch_stnd

yrep_sockeye_model3  <- fitted(sockeye_model3, scale = "response", summary = FALSE)

ppc_dens_overlay(y = y, yrep = yrep_sockeye_model3[sample(nrow(yrep_sockeye_model3), 25), ]) +
  ggtitle("sockeye_model3") # this looks better than model1

trace_plot(sockeye_model3$fit)


## 95% CI
ce1s_1 <- conditional_effects(sockeye_model3, effect = "far_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot <- ce1s_1$far_fac %>%
  select(far_fac, estimate__, lower__, upper__)

plot$far_fac <- reorder(plot$far_fac, desc(plot$far_fac))

g3 <- ggplot(plot, aes(far_fac, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Log (catch) anomaly") +
  xlab("FAR") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_x_discrete(labels=c(expression("<0.95"), expression("">=0.95))) +
  theme_bw()

print(g3)

ggsave("./CMIP6/figs/FAR_sockeye_categorical.png", width=1.5, height=2, units='in')

## expected return time ----------------------------------------

# load ersst anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

unique(ersst.anom$region)

ersst.anom <- ersst.anom %>%
  dplyr::filter(region == "Gulf_of_Alaska")

# plot distributions to check
ggplot(filter(ersst.anom, year %in% 1950:1999), aes(annual.anomaly.three.yr.running.mean)) +
  geom_density(fill = "grey")

# sst anomaly threshold predicted to correspond with highest 5 FAR values (2017 value is threshold) = 2.075
ersst.max <- 2.075

# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv")

# load estimated warming level timing for each model
timing <- read.csv("./CMIP6/summaries/model.north.pacific.warming.timing.csv")

# get vector of model names
models <- unique(cmip.anom$model)

# create df to catch outcomes for extreme runs
extreme.outcomes <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  pre.temp <- cmip.anom %>% 
    filter(experiment == "piControl",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "preindustrial",
                                       count = sum(pre.temp$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(pre.temp$annual.three.yr.running.mean))))
  
  
  
  
} # close i loop (models)


head(extreme.outcomes) 

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
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "1950_to_0.5",
                                       count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))) )
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "0.5_to_1.0",
                                       count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))))
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "1.0_to_1.5",
                                       count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "1.5_to_2.0",
                                       count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))))
  
  
  
} # close i loop (models)


# check
check <- extreme.outcomes %>%
  group_by(period) %>%
  summarise(count = sum(count),
            N = sum(N)) %>%
  mutate(prop = count/N)

View(check)

## fit brms model to estimate probabilities-----------

# model weights for extreme events in different periods - 
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

extremes <- left_join(extreme.outcomes, weights) %>%
  mutate(model_fac = as.factor(model))

## brms: setup ---------------------------------------------

form <-  bf(count | trials(N) + weights(total_weight, scale = TRUE) ~
              period + (1 | model_fac))

# fit model

extremes_brms <- brm(form,
                     data = extremes,
                     family = binomial(link = "logit"),
                     seed = 1234,
                     cores = 4, chains = 4, iter = 12000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.9, max_treedepth = 14))

saveRDS(extremes_brms, "./CMIP6/brms_output/extremes_binomial.rds")


# evaluate 

model <- readRDS("./CMIP6/brms_output/extremes_binomial.rds")

check_hmc_diagnostics(model$fit)
neff_lowest(model$fit) 
rhat_highest(model$fit)
summary(model)
bayes_R2(model) 
trace_plot(model$fit)

# and plot
new.dat <- data.frame(period = unique(extremes$period),
                      model = NA,
                      N = 1000) 

plot.dat <- data.frame()

probs <- posterior_epred(model, newdata = new.dat, re_formula = NA)/1000 # dive by N to get probability

plot.dat <- rbind(plot.dat,
                  data.frame(period = new.dat$period,
                             prob = apply(probs, 2, median),
                             lower = apply(probs, 2, quantile, probs = 0.025),
                             upper = apply(probs, 2, quantile, probs = 0.975)))


# calculate inverse to get expected return time
plot.dat[,c(2:4)] <- 1/plot.dat[,c(2:4)]

# and change values above 10^4 to 10^4

change <- plot.dat[,c(2:4)] > 10^4

plot.dat[,c(2:4)][change] <- 10^4

# set periods in order
period.order <- data.frame(period = unique(plot.dat$period),
                           period.order = 1:5)

plot.dat <- left_join(plot.dat, period.order) %>%
  mutate(period = reorder(period, period.order))


ggplot(plot.dat, aes(period, prob)) +
  geom_errorbar(aes(x = period, ymin = lower, ymax = upper), width = 0.3) +
  geom_point(color = "red", size = 4) +
  scale_y_continuous(breaks=c( 1,2,5,10,20,50,100,200,500,1000,2000,5000),
                     minor_breaks = c(2:9, 
                                      seq(20, 90, by = 10),
                                      seq(200, 900, by = 100),
                                      seq(2000, 9000, by = 1000))) +
  coord_trans(y = "pseudo_log") +
  ylab("Expected return time (years)") + 
  xlab("North Pacific warming") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

ggsave("./CMIP6/figs/sockeye_extreme_return_time.png", width = 3, height = 6, units = 'in')        

## different approach - plot hindcast and projected pdfs for sst anomalies --------------------------
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

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(resample.pdf, aes(period, anomaly)) +
  geom_violin(fill = cb[6], lty = 0, alpha = 0.5) +
  coord_flip() +
  xlab("North Pacific warming") +
  ylab("Anomaly (Std. Dev.)") +
  geom_hline(yintercept = ersst.max, lty = 2) 


ggsave("./CMIP6/figs/sockeye_anomaly_pdfs.png", width = 3, height = 3, units = 'in')        

