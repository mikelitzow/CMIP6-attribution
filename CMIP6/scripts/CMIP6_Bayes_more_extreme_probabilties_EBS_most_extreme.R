## estimate probability of more extreme annual sst anomalies for EBS - 
## greater than yet observed (5.07SD)

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

cb <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv") %>%
  filter(region == "Eastern_Bering_Sea")

# find max annual sst anomaly for each time series
ersst.max <- max(ersst.anom$annual.anomaly.unsmoothed) # 5.08


# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv") %>%
  filter(region == "Eastern_Bering_Sea")

# load CMIP6 model weights
model.weights <- read.csv("./CMIP6/summaries/normalized_CMIP6_weights.csv")  %>%
  filter(region == "Eastern_Bering_Sea")

# load estimated warming level timing for each model
timing <- read.csv("./CMIP6/summaries/CMIP6_brms_warming_timing_ssp585.csv")

# and round decimal years for timing to integer
timing <- timing %>%
  mutate(year = round(year))

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
               model == models[i])
      
      
      # record how many model years are more extreme
      extreme.outcomes <- rbind(extreme.outcomes,
                                data.frame(model = models[i],
                                           period = "preindustrial",
                                           count = sum(pre.temp$annual.unsmoothed >= 5.08, na.rm = T),
                                           N = length(!is.na(pre.temp$annual.unsmoothed))))
      
    
  } # close i loop (models)
  

## record outcomes using different warming levels from hist.585 ----------------

  # loop through each model
  for(i in 1:length(models)){ # start i loop (models)
    # i <- 1

      
      # separate model and region of interest
      hist.temp <- cmip.anom %>% 
        filter(experiment == "hist_ssp585",
               model == models[i])

      
      ## pull 1950 - 0.5 degrees warming
      
      use = 1950:timing$year[timing$model == models[i] & timing$warming == 0.5]
      
      # and limit hist.temp to these years
      hist.temp.use <- hist.temp %>%
        filter(year %in% use)
      
      # record how many model years are more extreme
      extreme.outcomes <- rbind(extreme.outcomes,
                                data.frame(model = models[i],
                                           period = "1950_to_0.5",
                                           count = sum(hist.temp.use$annual.unsmoothed >= 5.08, na.rm = T),
                                           N = length(!is.na(hist.temp.use$annual.unsmoothed))))
      
      ## pull 0.5 - 1.0 degrees warming
      
      use = timing$year[timing$model == models[i] & timing$warming == 0.5]:timing$year[timing$model == models[i] & timing$warming == 1.0]
      
      # and limit hist.temp to these years
      hist.temp.use <- hist.temp %>%
        filter(year %in% use)
      
      # record how many model years are more extreme
      extreme.outcomes <- rbind(extreme.outcomes,
                                data.frame(model = models[i],
                                            period = "0.5_to_1.0",
                                           count = sum(hist.temp.use$annual.unsmoothed >= 5.08, na.rm = T),
                                           N = length(!is.na(hist.temp.use$annual.unsmoothed))))
      
      ## pull 1.0 - 1.5 degrees warming
      
      use = timing$year[timing$model == models[i] & timing$warming == 1.0]:timing$year[timing$model == models[i] & timing$warming == 1.5]
      
      # and limit hist.temp to these years
      hist.temp.use <- hist.temp %>%
        filter(year %in% use)
      
      # record how many model years are more extreme
      extreme.outcomes <- rbind(extreme.outcomes,
                                data.frame(model = models[i],
                                           period = "1.0_to_1.5",
                                           count = sum(hist.temp.use$annual.unsmoothed >= 5.08, na.rm = T),
                                           N = length(!is.na(hist.temp.use$annual.unsmoothed))))
      
      
      ## pull 1.5 - 2.0 degrees warming
      
      use = timing$year[timing$model == models[i] & timing$warming == 1.5]:timing$year[timing$model == models[i] & timing$warming == 2.0]
      
      # and limit hist.temp to these years
      hist.temp.use <- hist.temp %>%
        filter(year %in% use)
      
      # record how many model years are more extreme
      extreme.outcomes <- rbind(extreme.outcomes,
                                data.frame(model = models[i],
                                           period = "1.5_to_2.0",
                                           count = sum(hist.temp.use$annual.unsmoothed >= 5.08, na.rm = T),
                                           N = length(!is.na(hist.temp.use$annual.unsmoothed))))
      
  
  } # close i loop (models)
  


# check
check <- extreme.outcomes %>%
  group_by(period) %>%
  summarise(count = sum(count),
            N = sum(N)) %>%
  mutate(prop = count/N)

View(check)



# model weights for extreme events in different periods - 
# product of regional weighting (based on ar(1), correlation, bias) and 
# prediction of observed N. Pac. weighting

# load regional model weights 
regional_weights <- read.csv("./CMIP6/summaries/normalized_CMIP6_weights.csv") %>%
  filter(region == "Eastern_Bering_Sea")

# calculate region-specific model warming weights (based on prediction of experienced warming)
# save these values for future analysis
ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv") %>%
  filter(region == "Eastern_Bering_Sea")

ersst <- ersst %>%
  select(year, annual.unsmoothed) %>%
  mutate(model = "ersst")

models <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv") %>%
  filter(region == "Eastern_Bering_Sea")


data <- models %>% 
  filter(experiment == "hist_ssp585") %>%
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



extreme.outcomes <- left_join(extreme.outcomes, regional_weights) %>%
  mutate(model_fac = as.factor(model))


## brms: setup ---------------------------------------------

form <-  bf(count | trials(N) + weights(normalized_weight, scale = TRUE) ~
              period + (1 | model_fac))


    
    extremes_brms_EBS <- brm(form,
                         data = extreme.outcomes,
                         family = binomial(link = "logit"),
                         seed = 1299,
                         cores = 4, chains = 4, iter = 20000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.9, max_treedepth = 20))
    
    saveRDS(extremes_brms_EBS, "./CMIP6/brms_outputEBS_5.07_SD_extremes_binomial.rds")
    
  
    check_hmc_diagnostics(extremes_brms_EBS$fit) # GOA SD4, BCC SD4  increase tree_depth
    neff_lowest(extremes_brms_EBS$fit) # GOA SD only ~100-200;   BC SD4 only 300-400
    rhat_highest(extremes_brms_EBS$fit)
    summary(extremes_brms_EBS)

# summarize
    new.dat <- data.frame(period = unique(extreme.outcomes$period),
                          model = NA,
                          N = 1000) 

      
      probs <- posterior_epred(model, newdata = new.dat, re_formula = NA)/1000 # dive by N to get probability
      
      plot.dat <- data.frame(
                                   period = new.dat$period,
                                   prob = apply(probs, 2, median),
                                   lower = apply(probs, 2, quantile, probs = 0.025),
                                   upper = apply(probs, 2, quantile, probs = 0.975))
    
      # calculate inverse to get expected return time
      plot.dat[,c(2:4)] <- 1/plot.dat[,c(2:4)]

      plot.dat

      change <- plot.dat[,c(2:4)] > 10^4
      
      plot.dat[,c(2:4)][change] <- 10^4
      
      period.order <- data.frame(period = unique(plot.dat$period),
                                 period.order = 1:5)
      
      plot.dat <- left_join(plot.dat, period.order) %>%
        mutate(period = reorder(period, period.order))
      
      pos_dodge = position_dodge(width = 0.2)
      
      ggplot(plot.dat, aes(period, prob)) +
        geom_errorbar(aes(x = period, ymin = lower, ymax = upper), width = 0.3, position = pos_dodge) +
        geom_point(size = 4, position = pos_dodge, color = "red") +
        scale_y_continuous(breaks=c(1,10,100,1000,10000),
                           labels = c(expression(10^0), 
                                      expression(10^1),
                                      expression(10^2),
                                      expression(10^3),
                                      expression(">"~10^4)),
                           # labels = c("1", "10", "100", "1000", ">10,000"),
                           minor_breaks = c(2:9, 
                                            seq(20, 90, by = 10),
                                            seq(200, 900, by = 100),
                                            seq(2000, 9000, by = 1000))) +
        scale_x_discrete(labels = c("Preindustrial", "1950 to 0.5°", "0.5° to 1.0°", "1.0° to 1.5°", "1.5° to 2.0°")) +
        coord_trans(y = "pseudo_log") +
        ylab("Expected return time (years)") + 
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45,
                                         hjust = 1))

      ggsave("./CMIP6/figs/EBS_5.07SD_extreme_return_time.png", width = 6, height = 8, units = 'in')             
      