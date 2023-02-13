# compare each model with every other model to tune sigma-d for model weighting

library(tidyverse)



# set theme
theme_set(theme_bw())

## load model time series ---------------------------
cmip <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv") %>%
  filter(experiment == "hist_ssp585") %>%
  select(region, model, year, annual.unsmoothed) %>%
  rename(annual.sst = annual.unsmoothed)

# load model similarities for weighting
sims <- read.csv("./CMIP6/summaries/CMIP6_time_series_model_similarities.csv") 

## model evaluation -------------------------------------------------

# get vector of models to roll through
models <- unique(cmip$model)

# get regions to roll through
regions <- unique(cmip$region)

# loop through each model as "truth" and for each region predict
# mean sst 1950:2014 and 2015:2044
# create object to catch results
perfect_model_prediction  <- data.frame()

for(m in 1:length(models)){
  # m <- 2 
  # create vector of models to compare with
  compare.cmip <- models[-m]
  
  for(r in 1:length(regions)){
    # r <- 1
    true.cmip <- cmip %>%
      filter(region == regions[r],
             model == models[m]) 
  
    # calculate "true" mean sst for 1950:2014 and 2015:2044
    # these climatologies are critical to FAR calculations (former)
    # and sst projections under warming (latter)
    true.1850.1949 <- mean(true.cmip$annual.sst[true.cmip$year %in% 1850:1949])
    true.1950.2014 <- mean(true.cmip$annual.sst[true.cmip$year %in% 1950:2014])
    true.2015.2044 <- mean(true.cmip$annual.sst[true.cmip$year %in% 2015:2044])
    
    # and calculate sd, AR(1), and trend as in the actual weighting
    true.sd <- sd(true.cmip$annual.sst[true.cmip$year %in% 1950:2014])
    true.ar <- ar(true.cmip$annual.sst[true.cmip$year %in% 1950:2014], order.max = 1, aic = F)$ar   
    true.trend <- summary(lm(annual.sst ~ year, 
                             data = true.cmip[true.cmip$year %in% 1973:2022,]))$coefficients[2,1]
    
    for(c in 1: length(compare.cmip)){
      # c <- 1
      comp.cmip <- cmip %>%
        filter(region == regions[r],
               model == compare.cmip[c]) 
      
      # calculate predicted values
      pred.1850.1949 <- mean(comp.cmip$annual.sst[true.cmip$year %in% 1850:1949])
      pred.1950.2014 <- mean(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014])
      pred.2015.2044 <- mean(comp.cmip$annual.sst[comp.cmip$year %in% 2015:2044])
      
      cmip.sd <- sd(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014])
      cmip.ar <- ar(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014], order.max = 1, aic = F)$ar
      cmip.trend <- summary(lm(annual.sst ~ year, 
                               data = comp.cmip[comp.cmip$year %in% 1973:2022,]))$coefficients[2,1]

      # add to output
      perfect_model_prediction <- rbind(perfect_model_prediction,
                      data.frame(true_model = models[m],
                                 region = regions[r],
                                 prediction_model = compare.cmip[c],
                                 climatology_diff = abs(true.1950.2014 - pred.1950.2014),
                                 sd_diff = abs(true.sd - cmip.sd),
                                 ar_diff = abs(true.ar - cmip.ar),
                                 trend_diff = abs(true.trend - cmip.trend),
                                 true.sd = true.sd,
                                 true.ar = true.ar,
                                 cmip.sd = cmip.sd,
                                 cmip.ar = cmip.ar,
                                 true.1850.1949 = true.1850.1949,
                                 true.1950.2014 = true.1950.2014,
                                 true.2015.2044 = true.2015.2044,
                                 pred.1850.1949 = pred.1850.1949,
                                 pred.1950.2014 = pred.1950.2014,
                                 pred.2015.2044 = pred.2015.2044))
      
    } # close c loop (comparison models)
  } # close m loop (models)
} # close r loop (regions) 

## loop through each model/region comparison and calculate prediction model for distance from "true" model

# set a range of sigma_s, centered around 0.675, to tune sigma_d (following Zhao et al. 2022)
sigma_s = seq(0.425, 0.925, by = 0.125)

# define a range of possible sigma_d values between 0.2 and 1.1
# (range used in Zhao et al. 2022 Earth's Future) 
sigma_d <- seq(0.2, 1.1, by = 0.005)

# vector of true models to predict
true_models <- unique(perfect_model_prediction$true_model)

# create object to summarize D and S for each true model prediction
predict_D_S <- data.frame()

# create object to collect tuning results
sigma_s_d_tuning <- data.frame()


# start by looping through sigma_s
for(sig.s in 1:length(sigm_s)){
  # sig.s <- 1
for(r in 1:length(regions)){
  # r <- 1

  for(t in 1:length(true_models)){
   # t <- 1
    
    temp_predict <- perfect_model_prediction %>%
      filter(region == regions[r],
             true_model == true_models[t])
      
    # scale each difference by the median, and average to get D
    temp_predict <- temp_predict %>%
      mutate(diff1 = climatology_diff / median(temp_predict$climatology_diff),
             diff2 = sd_diff / median(temp_predict$sd_diff),
             diff3 = ar_diff / median(temp_predict$ar_diff),
             diff4 = trend_diff / median(temp_predict$trend_diff)) 
    
    temp_predict$D = apply(temp_predict[,18:21], 1, mean)

    # simplify 
    temp_predict <- temp_predict %>%
      select(true_model, region, prediction_model, D, 
             true.1850.1949, true.1950.2014, true.2015.2044, pred.1850.1949, pred.1950.2014, pred.2015.2044, true.sd, true.ar, cmip.sd, cmip.ar)
    
    # loop through each prediction model and calculate independence weight
    
    # get vector of prediction models
    pred_mods <- unique(temp_predict$prediction_model)
    
    # create data frame to store similarities among prediction models
    all_similarities <- data.frame()
    
    for(p in 1:length(pred_mods)){
      # p <- 1
      
      temp_sim <- sims %>%
        filter(region == regions[r],
               model == pred_mods[p],
               comparison != temp_predict$true_model) %>% # do not include similarity with the "true" model we are predicting!
        mutate(sim1 = climatology_diff / median(temp_sim$climatology_diff),
              sim2 = sd_diff / median(temp_sim$sd_diff),
              sim3 = ar_diff / median(temp_sim$ar_diff),
              sim4 = trend_diff / median(temp_sim$trend_diff))
    
      temp_sim$S = apply(temp_sim[,8:11], 1, mean)
      
      # collate similarity for every prediction model combination
      all_similarities <- rbind(all_similarities,
                                data.frame(region = regions[r],
                                           true_model = true_models[t],
                                           prediction_model = temp_sim$model,
                                           sigma_s = sigma_s[sig.s],
                                           independence_skill = 1/(1+sum(exp(-(temp_sim$S^2/sigma_s[sig.s]^2)))))) %>%
        group_by(region, true_model, prediction_model) %>%
        summarise(sigma_s = sigma_s[sig.s],
                  independence_skill = mean(independence_skill))
    
      # # check calculation above
      # check <- NA
      # 
      # for(i in 1:nrow(temp_sim)){
      # check[i] <- exp(-(temp_sim$S[i]^2/sigma_s[sig.s]^2))
      # }
      # 
      # check_independence <- 1/(1+sum(check))
      # check_independence # correct!
      
    } # close p loop (prediction models)
      
    # join with temp_predict
    temp_predict <- left_join(temp_predict, all_similarities) 
    
    predict_D_S <- rbind(predict_D_S, temp_predict)
    
  } # close t loop (true models)
}  # close r loop (regions)

  
  
  
# loop through each region and true model, calculate and save weighted prediction for each target quantity

# create object for catching results
regional_prediction <- data.frame()

for(r in 1:length(regions)){
  # r <- 1
  
  for(t in 1:length(true_models)){
    # t <- 1
    
    # subset predict_D_S
    temp_dat <- predict_D_S %>%
      filter(region == regions[r],
             true_model == true_models[t])
    
    # loop through different values of sigma_D
    for(s in 1:length(sigma_d)){
      # s <- 1
      temp_dat <- temp_dat %>%
        mutate(W = exp(-(D^2/sigma_d[s]))*independence_skill,
               norm_W = W/mean(W)) # normalize by the mean value of W

      # temp_plot <- temp_dat %>%
      #   mutate(prediction_model = reorder(prediction_model, desc(norm_W)))
      # 
      # ggplot(temp_plot, aes(prediction_model, norm_W)) +
      #   geom_point() +
      #   geom_path()
      
      # calculate 10th and 90th quantile of weighted prediction
      predict_1850.1949 <- whdquantile(x = temp_dat$pred.1850.1949, weights = temp_dat$norm_W, probs = c(0.1, 0.9))      
      predict_1950.2014 <- whdquantile(x = temp_dat$pred.1950.2014, weights =  temp_dat$norm_W, probs = c(0.1, 0.9))
      predict_2015.2044 <- whdquantile(x = temp_dat$pred.2015.2044, weights =  temp_dat$norm_W, probs = c(0.1, 0.9))
      
      predict_sd <- whdquantile(x = temp_dat$cmip.sd, weights =  temp_dat$norm_W, probs = c(0.1, 0.9))
      predict_ar <- whdquantile(x = temp_dat$cmip.ar, weights =  temp_dat$norm_W, probs = c(0.1, 0.9))
      
      # save results
      regional_prediction <- rbind(regional_prediction,
                                   data.frame(region = regions[r],
                                              true_model = true_models[t],
                                              sigma_s = sigma_s[sig.s],
                                              sigma_d = sigma_d[s],
                                              variable = c("mean_1850-1949", "mean_1950-2014", "mean_2015-2044", "sd_1950-2014", "ar_1950-2014"),
                                              true = c(mean(temp_dat$true.1850.1949), mean(temp_dat$true.1950.2014), mean(temp_dat$true.2015.2044),
                                                       mean(temp_dat$true.sd), mean(temp_dat$true.ar)),
                                              predicted.10 = c(predict_1850.1949[1], predict_1950.2014[1], predict_2015.2044[1], predict_sd[1], predict_ar[1]),
                                              predicted.90 = c(predict_1850.1949[2], predict_1950.2014[2], predict_2015.2044[2], predict_sd[2], predict_ar[2])))
      
    } # close s loop (sigma_d)
    
  } # close t loop (true models)
  
} # close r loop (regions)

# add column to record if true value is within 10th - 90th quantiles

regional_prediction$between <- NA

for(i in 1:nrow(regional_prediction)){
  
regional_prediction$between[i] <- (between = if_else(between(regional_prediction$true[i],
                                                          regional_prediction$predicted.10[i], 
                                                          regional_prediction$predicted.90[i]), 1, 0))

}

# for each region, calculate the proportion of predictions in the 10th-90th prediction quantile for each value of sigma-d
# and calculate true-prediction correlation



variables <- unique(regional_prediction$variable)

# for(r in 1:length(regions)){
  # r <- 1
  
  for(v in 1:length(variables)){
  #   # v <- 1
    
    for(s in 1:length(sigma_d)){
    # s <- 1
      
    temp_dat <- regional_prediction %>%
      filter(variable == variables[v],
             sigma_d == sigma_d[s])
    
    sigma_d_tuning <- rbind(sigma_d_tuning,
                            data.frame(variable = variables[v],
                                       sigma_s = sigma_s[sig.s],
                                       sigma_d = sigma_d[s],
                                       proportion_between = sum(temp_dat$between)/
                                         length(temp_dat$between)))
    }
    }
  
} # close sig.s loop (sigma_s values)

ggplot(sigma_d_tuning, aes(sigma, proportion_between, color = variable)) +
  geom_line() +
  facet_wrap(~region)
