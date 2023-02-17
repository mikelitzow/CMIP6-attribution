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
# 
# ## model evaluation -------------------------------------------------
# 
# # get vector of models to roll through
# models <- unique(cmip$model)
# 
# # get regions to roll through
# regions <- unique(cmip$region)
# 
# # loop through each model as "truth" and for each region predict
# # climatology, trend, sd, and ar(1)
# # create object to catch results
# perfect_model_prediction  <- data.frame()
# 
# for(m in 1:length(models)){
#   # m <- 2 
#   # create vector of models to compare with
#   compare.cmip <- models[-m]
#   
#   for(r in 1:length(regions)){
#     # r <- 1
#     true.cmip <- cmip %>%
#       filter(region == regions[r],
#              model == models[m]) 
#   
#     # calculate "true" mean sst for 1950:2014 and 2015:2044
#     # these climatologies are critical to FAR calculations (former)
#     # and sst projections under warming (latter)
#     true.1950.2014 <- mean(true.cmip$annual.sst[true.cmip$year %in% 1950:2014])
#     
#     # and calculate sd, AR(1), and trend as in the actual weighting
#     true.sd <- sd(true.cmip$annual.sst[true.cmip$year %in% 1950:2014])
#     true.ar <- ar(true.cmip$annual.sst[true.cmip$year %in% 1950:2014], order.max = 1, aic = F)$ar   
#     true.trend <- summary(lm(annual.sst ~ year, 
#                              data = true.cmip[true.cmip$year %in% 1973:2022,]))$coefficients[2,1]
#     
#     for(c in 1: length(compare.cmip)){
#       # c <- 1
#       comp.cmip <- cmip %>%
#         filter(region == regions[r],
#                model == compare.cmip[c]) 
#       
#       # calculate predicted values
#       pred.1950.2014 <- mean(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014])
#       
#       cmip.sd <- sd(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014])
#       cmip.ar <- ar(comp.cmip$annual.sst[comp.cmip$year %in% 1950:2014], order.max = 1, aic = F)$ar
#       cmip.trend <- summary(lm(annual.sst ~ year, 
#                                data = comp.cmip[comp.cmip$year %in% 1973:2022,]))$coefficients[2,1]
# 
#       # add to output
#       perfect_model_prediction <- rbind(perfect_model_prediction,
#                       data.frame(true_model = models[m],
#                                  region = regions[r],
#                                  prediction_model = compare.cmip[c],
#                                  climatology_diff = abs(true.1950.2014 - pred.1950.2014),
#                                  sd_diff = abs(true.sd - cmip.sd),
#                                  ar_diff = abs(true.ar - cmip.ar),
#                                  trend_diff = abs(true.trend - cmip.trend),
#                                  true.climatology = true.1950.2014,
#                                  true.sd = true.sd,
#                                  true.ar = true.ar,
#                                  true.trend = true.trend,
#                                  pred.climatology = pred.1950.2014,
#                                  pred.sd = cmip.sd,
#                                  pred.ar = cmip.ar,
#                                  pred.trend = cmip.trend))
#       
#     } # close c loop (comparison models)
#   } # close m loop (models)
# } # close r loop (regions) 

## loop through each model/region comparison and calculate prediction model for distance from "true" model

# set a range of sigma_s guided by Zhao et al. 2022)
sigma_s = seq(0.2, 0.8, by = 0.1)

# define a range of possible sigma_d values between 0.2 and 1.1
# (range used in Zhao et al. 2022 Earth's Future) 
sigma_d <- seq(0.2, 2, by = 0.005)

# vector of true models to predict
true_models <- unique(perfect_model_prediction$true_model)

# create object to summarize D and S for each true model prediction
predict_D_S <- data.frame()

# create object to collect tuning results
sigma_s_d_tuning <- data.frame()

# start by looping through sigma_s
for(sig.s in 1:length(sigma_s)){
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
      mutate(diff1 = 2*climatology_diff / median(temp_predict$climatology_diff),
             diff2 = sd_diff / median(temp_predict$sd_diff),
             diff3 = ar_diff / median(temp_predict$ar_diff),
             diff4 = 2*trend_diff / median(temp_predict$trend_diff)) 
    
    temp_predict$D = apply(temp_predict[,16:19], 1, mean)

    # simplify 
    temp_predict <- temp_predict %>%
      select(true_model, region, prediction_model, D, 
             true.climatology, true.sd, true.ar, true.trend, 
             pred.climatology, pred.sd, pred.ar, pred.trend)
    
    # loop through each prediction model and calculate independence weight
    
    # get vector of prediction models
    pred_mods <- unique(temp_predict$prediction_model)
    
    # create data frame to store similarities among prediction models
    all_similarities <- data.frame()
    
    for(p in 1:length(pred_mods)){
      # p <- 1
     
      # do not include similarity with the "true" model we are predicting! 
      temp_sim <- sims %>%
        filter(region == regions[r],
               model == pred_mods[p],
               comparison != temp_predict$true_model) %>% 
        mutate(sim1 = 2*climatology_diff / median(temp_sim$climatology_diff),
              sim2 = sd_diff / median(temp_sim$sd_diff),
              sim3 = ar_diff / median(temp_sim$ar_diff),
              sim4 = 2*trend_diff / median(temp_sim$trend_diff))
    
      temp_sim$S = apply(temp_sim[,8:11], 1, mean)
      
      # collate similarity for every prediction model combination
      all_similarities <- rbind(all_similarities,
                                data.frame(region = regions[r],
                                           true_model = true_models[t],
                                           prediction_model = temp_sim$model,
                                           sigma_s = sigma_s[sig.s],
                                           independence_skill = 1/(1+sum(exp(-(temp_sim$S^2/sigma_s[sig.s]^2)))))) %>%
        group_by(region, true_model, prediction_model) %>%
        summarise(sigma_s = mean(sigma_s),
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
    for(sig.d in 1:length(sigma_d)){
      # sig.d <- 1
      temp_dat <- temp_dat %>%
        mutate(W = exp(-(D^2/sigma_d[sig.d]^2))*independence_skill,
               norm_W = W/mean(W)) # normalize by the mean value of W

      # temp_plot <- temp_dat %>%
      #   mutate(prediction_model = reorder(prediction_model, desc(norm_W)))
      # 
      # ggplot(temp_plot, aes(prediction_model, norm_W)) +
      #   geom_point() +
      #   geom_path()
      
      # calculate 10th, 50th, 90th quantile of weighted prediction
      predict_climatology <- whdquantile(x = temp_dat$pred.climatology, weights =  temp_dat$norm_W, probs = c(0.1, 0.5, 0.9))
      predict_trend <- whdquantile(x = temp_dat$pred.trend, weights =  temp_dat$norm_W, probs = c(0.1, 0.5, 0.9))
      predict_sd <- whdquantile(x = temp_dat$pred.sd, weights =  temp_dat$norm_W, probs = c(0.1, 0.5, 0.9))
      predict_ar <- whdquantile(x = temp_dat$pred.ar, weights =  temp_dat$norm_W, probs = c(0.1, 0.5, 0.9))
      
      # save results
      regional_prediction <- rbind(regional_prediction,
                                   data.frame(region = regions[r],
                                              true_model = true_models[t],
                                              sigma_s = sigma_s[sig.s],
                                              sigma_d = sigma_d[sig.d],
                                              variable = c("mean_1950-2014", "trend_1973-2022", "sd_1950-2014", "ar_1950-2014"),
                                              true = c(mean(temp_dat$true.climatology), mean(temp_dat$true.trend),
                                                       mean(temp_dat$true.sd), mean(temp_dat$true.ar)),
                                              predicted.10 = c(predict_climatology[1], predict_trend[1], predict_sd[1], predict_ar[1]),
                                              predicted.mean = c(predict_climatology[2], predict_trend[2], predict_sd[2], predict_ar[2]),
                                              predicted.90 = c(predict_climatology[3], predict_trend[3], predict_sd[3], predict_ar[3])))
      
    } # close sig.d loop (sigma_d)
    
  } # close t loop (true models)
  
} # close r loop (regions)

# add column to record if true value is within 10th - 90th quantiles

regional_prediction$between <- NA

for(i in 1:nrow(regional_prediction)){
  
regional_prediction$between[i] <- (between = if_else(between(regional_prediction$true[i],
                                                          regional_prediction$predicted.10[i], 
                                                          regional_prediction$predicted.90[i]), 1, 0))

} # 

# for each region, calculate the proportion of predictions in the 10th-90th prediction quantile for each value of sigma-d
# and calculate true-prediction correlation



variables <- unique(regional_prediction$variable)

for(r in 1:length(regions)){
# r <- 1
  
  for(v in 1:length(variables)){
  #   # v <- 1
    
    for(sig.d in 1:length(sigma_d)){
    # sig.d <- 1
      
    temp_dat <- regional_prediction %>%
      filter(region == regions[r],
              variable == variables[v],
             sigma_d == sigma_d[sig.d])
    
    sigma_s_d_tuning <- rbind(sigma_s_d_tuning,
                            data.frame(region = regions[r],
                                       variable = variables[v],
                                       sigma_s = sigma_s[sig.s],
                                       sigma_d = sigma_d[sig.d],
                                       true = temp_dat$true,
                                       predicted = temp_dat$predicted.mean,
                                       proportion_between = sum(temp_dat$between)/
                                         length(temp_dat$between)))
    } # close sig.d loop
    } # close variable loop
} # close region loop
} # close sig.s loop (sigma_s values)
# check <- read.csv("./CMIP6/summaries/tuning_results_sigma_perfect_model.csv")

# save results
write.csv(sigma_s_d_tuning, "./CMIP6/summaries/tuning_results_sigma_perfect_model.csv", row.names = F)

## process results-------------------

sigma_s_d_tuning <- read.csv("./CMIP6/summaries/tuning_results_sigma_perfect_model.csv")


# summarize
perfect_model_out <- sigma_s_d_tuning %>%
  group_by(region, sigma_s, sigma_d) %>%
  summarise(mean_proportion_between = mean(proportion_between)) %>%
  ungroup()

# this is the proportion of "true" values
# within 10th-90th quantiles of weighted 
# prediction across all 23 "true" models - 
# averaged b/c it is the same value in dataframe across
# all 23 "true" models

# check
sum(is.na(perfect_model_out$sigma_d))

# this is the correlation between "true" and weighted 
# prediction for each "true" model

# clean up to plot
perfect_model_plot <- perfect_model_out %>%
  mutate(sigma_s = as.factor(sigma_s))

# clean up and order region names
clean_names <- data.frame(region = unique(perfect_model_plot$region),
                   plot_region = as.factor(c("British Columbia Coast",
                                   "Eastern Bering Sea",
                                   "Gulf of Alaska",
                                   "North Pacific",
                                   "Northern California Current",
                                   "Southern California Current")),
                   plot_order = c(4,2,3,1,5,6))

perfect_model_plot <- left_join(perfect_model_plot, clean_names) %>%
  mutate(plot_region = reorder(plot_region, plot_order))

# custom colorblind colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(perfect_model_plot, aes(sigma_d, mean_proportion_between, color = as.factor(sigma_s))) +
  geom_line() +
  facet_wrap(~plot_region) +
  geom_hline(yintercept = 0.8, lty = 2) +
  guides(color = guide_legend(title = expression(sigma["s"]))) +
  scale_color_manual(values = cb[c(1:4, 6:8)]) +
  xlab(expression(sigma["d"])) +
  ylab(expression("Proportion between 10"^"th"~"and 90"^"th"~"quantiles"))

ggsave("./CMIP6/figs/combined_sigma_d_sigma_s_one_panel_per_region.png", width = 10, height = 6, units = 'in')

# look at effect of sigma-s on independence skill

look <- data.frame(sigma_s = seq(0.2, 0.8, by = 0.1),
                   S = 0.5) 
look$independence_skill = 1/(1+exp(-(look$S^2/look$sigma_s^2)))
look
# as sigma-s increases, penalty for non-independence increases

## get correlations
perfect_model_cor <- sigma_s_d_tuning %>%
  filter(sigma_s == 0.4) %>%
  group_by(region, variable, sigma_d) %>%
  summarize(correlation = cor(true, predicted)) %>%
  ungroup()

ggplot(perfect_model_cor, aes(sigma_d, correlation, color = variable)) +
  geom_line() +
  facet_wrap(~region)


# calculate sigma-d for each system
max_proportion <- perfect_model_out %>%
  filter(sigma_s == 0.4) %>%
  group_by(region) %>%
  summarise(max_prop = max(mean_proportion_between)) # all regions reach 0.8 except SCC

# query ssc separately
scc_sigma <- perfect_model_out %>%
  filter(region == "Southern_California_Current",
         sigma_s == 0.4,
         mean_proportion_between == max_proportion$max_prop[max_proportion$region == "Southern_California_Current"]) 

# now select sigma-d for other regions 
select_sigma <- perfect_model_out %>%
  filter(sigma_s == 0.4) %>%
  mutate(diff = mean_proportion_between - 0.8) %>%
  filter(diff >= 0) %>%
  group_by(region) %>%
  summarise(sigma_d = min(sigma_d)) %>%

select_sigma  <-  rbind(select_sigma, data.frame(region = "Southern_California_Current",
                                                sigma_d = min(scc_sigma$sigma_d)))
select_sigma <- select_sigma %>%
  mutate(sigma_s = 0.4)

# save
write.csv(select_sigma, "./CMIP6/summaries/model_weighting_sigma.csv", row.names = F)

plot_dat <- left_join(select_sigma, sigma_s_d_tuning)  

plot_dat <- left_join(plot_dat, clean_names) %>%
  mutate(plot_region = reorder(plot_region, plot_order))

# and scale to plot!
plot_scale <- plot_dat %>%
  group_by(region, variable) %>% 
  summarise(mean = mean(true),
            sd = sd(true))

plot_dat <- left_join(plot_dat, plot_scale)

# now scale
plot_dat <- plot_dat %>%
  mutate(sc_true = (true-mean)/sd,
         sc_predicted = (predicted-mean)/sd)

ggplot(plot_dat, aes(sc_predicted, sc_true)) +
  geom_point() +
  facet_grid(region~variable, scales = "free")
