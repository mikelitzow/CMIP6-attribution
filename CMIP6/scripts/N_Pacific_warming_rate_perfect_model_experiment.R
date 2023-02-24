# tune sigma-d for model weights for estimating N. Pacific warming rate 
# fixing sigma-s = 0.4 based on earlier experiments

library(tidyverse)

source("./CMIP6/scripts/weighted_quantiles.R")

# set theme
theme_set(theme_bw())

# custom colorblind colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load model weights 
# (RMSE for comparison between model and ERSST anomalies for the years 1973:2022 [wrt 1850/54 - 1949])
weights <- read.csv("./CMIP6/summaries/N_Pac_warming_model_weights.csv")

# load model anomaly time series (for prediction)
anomalies <- read.csv("./CMIP6/summaries/ne_pacific_annual_modeled_sst.csv", row.names = 1) %>%
  rename(model = name)

# compare with updated anomaly time series (that have been corrected by removing lat < 20N)
check <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv") %>%
  filter(region == "North_Pacific",
         experiment == "hist_ssp585") %>%
  select(model, year, annual.unsmoothed) %>%
  rename(sst.check = annual.unsmoothed)

compare <- left_join(anomalies, check)

ggplot(compare, aes(value, sst.check)) +
  geom_point() +
  facet_wrap(~model, scales = "free") 
## THIS IS THE INCORRECT COMPARISON - CMIP6.anomaly.time.series.csv IS IN UNITS OF SD; ne_pacific_annual_modeled_sst.csv IS IN UNITS OF DEGREES

# vector of model names
models <- weights$model

# loop through each model as "truth" and predict
# mean anomaly for 2001-2015 and 2031-2045

perfect_model_prediction  <- data.frame()

for(m in 1:length(models)){
  # m <- 2
  
  # create vector of models to compare with
  compare.cmip <- models[-m]
  
  # get anomaly time series for "true" model
    true.anom <- anomalies %>%
      filter(name == models[m])
    
    # calculate "true" mean anomaly for 2001-2015 and 2031-2045
    true.2001.2015 <- mean(true.anom$value[true.anom$year %in% 2001:2015])
    true.2031.2045 <- mean(true.anom$value[true.anom$year %in% 2031:2045])    
    
    # calculate predicted anomalies for each period for each prediction model
    for(c in 1: length(compare.cmip)){
      # c <- 1
      comp.anom <- anomalies %>%
        filter(name == compare.cmip[c])
      
      # calculate predicted values
      pred.2001.2015 <- mean(comp.anom$value[comp.anom$year %in% 2001:2015])
      pred.2031.2045 <- mean(comp.anom$value[comp.anom$year %in% 2031:2045])
      

      # add to output
      perfect_model_prediction <- rbind(perfect_model_prediction,
                                        data.frame(true_model = models[m],
                                                   prediction_model = compare.cmip[c],
                                                   true.2001.2015 = true.2001.2015,
                                                   true.2031.2045 = true.2031.2045,
                                                   pred.2001.2015 = pred.2001.2015,
                                                   pred.2031.2045 = pred.2031.2045))
      
    } # close c loop (comparison models)
  } # close m loop (models)

# loop through each model and calculate pairwise difference from every other model

pairwise_similarities <- data.frame()

for(m in 1:length(models)){
  # m <- 1
  
  # create vector of models to compare with
  compare.cmip <- models[-m]
  
  temp.sim <- weights %>%
    filter(model != models[m])
  
  pairwise_similarities <- rbind(pairwise_similarities,
                                 data.frame(model = models[m],
                                            similarity = abs(weights$RMSE[weights$model == models[m]]-
                                              temp.sim$RMSE)))
  
}

# scale by median
pairwise_similarities <- pairwise_similarities %>%
  mutate(sc_similarity = similarity / median(similarity))

# calculate independence weight for each model, using sigma-s = 0.4

independence <- data.frame()

for(m in 1:length(models)){
  # m <- 1
  
  temp_sim <- pairwise_similarities %>%
    filter(model == models[m])
  
  
  this_weight = 1/(1+sum(exp(-(temp_sim$sc_similarity^2/0.4^2))))
  
  independence <- rbind(independence,
                              data.frame(model = models[m],
                                         independence_weight = this_weight))
  
}

# scale RMSE by median to get scaled differences and join with independence weight
weights <- weights %>%
  mutate(D = RMSE / median(RMSE))

weights <-  left_join(weights, independence)

# loop through possible sigma-d and record proportion of true values in 10th-90th quantiles of prediction

sigma_d <- seq(0.2, 2.5, by = 0.005)

prediction_out <- data.frame()

# rename model column in weights
weights <- weights %>%
  rename(prediction_model = model)

for(s in 1:length(sigma_d)){
  # s <- 1
  
  for(m in 1:length(models)){
    # m <- 1
    
  temp_pred <- perfect_model_prediction %>%
    filter(true_model == models[m])
  
  temp_pred <- left_join(temp_pred, weights)
  
  # calculate combined skill and independence weight
  temp_pred <- temp_pred %>%
    mutate(skill_weight = (exp(-(D^2/sigma_d[s]^2))),
           combined_weight = skill_weight*independence_skill,
           normalized_weight = combined_weight/mean(combined_weight))
  
  # get weighted predictions for each period
  wpred_2001_2015 <- whdquantile(temp_pred$pred.2001.2015, weights = temp_pred$normalized_weight, probs = c(0.1, 0.5, 0.9))
  wpred_2031_2045 <- whdquantile(temp_pred$pred.2031.2045, weights = temp_pred$normalized_weight, probs = c(0.1, 0.5, 0.9))
  
  # and save results
  prediction_out <- rbind(prediction_out,
                          data.frame(sigma_d = sigma_d[s],
                                     true = c(mean(temp_pred$true.2001.2015), mean(temp_pred$true.2031.2045)),
                                     weighted_10 = c(wpred_2001_2015[1], wpred_2031_2045[1]),
                                     weighted_mean = c(wpred_2001_2015[2], wpred_2031_2045[2]),
                                     weighted_90 = c(wpred_2001_2015[3], wpred_2031_2045[3]))) 
  }
  
}

# save 
write.csv(prediction_out, "./CMIP6/summaries/N_Pac_warming_perfect_model_tuning.csv", row.names = F)


prediction_out <- read.csv("./CMIP6/summaries/N_Pac_warming_perfect_model_tuning.csv")

# determine whether each true value falls in the 10th-90th quantile band
prediction_out$between <- NA

for(i in 1:nrow(prediction_out)){
  
  prediction_out$between[i] <- if_else(between(prediction_out$true[i],
                                               prediction_out$weighted_10[i],
                                               prediction_out$weighted_90[i]), 1, 0)
  
}


# summarize results
summary_out <- prediction_out %>%
  group_by(sigma_d) %>%
  summarise(proportion = sum(between)/length(between))

ggplot(summary_out, aes(sigma_d, proportion)) +
  geom_line() +
  geom_hline(yintercept = 0.8, lty = 2) + xlab(expression(sigma["d"])) +
  ylab(expression("Proportion between 10"^"th"~"and 90"^"th"~"quantiles"))

ggsave("./CMIP6/figs/N_Pac_warming_sigma_d_tuning.png", width = 6, height = 4, units = 'in')

# find minimum value of sigma-d that produces proportion of 0.8
keep <- summary_out %>%
  filter(proportion >= 0.8)

# sigma_d = 0.82

# figure skill weight and combined weight with this value of sigma_d
weights <- weights %>%
  mutate(skill_weight = (exp(-(D^2/0.82^2))),
         combined_weight = skill_weight*independence_weight,
         normalized_weight = combined_weight/mean(combined_weight)) %>%
  arrange(desc(skill_weight)) %>%
  rename("Independence weight" = independence_weight,
         "Skill weight" = skill_weight,
         "Combined weight" = combined_weight)

# wrangle to plot
plot_weights <- weights %>%
  select(-RMSE, -D, -normalized_weight) %>% 
  mutate(model_number = 1:23) %>%
  pivot_longer(cols = c(-prediction_model, -model_number))
  
# order
plot_order <- data.frame(name = unique(plot_weights$name),
                         order = c(2, 1, 3))

plot_weights <- left_join(plot_weights, plot_order) %>%
  mutate(name = reorder(name, order))

ggplot(plot_weights, aes(model_number, value, color = name)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = cb[c(6,4,2)], labels = c("Skill weight", "Independence weight", "Combined weight")) +
  labs(x = "Model number",
       y = "Value") +
  theme(legend.title = element_blank())

ggsave("./CMIP6/figs/N_Pac_warming_skill_independence_combined_weight.png", width = 6, height = 3.25, units = 'in')

# plot normalized combined weights
plot_weights <- weights %>%
  arrange(desc(normalized_weight)) %>%
  mutate(model_number = 1:23)

ggplot(plot_weights, aes(model_number, normalized_weight)) +
  geom_point() +
  geom_line() +
  labs(x = "Model number",
       y = "Model weight")

# and save
write.csv(weights, "./CMIP6/summaries/N_Pac_warming_model_weights.csv", row.names = F)
