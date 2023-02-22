# calculate weights from ersst differences, model independence, and shaping parameters

library(tidyverse)

# set theme
theme_set(theme_bw())

# custom colorblind colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## load differences, similarities, and shape variables and calculate weight-----------------------

diff <- read.csv("./CMIP6/summaries/CMIP6_time_series_differences.csv")
sim <- read.csv("./CMIP6/summaries/CMIP6_time_series_model_similarities.csv")
shape <- read.csv("./CMIP6/summaries/model_weighting_sigma.csv")


# object to catch results
CMIP6_weights <- data.frame()

regions <- unique(diff$region)

# set sigma_s = 0.4
sigma_s = 0.4

for(r in 1:length(regions)){
  # r <- 1
  
  diff_temp <- diff %>% 
    filter(region == regions[r])

  shape_temp <- shape %>%
    filter(region == regions[r])
  
  # scale differences for this region by the median and get mean difference (D)
  diff_temp <- diff_temp %>%
    mutate(diff1 = climatology_diff/median(climatology_diff),
           diff2 = sd_diff/median(sd_diff),
           diff3 = ar_diff/median(ar_diff),
           diff4 = trend_diff/median(trend_diff))
  
  diff_temp$D = apply(diff_temp[,7:10], 1, mean)
  
  # loop through each model and calculate similarity to other models
  models <- diff_temp$model
  
  for(m in 1:length(models)){
    # m <- 1
    
    sim_temp <- sim %>%
      filter(region == regions[r],
             model == models[m]) %>% 
      mutate(sim1 = climatology_diff / median(climatology_diff),
             sim2 = sd_diff / median(sd_diff),
             sim3 = ar_diff / median(ar_diff),
             sim4 = trend_diff / median(trend_diff))
    
    sim_temp$S = apply(sim_temp[,8:11], 1, mean)
    
    # get independence weight
      independence_weight = 1/(1+sum(exp(-(sim_temp$S^2/sigma_s^2))))

    # # check calculation above
    # check <- NA
    # 
    # for(i in 1:nrow(sim_temp)){
    # check[i] <- exp(-(sim_temp$S[i]^2/sigma_s^2))
    # }
    # 
    # check_independence <- 1/(1+sum(check))
    # check_independence # correct!
    
    # save results
      
    CMIP6_weights <- rbind(CMIP6_weights,
                           data.frame(region = regions[r],
                                      model = models[m],
                                      skill_weight = exp(-(diff_temp$D[diff_temp$model== models[m]]^2/shape_temp$sigma_d^2)),
                                      independence_weight = independence_weight))
    
  }
  
}

# get combined weight
CMIP6_weights <- CMIP6_weights %>%
  mutate(combined_weight = skill_weight*independence_skill)

# plot non-normalized
plot_CMIP6_weights <- CMIP6_weights %>%
  arrange(region, desc(combined_weight)) %>%
  mutate(model_number = rep(1:23, length.out = 6*23)) %>%
  pivot_longer(cols = c(-region, -model, -model_number))

# clean up and order region names
clean_names <- data.frame(region = unique(plot_CMIP6_weights$region),
                          plot_region = as.factor(c("British Columbia Coast",
                                                    "Eastern Bering Sea",
                                                    "Gulf of Alaska",
                                                    "North Pacific",
                                                    "Northern California Current",
                                                    "Southern California Current")),
                          plot_order = c(4,2,3,1,5,6))

plot_CMIP6_weights <- left_join(plot_CMIP6_weights, clean_names) %>%
  mutate(plot_region = reorder(plot_region, plot_order)) %>%
  mutate(weight_order = case_when(name == "skill_weight" ~ 1,
                                  name == "independence_weight" ~ 2,
                                  name == "combined_weight" ~3),
         name = reorder(name, weight_order)) 


ggplot(plot_CMIP6_weights, aes(model_number, value, color = name)) +
  geom_point() +
  geom_line() + 
  facet_wrap(~plot_region, scales = "free_y") +
  scale_color_manual(values = cb[c(6,4,2)], labels = c("Skill weight", "Independence weight", "Combined weight")) +
  labs(x = "Model number",
       y = "Value") +
  theme(legend.title = element_blank())

ggsave("./CMIP6/figs/combined_weights_one_panel_per_region.png", width = 10, height = 5, units = 'in')

## plot normalized combined weight
CMIP6_normalized_weights <- data.frame()

# now normalize so that average weight = 1 in each region
for(r in 1:length(regions)){
  # r <- 1
  
  temp_weight <- CMIP6_weights %>%
    filter(region == regions[r])
  
  CMIP6_normalized_weights <- rbind(CMIP6_normalized_weights,
                                    data.frame(region = regions[r],
                                               model = temp_weight$model,
                                               combined_weight = temp_weight$combined_weight/mean(temp_weight$combined_weight)))
}

CMIP6_normalized_weights <- CMIP6_normalized_weights %>%
  arrange(region, desc(combined_weight)) %>%
  mutate(model_number = rep(1:23, length.out = 6*23)) 

# clean up and order names
CMIP6_normalized_weights <- left_join(CMIP6_normalized_weights, clean_names) %>%
  mutate(plot_region = reorder(plot_region, plot_order)) 


ggplot(CMIP6_normalized_weights, aes(model_number, combined_weight)) +
  geom_point() +
  geom_line() + 
  geom_hline(yintercept = 1, lty = 2) +
  facet_wrap(~plot_region) + 
  labs(x = "Model number",
       y = "Model weight")

ggsave("./CMIP6/figs/total_normalized_weight_one_panel_per_region.png", width = 10, height = 5, units = 'in')

ggplot(CMIP6_normalized_weights, aes(combined_weight)) +
  geom_histogram(bins = 15, fill = "grey", color = "black") +
  facet_wrap(~plot_region)



## plot sst time series relative to observations, color coding for model weight

cmip <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.csv") %>%
  filter(experiment == "hist_ssp585") %>%
  select(region, model, year, annual.unsmoothed) %>%
  rename(annual.sst = annual.unsmoothed)

ersst <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv") %>%
  select(region, year, annual.unsmoothed) %>%
  rename(annual.sst = annual.unsmoothed)

# data to plot
plot <- cmip %>%
  rename(model_sst = annual.sst) %>%
  left_join(.,ersst) %>%
  rename(observed_sst = annual.sst) 

# add model weights
plot_weights <- CMIP6_normalized_weights %>%
  select(-model_number) %>%
  mutate(proportional_weight = case_when(combined_weight >= 2 ~ "> 2",
                                  combined_weight >= 1.5 & combined_weight < 2 ~ "1.5 - 2",
                                  combined_weight >= 1.25 & combined_weight < 1.5 ~ "1.25 - 1.5",
                                  combined_weight >= 0.75 & combined_weight < 1.25 ~ "0.75 - 1.25",
                                  combined_weight >= 0.5 & combined_weight < 0.75 ~ "0.5 - 0.75",
                                  combined_weight < 0.5 ~ "< 0.5"
                                  ),
         weight_order = case_when(proportional_weight == "> 2" ~ 6,
                                  proportional_weight == "1.5 - 2" ~ 5,
                                  proportional_weight == "1.25 - 1.5" ~ 4,
                                  proportional_weight == "0.75 - 1.25" ~ 3,
                                  proportional_weight == "0.5 - 0.75" ~ 2,
                                  proportional_weight == "< 0.5" ~ 1),
         proportional_weight = reorder(proportional_weight, weight_order)) %>%
  ungroup()
 
plot <- left_join(plot, plot_weights) 

plot <- left_join(plot, clean_names)

# set colors
library(RColorBrewer)

display.brewer.pal(n = 9, "YlOrRd")
display.brewer.pal(n = 9, "Greys")
my.col1 <- brewer.pal(n = 9, "YlOrRd")
my.col2 <- brewer.pal(n = 9, "Greys")

my.col <- c(my.col2[c(3,5)], my.col1[c(3,5,7,9)])

ggplot(plot, aes(year, model_sst, group = interaction(model, proportional_weight), color = proportional_weight)) +
  geom_line(size = 0.3, alpha = 0.5) +
  facet_wrap(~plot_region, scales = "free_y") +
  scale_color_discrete(type = my.col) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  geom_line(aes(year, observed_sst), color = "black", size = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "SST (°C)",
       color = "Model weight")

ggsave("./CMIP6/figs/weighted_models_vs_ersst.png", width = 9, height = 4.5, units = 'in')
