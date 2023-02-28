## estimate preindustrial and historical 
## probabilities for given 3-yr rolling mean sst anomalies

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

# load region names
regions <- read.csv("./CMIP6/summaries/clean_region_names.csv")
regions <- regions[1:6,1]

for(i in 1:length(regions)){
    # i <-5
## load preindustrial and historical outcomes
preindustrial <- read.csv(paste("./CMIP6/summaries/", regions[i], "_preindustrial_outcomes.csv", sep = ""))

historical <- read.csv(paste("./CMIP6/summaries/", regions[i], "_historical_outcomes_rolling_window.csv", sep = ""))


# combine
outcomes <- rbind(preindustrial, historical)

# load model weights (based on ar(1), correlation, bias)
weights <- read.csv("./CMIP6/summaries/normalized_CMIP6_weights.csv")

weights <- weights %>%
   filter(region == regions[i]) %>%
   select(model, normalized_weight) %>%
   rename(model_weight = normalized_weight)

outcomes <- left_join(outcomes, weights) %>%
    select(model, model_weight, period, region, ersst.year, annual.anomaly.1yr, annual.1yr.events)

# # drop NAs due to smoothing
# drop <- is.na(outcomes$annual.anomaly.3yr)
# outcomes <- outcomes[!drop,]

## brms: setup ---------------------------------------------

library(data.table)
dt <- as.data.table(outcomes)
dt[model == "ACCESS-CM2" & period == "historical", .(N = .N), by = .(period, model, ersst.year)]
dt[model == "ACCESS-CM2" & period == "historical" & ersst.year == 2017, ]

prop <- dt[ , .(count = sum(annual.3yr.events, na.rm = T),
                N = .N,
                annual.anomaly.1yr = unique(annual.anomaly.1yr),
                model_weight = unique(model_weight)),
           by = .(period, model, ersst.year)]
prop[ , prop := count / N]
prop[ , model_fac := as.factor(model)]
prop[ , year_fac := as.factor(ersst.year)]


head(prop[period == "preindustrial", ])
head(prop[period == "historical", ])

g <- ggplot(prop) +
    geom_point(aes(x = year_fac, y = prop, color = model_fac)) +
    facet_wrap( ~ period) +
    theme(legend.position="none")
print(g)


form <-  bf(count | trials(N) + weights(model_weight, scale = TRUE) ~
            period + s(annual.anomaly.1yr, by = period, k = 6) +
            s(ersst.year, by = period, k = 6) + (1 | model_fac))

far_brms2 <- brm(form,
                 data = prop,
                 family = binomial(link = "logit"),
                 seed = 1234,
                 cores = 4, chains = 4, iter = 4000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 15))

# saveRDS(far_brms2, paste("./CMIP6/brms_output/", regions[i], "_binomial2.rds", sep = ""))

saveRDS(far_brms2, paste("./CMIP6/brms_output/", regions[i], "mean_annual_sst_rolling_window_binomial2.rds", sep = ""))

 }
## ~2.5 hours run time

far_brms2 <- readRDS("./CMIP6/brms_output/Gulf_of_Alaska_binomial2.rds")

ce <- conditional_effects(far_brms2)
plot(ce, ask = FALSE)

ce2 <- conditional_effects(far_brms2, effect = "annual.anomaly.1yr:period", re_formula = NA)
plot(ce2)
head(ce2[[1]])
tail(ce2[[1]])



## setup new data
x <- ce2[[1]]
nd <- x[ , c("period", "annual.anomaly.1yr", "N", "model_fac")]
nd$N <- 1000
nd_pre <- nd[nd$period == "preindustrial", ]
nd_his <- nd[nd$period == "historical", ]

## make predictions
## exclude random effects for model_fac
pre_pp <- posterior_epred(far_brms2, newdata = nd_pre, re_formula = NA)
his_pp <- posterior_epred(far_brms2, newdata = nd_his, re_formula = NA)

## Calc probabilities
## These are our posterior probabilities to use for FAR calculation
pre_prob <- pre_pp / unique(nd$N)
his_prob <- his_pp / unique(nd$N)

range(pre_prob)
range(his_prob)
# plot(as.vector(pre_prob))
# plot(as.vector(his_prob))



## Plot conditional effects
pre_pred <- data.frame(period = "preindustrial",
                       annual.anomaly.1yr = nd_pre$annual.anomaly.1yr,
                       prob = apply(pre_prob, 2, mean),
                       lower = apply(pre_prob, 2, quantile, probs = 0.025),
                       upper = apply(pre_prob, 2, quantile, probs = 0.975))
his_pred <- data.frame(period = "historical",
                       annual.anomaly.1yr = nd_his$annual.anomaly.1yr,
                       prob = apply(his_prob, 2, mean),
                       lower = apply(his_prob, 2, quantile, probs = 0.025),
                       upper = apply(his_prob, 2, quantile, probs = 0.975))
pred <- rbind(pre_pred, his_pred)

g <- ggplot(pred) +
    geom_line(aes(x = annual.anomaly.1yr, y = prob, color = period), size = 0.8) +
    geom_ribbon(aes(x = annual.anomaly.1yr, ymin = lower, ymax = upper, fill = period), alpha = 0.15) +
    scale_color_manual(values = c("tomato", "steelblue")) +
    scale_fill_manual(values = c("tomato", "steelblue"))
print(g)



## Calc FAR
far <- 1 - (pre_prob / his_prob)
range(far, na.rm = TRUE)


far_pred <- data.frame(annual.anomaly.1yr = nd_pre$annual.anomaly.1yr,
                       prob = apply(far, 2, mean),
                       lower = apply(far, 2, quantile, probs = 0.025),
                       upper = apply(far, 2, quantile, probs = 0.975))

g <- ggplot(far_pred) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_line(aes(x = annual.anomaly.1yr, y = prob), size = 0.8) +
    geom_ribbon(aes(x = annual.anomaly.1yr, ymin = lower, ymax = upper), alpha = 0.15)
print(g)
