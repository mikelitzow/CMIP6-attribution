## experimental - logistic model to estimate preindustrial and historical 
## probabilities for given sst anomalies

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## load preindustrial and historical outcome for the GOA
preindustrial <- read.csv("./CMIP6/summaries/Gulf_of_Alaska_preindustrial_outcomes.csv")
historical <- read.csv("./CMIP6/summaries/Gulf_of_Alaska_historical_outcomes.csv")

# combine
outcomes <- rbind(preindustrial, historical)

# load model weights (based on ar(1), correlation, bias)
weights <- read.csv("./CMIP6/summaries/CMIP6_model_weights_by_region_window.csv")

weights <- weights %>%
   filter(window == "annual",
          region == "Gulf_of_Alaska") %>%
   select(model, scaled.total.weight) %>%
   rename(model_weight = scaled.total.weight)

outcomes <- left_join(outcomes, weights)




##################################

library(data.table)
dt <- as.data.table(outcomes)
dt[model == "ACCESS-CM2" & period == "historical", .(N = .N), by = .(period, model, ersst.year)]
dt[model == "ACCESS-CM2" & period == "historical" & ersst.year == 2017, ]

prop <- dt[ , .(count = sum(annual.1yr.events),
                N = .N,
                annual.anomaly.1yr = unique(annual.anomaly.1yr),
                model_weight = unique(model_weight)),
           by = .(period, model, ersst.year)]
prop[ , prop := count / N]
prop[ , model_fac := as.factor(model)]

head(prop[period == "preindustrial", ])
head(prop[period == "historical", ])

g <- ggplot(prop) +
    geom_point(aes(x = annual.anomaly.1yr, y = prop, color = model_fac)) +
    facet_wrap( ~ period) +
    theme(legend.position="none")
print(g)


## What is the response here?
form <-  bf(count | trials(N) + weights(model_weight, scale = TRUE) ~
            period + s(annual.anomaly.1yr, by = period, k = 6) + (1 | model_fac))

far_brms2 <- brm(form,
                 data = prop,
                 family = binomial(link = "logit"),
                 seed = 1234,
                 cores = 4, chains = 4, iter = 4000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.9, max_treedepth = 13))
saveRDS(far_brms2, "./CMIP6/brms_output/Gulf_of_Alaska_binomial2.rds")
## ~2.5 hours run time

# far_brms2 <- readRDS("./CMIP6/brms_output/Gulf_of_Alaska_binomial2.rds")

ce <- conditional_effects(far_brms2)
plot(ce, ask = FALSE)

ce2 <- conditional_effects(far_brms2, effect = "annual.anomaly.1yr:period", re_formula = NULL)
plot(ce2)
head(ce2[[1]])



## setup new data
x <- ce2[[1]]
nd <- x[ , c("period", "annual.anomaly.1yr", "N", "model_fac")]
nd$N <- 100
nd$count <- 1
nd_pre <- nd[nd$period == "preindustrial", ]
nd_his <- nd[nd$period == "historical", ]

## make predictions
pre_pp <- posterior_predict(far_brms2, newdata = nd_pre, re_formula = NULL)
his_pp <- posterior_predict(far_brms2, newdata = nd_his, re_formula = NULL)

## Calc probabilities
## These are our posterior probabilities to use for FAR calculation
pre_prob <- pre_pp / 100
his_prob <- his_pp / 100


sst_obs <- nd_pre$annual.anomaly.1yr

pre_prob_mean <- apply(pre_prob, 2, mean)
his_prob_mean <- apply(his_prob, 2, mean)
plot(x = sst_obs, y = pre_prob_mean, type = "l")
lines(x = sst_obs, y = his_prob_mean, col = "red3")



##################################


## brms: setup ---------------------------------------------

# setup variables - model as factor
outcomes$model_fac <- as.factor(outcomes$model)

## fit: brms --------------------------------------

# define model formula

far_formula <-  bf(annual.1yr.events | weights(model_weight, scale = TRUE) ~
                      s(annual.anomaly.1yr, by = period, k = 6) + period + (1 | model_fac))

# run with default priors

far_brms <- brm(far_formula,
                     data = outcomes,
                     family = bernoulli(link = "logit"),
                     cores = 4, chains = 4, iter = 3000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 15))

saveRDS(far_brms, "./CMIP6/brms_output/Gulf_of_Alaska_binomial.rds")
