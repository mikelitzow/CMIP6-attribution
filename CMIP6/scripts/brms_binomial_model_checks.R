# check binomial models used to estimate FAR for different "presents"
library(tidyverse)
library(brms)
library(rstan)

source("./CMIP6/scripts/stan_utils.R")


model <- readRDS("./CMIP6/brms_output/North_Pacific_rolling_window_binomial2.rds") # checks out

check_hmc_diagnostics(model$fit)
neff_lowest(model$fit) 
rhat_highest(model$fit)
summary(model)
bayes_R2(model) # wow - 0.99!
trace_plot(model$fit)

##
model <- readRDS("./CMIP6/brms_output/Eastern_Bering_Sea_rolling_window_binomial2.rds") # checks out

check_hmc_diagnostics(model$fit)
neff_lowest(model$fit) 
rhat_highest(model$fit)
summary(model)
bayes_R2(model) # wow - 0.99!
trace_plot(model$fit)

##
model <- readRDS("./CMIP6/brms_output/Gulf_of_Alaska_rolling_window_binomial2.rds") # checks out

check_hmc_diagnostics(model$fit)
neff_lowest(model$fit) 
rhat_highest(model$fit)
summary(model)
bayes_R2(model) # 0.99
trace_plot(model$fit)

##
model <- readRDS("./CMIP6/brms_output/British_Columbia_Coast_rolling_window_binomial2.rds") # checks out

check_hmc_diagnostics(model$fit) # 2 divergent transitions
neff_lowest(model$fit) 
rhat_highest(model$fit)
summary(model)
bayes_R2(model)
trace_plot(model$fit)

##
model <- readRDS("./CMIP6/brms_output/Northern_California_Current_rolling_window_binomial2.rds") # checks out

check_hmc_diagnostics(model$fit) # two divergent transitions
neff_lowest(model$fit) 
rhat_highest(model$fit)
summary(model)
bayes_R2(model)
trace_plot(model$fit)


##
model <- readRDS("./CMIP6/brms_output/Southern_California_Current_rolling_window_binomial2.rds") # checks out

check_hmc_diagnostics(model$fit)
neff_lowest(model$fit) # a little low
rhat_highest(model$fit)
summary(model)
bayes_R2(model)
trace_plot(model$fit)

