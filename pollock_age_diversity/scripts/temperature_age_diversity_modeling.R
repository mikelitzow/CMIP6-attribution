## effect of temperature on age diversity

library(ggplot2)
library(plyr)
library(dplyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./pollock_age_diversity/scripts/stan_utils.R")


## Read in data --------------------------------------------
data <- read.csv("./pollock_age_diversity/data/shannon.age.diversity.acoustic.trawl.data.csv",
                 row.names = 1)


# add sst data
sst <- read.csv("./pollock_age_diversity/data/goa.sst.csv",
                row.names = 1)

data <- left_join(data, sst)

# drop missing years 
data <- na.omit(data)
  

## brms: setup ---------------------------------------------

## Define model formulas - 
## comparing three levels of smoothing in sst data
## (annual means, 2-yr running means, 3-yr running means)
age1_formula <-  bf(shannon ~ s(annual.sst, k = 4) + ar())

age2_formula <-  bf(shannon ~ s(two.yr.sst, k = 4) + ar())

age3_formula <-  bf(shannon ~ s(three.yr.sst, k = 4) + ar())


## Show default priors
get_prior(age1_formula, data)


## Set priors
age_priors <- c(set_prior("normal(0.5, 3)", class = "ar"), # mean ar = 0.5
                set_prior("normal(0, 3)", class = "b"),
                set_prior("student_t(3, 0, 3)", class = "Intercept"),
                set_prior("student_t(3, 0, 3)", class = "sds"),
                set_prior("student_t(3, 0, 3)", class = "sigma"))

## fit models --------------------------------------
age_sst1 <- brm(age1_formula,
                 data = data,
                 prior = age_priors,
                 seed = 1234,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
saveRDS(age_sst1, file = "./pollock_age_diversity/output/age_sst1.rds")

# age_sst1 <- add_criterion(age_sst1, "loo",
#                                moment_match = TRUE)

saveRDS(age_sst1, file = "./pollock_age_diversity/output/age_sst1.rds")

age_sst1 <- readRDS("./pollock_age_diversity/output/age_sst1.rds")
check_hmc_diagnostics(age_sst1$fit)
neff_lowest(age_sst1$fit)
rhat_highest(age_sst1$fit)
summary(age_sst1)
bayes_R2(age_sst1)
plot(age_sst1$criteria$loo, "k")
plot(conditional_smooths(age_sst1), ask = FALSE)
y <- data$shannon
yrep_age_sst1  <- fitted(age_sst1, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_age_sst1[sample(nrow(yrep_age_sst1), 25), ]) +
    xlim(0, 3) +
    ggtitle("age_sst1")
pdf("./pollock_age_diversity/figs/trace_age_sst1.pdf", width = 6, height = 4)
    trace_plot(age_sst1$fit)
dev.off()


age_sst2 <- brm(age2_formula,
                data = data,
                prior = age_priors,
                seed = 1234,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99, max_treedepth = 10))
saveRDS(age_sst2, file = "./output/age_sst2.rds")
# 
# age_sst2 <- add_criterion(age_sst2, "loo",
#                           moment_match = TRUE)

saveRDS(age_sst2, file = "./pollock_age_diversity/output/age_sst2.rds")

age_sst2 <- readRDS("./pollock_age_diversity/output/age_sst2.rds")
check_hmc_diagnostics(age_sst2$fit)
neff_lowest(age_sst2$fit)
rhat_highest(age_sst2$fit)
summary(age_sst2)
bayes_R2(age_sst2)
plot(age_sst2$criteria$loo, "k")
plot(conditional_smooths(age_sst2), ask = FALSE)
y <- data$shannon
yrep_age_sst2  <- fitted(age_sst2, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_age_sst2[sample(nrow(yrep_age_sst2), 25), ]) +
  xlim(0, 3) +
  ggtitle("age_sst2")
pdf("./pollock_age_diversity/figs/trace_age_sst2.pdf", width = 6, height = 4)
trace_plot(age_sst2$fit)
dev.off()


age_sst3 <- brm(age3_formula,
                data = data,
                prior = age_priors,
                seed = 1234,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.9999, max_treedepth = 10))
saveRDS(age_sst3, file = "./pollock_age_diversity/output/age_sst3.rds")

# age_sst3 <- add_criterion(age_sst3, "loo",
#                           moment_match = TRUE)

saveRDS(age_sst3, file = "./pollock_age_diversity/output/age_sst3.rds")

age_sst3 <- readRDS("./pollock_age_diversity/output/age_sst3.rds")
check_hmc_diagnostics(age_sst3$fit)
neff_lowest(age_sst3$fit)
rhat_highest(age_sst3$fit)
summary(age_sst3)
bayes_R2(age_sst3)
plot(age_sst3$criteria$loo, "k")
plot(conditional_smooths(age_sst3), ask = FALSE)
y <- data$shannon
yrep_age_sst3  <- fitted(age_sst3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_age_sst3[sample(nrow(yrep_age_sst3), 25), ]) +
  xlim(0, 3) +
  ggtitle("age_sst3")
pdf("./pollock_age_diversity/figs/trace_age_sst3.pdf", width = 6, height = 4)
trace_plot(age_sst3$fit)
dev.off()

## Model selection -----------------------------------------
age_sst1   <- readRDS("./pollock_age_diversity/output/age_sst1.rds")
age_sst2   <- readRDS("./pollock_age_diversity/output/age_sst2.rds")
age_sst3   <- readRDS("./pollock_age_diversity/output/age_sst3.rds")

# loo(age_sst1, age_sst2, age_sst3, moment_match = T, reloo = T)


## Predicted effects ---------------------------------------
age_sst3   <- readRDS("./pollock_age_diversity/output/age_sst3.rds")

## three-year running mean sst predictions ##

## 95% CI
ce1s_1 <- conditional_effects(age_sst3, effect = "three.yr.sst", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(age_sst3, effect = "three.yr.sst", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(age_sst3, effect = "three.yr.sst", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$three.yr.sst
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$three.yr.sst[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$three.yr.sst[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$three.yr.sst[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$three.yr.sst[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(data$three.yr.sst), amount = 0.1),
                          rep(NA, 100-length(unique(data$three.yr.sst))))

ggplot(dat_ce) +
    aes(x = effect1__, y = estimate__) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
    geom_line(size = 1, color = "red3") +
    labs(x = "Three-year running mean SST (anomaly)", y = "Age diversity") +
    theme_bw()+
    geom_rug(aes(x=rug.anom, y=NULL))


ggsave("./pollock_age_diversity/figs/three.yr_sst_predicted_age_diversity.png", width = 6, height = 4)



## refit best model with FAR values for 3-r running mean SST
## Read in data --------------------------------------------
data <- read.csv("./pollock_age_diversity/data/shannon.age.diversity.acoustic.trawl.data.csv",
                 row.names = 1)


# add sst.3 FAR estimates
far <- read.csv("./CMIP6/summaries/observed_3yr_running_mean_sst_FAR_bayes_estimates.csv",
                row.names = 1)

# just need year and FAR estimate
far <- far %>%
  select(year, estimate__) %>%
  rename(far = estimate__)

data <- left_join(data, far)

# drop missing years 
data <- na.omit(data)


## brms: setup ---------------------------------------------

## Define model formulas

age3_far_formula <-  bf(shannon ~ s(far, k = 4) + ar())


## Show default priors
get_prior(age3_far_formula, data)


## Set priors
age_priors <- c(set_prior("normal(0.5, 3)", class = "ar"), # mean ar = 0.5
                set_prior("normal(0, 3)", class = "b"),
                set_prior("student_t(3, 0, 3)", class = "Intercept"),
                set_prior("student_t(3, 0, 3)", class = "sds"),
                set_prior("student_t(3, 0, 3)", class = "sigma"))

## fit model --------------------------------------
age_far3 <- brm(age3_far_formula,
                data = data,
                prior = age_priors,
                seed = 1234,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.9999, max_treedepth = 10))
saveRDS(age_far3, file = "./pollock_age_diversity/output/age_far3.rds")

check_hmc_diagnostics(age_far3$fit)
neff_lowest(age_far3$fit)
rhat_highest(age_far3$fit)
summary(age_far3)
bayes_R2(age_far3)
plot(conditional_smooths(age_far3), ask = FALSE)
y <- data$shannon
yrep_age_far3  <- fitted(age_far3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_age_far3[sample(nrow(yrep_age_far3), 25), ]) +
  xlim(0, 3) +
  ggtitle("age_far3")
pdf("./pollock_age_diversity/figs/trace_age_far3.pdf", width = 6, height = 4)
trace_plot(age_far3$fit)
dev.off()


## Predicted effects ---------------------------------------
age_far3   <- readRDS("./pollock_age_diversity/output/age_far3.rds")

## three-year running mean sst predictions ##

## 95% CI
ce1s_1 <- conditional_effects(age_far3, effect = "far", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(age_far3, effect = "far", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(age_far3, effect = "far", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$far
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$far[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$far[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$far[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$far[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(data$far), amount = 0.01),
                          rep(NA, 100-length(unique(data$far))))

ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "FAR: three-year running mean SST", y = "Age diversity") +
  theme_bw()+
  geom_rug(aes(x=rug.anom, y=NULL))


ggsave("./pollock_age_diversity/figs/three.yr_sst_predicted_age_diversity.png", width = 6, height = 4)
