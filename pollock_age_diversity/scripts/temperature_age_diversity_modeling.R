## effect of temperature on age diversity

library(ggplot2)
library(plyr)
library(dplyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------
data <- read.csv("./data/shannon.age.diversity.acoustic.trawl.data.csv",
                 row.names = 1)


# add sst data
sst <- read.csv("./data/goa.sst.csv",
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
saveRDS(age_sst1, file = "./output/age_sst1.rds")

# age_sst1 <- add_criterion(age_sst1, "loo",
#                                moment_match = TRUE)

saveRDS(age_sst1, file = "./output/age_sst1.rds")

age_sst1 <- readRDS("./output/age_sst1.rds")
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
pdf("./figs/trace_age_sst1.pdf", width = 6, height = 4)
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

saveRDS(age_sst2, file = "./output/age_sst2.rds")

age_sst2 <- readRDS("./output/age_sst2.rds")
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
pdf("./figs/trace_age_sst2.pdf", width = 6, height = 4)
trace_plot(age_sst2$fit)
dev.off()


age_sst3 <- brm(age3_formula,
                data = data,
                prior = age_priors,
                seed = 1234,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.9999, max_treedepth = 10))
saveRDS(age_sst3, file = "./output/age_sst3.rds")

# age_sst3 <- add_criterion(age_sst3, "loo",
#                           moment_match = TRUE)

saveRDS(age_sst3, file = "./output/age_sst3.rds")

age_sst3 <- readRDS("./output/age_sst3.rds")
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
pdf("./figs/trace_age_sst3.pdf", width = 6, height = 4)
trace_plot(age_sst3$fit)
dev.off()

## Model selection -----------------------------------------
age_sst1   <- readRDS("./output/age_sst1.rds")
age_sst2   <- readRDS("./output/age_sst2.rds")
age_sst3   <- readRDS("./output/age_sst3.rds")

# loo(age_sst1, age_sst2, age_sst3, moment_match = T, reloo = T)


## Predicted effects ---------------------------------------
age_sst3   <- readRDS("./output/age_sst3.rds")

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


ggsave("./figs/temp.anom_predicted_effect_cod2sg_zinb_k3.png", width = 3, height = 2)
## this is Fig. 2a in the draft


## Julian predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod2sg_zinb_k3, effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod2sg_zinb_k3, effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod2sg_zinb_k3, effect = "julian", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$julian
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$julian[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$julian[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$julian[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$julian[["lower__"]]

g <- ggplot(dat_ce) +
    aes(x = effect1__, y = estimate__) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
    geom_line(size = 1.5, color = "red3") +
    coord_trans(y = "pseudo_log") +
    labs(x = "Day of year", y = "Cod abundance") +
    theme_bw()
print(g)
ggsave("./figs/julian_predicted_effect_cod2sg_zinb_k3.png", width = 5, height = 4)


## SSB predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod2sg_zinb_k3, effect = "ssb", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod2sg_zinb_k3, effect = "ssb", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod2sg_zinb_k3, effect = "ssb", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$ssb
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$ssb[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$ssb[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$ssb[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$ssb[["lower__"]]

g <- ggplot(dat_ce) +
    aes(x = effect1__, y = estimate__) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
    geom_line(size = 1.5, color = "red3") +
    labs(x = "SSB", y = "Cod abundance") +
    coord_trans(y = "pseudo_log") +
    theme_bw()
print(g)
ggsave("./figs/SSB_predicted_effect_cod2sg3_zinb_k3.png", width = 5, height = 4)
