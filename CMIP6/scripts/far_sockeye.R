## SST-catch and FAR-catch relationships for GOA sockeye

library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

# set palette for plotting
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## Read in far,sst, and catch data --------------------------------------------
data <- read.csv("./CMIP6/data/GOA_sockeye_catch_sst_far_no_s.peninsula.csv")


## Analyze sst-catch relationship ---------------------------------------------

# plot time series and experienced sst
ggplot(data, aes(annual_sst_3, log_catch)) +
  geom_point()

# fit brms model
sst_catch_formula <- bf(log_catch_stnd ~ s(annual_sst_3) +  year)

sst_catch_brm <- brm(sst_catch_formula,
                      data = data,
                     seed = 989,
                      cores = 4, chains = 4, iter = 2000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))


saveRDS(sst_catch_brm, file = "./CMIP6/brms_output/sst_catch.rds")

sst_catch_brm <- readRDS("./CMIP6/brms_output/sst_catch.rds")
check_hmc_diagnostics(sst_catch_brm$fit)
neff_lowest(sst_catch_brm$fit)
rhat_highest(sst_catch_brm$fit)
summary(sst_catch_brm)
bayes_R2(sst_catch_brm)

conditional_effects(sst_catch_brm)

y <- data$log_catch_stnd
yrep_sst_catch_brm  <- fitted(sst_catch_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_sst_catch_brm[sample(nrow(yrep_sst_catch_brm), 25), ]) +
  xlim(-6, 6) 

trace_plot(sst_catch_brm$fit)



##############

# plot time series and experienced FAR
ggplot(data, aes(annual_far_3, log_catch)) +
  geom_point()


# fit brms model
far_catch_formula <- bf(log_catch_stnd ~ s(annual_far_3) + year)

far_catch_brm <- brm(far_catch_formula,
                     data = data,
                     cores = 4, chains = 4, iter = 2000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 10))


saveRDS(far_catch_brm, file = "./CMIP6/brms_output/far_catch.rds")

far_catch_brm <- readRDS("./CMIP6/brms_output/far_catch.rds")
check_hmc_diagnostics(far_catch_brm$fit)
neff_lowest(far_catch_brm$fit)
rhat_highest(far_catch_brm$fit)
summary(far_catch_brm)
bayes_R2(far_catch_brm)

conditional_effects(far_catch_brm)

y <- data$log_catch_stnd
yrep_far_catch_brm  <- fitted(far_catch_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_far_catch_brm[sample(nrow(yrep_far_catch_brm), 25), ]) +
  xlim(-6, 6) 

trace_plot(far_catch_brm$fit)


## plot far and sst brms models -----------------------

## first, sst
## 95% CI
ce1s_1 <- conditional_effects(sst_catch_brm, effect = "annual_sst_3", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(sst_catch_brm, effect = "annual_sst_3", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(sst_catch_brm, effect = "annual_sst_3", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$annual_sst_3
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$annual_sst_3[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$annual_sst_3[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$annual_sst_3[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$annual_sst_3[["lower__"]]

sst.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Sea surface temperature anomaly", y = "Log catch anomaly", tag = "A") +
  geom_text(data=data, aes(annual_sst_3, log_catch_stnd, label = year), size=2.5) + ## TODO is this right?
  theme_bw()

sst.plot

## plot far
## 95% CI
ce1s_1 <- conditional_effects(far_catch_brm, effect = "annual_far_3", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(far_catch_brm, effect = "annual_far_3", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(far_catch_brm, effect = "annual_far_3", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$annual_far_3
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$annual_far_3[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$annual_far_3[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$annual_far_3[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$annual_far_3[["lower__"]]

far.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Fraction of attributable risk", y = "Log catch anomaly", tag = "B") +
  geom_text(data=data, aes(annual_far_3, log_catch_stnd, label = year), size=2.5) + ## TODO is this right?
  theme_bw()

far.plot


plot.dat <- data %>%
  select(year, annual_far_3, log_catch) %>%
  rename('a) Ocean entry FAR' = annual_far_3,
         'b) Log(catch)' = log_catch) %>%
  pivot_longer(cols = -year) %>%
  mutate(order = if_else(name == "b) Log(catch)", 2, 1),
         name = reorder(name, order))

far_catch <- ggplot(plot.dat, aes(year, value)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~name, scale = "free_y", ncol = 1) +
  theme(axis.title.x = element_blank()) +
  ylab("Value") 

far_catch

# separate plots for FAR and catch

far_plot <- ggplot(data, aes(year, annual_far_3)) +
  geom_point(size = 2) +
  geom_line() +  
  theme(axis.title.x = element_blank()) +
  labs(y = "Ocean entry FAR", tag = "A")

far_plot

# alternate plot for stakeholder talk
alt_plot <- ggplot(data, aes(year, annual_far_3)) +
  geom_point(size = 2) +
  geom_line() +  
  theme(axis.title.x = element_blank()) +
  labs(y = "% human risk") +
  scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1.0),
                     labels = c("40%", "60%", "80%", "100%"))

alt_plot

ggsave("./CMIP6/figs/alternate_FAR_plot_stakeholders.png", width = 6, height = 4, units = 'in')

catch_plot <- ggplot(data, aes(year, log_catch)) +
  geom_point(size = 2) +
  geom_line() +  
  theme(axis.title.x = element_blank()) +
  labs(y = "Ln(sockeye catch)", tag = "B")

catch_plot

# alternate for stakeholders
alt_catch_plot <- ggplot(data, aes(year, log_catch)) +
  geom_point(size = 2) +
  geom_line() +  
  theme(axis.title.x = element_blank()) +
  labs(y = "Millions of fish") +
  scale_y_continuous(breaks = c(15.5, 16, 16.5),
                     labels = c(round(exp(15.5)/1e6, 1), 
                                round(exp(16)/1e6, 1),
                                round(exp(16.5)/1e6, 1)))

alt_catch_plot

# save combined
png("./CMIP6/figs/combined_stakeholder_sockeye_plot.png", width = 4, height = 4, units = 'in', res = 300)

ggpubr::ggarrange(alt_plot, alt_catch_plot, ncol = 1)

dev.off()


## brms: setup ---------------------------------------------

sockeye_form <- bf(log_catch_stnd ~ 1 + s(annual_far_3) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term --------------------------------------
sockeye_model1 <- brm(sockeye_form,
                   data = data,
                   prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(sockeye_model1, file = "./CMIP6/brms_output/sockeye_far.rds")

sockeye_model1 <- readRDS("./CMIP6/brms_output/sockeye_far.rds")

check_hmc_diagnostics(sockeye_model1$fit)

neff_lowest(sockeye_model1$fit)

rhat_highest(sockeye_model1$fit)

summary(sockeye_model1)

bayes_R2(sockeye_model1)

plot(conditional_smooths(sockeye_model1), ask = FALSE)

plot(sockeye_model1, ask = FALSE)

y <- data$log_catch_stnd

yrep_sockeye_model1  <- fitted(sockeye_model1, scale = "response", summary = FALSE)

ppc_dens_overlay(y = y, yrep = yrep_sockeye_model1[sample(nrow(yrep_sockeye_model1), 25), ]) +
  ggtitle("sockeye_model1")

trace_plot(sockeye_model1$fit)

## refit with cor, not cov in ar() ---------------------

sockeye_form2 <- bf(log_catch_stnd ~ 1 + s(annual_far_3) + ar(time = year, p = 1))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit --------------------------------------
sockeye_model2 <- brm(sockeye_form2,
                     data = data,
                     prior = priors,
                     cores = 4, chains = 4, iter = 3000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(sockeye_model2, file = "./CMIP6/brms_output/sockeye_far2.rds")

sockeye_model2 <- readRDS("./CMIP6/brms_output/sockeye_far2.rds")

check_hmc_diagnostics(sockeye_model2$fit)

neff_lowest(sockeye_model2$fit)

rhat_highest(sockeye_model2$fit)

summary(sockeye_model2)

bayes_R2(sockeye_model2)

plot(conditional_smooths(sockeye_model2), ask = FALSE)

plot(sockeye_model2, ask = FALSE)

y <- data$log_catch_stnd

yrep_sockeye_model2  <- fitted(sockeye_model2, scale = "response", summary = FALSE)

ppc_dens_overlay(y = y, yrep = yrep_sockeye_model2[sample(nrow(yrep_sockeye_model2), 25), ]) +
  ggtitle("sockeye_model2") # this looks better than model1

trace_plot(sockeye_model2$fit)

## plot model2 ------------------------

## 95% CI
ce1s_1 <- conditional_effects(sockeye_model2, effect = "annual_far_3", re_formula = NA,
                              prob = 0.95)
## 90% CI
ce1s_2 <- conditional_effects(sockeye_model2, effect = "annual_far_3", re_formula = NA,
                              prob = 0.9)
## 80% CI
ce1s_3 <- conditional_effects(sockeye_model2, effect = "annual_far_3", re_formula = NA,
                              prob = 0.8)
dat_ce <- ce1s_1$annual_far_3

#########################

dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$annual_far_3[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$annual_far_3[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$annual_far_3[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$annual_far_3[["lower__"]]
# dat_ce[["rug.anom"]] <- c(jitter(unique(data$annual_far_3), amount = 0.0051),
                          # rep(NA, 100-length(unique(data$annual_far_3))))


g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Fraction of Attributable Risk", y = "log(catch)") +
  theme_bw() +
  geom_text(data=data, aes(annual_far_3, log_catch_stnd, label = year), size=2.5) +
  geom_hline(yintercept = 0, lty = 2) 
  


+
  geom_rug(aes(x=rug.anom, y=NULL)) 

print(g2)

## fit third model with categorical FAR -----------------------------


# add categorical variable - is FAR above or below 0.98?

# examine FAR distribution
ggplot(data, aes(annual_far_3)) +
  geom_histogram(fill = "grey", color = "black", bins = 14)

far_check <- data %>%
  arrange(desc(annual_far_3))

far_check # supports 0.98 as a cutoff


  data <- data %>%
  mutate(far_fac = as.factor(if_else(annual_far_3 > 0.98, "above", "below")))

sockeye_form3 <- bf(log_catch_stnd ~ 1 + far_fac + ar(time = year, p = 1))


get_prior(sockeye_form3, data = data)

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit
sockeye_model3 <- brm(sockeye_form3,
                      data = data,
                      prior = priors,
                      cores = 4, chains = 4, iter = 2000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(sockeye_model3, file = "./CMIP6/brms_output/sockeye_far3.rds")

sockeye_model3 <- readRDS("./CMIP6/brms_output/sockeye_far3.rds")

check_hmc_diagnostics(sockeye_model3$fit)

neff_lowest(sockeye_model3$fit)

rhat_highest(sockeye_model3$fit)

summary(sockeye_model3)

bayes_R2(sockeye_model3)

plot(sockeye_model3, ask = FALSE)

y <- data$log_catch_stnd

yrep_sockeye_model3  <- fitted(sockeye_model3, scale = "response", summary = FALSE)

ppc_dens_overlay(y = y, yrep = yrep_sockeye_model3[sample(nrow(yrep_sockeye_model3), 25), ]) +
  ggtitle("sockeye_model3") # this looks better than model1

trace_plot(sockeye_model3$fit)


## 95% CI
ce1s_1 <- conditional_effects(sockeye_model3, effect = "far_fac", re_formula = NA,
                              prob = 0.95)  

plot <- ce1s_1$far_fac %>%
  select(far_fac, estimate__, lower__, upper__)

plot$far_fac <- reorder(plot$far_fac, desc(plot$far_fac))

g3 <- ggplot(plot, aes(far_fac, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Ln(catch anomaly)") +
  xlab("Fraction of Attibutable Risk (FAR)") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_x_discrete(labels=c(expression(""<=0.91), expression("">=0.98))) +
  theme_bw() + 
  labs(tag = "C")

print(g3)

ggsave("./CMIP6/figs/FAR_sockeye_categorical.png", width=1.5, height=2, units='in')


# plot posterior distributions

posteriors <- as.data.frame(posterior_epred(sockeye_model3))

names(posteriors) <- data$far_fac

posteriors <- posteriors %>%
  pivot_longer(everything(), names_to = "FAR", values_to = "Catch")

  
categorical.plot <- ggplot(posteriors, aes(Catch, fill = FAR)) +
  geom_density(color = NA,  alpha = 0.7) +
  scale_fill_manual(values = cb[c(8,4)], labels = c("\u2265 0.98", "< 0.91")) +
  labs(x = "Log catch anomaly", y = "Density", tag = "C") +
  theme(legend.position = c(0.2, 0.8)) +
  geom_hline(yintercept  = 0)

categorical.plot


# alternate plot for stakeholders
mu <- 16.19
sd <- 0.345


alt_g3 <- ggplot(plot, aes(far_fac, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Millions of fish") +
  xlab("% of human risk") +
  scale_y_continuous(breaks = c(0, -1, -2),
                     labels = c(round(exp(mu)/1e6, 1),
                                round(exp(mu-sd)/1e6, 1),
                                round(exp(mu-2*sd)/1e6, 1))) +
scale_x_discrete(labels=c("< 91%", "> 98%")) +
  theme_bw() 

print(alt_g3)

ggsave("./CMIP6/figs/alt_FAR_sockeye_categorical.png", width=3, height=2.5, units='in')


## expected return time ----------------------------------------





# load ersst anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

unique(ersst.anom$region)

ersst.anom <- ersst.anom %>%
  dplyr::filter(region == "Gulf_of_Alaska")

# plot distributions to check
ggplot(filter(ersst.anom, year %in% 1950:1999), aes(annual.anomaly.three.yr.running.mean)) +
  geom_density(fill = "grey")


# sst anomaly threshold predicted to correspond with highest 5 FAR values (2017 value is threshold) = 2.075
ersst.max <- 2.075

# load CMIP6 anomalies
cmip.anom <- read.csv("./CMIP6/summaries/CMIP6.anomaly.time.series.csv")

# load estimated warming level timing for each model
timing <- read.csv("./CMIP6/summaries/model.north.pacific.warming.timing.csv")

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
           model == models[i],
           region == "Gulf_of_Alaska")
  
  
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "preindustrial",
                                       count = sum(pre.temp$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(pre.temp$annual.three.yr.running.mean))))
  
  
  
  
} # close i loop (models)


head(extreme.outcomes) 

## record outcomes using different warming levels from hist.585 ----------------


# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  hist.temp <- cmip.anom %>% 
    filter(experiment == "hist_ssp585",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  # # separate this region from ersst.max
  # ersst.temp <- ersst.max %>%
  #   filter(region == regions[j])
  
  ## pull 1950 - 0.5 degrees warming
  
  use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "1950_to_0.5",
                                       count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))) )
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "0.5_to_1.0",
                                       count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))))
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "1.0_to_1.5",
                                       count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record how many model years are more extreme
  extreme.outcomes <- rbind(extreme.outcomes,
                            data.frame(model = models[i],
                                       period = "1.5_to_2.0",
                                       count = sum(hist.temp.use$annual.three.yr.running.mean >= ersst.max, na.rm = T),
                                       N = length(na.omit(hist.temp.use$annual.three.yr.running.mean))))
  
  
  
} # close i loop (models)


# check
check <- extreme.outcomes %>%
  group_by(period) %>%
  summarise(count = sum(count),
            N = sum(N)) %>%
  mutate(prop = count/N)

View(check)




### plot hindcast and projected pdfs for sst anomalies --------------------------
# create df to catch outcomes for extreme runs
anomaly.pdfs <- data.frame()

# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  pre.temp <- cmip.anom %>% 
    filter(experiment == "piControl",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                            data.frame(model = models[i],
                                       period = "preindustrial",
                                       anomaly = na.omit(pre.temp$annual.three.yr.running.mean)))
  
  
  
} # close i loop (models)


## record outcomes using different warming levels from hist.585 ----------------


# loop through each model
for(i in 1:length(models)){ # start i loop (models)
  # i <- 1
  
  # separate model and region of interest
  hist.temp <- cmip.anom %>% 
    filter(experiment == "hist_ssp585",
           model == models[i],
           region == "Gulf_of_Alaska")
  
  # # separate this region from ersst.max
  # ersst.temp <- ersst.max %>%
  #   filter(region == regions[j])
  
  ## pull 1950 - 0.5 degrees warming
  
  use = 1950:timing$year[timing$model == models[i] & timing$level == 0.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1950_to_0.5",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  ## pull 0.5 - 1.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 0.5]:timing$year[timing$model == models[i] & timing$level == 1.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "0.5_to_1.0",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  ## pull 1.0 - 1.5 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.0]:timing$year[timing$model == models[i] & timing$level == 1.5]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.0_to_1.5",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  ## pull 1.5 - 2.0 degrees warming
  
  use = timing$year[timing$model == models[i] & timing$level == 1.5]:timing$year[timing$model == models[i] & timing$level == 2.0]
  
  # and limit hist.temp to these years
  hist.temp.use <- hist.temp %>%
    filter(year %in% use)
  
  # record anomalies
  anomaly.pdfs <- rbind(anomaly.pdfs,
                        data.frame(model = models[i],
                                   period = "1.5_to_2.0",
                                   anomaly = na.omit(hist.temp.use$annual.three.yr.running.mean)))
  
  
  
} # close i loop (models)

# load CMIP6 model weights
model.weights <- read.csv("./CMIP6/summaries/normalized_CMIP6_weights.csv") 

# clean up model weights 
model.weights <- model.weights %>%
  filter(region == "Gulf_of_Alaska") %>%
  select(model, normalized_weight) 


# join anomalies with model weights
anomaly.pdfs <- left_join(anomaly.pdfs, model.weights) 


# resample to weight models
resample.pdf <- data.frame()

periods <- unique(anomaly.pdfs$period)

for(i in 1:length(periods)){
  # i <- 1

  temp <- anomaly.pdfs[anomaly.pdfs$period == periods[i],]
  set.seed(999)
  resample.pdf <- rbind(resample.pdf,
                        data.frame(period = periods[i],
                                   anomaly = sample(temp$anomaly, 500, replace = T, prob = temp$total_weight)))
  
  
  
}

# reorder
plot.order <- data.frame(period = unique(resample.pdf$period),
                         order = 1:5)



resample.pdf <- left_join(resample.pdf, plot.order) %>%
  mutate(period =  reorder(period, order))

# and plot

# get good labels for each warming period
labs <- data.frame(period = unique(resample.pdf$period),
                   plot_period = c("Preindustrial",
                                   "1950 to 0.5°",
                                   "0.5° to 1.0°",
                                   "1.0° to 1.5°",
                                   "1.5° to 2.0°"),
                   plot_order = 1:5)

resample.pdf <- left_join(resample.pdf, labs) %>%
  mutate(plot_period = reorder(plot_period, plot_order))


cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf_plot <- ggplot(resample.pdf, aes(plot_period, anomaly)) +
  geom_violin(fill = cb[6], lty = 0, alpha = 0.5) +
  coord_flip() +
  xlab("North Pacific warming") +
  ylab("SST anomaly (Std. Dev.)") +
  geom_hline(yintercept = ersst.max, lty = 2) +
  labs(tag = "D")

pdf_plot

# and get summary statistics (proportion above ersst.max threshold)

summary_pdf <- resample.pdf %>%
  group_by(period) %>%
  summarise(proportion = sum(anomaly >= ersst.max) / length(anomaly))

summary_pdf


ggsave("./CMIP6/figs/sockeye_anomaly_pdfs.png", width = 3, height = 3, units = 'in')        

# combine plots

void_plot <- ggplot() + theme_void()

png("./CMIP6/figs/combined_sockeye_FAR_plot.png", width = 9, height = 7, units = 'in', res = 300)

ggpubr::ggarrange(ggpubr::ggarrange(far_plot,  catch_plot, ggpubr::ggarrange(void_plot, g3, ncol = 2, widths = c(0.3, 0.7)), ncol = 1),
pdf_plot, ncol = 2, widths = c(0.45, 0.55))

dev.off()

