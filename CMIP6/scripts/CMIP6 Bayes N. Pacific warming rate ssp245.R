## estimate projected warming rate for the N. Pacific from CMIP6
## weighting models by their performance wrt 1972-2021 observed warming rate
## for ssp245

library(tidyverse)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./CMIP6/scripts/stan_utils.R")

theme_set(theme_bw())

## load model sst and create warming time series wrt 1850-1949 from historical/ssp245 runs --------
sst.245 <- read.csv("./CMIP6/summaries/CMIP6.sst.time.series.ssp245.csv")

# rearrange
sst.245 <- sst.245 %>%
  filter(region == "North_Pacific") %>%
  select(year, model, annual.unsmoothed) %>%
  pivot_wider(names_from = model, values_from = annual.unsmoothed)

model.warming.rate <- sst.245

# now get anomaly wrt 1850:1949
for(j in 2:ncol(sst.245)){
  
  model.warming.rate[,j] <- sst.245[,j] - colMeans(sst.245[1:100,])[j]
  
}

model.warming.rate <- model.warming.rate %>%
  pivot_longer(cols = -year)

# load model weights
model.weights <- read.csv("./CMIP6/summaries/N_Pac_warming_model_weights_ssp245.csv")

# wrangle data

# simplify weights
weights <- model.weights  %>%
  mutate(model_fac = as.factor(model)) %>%
  select(model_fac, weight)

dat <- model.warming.rate %>%
  mutate(model_fac = as.factor(name)) %>%
  select(-name) %>%
  rename(warming = value) 

levels(weights$model_fac); levels(dat$model_fac)

dat <- left_join(dat, weights)

# check for NA
sum(is.na(dat))

## run brms ---------------------------------

warming_formula <-  bf(warming | weights(weight) ~ s(year) + (1 | model_fac))

## Show default priors
get_prior(warming_formula, dat)

warming_brm <- brm(warming_formula,
                    data = dat,
                    cores = 4, chains = 4, iter = 5000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 16))

saveRDS(warming_brm, file = "./CMIP6/brms_output/warming_brm.rds")

warming_brm <- readRDS("./CMIP6/brms_output/warming_brm.rds")

check_hmc_diagnostics(warming_brm$fit)
neff_lowest(warming_brm$fit) # ??
rhat_highest(warming_brm$fit)
summary(warming_brm)
bayes_R2(warming_brm)
plot(warming_brm$criteria$loo, "k")
plot(conditional_effects(warming_brm), ask = FALSE)
y <- dat$warming
yrep_warming_brm  <- fitted(warming_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_warming_brm[sample(nrow(yrep_warming_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("warming_brm")
pdf("./figs/trace_warming_brm.pdf", width = 6, height = 4)
trace_plot(warming_brm$fit)
dev.off()


## SST anomaly predictions #### 95% CI
ce1s_1 <- conditional_effects(warming_brm, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975))

dat_ce <- ce1s_1$year %>%
  select(year, estimate__, lower__, upper__) %>%
  mutate(Source = "Modeled warming")

# add observations
ersst <- read.csv("./CMIP6/summaries/North_Pacific_ersst_1854-2021.csv")
ersst$warming.wrt.1854.1949 <- ersst$annual.unsmoothed - mean(ersst$annual.unsmoothed[ersst$year %in% 1854:1949])

obs <- ersst %>%
  select(year, warming.wrt.1854.1949) %>%
  rename(estimate__ = warming.wrt.1854.1949) %>%
  mutate(Source = "Actual temperature",
         lower__ = NA,
         upper__ = NA)

dat_ce <- rbind(dat_ce, obs)

dat_ce$order <- if_else(dat_ce$Source == "Actual temperature", 2, 1)

dat_ce$Source <- reorder(dat_ce$Source, dat_ce$order)

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(dat_ce) +
  aes(x = year, y = estimate__, color = Source, fill = Source) +
  geom_line() +
  geom_ribbon( 
              aes(x = year, ymin = lower__, ymax = upper__, fill = Source), alpha = 0.15, lty = 0) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  scale_color_manual(values = cb[c(2,4)]) +
  scale_fill_manual(values = cb[c(2,4)])


ggsave("./CMIP6/figs/modeled_observed_warming.png", width = 5, height = 3)


# get predicted warming for each year
new.dat <- data.frame(year = 1940:2100,
                      model = NA)

pred <- posterior_epred(warming_brm, newdata = new.dat, re_formula = NA)

new.dat$pred_mean <- colMeans(pred) 
new.dat$lower = apply(pred, 2, quantile, probs = 0.025)
new.dat$upper = apply(pred, 2, quantile, probs = 0.975)

ggplot(new.dat, aes(year, pred_mean)) +
  geom_line()

View(new.dat)

# save 
write.csv(new.dat, "./CMIP6/summaries/brms_predicted_North_Pac_warming.csv", row.names = F)

## warming timings under ssp585
# 2003 is the best timing for 0.5 degrees warming
# 2019 is the best timing for 1.0 degrees warming
# 2021 warming is c. 1.09 degrees
# 2033 warming is 1.5 degrees
# 2045 warming is 2.0 degrees

## fit inverse model - year as a function of warming -----------------------------

inverse_formula <-  bf(year | weights(weight) ~ s(warming) + (1 | model_fac))

## Show default priors
get_prior(inverse_formula, dat)

inverse_warming_brm <- brm(inverse_formula,
                           data = filter(dat, year >= 1972), # limit to 1972-on to ease fitting 
                           cores = 4, chains = 4, iter = 5000,
                           save_pars = save_pars(all = TRUE),
                           control = list(adapt_delta = 0.99, max_treedepth = 16))

saveRDS(inverse_warming_brm, file = "./CMIP6/brms_output/inverse_warming_brm_ssp245.rds")

inverse_warming_brm <- readRDS("./CMIP6/brms_output/inverse_warming_brm_ssp245.rds")

check_hmc_diagnostics(inverse_warming_brm$fit)
neff_lowest(inverse_warming_brm$fit) # ??
rhat_highest(inverse_warming_brm$fit)
summary(inverse_warming_brm)
bayes_R2(inverse_warming_brm)
plot(inverse_warming_brm$criteria$loo, "k")
plot(conditional_effects(inverse_warming_brm), ask = FALSE)

trace_plot(inverse_warming_brm$fit)

# predict for 0.5, 1.0, 1.5, 2.0 degrees

new.dat <- data.frame(warming = c(0.5, 1.0, 1.5, 2.0),
                      model_fac = NA, weight = 1)


pred <- posterior_epred(inverse_warming_brm, newdata = new.dat)

uci <- function(x) quantile(x, 0.975)

pred.summary <- data.frame(warming = new.dat$warming,
                           year = colMeans(pred),
                           UCI = apply(pred, 2, uci))


## SST anomaly predictions #### 95% CI
ce1s_1 <- conditional_effects(inverse_warming_brm, effect = "warming", re_formula = NA,
                              probs = c(0.025, 0.975), resolution = 10000)

index <- ce1s_1$warming$warming

choose <- c(which.min(abs(index - 0.5)),
            which.min(abs(index - 1.0)),
            which.min(abs(index - 1.5)),
            which.min(abs(index - 2.0)))

pred.plot <- data.frame(warming = c(0.5, 1.0, 1.5, 2.0),
                        year = ce1s_1$warming$estimate__[choose],
                        UCI = ce1s_1$warming$upper__[choose],
                        LCI = ce1s_1$warming$lower__[choose])


ggplot(pred.plot, aes(warming, year)) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
  geom_point(color = "red", size = 4) +
  labs(x = "North Pacific warming (°C)",
       y = "Year reached")

ggsave("./CMIP6/figs/Bayes_estimated_warming_timing.png", width = 3, height = 4, units = 'in') 
