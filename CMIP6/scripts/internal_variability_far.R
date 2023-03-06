# time-varying attribution of PDO/NPGO events

library(tidyverse)
library(brms)

theme_set(theme_bw())

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load data

pdo <- read.csv("./CMIP6/data/pdo.csv")

pdo <- pdo %>%
  pivot_longer(cols = -Year) %>%
  rename(year = Year) %>%
  group_by(year) %>%
  filter(year %in% 1950:2022) %>%
  summarise(pdo = mean(value))



  mutate(month = as.numeric(months(date)),
         year = as.numeric(as.character(chron::years(date)))) %>%
  select(-date)


# pdo$year <- as.numeric(pdo$year)

npgo <- read.csv("./CMIP6/data/npgo.csv")

str(pdo)

dat <- left_join(pdo, npgo) %>%
  filter(year %in% 1950:2021) 


plot.dat <- dat %>%
  pivot_longer(cols = c(-month, -year)) %>%
  mutate(dec.yr = year + (month - 0.5)/12)


ggplot(plot.dat, aes(dec.yr, value, color = name)) +
  geom_line()

# # group by winter year
# dat <- dat %>%
#   filter(month %in% c(11,12,1:3)) %>%
#   mutate(winter.year = if_else(month %in% 11:12, year + 1, year)) %>%
#   filter(winter.year %in% 1951:2021) %>%
#   group_by(winter.year) %>%
#   summarise(pdo = mean(pdo),
#             npgo = mean(npgo))
# 
# plot.dat <- dat %>%
#   pivot_longer(cols = c(-winter.year)) 
# 
# ggplot(plot.dat, aes(winter.year, value, color = name)) +
#   geom_line()

# group by year

dat <- dat %>% 
  group_by(year) %>%
  summarise(pdo = mean(pdo),
            npgo = mean(npgo))

## load far values----------------

# load region names
regions <- read.csv("./CMIP6/summaries/clean_region_names.csv")
regions <- regions[1:6,1]

# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# limit to 1950:2021
ersst.anom <- ersst.anom %>%
  filter(year %in% 1950:2021)

far_pred_annual <- data.frame()

for(i in 2:length(regions)){ # loop through regions, excluding N. Pacific
  
  # i <- 1
  
  # subset ersst.anom
  ersst.temp <- ersst.anom %>%
    filter(region == regions[i]) %>%
    select(year, annual.anomaly.unsmoothed) %>%
    rename(annual.anomaly.1yr = annual.anomaly.unsmoothed)
  
  # load regional model
  mod <- readRDS(paste("./CMIP6/brms_output/", regions[i], "_rolling_window_binomial2.rds", sep = ""))
  
  ## setup new data
  nd <- data.frame(period = c("historical", "preindustrial"),
                   ersst.year = rep(ersst.temp$year, each = 2),
                   annual.anomaly.1yr = rep(ersst.temp$annual.anomaly.1yr, each = 2),
                   N = 1000,
                   model_fac = NA)
  
  nd_pre <- nd[nd$period == "preindustrial", ]
  nd_his <- nd[nd$period == "historical", ]
  
  ## make predictions
  ## exclude random effects for model_fac
  pre_pp <- posterior_epred(mod, newdata = nd_pre, re_formula = NA)
  his_pp <- posterior_epred(mod, newdata = nd_his, re_formula = NA)
  
  ## Calc probabilities
  ## These are our posterior probabilities to use for FAR calculation
  pre_prob <- pre_pp / unique(nd$N)
  his_prob <- his_pp / unique(nd$N)
  
  
  ## Calc FAR
  far <- 1 - (pre_prob / his_prob)
  range(far, na.rm = TRUE)
  
  
  far_pred_annual <- rbind(far_pred_annual,
                           data.frame(region = regions[i],
                                      year = nd_pre$ersst.year,
                                      prob = apply(far, 2, mean),
                                      lower = apply(far, 2, quantile, probs = 0.025),
                                      upper = apply(far, 2, quantile, probs = 0.975)))
  
  
} # close i loop

dat <- left_join(far_pred_annual, dat)

## loop through 15 year window and calculate correlations and regression coefficients for far-pdo and far-npgo

regions <- unique(dat$region)

output <- data.frame()

for(i in 1969:2021){
  # i <- 1964
  
  temp <- dat %>%
    filter(year %in% (i-15):i)
  
  for(j in 1:length(regions)){
   # j <- 1 
  
   temp.region <- temp %>%
    filter(region == regions[j])
  
   mod1 <- lm(prob ~ pdo, data = temp.region)  
   mod2 <- lm(prob ~ npgo, data = temp.region)
   
   output <- rbind(output,
                   data.frame(region = regions[j],
                              year = i,
                              pdo.coeff = coefficients(mod1)[2],
                              npgo.coeff = coefficients(mod2)[2],
                              pdo.cor = cor(temp.region$pdo, temp.region$prob),
                              npgo.cor = cor(temp.region$npgo, temp.region$prob)))
    
  }
}

plot.out <- output %>%
  pivot_longer(cols = c(-region, -year))


ggplot(plot.out, aes(year, value, color = region)) +
  geom_line() +
  facet_wrap(~name, scale = "free_y") +
  scale_color_manual(values = cb[c(2:6)])

# get Kendall's tau for each and plot only correlation time series

plot.out <- output %>%
  select(-pdo.coeff, -npgo.coeff) %>%
  rename(PDO = pdo.cor,
         NPGO = npgo.cor)

f.p <- function(x) cor.test(x,1969:2021, method="kendall")$p.value
f.tau <- function(x) cor.test(x,1969:2021, method="kendall")$statistic


test <- plot.out %>%
  group_by(region) %>%
  summarise(PDO_p = f.p(PDO),
            PDO_tau = f.tau(PDO),
            NPGO_p = f.p(NPGO),
            NPGO_tau = f.tau(NPGO))

test  # all of them are getting stonger!

# save test output
write.csv(test, "./CMIP6/summaries/Kendalls_tau_PDO_NPGO_FAR.csv", row.names = F)

# and plot 

plot.out <- plot.out %>%
  pivot_longer(cols = c(-region, -year), values_to = "Correlation")

plot.order <- data.frame(name = c("PDO", "NPGO"),
                         order = c(1,2))

plot.out <- left_join(plot.out, plot.order)

plot.out$name <- reorder(plot.out$name, plot.out$order)

ggplot(plot.out, aes(year, Correlation, color = region)) +
  geom_line() +
  facet_wrap(~name, scale = "free_y") +
  scale_color_manual(values = cb[c(2:6)]) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  geom_hline(lty = 2, yintercept = 0)

# and save 
ggsave("./CMIP6/figs/internal_variability_far_correlations.png", width = 10, height = 3.5, units = 'in')

## EBS and GOA - RR for 76/77 PDO event-----------------------

# load RR 
RR <- read.csv("./CMIP6/summaries/complete_FAR_RR_time_series_with_uncertainty.csv")

# clean up 
dat <- RR %>%
  filter(region %in% c("Eastern_Bering_Sea", "Gulf_of_Alaska"),
         window_plot == 'Annual SST', 
         year %in% 1967:1986) %>%
  select(region, year, FAR, RR)

pdo <- read.csv("./CMIP6/data/pdo.csv")

pdo <- pdo %>%
  mutate(year = as.numeric(as.character(chron::years(date)))) %>%
  group_by(year) %>%
  summarize(PDO = mean(pdo))

dat <- left_join(dat, pdo)

ggplot(dat, aes(PDO, FAR)) +
  geom_point() +
  geom_smooth(method = "gam", se = F) +
  facet_wrap(~region, scales = "free_y", ncol = 1)

ggplot(dat, aes(PDO, RR)) +
  geom_point() +
  geom_smooth(method = "gam", se = F) +
  facet_wrap(~region, scales = "free_y", ncol = 1)
