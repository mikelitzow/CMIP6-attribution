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



npgo <- read.csv("./CMIP6/data/npgo.csv")

npgo <- npgo %>%
  filter(year %in% 1950:2021) %>% # 2022 not complete as of 3/6/23, use 1950-2021 only
  group_by(year) %>%
  summarize(npgo = mean(npgo))


dat <- left_join(pdo, npgo)

## load far values----------------

# load region names
regions <- read.csv("./CMIP6/summaries/clean_region_names.csv")
regions <- regions[1:6,1]

# load ERSST anomalies
ersst.anom <- read.csv("./CMIP6/summaries/regional_north_pacific_ersst_anomaly_time_series.csv")

# limit to 1950:2022
ersst.anom <- ersst.anom %>%
  filter(year %in% 1950:2022)

far_pred_annual <- data.frame()

for(i in 2:length(regions)){ # loop through regions, excluding N. Pacific
  
  # i <- 1
  
  # subset ersst.anom
  ersst.temp <- ersst.anom %>%
    filter(region == regions[i]) %>%
    select(year, annual.anomaly.unsmoothed) %>%
    rename(annual.anomaly.1yr = annual.anomaly.unsmoothed)
  
  # load regional model
  mod <- readRDS(paste("./CMIP6/brms_output/", regions[i], "annual_sst_rolling_window_binomial2.rds", sep = ""))
  
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

for(i in 1974:2022){
  # i <- 1964
  
  temp <- dat %>%
    filter(year %in% (i-14):i)
  
  for(j in 1:length(regions)){
   # j <- 2
  
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

f.p <- function(x) cor.test(x,1974:2022, method="kendall", na.action = "na.omit")$p.value
f.tau <- function(x) cor.test(x,1974:2022, method="kendall", na.action = "na.omit")$statistic


test <- plot.out %>%
  group_by(region) %>%
  summarise(PDO_p = f.p(PDO),
            PDO_tau = f.tau(PDO),
            NPGO_p = f.p(NPGO),
            NPGO_tau = f.tau(NPGO))

test  # all of them are getting stronger except GOA!

# save test output
write.csv(test, "./CMIP6/summaries/Kendalls_tau_PDO_NPGO_FAR.csv", row.names = F)

# and plot 

plot.out <- plot.out %>%
  pivot_longer(cols = c(-region, -year), values_to = "Correlation")

plot.order <- data.frame(name = c("PDO", "NPGO"),
                         order = c(1,2))

plot.out <- left_join(plot.out, plot.order)

plot.out$name <- reorder(plot.out$name, plot.out$order)

# clean region names
region.order <- data.frame(region = regions,
                           region_plot = str_replace_all(regions, "_", " "),
                           region_order = 1:5)

plot.out <- left_join(plot.out, region.order)

plot.out$region_plot <- reorder(plot.out$region_plot, plot.out$region_order)

ggplot(plot.out, aes(year, Correlation, color = region_plot)) +
  geom_line() +
  facet_wrap(~name) +
  scale_color_manual(values = cb[c(2:6)]) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  geom_hline(lty = 2, yintercept = 0)

# and save 
ggsave("./CMIP6/figs/internal_variability_far_correlations.png", width = 10, height = 3.5, units = 'in')
