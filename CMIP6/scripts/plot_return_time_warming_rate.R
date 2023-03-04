# plot combined predicted warming rate ssp245, ssp585, ersst
# and combine with expected extreme return times

library(tidyverse)
library(brms)

# set palette and theme

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_set(theme_bw())

# load model objects

warming.245 <- readRDS("./CMIP6/brms_output/inverse_warming_brm_ssp245.rds")

warming.585 <- readRDS("./CMIP6/brms_output/inverse_warming_brm.rds")

warming.ersst<- readRDS("./CMIP6/brms_output/inverse_warming_brm_ersst.rds")


# plot

new.dat <- data.frame(warming = c(0.5, 1.0, 1.5, 2.0),
                      model_fac = NA, weight = 1)

pred.245 <- posterior_epred(warming.245, newdata = new.dat)

pred.585 <- posterior_epred(warming.585, newdata = new.dat)

new.dat <- data.frame(ersst.warming = c(0.5, 1.0))

pred.ersst <- posterior_epred(warming.ersst, newdata = new.dat) 


## SST anomaly predictions #### 95% CI

# 245

ce1s_245 <- conditional_effects(warming.245, effect = "warming", re_formula = NA,
                              probs = c(0.025, 0.975), resolution = 10000)

index <- ce1s_245$warming$warming

choose <- c(which.min(abs(index - 0.5)),
            which.min(abs(index - 1.0)),
            which.min(abs(index - 1.5)),
            which.min(abs(index - 2.0)))

pred.plot.245 <- data.frame(source = "SSP 245",
                            warming = c(0.5, 1.0, 1.5, 2.0),
                        year = ce1s_245$warming$estimate__[choose],
                        UCI = ce1s_245$warming$upper__[choose],
                        LCI = ce1s_245$warming$lower__[choose])


# 585 

ce1s_585 <- conditional_effects(warming.585, effect = "warming", re_formula = NA,
                                probs = c(0.025, 0.975), resolution = 10000)

index <- ce1s_585$warming$warming

choose <- c(which.min(abs(index - 0.5)),
            which.min(abs(index - 1.0)),
            which.min(abs(index - 1.5)),
            which.min(abs(index - 2.0)))

pred.plot.585 <- data.frame(source = "SSP 585",
                            warming = c(0.5, 1.0, 1.5, 2.0),
                            year = ce1s_585$warming$estimate__[choose],
                            UCI = ce1s_585$warming$upper__[choose],
                            LCI = ce1s_585$warming$lower__[choose])


# ersst
ce1s_ERSST <- conditional_effects(warming.ersst, effect = "ersst.warming", re_formula = NA,
                                probs = c(0.025, 0.975), resolution = 10000)

index <- ce1s_ERSST$ersst.warming$ersst.warming

choose <- c(which.min(abs(index - 0.5)),
            which.min(abs(index - 1.0)))

pred.plot.ERSST <- data.frame(source = "ERSST",
                            warming = c(0.5, 1.0),
                            year = ce1s_ERSST$ersst.warming$estimate__[choose],
                            UCI = ce1s_ERSST$ersst.warming$upper__[choose],
                            LCI = ce1s_ERSST$ersst.warming$lower__[choose])

pred.plot <- rbind(pred.plot.245,
                   pred.plot.585,
                   pred.plot.ERSST)
pred.plot$warming <- as.factor(pred.plot$warming)

pos_dodge = position_dodge(width = 0.2)



time.plot <- ggplot(pred.plot, aes(warming, year, color = source)) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2, position = pos_dodge) +
  geom_point(size = 4, position = pos_dodge) +
  labs(x = "North Pacific warming (°C)",
       y = "Year reached") +
  scale_color_manual(values = cb[c(2,6,7)]) +
  theme(legend.position = c(0.75, 0.15),
        legend.title = element_blank())

# get extremes.plot from CMIP6_Bayes_more_extreme_probabilities.R

# and combine

# get blank plot
blank <- ggplot + theme_void()


png("./CMIP6/figs/extremes_probability_warming_timing.png", width = 10, height = 5, units = 'in', res = 300) 

ggpubr::ggarrange(extremes.plot, 
                    ggpubr::ggarrange(time.plot,
                                      blank,
                                      nrow = 2,
                                      heights = c(0.8, 0.2)),
                  ncol = 2,
                  widths = c(0.7, 0.3),
                  labels = c("a", "b"))

dev.off()
