library(tidyverse)

theme_set(theme_bw())
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


dat <- read.csv("./data/pollock_assessment_age_structure.csv")

head(dat)

ff <- function(x) (x-mean(x)) / sd(x)


dat[,2:11] <- apply(log(dat[,2:11]+0.01), 2, ff)

colMeans(dat[,2:11])


dat_order <- data.frame(name = colnames(dat)[2:11],
                        order = 1:10)

dat <- dat %>%
  pivot_longer(cols = -year)

dat <- left_join(dat, dat_order)

dat$name <- reorder(dat$name, dat$order)

dat <- dat %>%
  mutate(year_class = year-order, # since order = age!
         color = if_else(year_class == 2012, cb[7], "black")) 

ggplot(dat, aes(year, value)) +
  geom_hline(yintercept = 0, color = "dark grey") +
  geom_line() +
  geom_point(color = dat$color) +
  facet_wrap(~name, ncol = 5) +
  geom_vline(xintercept = c(2014, 2019), lty = 2, color = cb[8]) +
  theme(axis.title.x = element_blank()) +
  ylab("Scaled log(abundance)")

ggsave("./figs/assessment_model_scaled_abundance_by_age_time_series.png", width = 12, height = 4, units = "in")


ggplot(filter(dat, order %in% 4:9), aes(year, value, color = as.factor(order))) +
  geom_hline(yintercept = 0, color = "dark grey") +
  geom_line() +
  geom_vline(xintercept = c(2014, 2019), lty = 2, color = cb[8]) +
  theme(axis.title.x = element_blank())

ggplot(dat, aes(as.factor(order), value)) + # remember, order is the same as age!
  geom_col(fill = "grey90", color="black", position = "dodge", width = 0.5) +
  facet_wrap(~year, ncol = 5) 

dat <- dat %>%
  mutate(year_class = year-order, # since order = age!
         fill = if_else(year_class == 2012, cb[7], "grey90")) 

ggplot(filter(dat, order %in% 4:9, year >= 1983), aes(as.factor(order), value)) + # remember, order is the same as age!
  geom_col(aes(fill = fill), position = "dodge", width = 0.5) +
  facet_wrap(~year, ncol = 5) 

# alternate version - acoustic data from survey ---------------------

dat <- read.csv("./data/acoustic_age_abundance.csv")

dat[,2:16] <- apply(log(dat[,2:16]+0.01), 2, ff)

colMeans(dat[,2:16])


dat_order <- data.frame(name = colnames(dat)[2:16],
                        order = 1:15)

dat <- dat %>%
  pivot_longer(cols = -year)

dat <- left_join(dat, dat_order)

dat$name <- reorder(dat$name, dat$order)

dat <- dat %>%
  mutate(year_class = year-order, # since order = age!
         color = if_else(year_class == 2012, cb[7], "black")) 

ggplot(dat, aes(year, value)) +
  geom_hline(yintercept = 0, color = "dark grey") +
  geom_line() +
  geom_point(color = dat$color) +
  facet_wrap(~name, ncol = 5) +
  geom_vline(xintercept = c(2014, 2019), lty = 2, color = cb[8]) +
  theme(axis.title.x = element_blank())
