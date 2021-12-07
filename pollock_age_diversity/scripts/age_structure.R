# evaluate the effect of SST on pollock age structure as estimated by acoustic trawl survey

library(tidyverse)
theme_set(theme_bw())

# load data
dat <- read.csv("./pollock_age_diversity/data/estimated_millions_at_age_tab1.10_GOA_2021assessment.csv")

str(dat)

names(dat) <- c("year", "age1", "age2", "age3", "age4", "age5", "age6", "age7",
                "age8", "age9", "age10", "age11", "age12", "age13", "age14", "age15", "total")

# calculate age proportions
age.prop <- dat %>%
  select(-year, -total)

# remove age 1-3 (not mature!) from dat$total
for(i in 1:nrow(dat)){
dat$total[i] <- sum(dat[i,(5:16)])
}

for(i in 1:nrow(age.prop)){
age.prop[i,] <- age.prop[i,] / dat$total[i]
}

# limit to mature animals, age 3+
age.prop <- age.prop %>%
  select(-age1, -age2, -age3)

# check
rowSums(age.prop) # looks right

# plot to check
plot <- age.prop %>%
  mutate(year = dat$year) %>%
  pivot_longer(cols = -year)


ggplot(plot, aes(name, value)) +
  geom_bar(stat = "identity") +
  facet_wrap(~year)


## use Shannon-Weiner to calculate age diversity

# drop estimates of 0
plot.trimmed <- plot %>%
  dplyr::filter(value > 0)

plot.trimmed$product <- plot.trimmed$value * log(plot.trimmed$value)

shannon <- plot.trimmed %>%
  group_by(year) %>%
  summarise(shannon = -sum(product))

# plot to check
ggplot(shannon, aes(year, shannon)) +
  geom_line() +
  geom_point()





year.total <- age.diversity %>%
  group_by(year) %>%
  summarise(year.total = sum(count))

age.diversity <- left_join(age.diversity, year.total) %>%
  mutate(proportion = count / year.total,
         ln.proportion = log(proportion),
         product = proportion*ln.proportion)

age.shannon <- age.diversity %>%
  group_by(year) %>%
  summarise(shannon = -sum(product))

plot.shannon <- data.frame(year = 1986:2020,
                           shannon = NA)

plot.shannon$shannon <- age.shannon$shannon[match(plot.shannon$year, age.shannon$year)]

ggplot(plot.shannon, aes(year, shannon)) +
  geom_line() +
  geom_point()


# load sst as a covariate
sst <- read.csv("./data/monthly.western.GOA.SST.anomalies.wrt.1980-2020.csv")

# now cross-correlate sst-diversity

sst$shannon <- ifelse(sst$month == 3, age.shannon$shannon[match(sst$year, age.shannon$year)], NA)

sst$lag.neg.36 <- lag(sst$anom, 36)
sst$lag.neg.35 <- lag(sst$anom, 35)
sst$lag.neg.34 <- lag(sst$anom, 34)
sst$lag.neg.33 <- lag(sst$anom, 33)
sst$lag.neg.32 <- lag(sst$anom, 32)
sst$lag.neg.31 <- lag(sst$anom, 31)
sst$lag.neg.30 <- lag(sst$anom, 30)
sst$lag.neg.29 <- lag(sst$anom, 29)
sst$lag.neg.28 <- lag(sst$anom, 28)
sst$lag.neg.27 <- lag(sst$anom, 27)
sst$lag.neg.26 <- lag(sst$anom, 26)
sst$lag.neg.25 <- lag(sst$anom, 25)

sst$lag.neg.24 <- lag(sst$anom, 24)
sst$lag.neg.23 <- lag(sst$anom, 23)
sst$lag.neg.22 <- lag(sst$anom, 22)
sst$lag.neg.21 <- lag(sst$anom, 21)
sst$lag.neg.20 <- lag(sst$anom, 20)
sst$lag.neg.19 <- lag(sst$anom, 19)
sst$lag.neg.18 <- lag(sst$anom, 18)
sst$lag.neg.17 <- lag(sst$anom, 17)
sst$lag.neg.16 <- lag(sst$anom, 16)
sst$lag.neg.15 <- lag(sst$anom, 15)
sst$lag.neg.14 <- lag(sst$anom, 14)
sst$lag.neg.13 <- lag(sst$anom, 13)

sst$lag.neg.12 <- lag(sst$anom, 12)
sst$lag.neg.11 <- lag(sst$anom, 11)
sst$lag.neg.10 <- lag(sst$anom, 10)
sst$lag.neg.9 <- lag(sst$anom, 9)
sst$lag.neg.8 <- lag(sst$anom, 8)
sst$lag.neg.7 <- lag(sst$anom, 7)
sst$lag.neg.6 <- lag(sst$anom, 6)
sst$lag.neg.5 <- lag(sst$anom, 5)
sst$lag.neg.4 <- lag(sst$anom, 4)
sst$lag.neg.3 <- lag(sst$anom, 3)
sst$lag.neg.2 <- lag(sst$anom, 2)
sst$lag.neg.1 <- lag(sst$anom, 1)

sst$lag.0 <- sst$anom

sst$lag.pos.1 <- lead(sst$anom, 1)
sst$lag.pos.2 <- lead(sst$anom, 2)
sst$lag.pos.3 <- lead(sst$anom, 3)
sst$lag.pos.4 <- lead(sst$anom, 4)
sst$lag.pos.5 <- lead(sst$anom, 5)
sst$lag.pos.6 <- lead(sst$anom, 6)
sst$lag.pos.7 <- lead(sst$anom, 7)
sst$lag.pos.8 <- lead(sst$anom, 8)
sst$lag.pos.9 <- lead(sst$anom, 9)
sst$lag.pos.10 <- lead(sst$anom, 10)
sst$lag.pos.11 <- lead(sst$anom, 11)
sst$lag.pos.12 <- lead(sst$anom, 12)

sst$lag.pos.13 <- lead(sst$anom, 13)
sst$lag.pos.14 <- lead(sst$anom, 14)
sst$lag.pos.15 <- lead(sst$anom, 15)
sst$lag.pos.16 <- lead(sst$anom, 16)
sst$lag.pos.17 <- lead(sst$anom, 17)
sst$lag.pos.18 <- lead(sst$anom, 18)
sst$lag.pos.19 <- lead(sst$anom, 19)
sst$lag.pos.20 <- lead(sst$anom, 20)
sst$lag.pos.21 <- lead(sst$anom, 21)
sst$lag.pos.22 <- lead(sst$anom, 22)
sst$lag.pos.23 <- lead(sst$anom, 23)
sst$lag.pos.24 <- lead(sst$anom, 24)

cor.dat <- sst[,6:67]

cors <- cor(cor.dat, use = "p")

dim(cors)

cor.plot <- data.frame(lag = -36:24,
                       cor = cors[1,2:62])


ggplot(cor.plot, aes(lag, cor)) +
  geom_bar(stat = "identity", fill = "grey", color = "dark grey") +
  scale_x_continuous(breaks = seq(-36, 24, by = 6)) +
  labs(x = "lag (months; 0 = March of sampling year",
       y = "SST - Shannon correlation")

## average for the two years prior as covariate for age diversity

sst.annual <- tapply(sst$anom, sst$year, mean)

sst.annual <- data.frame(year = as.numeric(names(sst.annual)),
                        sst = sst.annual,
                        sst.2 = zoo::rollmean(sst.annual, 2, align = "right", fill = NA),
                        sst.3 = zoo::rollmean(sst.annual, 3, align = "right", fill = NA))


ggplot(sst.annual, aes(year, sst)) +
  geom_line() + 
  geom_point() +
  geom_line(aes(year, sst.2), color = "red")


sst.annual <- left_join(sst.annual, age.shannon)

sst.annual <- sst.annual %>%
  pivot_longer(cols = c(-year, -shannon))

ggplot(sst.annual, aes(value, shannon)) +
  geom_point() +
  facet_wrap(~name)

# compare predictive skill / model fit
library(mgcv)

sst.annual <- sst.annual %>%
  pivot_wider(names_from = name)

mod1 <- gam(shannon ~ s(sst, k = 4), data = sst.annual)
mod2 <- gam(shannon ~ s(sst.2, k = 4), data = sst.annual)
mod3 <- gam(shannon ~ s(sst.3, k = 4), data = sst.annual)

MuMIn::AICc(mod1, mod2, mod3) # ain't no denying - mod3 best by far
summary(mod3)
plot(mod3, residuals = T, pch = 19)

# now need to scale weights by age and maturity stage

# first, clean up - only known sex, known age, age 4-10,
# stages developing - spent

mature.weights <- dat %>%
  filter(maturity_table_3 %in% 2:5,
         Age %in% 4:10,
         sex.code %in% 1:2)

# check weight distribution
ggplot(mature.weights, aes(weight)) +
  geom_histogram(fill = "dark grey", color = "black") +
  facet_grid(maturity_table_3~Age, scales = "free")

# few individuals at some age-maturity combinations - this may be 
# important to understanding some of the age class dynamics
# many developing age 4, few spawning
# few developing age 6+

# log transform weight
mature.weights$log.weight <- log(mature.weights$weight)

# check for NAs (age 7 below is plotting as no data)
check <- is.na(mature.weights$log.weight)
sum(check)

filter(mature.weights, Age == 7)
mature.weights[which.min(mature.weights$log.weight),] # one value of weight = 0!!

mature.weights <- mature.weights %>%
  filter(weight > 0)

# need to scale weight by age
mature.weights <- plyr::ddply(mature.weights, "Age", transform, sc.weight = scale(log.weight))

# check
ggplot(mature.weights, aes(log.weight, sc.weight)) +
  geom_point() +
  facet_wrap(~Age, scales = "free")

#  model the progression of weight for each year class through time

# assign year class
levels(mature.weights$year)
mature.weights$year <- as.numeric(as.character(mature.weights$year))

mature.weights$year.class <- mature.weights$year - mature.weights$Age

# plot to check
ggplot(mature.weights, aes(Age)) +
  geom_histogram() +
  facet_wrap(~year.class)

# begin with an analysis of SST and density-dependent effects on weight-at-age
# not using a cohort approach for now!

# load covariates

d1 <- read.csv("./data/western.goa.sst.csv")
d2 <- read.csv("./data/pollock_SSB_2020_SAFE.csv")

# add to mature weights (and simplify)

all.dat <- mature.weights %>%
  select(-DateTime, -Sex)

all.dat <- left_join(all.dat, d1)
all.dat <- left_join(all.dat, d2)
str(all.dat)

library(fRegression)

all.dat$maturity <- as.factor(as.character(all.dat$maturity_table_3))
all.dat$sex.code <- as.factor(all.dat$sex.code)

ggplot(all.dat, aes(nov.feb.wSST, sc.weight)) +
  geom_point() +
  facet_wrap(~Age) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4))


ggplot(all.dat, aes(year, sc.weight)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4))

ggplot(all.dat, aes(Age, sc.weight)) +
  geom_point() +
  facet_wrap(~year.class) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4))

ggplot(all.dat, aes(nov.feb.wSST, sc.weight)) +
  geom_point() +
  facet_wrap(~year.class) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4)) +
  geom_hline(yintercept = 0)

ggplot(dplyr::filter(all.dat, Age %in% 4:10, year.class %in% 1979:2013),
       aes(Age, weight, color = as.factor(as.character(year.class)))) +
  geom_point(stroke = 0, size=0) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), se = F) +
  coord_cartesian(ylim = c(200, 2200)) + 
  scale_color_viridis_d()

m1 <- gamm(weight ~  maturity + sex.code + Age +
             s(nov.feb.wSST, k = 4) + s(poll.SSB.2020, k = 4),
           random = list(year.class=~1), data=all.dat) 

summary(m1$gam)
plot(m1$gam)


all.dat$era <- as.factor(if_else(all.dat$year < 2005, 1,
                       if_else(all.dat$year %in% 2006:2015, 2, 3)))

m2 <- gamm(weight ~  maturity + sex.code + Age +
             s(nov.feb.wSST, k = 4, by = era) + s(poll.SSB.2020, k = 4),
           random = list(year.class=~1), data=all.dat) 

summary(m2$gam)
plot(m2$gam)


m3 <- gamm(weight ~  maturity + sex.code + Age +
             s(nov.feb.wSST, k = 4) + s(poll.SSB.2020, k = 4, by = era),
           random = list(year.class=~1), data=all.dat) 

summary(m3$gam)
plot(m3$gam)

mod.out <- MuMIn::AICc(m1, m2, m3)

mod.out$delta.AICc <- mod.out$AICc - min(mod.out$AICc)

mod.out

# add lagged spring/summer temps

d1 <- read.csv("./data/western.goa.sst.csv")
head(d1)

lagged.temp <- d1 %>%
  select(-nov.feb.wSST) %>%
  mutate(year=year+1)

names(lagged.temp)[2] <- "apr.jul.wSST.lag1"

all.dat <- left_join(all.dat, lagged.temp)

# repeat models m1-m3 with lagged spring/summer temp instead of winter temp
m4 <- gamm(weight ~  maturity + sex.code + Age +
             s(apr.jul.wSST.lag1, k = 4) + s(poll.SSB.2020, k = 4),
           random = list(year.class=~1), data=all.dat) 

summary(m4$gam)
plot(m4$gam)

m5 <- gamm(weight ~  maturity + sex.code + Age +
             s(apr.jul.wSST.lag1, k = 4, by = era) + s(poll.SSB.2020, k = 4),
           random = list(year.class=~1), data=all.dat) 

summary(m5$gam)
plot(m5$gam)

m6 <- gamm(weight ~  maturity + sex.code + Age +
             s(apr.jul.wSST.lag1, k = 4) + s(poll.SSB.2020, k = 4, by = era),
           random = list(year.class=~1), data=all.dat) 

summary(m6$gam)
plot(m6$gam)

mod.out <- MuMIn::AICc(m1, m2, m3, m4, m5, m6)

mod.out$delta.AICc <- mod.out$AICc - min(mod.out$AICc)

mod.out # lagged spring-summer temperatures are much better!

# does the class of temperature drive density dependence?

ggplot(all.dat, aes(year, apr.jul.wSST.lag1)) +
  geom_point() +
  geom_line()

hist(unique(all.dat$apr.jul.wSST.lag1))


all.dat$sst.class <- as.factor(if_else(all.dat$apr.jul.wSST.lag1 < 7, 1,
                             if_else(all.dat$apr.jul.wSST.lag1 > 8, 3, 2)))


m7 <- gamm(weight ~  maturity + sex.code + Age + s(poll.SSB.2020, k = 4, by = sst.class),
           random = list(year.class=~1), data=all.dat) 

summary(m7$gam)
plot(m7$gam)

mod.out <- MuMIn::AICc(m1, m2, m3, m4, m5, m6, m7)

mod.out$delta.AICc <- mod.out$AICc - min(mod.out$AICc)

mod.out # m7 is not very good!

m8 <- gamm(weight ~  maturity + sex.code + Age +
             s(apr.jul.wSST.lag1, k = 4, by = era) + s(poll.SSB.2020, k = 4, by = era),
           random = list(year.class=~1), data=all.dat) 

summary(m8$gam)
plot(m8$gam)

mod.out <- MuMIn::AICc(m1, m2, m3, m4, m5, m6, m7, m8)

mod.out$delta.AICc <- mod.out$AICc - min(mod.out$AICc)

mod.out # mod 8 is best by a mile

ggplot(all.dat, aes(poll.SSB.2020, apr.jul.wSST.lag1)) +
  geom_point() +
  facet_wrap(~era)

m9 <- gamm(weight ~  maturity + sex.code + Age +
             te(apr.jul.wSST.lag1, poll.SSB.2020, k = 4),
           random = list(year.class=~1), data=all.dat) 

summary(m9$gam)
plot(m9$gam)

mod.out <- MuMIn::AICc(m1, m2, m3, m4, m5, m6, m7, m8, m9)

mod.out$delta.AICc <- mod.out$AICc - min(mod.out$AICc)

mod.out # m9 isn't so good

vis.gam(m9$gam, view=c("apr.jul.wSST.lag1", "poll.SSB.2020"),
        plot.type='contour', color='topo')


ggplot(dplyr::filter(all.dat, Age %in% 4:10, year.class %in% 1979:2013),
       aes(poll.SSB.2020, sc.weight)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~year.class)

  