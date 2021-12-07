library(tidyverse)
theme_set(theme_bw())

# first, combine data for different periods into a single df

d1 <- read.csv("./data/ShelikofHistorical_86_94_v2.csv")

d2 <- read.csv("./data/ShelikofHistorical_1995_2002.csv")

d3 <- read.csv("./data/ShelikofHistorical_2003_2020.csv")

d4 <- read.csv("./data/ShelikofHistorical_0704.csv")

head(d1); head(d2); head(d3); head(d4)

# check dates on d4 / d3
unique(d4$DateTime) # March 13-28 2007
unique(d3$DateTime)

check <- grep("2007", d3$DateTime)
d3[check,] # no matches - that's good!


# make d1/d2 names compatible
names(d1)[11] <- "Sex"

identical(names(d1), names(d2))

identical(names(d2), names(d3))

names(d2); names(d3)

# examine maturity table
mat <- read.csv("./data/maturity.csv")
mat

# see how these are used in the different data sets
unique(d1$Maturity)
unique(d1$MaturityTable)

unique(d2$Maturity)
unique(d2$MaturityTable)

unique(d3$Maturity)
unique(d3$MaturityIndex)

unique(d4$Maturity)
unique(d4$MaturityIndex)

# for an initial look, I will lump into the original maturity schedule (Table 3)
d1$maturity_table_3 <- d1$Maturity

d2$maturity_table_3 <- ifelse(d2$MaturityTable == 3, d2$Maturity,
                               ifelse(d2$MaturityTable == 11 & d2$Maturity == 1, 1,
                                       ifelse(d2$MaturityTable == 11 & d2$Maturity %in% c(2,3), 2,
                                               ifelse(d2$MaturityTable == 11 & d2$Maturity %in% c(4,5), 3,
                                                       ifelse(d2$MaturityTable == 11 & d2$Maturity == 6, 4, 5)))))
d3$maturity_table_3 <- case_when(
  d3$Maturity == "Immature" ~ 1,
  d3$Maturity %in% c("Developing I", "Developing II", "Developing") ~ 2,
  d3$Maturity %in% c("Prespawning I", "Prespawning II", "Prespawning") ~ 3,
  d3$Maturity == "Spawning" ~ 4,
  d3$Maturity %in% c("Spent I", "Spent II", "Spent") ~ 5)

d4$maturity_table_3 <- case_when(
  d4$Maturity == "Immature" ~ 1,
  d4$Maturity %in% c("Developing I", "Developing II", "Developing") ~ 2,
  d4$Maturity %in% c("Prespawning I", "Prespawning II", "Prespawning") ~ 3,
  d4$Maturity == "Spawning" ~ 4,
  d4$Maturity %in% c("Spent I", "Spent II", "Spent") ~ 5)

# clean up names
names(d1)[6:7] <- names(d2)[6:7] <- names(d3)[6:7] <- names(d4)[6:7] <- c("length", "weight")

# combine
d1 <- d1 %>% 
  select(Survey, Haul, Latitude, Longitude, DateTime, length, weight, Age, Sex, maturity_table_3)

d2 <- d2 %>% 
  select(Survey, Haul, Latitude, Longitude, DateTime, length, weight, Age, Sex, maturity_table_3)

d3 <- d3 %>% 
  select(Survey, Haul, Latitude, Longitude, DateTime, length, weight, Age, Sex, maturity_table_3)

d4 <- d4 %>% 
  select(Survey, Haul, Latitude, Longitude, DateTime, length, weight, Age, Sex, maturity_table_3)

# change d3/d4 weight to g
d3$weight <- 1000*d3$weight
d4$weight <- 1000*d4$weight

dat <- rbind(d1, d2, d3, d4)

# # limit to known-age
# 
# drop <- is.na(dat$Age)
# 
# dat <- dat[!drop,]

# get julian day
temp <- dat %>%
  select(DateTime) 

temp.more <- separate(temp, DateTime, into = c("date", "time"), sep = "\\s")

temp.more$date <- chron::dates(temp.more$date)

dat$julian <- lubridate::yday(temp.more$date)
dat$year <- lubridate::year(temp.more$date)
unique(dat$year)

hist(dat$julian) # looks right!

# # drop age -0
# dat <- dat %>%
#   filter(Age > 0)

ggplot(dat, aes(Age, weight)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), se = F) +
  facet_wrap(~year)

# becomes very age-truncated in recent years!
# examine just age 1-10
ggplot(dplyr::filter(dat, Age %in% 1:10), aes(Age, weight)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), se = F) +
  facet_wrap(~year)

ggplot(dplyr::filter(dat, Age %in% 4:10), aes(Age, weight, color = as.factor(as.character(year)))) +
  geom_point(stroke = 0, size=0) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), se = F) +
  coord_cartesian(ylim = c(200, 2200)) + 
  scale_color_viridis_d()

ggplot(dat, aes(Age)) +
  geom_histogram(bins = 23) +
  facet_wrap(~year, scale = "free_y")

# age structure collapsed (?!)

# examine just spawning fish < 15 years old
ggplot(filter(dat, maturity_table_3 == 4, Age < 15), aes(Age)) +
  geom_histogram(bins = 23) +
  facet_wrap(~year, scale = "free_y")

# and weight of spawning fish
ggplot(filter(dat, maturity_table_3 == 4), aes(weight)) +
  geom_histogram() +
  facet_wrap(~year, scale = "free_y")

spawn.age <- dat %>%
  filter(Age %in% 1:15, maturity_table_3 == 4)

ggplot(spawn.age, aes(as.factor(Age), weight)) +
  geom_violin() +
  facet_wrap(~year)

# mean and sd of spawning age by year
spawn.summary <- spawn.age %>%
  group_by(year) %>%
  summarize(mean.age = mean(Age), sd.age = sd(Age)) %>%
  pivot_longer(cols = - year) %>%
  mutate(year = as.numeric(as.character(year)))

ggplot(spawn.summary, aes(year, value)) +
  geom_line() +
  geom_point() + 
  facet_wrap(~name, scales = "free_y", ncol=1)

# spawn date
spawn.all <- dat %>%
  filter(maturity_table_3 == 4)

ggplot(spawn.all, aes(julian)) +
  geom_histogram(bins=20, fill = "dark grey", color = "black") + 
  facet_wrap(~year, scales = "free_y") +
  coord_cartesian(xlim = c(60, 90))

# add sex 
dat$sex.code <- if_else(dat$Sex == 1, 1,
                        if_else(dat$Sex == 2, 2,
                                if_else(dat$Sex == "Male", 1, 
                                        ifelse(dat$Sex == "Female", 2, NA))))

# fit an exploratory model to weight of spawners with age, sex, and year effects
library(mgcv)

spawn.age <- dat %>%
  filter(Age %in% 1:20, maturity_table_3 == 4)

spawn.age$year <- as.factor(as.character(spawn.age$year))
spawn.age$sex.code <- as.factor(spawn.age$sex.code)

mod1 <- gam(weight ~ s(Age, k = 4) + sex.code + year, na.action = "na.exclude",
           data = spawn.age)

summary(mod1)

plot(mod1)

mod2 <- gam(weight ~ s(Age, by = year, k = 4) + sex.code, na.action = "na.exclude",
           data = spawn.age)

summary(mod2)

# save predicted values to plot
age.range <- spawn.age %>%
  group_by(year) %>%
  summarise(min.age = min(Age),
            max.age = max(Age))

plot.pred <- data.frame()

for(i in 1:nrow(age.range)){
  temp.age <- seq(age.range$min.age[i],
                   age.range$max.age[i],
                   length.out = 100)
  
  temp <- data.frame(
                     Age = temp.age,
                     year = age.range$year[i],
                     sex.code = 2)
  plot.pred <- rbind(plot.pred, temp)
  
}

plot.pred$pred.weight <- predict(mod2, type = "response", newdata = plot.pred)

ggplot(plot.pred, aes(Age, pred.weight, color = year)) +
  geom_line() +
  scale_color_viridis_d()


## look at weight distributions by age;
## include spawning and pre-spawning for now

spawn.prespawn  <- dat %>%
  filter(maturity_table_3 %in% 3:4, 
         Age %in% 1:20)

ggplot(spawn.prespawn, aes(weight)) +
  geom_histogram(fill = "grey90", color="black") + 
  facet_wrap(~Age, scale = "free_y")

# so, few that are under 4 or older than 10
# and the distributions are definitely skewed

# I'll look at the numbers for developing, spawning, prespawning, spent separately
develop.spawn.prespawn.spent  <- dat %>%
  filter(maturity_table_3 %in% 2:5, 
         Age %in% 1:20)

ggplot(filter(develop.spawn.prespawn.spent, maturity_table_3 == 2), aes(as.factor(Age))) +
  geom_histogram(fill = "grey90", color="black", stat = "count") + 
  facet_wrap(~year, scale = "free_y")

ggplot(filter(develop.spawn.prespawn.spent, maturity_table_3 == 3), aes(as.factor(Age))) +
  geom_histogram(fill = "grey90", color="black", stat = "count") + 
  facet_wrap(~year, scale = "free_y")

ggplot(filter(develop.spawn.prespawn.spent, maturity_table_3 == 4), aes(as.factor(Age))) +
  geom_histogram(fill = "grey90", color="black", stat = "count") + 
  facet_wrap(~year, scale = "free_y")

ggplot(filter(develop.spawn.prespawn.spent, maturity_table_3 == 2), aes(as.factor(Age))) +
  geom_histogram(fill = "grey90", color="black", stat = "count") + 
  facet_wrap(~year, scale = "free_y")

## so all 4 classes look reasonable to include!


# 4-8 might capture most of the population

# look at counts per year
ggplot(filter(spawn.prespawn, maturity_table_3 == 3,
              Age %in% 4:10), aes(Age)) +
  geom_histogram(fill = "grey90", color="black") + 
  facet_wrap(~year, scale = "free_y")

ggplot(filter(spawn.prespawn, maturity_table_3 == 4,
              Age %in% 4:10), aes(Age)) +
  geom_histogram(fill = "grey90", color="black") + 
  facet_wrap(~year, scale = "free_y")

# not many spawner samples in some years - 2000-2003, 2009, 2012

# calculate the proportion of developing / pre-spawning / spawning / spent 
# that are in the dominant year class each year (and sample size!)

mature <- dat %>%
  filter(maturity_table_3 %in% 2:5,
         Age %in% 4:10) %>%
  group_by(year, Age) %>%
  summarize(count = n())

ggplot(mature, aes(Age, count)) +
  geom_col(fill = "grey90", color="black", position = "dodge", width = 0.5) + 
  facet_wrap(~year, scale = "free_y")+
  geom_hline(yintercept = 0)

mature <- mature %>%
  pivot_wider(names_from = Age, values_from = count)

mature$n <- apply(mature[,2:8], 1, sum, na.rm=T)

mature$max <- apply(mature[,2:8], 1, max, na.rm=T)

mature$max.prop <- mature$max/mature$n

plot.mature <- mature %>% 
  select(year, n, max.prop) %>%
  pivot_longer(cols = -year)

ggplot(plot.mature, aes(as.numeric(as.character(year)), value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 1)

## addition - use Shannon-Weiner to calculate age diversity of 
# developing / pre-spawning / spawning / spent fish
age.diversity <- dat %>%
  filter(maturity_table_3 %in% 2:5,
         Age > 1) %>%
  group_by(year, Age) %>%
  summarize(count = n()) 

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

# save 
write.csv(plot.shannon, "./data/shannon.age.diversity.acoustic.trawl.data.csv")



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

  