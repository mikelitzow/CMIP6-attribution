# summarize annual GOA sst with and without smoothing
# and calculate corresponding FAR values - 
# both from observations and from models as current probability

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)

## annual SST means --------------------

# load
nc <- nc_open("./CMIP6/data/nceiErsstv5_5d38_fa81_4a70.nc")

# process
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- years(d)

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

SST <- aperm(SST, 3:1)  

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check

# need to drop Bristol Bay cells
BB <- c("N58E200", "N58E202", "N56E200")
SST[,BB] <- NA

# trim spatial area

drop <-  lat < 56 
SST[,drop] <- NA

drop <-  lat > 55 & lon == 198
SST[,drop] <- NA

# and check
temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(195,230), ylim=c(52,62))
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")

# calculate monthly mean
obs.sst <- rowMeans(SST, na.rm = T)

# and annual observed means
ann.sst <- tapply(obs.sst, as.numeric(as.character(yr)), mean)

# combine with 2-year and 3-yr rolling means

# first, limit to 1950-2021
obs.sst <- data.frame(year = 1950:2021,
                      sst = ann.sst[names(ann.sst) %in% 1950:2021])

obs.sst$sst2 <- rollmean(obs.sst$sst, 2, fill = NA, align = "right") # corresponds to year before and year of ocean entry!

obs.sst$sst3 <- rollmean(obs.sst$sst, 3, fill = NA) # year before, year of, year after ocean entry

## calculate FAR -----------------------

# load monthly sst model runs

mod.dat <- read.csv("./CMIP6/summaries/GOA_monthly_sst_piControl_hist585.csv")

# get annual means
yr <- as.numeric(as.character(years(mod.dat$date)))

# fix the bad years before 1969
yr[1:228] <- yr[1:228] - 100


# use 1981-2020 for calculating present-day probability (period with data for fishery)

# separate preindustrial runs; use the entire time series for preindustrial probability calculations

preindustrial <- plot.model %>%
  filter(name == "piControl")

hist.1981.2020 <- plot.model %>%
  filter(name == "hist_ssp585", year %in% 1981:2020)


far.plot <- ersst.plot %>%
  filter(year >= 1981) %>%
  mutate(obs.FAR = NA)

for(i in 1:nrow(far.plot)){
  
  temp <- far.plot$anomaly[i]
  far.plot$obs.FAR[i] <- 1 - ((sum(preindustrial$anomaly >= temp) / length(preindustrial$anomaly)) /
                                (sum(hist.1981.2020$anomaly >= temp) / length(hist.1981.2020$anomaly)))
  
}


ggplot(far.plot, aes(year, obs.FAR)) +
  geom_line()

# now FAR for mean model projections 2015 - 2040
mean.proj.2015.2040 <- plot.model %>%
  filter(name == "hist_ssp585", year %in% 2015:2040) %>%
  group_by(year) %>%
  summarise(anomaly = mean(anomaly)) %>%
  mutate(proj.FAR = NA)

proj.2015.2040 <- plot.model %>%
  filter(name == "hist_ssp585", year %in% 2015:2040)



for(i in 2015:2040){
  # i <- 2015
  temp <- mean.proj.2015.2040$anomaly[mean.proj.2015.2040$year == i]
  mean.proj.2015.2040$proj.FAR[mean.proj.2015.2040$year == i] <- 1 - ((sum(preindustrial$anomaly >= temp) / length(preindustrial$anomaly)) /
                                                                        (sum(hist.1981.2020$anomaly >= temp) / length(hist.1981.2020$anomaly)))
  
}

plot <- data.frame(year = 1981:2040,
                   observed = c(far.plot$obs.FAR, rep(NA, length(2021:2040))),
                   mean.projected = c(rep(NA, length(1981:2014)), mean.proj.2015.2040$proj.FAR)) %>%
  pivot_longer(cols = -year, values_to = "FAR")


ggplot(plot, aes(year, FAR, color = name)) +
  geom_line() +
  scale_color_manual(values = c("red", "black"))

ggsave("./figs/observed_and_mean_projected_FAR.png", width = 6, height = 3, units = 'in')


# next steps: fite Bayesian models to both observed (1981-2020) and hist/projected (1981-2040) TS
# fit each model separately (historical and preindustrial probability from only that model for both obs and hist/proj)
# then fit brms model with year as a factor and model ID as a group-level term; plot 90% CI
