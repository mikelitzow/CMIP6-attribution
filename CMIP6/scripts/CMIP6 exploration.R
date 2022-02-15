# explore CMIP6 runs from Trond

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)


# set palette
new.col <- oceColorsPalette(64)

# set theme
theme_set(theme_bw())

## process model outputs ---------------------

# compare original CMIP filenames (1900-2099)
# with new filenames (1850-2099)

files.orig <- list.files("./CMIP6/CMIP6_outputs/1900-2099_runs")

# save as a table for write-up
drop <- grep("245", files.orig)

files.orig <- files.orig[-drop] %>%
  str_remove(".nc")

write.csv(files.orig, "./CMIP6/figs/files.table.1900-2099.csv")

files.new <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585")

files.new <- files.new %>%
  str_remove(".nc")

write.csv(files.new, "./CMIP6/figs/files.table.1900-2099.csv")

# merge
new <- data.frame(models = files.new,
                  'run_1850-2099' = "yes")

orig <- data.frame(models = files.orig,
                   'run_1900-2099' = "yes")

both <- full_join(orig, new)

write.csv(both, "./CMIP6/figs/files.table.both.runs.csv", row.names = F)

# loop through each file, save the file name and experiment,
# capture the time series of raw temps for area of interest for each file-experiment comparison

# objects for saving dates and temps and file-experiment ID from each realization
dates <- temps  <- matrix()
experiment.file <- data.frame()
files.new <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585")

for(i in 1:length(files.new)){

  # i <- 1
  
  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.new[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # nc
  
  # extract one experiment at a time and save with file name
  
  # get list of experiments
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  
  for(j in 1:length(experiments)){
  
  # j <- 1
  
  this.one <- data.frame(file = files.new[i], 
                         experiment = ncvar_get(nc, "experiment", start = j, count = 1))
  
  experiment.file <- rbind(experiment.file,
                           this.one)
  
  # extract dates

  d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
  
  this.name <- paste(files.new[i], ncvar_get(nc, "experiment", start = j, count = 1), sep = "_")

  these.dates <- data.frame(this.name = as.character(d))
  dates <- cbind(dates, these.dates)

  # extract spatial area
  # 20-68 deg. N, 120-250 deg. E
  x <- ncvar_get(nc, "lon")
  y <- ncvar_get(nc, "lat")
  

  SST <- ncvar_get(nc, "tos", verbose = F, start = c(j,1,1,1), count = c(1,-1,-1,-1))


  # Change data from a 3-D array to a matrix of monthly data by grid point:
  # First, reverse order of dimensions ("transpose" array)
  SST <- aperm(SST, 3:1)

  # Change to matrix with column for each grid point, rows for monthly means
  SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

  # Keep track of corresponding latitudes and longitudes of each column:
  lat <- rep(y, length(x))   
  lon <- rep(x, each = length(y))   
  dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))


  drop <- lon < 198 | lat < 55 
  SST[,drop] <- NA
  
  drop <- lon == 198.5 & lat > 56
  SST[,drop] <- NA
  
  drop <- lon == 199.5 & lat > 56
  SST[,drop] <- NA
  
  drop <- lon == 200.5 & lat > 56
  SST[,drop] <- NA

  drop <- lon == 201.5 & lat > 57
  SST[,drop] <- NA  
  
  drop <- lon == 202.5 & lat > 58
  SST[,drop] <- NA  
  # and plot 
# 
  # # png("./figs/CMIP6_spatial_domain.png", 6, 4.5, units = "in", res = 300)
  # SST.mean <- colMeans(SST)
  # z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  # image(x,y,z, col=new.col, ylim = c(40, 64), xlim = c(140, 240))
  # contour(x, y, z, add=T, col="white")
  # map('world2Hires',fill=F,add=T, lwd=2)
  # # 
  # # dev.off()

  
  # get average temperature
  these.temps <- data.frame(this.name = rowMeans(SST, na.rm = T))
  names(these.temps) <- this.name
  
  temps <- cbind(temps, these.temps)
  
  }
}

# check that all dates are the same
View(dates)

# just checked visually, but yes, they appear to be!


## load and process ERSSTv5 observations ----------------------------------

# load ERSST for comparison - 1900 through 2020

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(2021-12-01T00:00:00Z)][(0.0):1:(0.0)][(52):1:(62)][(198):1:(226)]", "~temp")


# load and process SST data
# nc <- nc_open("~temp")

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

# and check
temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(180,240), ylim=c(40,62))
contour(x, y, z, add=T)  
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")

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

## model evaluation --------------------------------------------

# evaluate behavior of ERSST - in the past we've seen implausible
# monthly values in the 1st half of the 20th century

# first, calculate monthly anomalies

f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
mu <- apply(SST, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location

sst.anom <- SST - mu   # compute matrix of anomalies

obs.plot <- data.frame(monthly.anom = rowMeans(sst.anom, na.rm = T),
                       date = lubridate::parse_date_time(x = paste(yr,as.numeric(m),"01"),orders="ymd",tz="America/Anchorage"))

ggplot(obs.plot, aes(date, monthly.anom)) +
  geom_line(lwd = 0.02)

ggsave("./CMIP6/figs/GOA_ERSST_monthly_amonalies_1854-2021.png", width = 7, height = 4, units = 'in')

# Strong pattern - large monthly anomalies before 1950
# I think we should limit comparison to 1950:2014
# combine observations and historical simulations in one df

compare.sst <- data.frame(year = as.numeric(names(ann.sst)),
                          ERSST = ann.sst) %>%
  filter(year %in% 1950:2014)

ggplot(compare.sst, aes(year, ERSST)) +
  geom_line()

# now get annual means of historical experiments!

# start with vector for selecting historical simulations

hist <- grep("hist", experiment.file$experiment)

# drop first column of temps
temps.use <- temps[,2:ncol(temps)]

# check we have the correct number of temperature outputs
identical(ncol(temps.use), nrow(experiment.file)) # looks good

historical.temps <- temps.use[, hist]
names(historical.temps) <- experiment.file$file[hist]

# get annual means for each

# get years for one of the runs (they're all the same!)
path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.new[i], sep="")

nc <- nc_open(path)

d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))

years <- as.numeric(as.character(years(d)))


f <- function(x) tapply(x, years, mean)

ann.hist.temps <- as.data.frame(apply(historical.temps, 2, f))

ann.hist.temps$year <- 1850:2099

# this left_join restricts to the reference period 
# defined by the overlap of trustworthy ERSST  for GOA regions (1950-2012)
# and historical simulation period for CMIP6 (1850-2014)

compare.sst <- left_join(compare.sst, ann.hist.temps) 

# seems to be an offset in units??

# drop year
compare.sst <- compare.sst[,2:ncol(compare.sst)]

mean.compare <- data.frame(names = as.factor(names(compare.sst)),
  mean = colMeans(compare.sst))

# looks like some are in K - change
change <- mean.compare$mean > 200

compare.sst[,change] <- compare.sst[,change] - 273.15


mean.compare <- data.frame(names = as.factor(names(compare.sst)),
                           mean = colMeans(compare.sst))

ggplot(mean.compare, aes(names, mean)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))


# compare bias in climatology and correlation with historical

ersst <- compare.sst[,1]
historical.runs <- compare.sst[,2:ncol(compare.sst)]

bias <- colMeans(historical.runs) - mean(ersst)


plot.bias <- data.frame(name = names(bias),
                        bias = abs(bias)) 

# clean up model names
plot.bias$name <- str_remove(plot.bias$name, ".nc")

plot.bias$name <- reorder(plot.bias$name, plot.bias$bias)

theme_set(theme_bw())

ggplot(plot.bias, aes(name, bias)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  labs(title = "Bias, 1950-2014", y = "Absolute value of bias (°C)")

ggsave("./CMIP6/figs/model_bias.png", width = 7, height = 5, units = 'in')

# calculate correlation for low-frequency variability (smoothed with 10-yr running means)
ff <- function(x) zoo::rollmean(x, 10, fill = NA)

historical.smoothed <- apply(historical.runs, 2, ff)
ersst.smoothed <- zoo::rollmean(ersst, 10, fill = NA)

correlation <- cor(historical.smoothed, ersst.smoothed, use = "p")

plot.correlation <- data.frame(name = rownames(correlation),
                        correlation = correlation) 

plot.correlation$name <- str_remove(plot.correlation$name, ".nc")
plot.correlation$name <- reorder(plot.correlation$name, plot.correlation$correlation)

arrange(plot.correlation, desc(correlation))

ggplot(plot.correlation, aes(name, correlation)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) 

ggsave("./CMIP6/figs/model_correlation.png", width = 7, height = 5, units = 'in')

# combine a plot of bias and correlation
bias.corr <- left_join(plot.bias, plot.correlation)

ggplot(bias.corr, aes(bias, correlation, label = name)) +
  geom_text() 

# suggests a cluster of high correlation, low bias models

ggsave("./CMIP6/figs/model_bias_vs_correlation.png", width = 8, height = 6, units = 'in')

# finally, ar(1) values for ersst and models
ff <- function(x) ar(x, order.max = 1)$ar

historical.ar <- apply(historical.runs, 2, ff)
ersst.ar <- ar(ersst, order.max = 1)$ar

# plot difference from actual ar vs bias
plot.ar.difference <- data.frame(name = str_remove(names(historical.ar), ".nc"),
                                       ar.difference = abs(historical.ar - ersst.ar))

bias.ar <- left_join(plot.bias, plot.ar.difference)

ggplot(bias.ar, aes(bias, ar.difference, label = name)) +
  geom_text()

# suggest we weight models based on ar(1), bias, and correlation
weights <- left_join(bias.ar, plot.correlation)

# scale all three
ff <- function(x) as.vector(scale(x))

scaled.weights <- weights
scaled.weights[,2:4] <- apply(scaled.weights[,2:4], 2, ff)

# now change sign of weights - i.e., bias and ar.difference are better when low

signed.weights <- scaled.weights %>% 
  mutate(bias = -bias,
         ar.difference = -ar.difference)

# get scaling factor - set smallest value for each = 1

ff <- function(x) 1 - min(x)
  
scaling.factor <- apply(signed.weights[,2:4], 2, ff)

model.weights <- signed.weights

for(i in 2:4){

  model.weights[,i] <- model.weights[,i] + scaling.factor[(i-1)]

}


# evaluate weights squared
model.weights[,2:4] <- model.weights[,2:4]^2

model.weights$total <- apply(model.weights[,2:4], 1, sum)
max(model.weights$total)/min(model.weights$total)

plot.weights <- model.weights %>%
  rename(model = name) %>%
  pivot_longer(cols = -model)

ggplot(plot.weights, aes(value)) +
  geom_histogram(bins = 5, fill = "dark grey", color = "black") + 
  facet_wrap(~name, scales = "free_x")

# plot by model name
plot.weights$model <- reorder(plot.weights$model, plot.weights$value)

ggplot(plot.weights, aes(model, value)) +
  geom_bar(stat = "identity", fill = "dark grey", color = "black") + 
  facet_wrap(~name, scales = "free_y", ncol = 1) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# this looks reasonable to me! 
model.weights$total/min(model.weights$total)

# seven models with weights > 4 times greater than worst model
# so this will have the result of downweighting historically poor models
# and treating the generally-similar models as a body


# select all the models
best.model.smoothed <- as.data.frame(historical.smoothed)

names(best.model.smoothed) <- str_remove(names(best.model.smoothed), ".nc")

# and for bias
names(bias) <- str_remove(names(bias), ".nc")

best.model.bias <- bias

# check that names line up
identical(names(best.model.bias), names(best.model.smoothed)) # yes!

# bias-correct smoothed time series for best models -
# subtract bias from each
best.model.smoothed.bias.corrected <- best.model.smoothed - best.model.bias


# and scale
ff <- function(x) as.vector(scale(x))

best.model.smoothed.anomaly <- apply(best.model.smoothed, 2, ff)

best.model.smoothed.anomaly <- as.data.frame(best.model.smoothed.anomaly) %>% 
  mutate(year = 1950:2014) %>%
  pivot_longer(cols = -year, values_to = "anomaly") %>%
  rename(model = name)

ersst.smoothed.anomaly <- data.frame(year = 1950:2014,
                                     anomaly = as.vector(scale(ersst.smoothed)))


ggplot(best.model.smoothed.anomaly, aes(year, anomaly, color = model)) +
  geom_line() +
  geom_line(data = ersst.smoothed.anomaly, aes(year, anomaly, color = model), color = "black", lwd = 1) +
  ggtitle("10-yr running mean SST 1900-2020 (ERSST observations in black)") +
  ylab("SST (anomaly)") +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, color = "dark grey")


ggsave("./CMIP6/figs/ersst_vs_best_bias_models_smoothed.png", width = 8, height = 5, units = 'in')


## evaluate North Pacific-wide warming rate for the same models -------------------

# compare historical / ssp585 with pre-industrial 

# start from the beginning - get list of file names
files <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/")

# loop through each file, save the file name and experiment,
# capture the time series of raw temps for area of interest for each file-experiment comparison


# objects for saving dates and temps and file-experiment ID from each realization
piControl.temps <- hist.585.temps  <- matrix()

for(i in 1:length(files)){
  
  # i <- 1
  
  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # nc
  
  # extract one experiment at a time and save with file name
  
  # get list of experiments
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  
  for(j in 1:length(experiments)){
    
    # j <- 2
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "hist_ssp585"){
      
    temp <- data.frame(model =  files[i],
                       experiment = ncvar_get(nc, "experiment", start = j, count = 1))
    

    
    # extract dates
    
    d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
    yr <- as.numeric(as.character(years(d)))
    
    # extract spatial area
    x <- ncvar_get(nc, "lon")
    y <- ncvar_get(nc, "lat")
    
    SST <- ncvar_get(nc, "tos", verbose = F, start = c(j,1,1,1), count = c(1,-1,-1,-1))
    
    # Change data from a 3-D array to a matrix of monthly data by grid point:
    # First, reverse order of dimensions ("transpose" array)
    SST <- aperm(SST, 3:1)
    
    # Change to matrix with column for each grid point, rows for monthly means
    SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3])) 
    
    # remove  values south of 20N
    drop <- lat < 20
    
    SST[,drop] <- NA
    
    # plot to check 
    png(paste("./CMIP6/figs/", files[i], "_hist_spatial_means_1850-1949.png", sep = ""), 6, 4.5, units = "in", res = 300)
    SST.mean <- colMeans(SST[yr < 1950,])
    z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
    image(x,y,z, col=new.col)
    contour(x, y, z, add=T, col="white")
    map('world2Hires',fill=F,add=T, lwd=2)
    
    dev.off()
    
    these.temps <- data.frame(xx = rowMeans(SST, na.rm = T))
    names(these.temps)<-  str_remove(files[i], ".nc")
    
    hist.585.temps <- cbind(hist.585.temps, these.temps)
    
    } # close if 
    
    
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "piControl"){
      

      # extract spatial area
      x <- ncvar_get(nc, "lon")
      y <- ncvar_get(nc, "lat")
      
      # Keep track of corresponding latitudes and longitudes of each column:
      lat <- rep(y, length(x))   
      lon <- rep(x, each = length(y))   
      dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
      
      
      SST <- ncvar_get(nc, "tos", verbose = F, start = c(j,1,1,1), count = c(1,-1,-1,-1))
      
      
      # Change data from a 3-D array to a matrix of monthly data by grid point:
      # First, reverse order of dimensions ("transpose" array)
      SST <- aperm(SST, 3:1)
      
      # Change to matrix with column for each grid point, rows for monthly means
      SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3])) 
      
      # remove  values south of 20N
      drop <- lat < 20
      
      SST[,drop] <- NA
      
      # plot to check 
      png(paste("./CMIP6/figs/", files[i], "_piControl_spatial_means.png", sep = ""), 6, 4.5, units = "in", res = 300)
      SST.mean <- colMeans(SST)
      z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
      image(x,y,z, col=new.col)
      contour(x, y, z, add=T, col="white")
      map('world2Hires',fill=F,add=T, lwd=2)
      
      dev.off()
      

      these.temps <- data.frame(xx = rowMeans(SST, na.rm = T))
      names(these.temps)<-  str_remove(files[i], ".nc")
      
      piControl.temps <- cbind(piControl.temps, these.temps)
      
    } # close if 
    
  } # close j
} # close i

# remove leading column of NAs
piControl.temps <- piControl.temps[,-1]
hist.585.temps <- hist.585.temps[,-1]

# get annual means
ff <- function(x) tapply(x, as.numeric(as.character(years((d)))), mean)

hist.585.annual <- apply(hist.585.temps, 2, ff)
piControl.annual <- apply(piControl.temps, 2, ff)



# and ersst for the same area

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(2021-12-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(66)][(110):1:(250)]", "~temp")

# load and process SST data
# nc <- nc_open("~temp")

nc <- nc_open("./CMIP6/data/nceiErsstv5_c5fc_6a40_5e5b.nc")

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

temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

# calculate monthly mean
obs.n.pac.sst <- rowMeans(SST, na.rm = T)

# and annual observed means
n.pac.ann.obs.sst <- data.frame(year = 1854:2021,
                             ersst = tapply(obs.n.pac.sst, as.numeric(as.character(yr)), mean),
                             observations = "ERSST")

# and observed warming rate
n.pac.obs.warming <- data.frame(year = 1854:2021,
                                ersst.warming = n.pac.ann.obs.sst$ersst - 
                                  mean(n.pac.ann.obs.sst$ersst[n.pac.ann.obs.sst$year %in% 1854:1949]))
# switch to C from K!
K_C <- colMeans(hist.585.annual) > 200

hist.585.annual[,K_C] <- hist.585.annual[,K_C] - 273.15
piControl.annual[,K_C] <- piControl.annual[,K_C] - 273.15

# compare time series!

check.hist <- as.data.frame(hist.585.annual) %>%
  mutate(experiment = "hist.ssp585",
         year = 1850:2099) %>%
  pivot_longer(cols = c(-year, -experiment))

check.control <- as.data.frame(piControl.annual) %>%
  mutate(experiment = "piControl",
         year = 1850:2099) %>%
  pivot_longer(cols = c(-year, -experiment))

check <- rbind(check.control, check.hist)

ggplot(check, aes(year, value, color = experiment)) +
  geom_line() +
  geom_line(data = n.pac.ann.obs.sst, aes(year, ersst, color = observations), lwd = 1) +
  facet_wrap(~name)

ggsave("./CMIP6/figs/control_hist585_comparison.png", width = 8, height = 6, units = 'in')

# now calculate warming rate
warming.rate <- hist.585.annual

for(j in 1:ncol(hist.585.annual)){
  
  warming.rate[,j] <- hist.585.annual[,j] - colMeans(hist.585.annual[1:100,])[j]
  
}

colMeans(hist.585.annual[1:100,])

warming.rate <- as.data.frame(warming.rate) %>%
  mutate(year = 1850:2099) %>%
  pivot_longer(cols = -year)


ggplot(warming.rate, aes(year, value, color = name)) +
  geom_line() +
  geom_line(data = n.pac.obs.warming, aes(year, ersst.warming), color = "black", lwd = 1) +
  ggtitle("Annual SST 1900-2099 (ERSST observations in black)") +
  ylab("SST (°C)") +
  theme(axis.title.x = element_blank())

ggsave("./CMIP6/figs/n_pac_model_estimated_warming_rate.png", width = 7, height = 5, units = 'in')

# save these values for future analysis
write.csv(warming.rate, "./CMIP6/summaries/ne_pacific_annual_modeled_sst.csv")
write.csv(n.pac.obs.warming, "./CMIP6/summaries/ne_pacific_annual_observed_sst.csv")

# need to:

# 1) fit loess to full time series for each model and get estimated time
# for 0.5, 1, 1.5, 2 degrees of warming wrt 1850-1949

# 2) fit loess to ersst to calculate warming wrt 1854-1949 through 2021

# 3) use differences between model and observed to weight each model

# 4) fit Bayesian regression model to estimate overall warming event

# 1) fit loess smoothers to save predicted values

warming.rate <- warming.rate %>%
  pivot_wider(names_from = name, values_from = value)
warming.rate <- as.matrix(warming.rate)

# create an object to catch
n.pac.warming.timing <- warming.rate

for(j in 2:ncol(warming.rate)){
  
  mod <- loess(warming.rate[,j] ~ warming.rate[,1])
  n.pac.warming.timing[,j] <- predict(mod)
  
} 

# plot to check
check.plot <- as.data.frame(n.pac.warming.timing) %>%
  pivot_longer(cols = -year)

ggplot(check.plot, aes(year, value)) +
  geom_line() + 
  facet_wrap(~name) +
  geom_hline(yintercept = c(0.5, 1, 1.5, 2), lty = 2, color = "red") +
  coord_cartesian(xlim = c(2000,2050), ylim = c(0,4))

# looks good!

# get the year that each warming level is reached for each model
levels <- c(0.5, 1, 1.5, 2)

timing <- data.frame()
n.pac.warming.timing <- as.data.frame(n.pac.warming.timing)

for(j in 2:ncol(n.pac.warming.timing)){
  # j <- 2
  temp <- NA
  
  for(i in 1:length(levels)){
    # i <- 1
    
    temp[i] <- min(n.pac.warming.timing$year[n.pac.warming.timing[,j] >= levels[i]])
    
  }
  
  timing <- rbind(timing,
                  data.frame(
                    model = colnames(n.pac.warming.timing[j]),
                    level = levels,
                    year = temp))
}

# and plot
ggplot(timing, aes(year)) +
  geom_histogram(bins = 6) +
  facet_wrap(~level, scales = "free_x")

ggplot(timing, aes(level, year, color = model)) +
  geom_point() +
  labs(x = "N. Pacific warming wrt 1850-1949 (°C)",
       y = "Year first reached") +
  scale_y_continuous(breaks = seq(1960, 2090, by = 10))

timing$level <- as.factor(timing$level)

obs.timing <- data.frame()


# fit loess to ersst to get observed warming rate
mod <- loess(n.pac.obs.warming[,2] ~ n.pac.obs.warming[,1])
n.pac.obs.warming$trend <- predict(mod)

# plot to check
ggplot(n.pac.obs.warming, aes(year, trend)) +
  geom_line()

# looks right


# now get observed timing of different warming levels - 
# we've only reached 1, which is 0.5 

  i <- 1
  temp <- min(n.pac.obs.warming$year[n.pac.obs.warming$trend >= levels[i]]) 

  obs.timing <- data.frame(level = levels,
                           timing = c(2003, NA, NA, NA))
  
  # 2003
  
  
  
  
obs.timing$level <- as.factor(obs.timing$level)

ggplot(timing, aes(y = year, level)) +
  geom_boxplot() +
  labs(x = "N Pacific warming wrt 1850-1949 (°C)",
       y = "Year first reached",
       title = "Red dot = ERSSTv5 warming") +
  geom_point(data = obs.timing, aes(y = timing, as.factor(level)), color = "red") 

ggsave("./CMIP6/figs/ne_pacific_warming_rate_models_obs.png", width = 5, height = 4, units = 'in')



# save timing
write.csv(timing, "./CMIP6/summaries/model.ne.pacific.warming.timing.csv")


# evaluate ability of different models to predict warming for 1950-2021

predict.timing <- n.pac.warming.timing %>%
 filter(year %in% 1972:2021) %>%
  pivot_longer(cols = -year)

response.timing <- n.pac.obs.warming %>%
  filter(year %in% 1972:2021) %>%
  select(year, trend) %>%
  rename(ersst = trend)

predict.timing <- left_join(predict.timing, response.timing)

# plot to check
ggplot(predict.timing, aes(value, ersst)) +
  geom_point() +
  facet_wrap(~name) # avoids a 'fishook' around declining rate of 1950s

model.warming.evaluation <- data.frame()

models <- unique(predict.timing$name)

for(i in 1:length(models)){
  
  # i <- 1
  
  temp <- predict.timing %>%
    filter(name == models[i])
  
  linear.fit <- lm(ersst ~ value, data = temp)
  
  RSS <- c(crossprod(linear.fit$residuals))
  Pearson.resid <- RSS / linear.fit$df.residual
  
  
  model.warming.evaluation <- rbind(model.warming.evaluation,
                                    data.frame(model = models[i],
                                               coeff = coefficients(linear.fit)[2],
                                               Pearson.resid = Pearson.resid))
  
}


model.warming.evaluation # should use coefficients! (inverse of difference from 1)

model.warming.evaluation$coeff.from.one <- 1-model.warming.evaluation$coeff
model.warming.evaluation$weight <- abs(1/model.warming.evaluation$coeff.from.one)

ggplot(model.warming.evaluation, aes(weight)) +
  geom_histogram(bins = 8, fill = "grey", color = "black") # seems right!


ggplot(model.warming.evaluation, aes(abs(coeff.from.one), weight)) +
  geom_point()

# save 
write.csv(model.warming.evaluation, "./CMIP6/summaries/N_Pac_warming_model_weights.csv", row.names = F)


# next, brms estimates of warming timing


#######################################
## below is old
#######################################

# calculate warming rate wrt preindustrial runs
ne.pac.rate <- as.data.frame(warming.rate) %>%
  mutate(year = 1900:2099) %>%
  filter(year %in% 1980:2014)

ne.pac.rate <- left_join(ne.pac.rate, ne.pac.ann.obs.sst) %>%
  select(-year)

yrs <- 1980:2014

ff <- function(x) summary(lm(x ~ yrs))$coef[2,1]*10
rate <- apply(ne.pac.rate, 2, ff)

coefs <- data.frame(name = names(ne.pac.rate),
                    warming.rate = rate) 

coefs$name <- reorder(coefs$name, coefs$warming.rate)

ggplot(coefs, aes(name, warming.rate)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  ggtitle("NE Pacific warming rate, 1981-2014") +
  ylab("Warming rate (°C / decade)")

ggsave("./figs/model_warming_rate_NE_Pacific_1981-2014.png", width = 7, height = 5, units = 'in')


# now plot annual (non-smoothed) anomalies and also compare AR(1) values with observations

# clean up names from annual runs
names(historical.runs) <- str_remove(names(historical.runs), ".nc")

best.model.annual <- historical.runs[names(historical.runs) %in% use] 

# calculate anomalies
best.model.annual.anomaly <- as.data.frame(apply(best.model.annual, 2, ff)) %>%
  mutate(year = 1900:2020) %>%
  pivot_longer(cols = -year, values_to = "anomaly") %>%
  rename(model = name)

ersst.annual.anomaly <- data.frame(year = 1900:2020,
                                   anomaly = as.vector(scale(ersst)))


ggplot(best.model.annual.anomaly, aes(year, anomaly, color = model)) +
  geom_line() +
  geom_line(data = ersst.annual.anomaly, aes(year, anomaly, color = model), color = "black", lwd = 1) +
  ggtitle("Annual SST 1900-2020 (ERSST observations in black)") +
  ylab("SST (anomaly)") +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, color = "dark grey")

ggsave("./figs/ersst_vs_best_bias_models_annual.png", width = 8, height = 5, units = 'in')






## now plot full time series for selected models under both SSPs, along with observations

# identify the models listed in "use"; subset these out of the temps matrix using experiment.file
# clean up file names
experiment.file$file <- str_remove(experiment.file$file, ".nc")
experiment.file$file <- str_remove(experiment.file$file, "_245")

keep <- vector()

for(i in 1:length(use)){
  
  temp <- grep(use[i], experiment.file$file)
  keep <- c(keep, temp)
  
}

use.file <- experiment.file[keep,]

use.temps <- temps[, keep]

# plot annual values under ssp.585, ssp.245, pi.control

# trying this again, more deliberately!

# get annual means for temps we're going to use
years <- as.numeric(as.character(years(d)))

ff <- function(x) tapply(x, years, mean)

temps.annual <- as.data.frame(apply(temps, 2, ff))

# pick out models one at a time and plot

use.models <- unique(use.file$file)
model.output <- data.frame()

for(i in 1:length(use.models)){
  
# i <- 7

pick.mod <- grep(use.models[i], names(temps.annual)) 

pick.temps <- temps.annual[,pick.mod]

# change from Kelvin if need be
if(mean(pick.temps[,1]) > 200) {pick.temps <- pick.temps-273.15}

names(pick.temps) <- experiment.file$experiment[pick.mod]

plot.pick <- pick.temps %>%
  mutate(year = 1900:2099) %>%
  pivot_longer(cols = -year)

ggplot(plot.pick, aes(year, value, color = name)) +
  geom_line()

# convert to anomaly relative to 20th century climatology
# for(j in 1:ncol(pick.temps)){
#   
#   j <- 1
#   for(ii in 1:nrow(pick.temps)){
#   pick.temps[ii,j] <- (pick.temps[ii,j] - mean(pick.temps[1:100,j])) / sd(pick.temps[1:100,j])
#   
# }}


mean.pick <- colMeans(pick.temps[1:100,])
mean.pick <- matrix(mean.pick, nrow = length(1900:2099), ncol = 3, byrow = T)

pick.temps <- pick.temps - mean.pick


sd.pick <- apply(pick.temps[1:100,], 2, sd)
sd.pick <- matrix(sd.pick, nrow = length(1900:2099), ncol = 3, byrow = T)

pick.temps <- pick.temps / sd.pick

plot.pick <- pick.temps %>%
  mutate(year = 1900:2099) %>%
  pivot_longer(cols = -year)

ggplot(plot.pick, aes(year, value, color = name)) +
  geom_line() +
  ggtitle(use.models[i])

pick.temps$year <- 1900:2099
pick.temps$model <- use.models[i]

pick.temps <- pick.temps %>%
  pivot_longer(cols = c(-year, -model))

model.output <- rbind(model.output, pick.temps)

}

ggplot(model.output, aes(year, value, color = name)) +
  geom_line() +
  facet_wrap(~model)

ggplot(model.output, aes(year, value, color = model)) +
  geom_line() +
  facet_wrap(~name)

# limit to 585 and control
plot.model <- model.output %>%
  filter(name != "hist_ssp245") %>%
  rename(anomaly = value)


# calculate model mean for each experiment
mod.mean <- plot.model %>%
  group_by(year, name) %>% 
  summarise(model.mean = mean(anomaly)) 

# and get ersst anomaly wrt 20th century

ersst.plot <- data.frame(year = 1900:2020,
                         anomaly = (ersst - mean(ersst[1:100])) / sd(ersst[1:100]))

ggplot(plot.model, aes(year, anomaly, color = model)) +
  geom_line() +
  facet_wrap(~name, ncol = 1) +
  geom_line(data = mod.mean, aes(year, model.mean), color = "black", lwd = 1) +
  geom_line(data = ersst.plot, aes(year, anomaly), color = "red", lwd = 1) +
  labs(y = "SST anomaly", title = "ERSST in red, multi-model mean in black") +
  theme(axis.title.x = element_blank()) +
  scale_x_continuous(breaks = seq(1900, 2100, 20)) +
  coord_cartesian(xlim = c(1900, 2100))

ggsave("./figs/ssp585_picontrol_anomaly_time_series.png", width = 9, height = 7, units = 'in')

# save for analysis elsewhere!
write.csv(plot.model, "./summaries/selected_model_scaled_GOA_output.csv")
write.csv(ersst.plot, "./summaries/ersst_scaled_GOA_output.csv")

## get estimate of timing for 0.5, 1. 1.5, 2 degree warming for N. Pacific

# load time series
# save these values for future analysis
ne.pac.annual <- read.csv("./summaries/ne_pacific_annual_modeled_sst.csv", row.names = 1)
ne.pac.annual$model <- str_remove(ne.pac.annual$model, ".nc")

ne.pac.obs <- read.csv("./summaries/ne_pacific_annual_observed_sst.csv", row.names = 1)

# get observed warming wrt 1900-1950
warming <- vector()

for(i in 1:nrow(ne.pac.obs)){
  
  warming[i] <- (ne.pac.obs$ersst[i] - mean(ne.pac.obs$ersst[ne.pac.obs$year %in% 1900:1950])) /
                   sd(ne.pac.obs$ersst[ne.pac.obs$year %in% 1900:1950])
  
  
}

obs.warming <- data.frame(
  year = 1900:2020,
  warming = warming)

ggplot(obs.warming, aes(year, warming)) +
  geom_point() +
  geom_smooth() 

# note that this shows more rapid warming than any of the models
# it appears this is because of the pre-1950 trend
# might be better to include data back to the 1850s (for models!)

mod <- loess(obs.warming$warming ~ obs.warming$year)

obs.warming$trend <- predict(mod)

# for each model, get temperature relative to 1900-1950 mean
ne.pac.warming <- ne.pac.annual %>%
  pivot_wider(names_from = model, values_from = anomaly)

ff <- function(x) x - mean(x[1:51])

ne.pac.warming[,2:ncol(ne.pac.warming)] <- apply(ne.pac.warming[,2:ncol(ne.pac.warming)], 2, ff)


plot.warming <- ne.pac.warming %>%
  pivot_longer(cols = -year)


ggplot(plot.warming, aes(year, value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~name)

ggsave("./figs/ne_pacific_warming_rate.png", width = 9, height = 6, units = 'in')




### following is old, can be transferred to separate FAR script

## calculate FAR -----------------------

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
