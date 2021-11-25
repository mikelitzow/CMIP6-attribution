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


# get list of file names
files <- list.files("./model_outputs/")

# loop through each file, save the file name and experiment,
# capture the time series of raw temps for area of interest for each file-experiment comparison

# objects for saving dates and temps and file-experiment ID from each realization
dates <- temps  <- matrix()
experiment.file <- data.frame()

for(i in 1:length(files)){

  # i <- 1
  
  path <- paste("./model_outputs/", files[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # nc
  
  # extract one experiment at a time and save with file name
  
  # get list of experiments
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  
  for(j in 1:length(experiments)){
  
  # j <- 1
  
  this.one <- data.frame(file = files[i], 
                         experiment = ncvar_get(nc, "experiment", start = j, count = 1))
  
  experiment.file <- rbind(experiment.file,
                           this.one)
  
  # extract dates

  d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
  
  this.name <- paste(files[i], ncvar_get(nc, "experiment", start = j, count = 1), sep = "_")

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


  drop <- lon > 224 | lon < 202 | lat < 55 | lat > 60
  SST[,drop] <- NA
  
  drop <- lon == 202.5 & lat == 58.5
  SST[,drop] <- NA
  
  # and plot 
  # SST.mean <- colMeans(SST)
  # z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  # image(x,y,z, col=new.col, ylim = c(53, 64), xlim = c(194, 228))
  # contour(x, y, z, add=T, col="white")
  # map('world2Hires',fill=F,add=T, lwd=2)

  # get average temperature
  these.temps <- data.frame(this.name = rowMeans(SST, na.rm = T))
  
  temps <- cbind(temps, these.temps)
  
  }
}

# check that all dates are the same
View(dates)

# just checked visually, but yes, they appear to be!

# load ERSST for comparison - 1950 through 2014

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1950-01-01):1:(2014-12-01T00:00:00Z)][(0.0):1:(0.0)][(54):1:(62)][(200):1:(226)]", "~temp")


# load and process SST data
# nc <- nc_open("~temp")

nc <- nc_open("./data/nceiErsstv5_c83f_8d23_4777.nc")

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
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(195,230), ylim=c(53,62))
contour(x, y, z, add=T)  
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")

# WGOA for cod/pollock paper

drop <- lon > 210 | lon < 202 | lat < 55 | lat > 59
wSST <- SST
wSST[,drop] <- NA


# and check
temp.mean <- colMeans(wSST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(195,230), ylim=c(53,62))
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")

# calculate monthly mean
obs.sst <- rowMeans(wSST, na.rm = T)

# and annual observed means
ann.sst <- tapply(obs.sst, as.numeric(as.character(yr)), mean)

# combine observations and historical simulations in one df
compare.sst <- data.frame(year = as.numeric(names(ann.sst)),
                          historical = ann.sst)


ggplot(compare.sst, aes(year, historical)) +
  geom_line()

# now get annual means of historical experiments!

# start with vector for selecting historical simulations

hist <- grep("hist", experiment.file$experiment)

# drop first column of temps
temps <- temps[,2:ncol(temps)]

# check we have the correct number of temperature outputs
identical(ncol(temps), nrow(experiment.file)) # looks good

historical.temps <- temps[, hist]
names(historical.temps) <- experiment.file$file[hist]

# get annual means for each

# get years for one of the runs (they're all the same!)
path <- paste("./model_outputs/", files[i], sep="")

nc <- nc_open(path)

d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))

years <- as.numeric(as.character(years(d)))


f <- function(x) tapply(x, years, mean)

ann.hist.temps <- as.data.frame(apply(historical.temps, 2, f))

ann.hist.temps$year <- 1900:2099

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
  geom_bar(stat = "identity")


# compare bias in climatology and correlation with historical

ersst <- compare.sst[,1]
historical.runs <- compare.sst[,2:ncol(compare.sst)]

bias <- colMeans(historical.runs) - mean(ersst)

correlation <- cor(historical.runs, ersst)

plot.bias <- data.frame(name = names(bias),
                        bias = abs(bias)) 

# drop _245; these are repeats I assume with different scenarios

drop <- grep("_245", plot.bias$name)
plot.bias <- plot.bias[-drop,]

plot.bias$name <- reorder(plot.bias$name, plot.bias$bias)

theme_set(theme_bw())

ggplot(plot.bias, aes(name, bias)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())

ggsave("./figs/model_bias.png", width = 7, height = 5, units = 'in')


plot.correlation <- data.frame(name = rownames(correlation),
                        correlation = correlation) 

drop <- grep("_245", plot.correlation$name)
plot.correlation <- plot.correlation[-drop,]

plot.correlation$name <- reorder(plot.correlation$name, plot.correlation$correlation)


ggplot(plot.correlation, aes(name, correlation)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())


ggsave("./figs/model_correlation.png", width = 7, height = 5, units = 'in')


# select the models with lowest bias and plot against ERSST observations
plot.bias

use <- as.vector(plot.bias$name[plot.bias$bias < 0.5])

best.model.hist <- historical.runs[, use]

best.model.hist$year <- 1950:2014

best.model.hist <-best.model.hist %>%
  pivot_longer(cols = -year)

plot.ersst <- data.frame(year = 1950:2014,
                         ersst = ersst) %>%
  pivot_longer(cols = -year)

ggplot(best.model.hist, aes(year, value, color = name)) +
  geom_line() +
  geom_line(data = plot.ersst, aes(year, value, color = name), color = "black", lwd = 1) +
  ggtitle("ERSST observations in black") +
  ylab("SST (°C)") +
  theme(axis.title.x = element_blank())

ggsave("./figs/ersst_vs_low_bias_models.png", width = 7, height = 5, units = 'in')
