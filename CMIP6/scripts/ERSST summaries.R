# produce time series for annual and winter mean SST for 
# N. Pacific, Bering, GOA, BC coast, and CCE

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)

# set palettes
new.col <- oceColorsPalette(64)
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set theme
theme_set(theme_bw())


## load and process ERSST ------------------------

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
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "")
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

# extract study area
# 54-66 deg. N, 188-202 deg. E
ebs.x <- c(183, 183, 203, 203, 191) 
ebs.y <- c(53, 65, 65, 57.5, 53)

polygon(ebs.x, ebs.y, border = "red", lwd = 2)

# GOA polygon
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)

polygon(goa.x, goa.y, border = "red", lwd = 2)

# BC polygon
bc.x <- c(231, 238, 228, 225, 225)
bc.y <- c(55, 49, 49, 53, 55)

polygon(bc.x, bc.y, border = "red", lwd = 2)

# northern CCE polygon
ncc.x <- c(238, 238, 233, 233, 238)
ncc.y <- c(49, 41, 41, 49, 49)

polygon(ncc.x, ncc.y, border = "red", lwd = 2)

# southern CCE polygon
scc.x <- c(238, 243, 239, 235, 235, 238)
scc.y <- c(41, 33, 33, 39, 41, 41)

polygon(scc.x, scc.y, border = "red", lwd = 2)

# those  areas look fine

# define cells within each polygon and plot to check

#first, ebs
ebs.sst <- as.data.frame(SST)

xp <- cbind(ebs.x, ebs.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ebs.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(ebs.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# now, goa

goa.sst <- as.data.frame(SST)

xp <- cbind(goa.x, goa.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

goa.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(goa.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# BC coast

bc.sst <- as.data.frame(SST)

xp <- cbind(bc.x, bc.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

bc.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(bc.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")


# NCC 

ncc.sst <- as.data.frame(SST)

xp <- cbind(ncc.x, ncc.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ncc.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(ncc.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(210, 245), ylim = c(30, 52))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# SCC
scc.sst <- as.data.frame(SST)

xp <- cbind(scc.x, scc.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)
scc.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(scc.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(210, 245), ylim = c(30, 52))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")


## summarize time series for each region ------------------------------

# first, make a list of data frames - THIS IS CLUNKY PROGRAMMING, COULD BE IMPROVED!

# npac is the full grid

npac.sst <- as.data.frame(SST)

sst.data.frames <- list()
sst.data.frames[[1]] <- npac.sst
sst.data.frames[[2]] <- ebs.sst
sst.data.frames[[3]] <- goa.sst
sst.data.frames[[4]] <- bc.sst
sst.data.frames[[5]] <- ncc.sst
sst.data.frames[[6]] <- scc.sst

# save names of each df for processing

sst.data.names <- c("npac_sst",
               "ebs_sst",
               "goa_sst",
               "bc_sst",
               "ncc_sst",
               "scc_sst")

# and a vector of clean names
sst.clean.names <- c("North_Pacific",
                     "Eastern_Bering_Sea",
                     "Gulf_of_Alaska",
                     "British_Columbia_Coast",
                     "Northern_California_Current",
                     "Southern_California_Current")

# loop through each df, process, summarize, combine, save

# create vector of latitudes to weight mean sst by cell area

cell.weight <- sqrt(cos(lat*pi/180))

# plot to check
hist(cell.weight, breaks = 50)
unique(cell.weight) # looks right

# create a weighted mean function
weighted.cell.mean <- function(x) weighted.mean(x, w = cell.weight, na.rm = T)

# create a function to compute monthly anomalies
monthly.anomalies <- function(x) tapply(x, m, mean) 

# short year and month vectors for 1950-2021
yr.1950.2021 <- yr[yr %in% 1950:2021]
m.1950.2021 <- m[yr %in% 1950:2021]

# and define winter year
winter.year <- if_else(m.1950.2021  %in% c("Nov", "Dec"), as.numeric(as.character(yr.1950.2021)) +1, as.numeric(as.character(yr.1950.2021)))
winter.year <- winter.year[m.1950.2021  %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

# create blank data frame for catching results
temp.time.series <- data.frame()

# processing loop

for(i in 1: length(sst.data.names)){
  
  # pull out sst dataframe of interest
  temp.dat <- sst.data.frames[[i]]
  
  # first, calculate monthly anomalies
  mu <- apply(temp.dat, 2, monthly.anomalies)	# compute monthly means for each time series (cell)
  mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location

  temp.anom <- temp.dat - mu   # compute matrix of anomalies
  
  # calculate weighted monthly means
  temp.monthly.anom <- apply(temp.anom, 1, weighted.cell.mean)
    
  # make data frame for plotting
  obs.plot <- data.frame(monthly.anom = temp.monthly.anom,
                         date = lubridate::parse_date_time(x = paste(yr,as.numeric(m),"01"),orders="ymd",tz="America/Anchorage"))

  ggplot(obs.plot, aes(date, monthly.anom)) +
    geom_line(lwd = 0.02) +
    ylab("Monthly anomaly") +
    theme(axis.title.x = element_blank())

  ggsave(paste("./CMIP6/figs/", sst.data.names[i], "_monthly_amonalies_1854-2021.png", sep = ""), width = 7, height = 4, units = 'in')
  
  
  ## now create time series of annual means
  ## unsmoothed, and two- and three-year running means
 
  # calculate monthly mean temp weighted by area for 1950-2021
  # this time period motivated by assessment of period of trustworthy monthly data in GOA 
  temp.monthly.sst <- apply(temp.dat[yr %in% 1950:2021,], 1, weighted.cell.mean) 
  
  # calculate annual means
  
  temp.annual <- tapply(temp.monthly.sst, yr[yr %in% 1950:2021], mean) 
  
  temp.2yr <- rollmean(temp.annual, 2, fill = NA, align = "right") # for salmon - year before and year of ocean entry
  
  temp.3yr <- rollmean(temp.annual, 3, fill = NA, align = "center") # for salmon - year before, year of, and year after ocean entry

  # calculate winter means
  temp.winter.monthly.sst <- temp.monthly.sst[m.1950.2021 %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

  temp.winter <- tapply(temp.winter.monthly.sst, winter.year, mean)
  
  # change leading and trailing years to NA b/c these are incomplete
  leading <- min(names(temp.winter))
  trailing <- max(names(temp.winter))
  
  temp.winter[names(temp.winter) %in% c(leading, trailing)] <- NA
  
  temp.winter.2yr <- rollmean(temp.winter, 2, fill = NA, align = "right")
  
  temp.winter.3yr <- rollmean(temp.winter, 3, fill = NA, align = "center")
  
  
  # combine into data frame of time series by region
  temp.time.series <- rbind(temp.time.series,
                            data.frame(region = sst.clean.names[i],
                                       year = 1950:2021,
                                       annual.unsmoothed = temp.annual[names(temp.annual) %in% 1950:2021],
                                       annual.two.yr.running.mean = temp.2yr[names(temp.annual) %in% 1950:2021],
                                       annual.three.yr.running.mean = temp.3yr[names(temp.annual) %in% 1950:2021],
                                       winter.unsmoothed = temp.winter[names(temp.winter) %in% 1950:2021],
                                       winter.two.yr.running.mean = temp.winter.2yr[names(temp.winter) %in% 1950:2021],
                                       winter.three.yr.running.mean = temp.winter.3yr[names(temp.winter) %in% 1950:2021]))
  
  }


# plot

plot.dat <- temp.time.series %>%
  select(region,
         year,
         annual.unsmoothed,
         winter.unsmoothed) %>%
  pivot_longer(cols = c(-region, -year))

ggplot(plot.dat, aes(year, value, color = region)) +
  geom_line() +
  facet_wrap(~name, scale = "free_y", ncol = 1) +
  scale_color_manual(values = cb[c(2:7)])
  
  
# and save 
write.csv(temp.time.series, "./CMIP6/summaries/regional_north_pacific_ersst_time_series.csv", row.names = F)