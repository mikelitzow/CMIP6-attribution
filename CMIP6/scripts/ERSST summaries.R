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

# set palette
new.col <- oceColorsPalette(64)

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
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(45, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

# extract study area
# 54-66 deg. N, 188-202 deg. E
ebs.x <- c(183, 183, 203, 203, 191) 
ebs.y <- c(53, 65, 65, 57.5, 53)

polygon(ebs.x, ebs.y, border = "red", lwd = 2)
# ebs.slp <- slp
# 
# xp <- cbind(poly.x, poly.y)
# loc=cbind(lon, lat)
# check <- in.poly(loc, xp=xp)
# 
# ebs.slp[,!check] <- NA

# GOA polygon
goa.x <- c(199, 199, 200, 205, 208, 225, 231, 199)
goa.y <- c(55, 55.6, 56.1, 59, 61, 61, 55, 55)

polygon(goa.x, goa.y, border = "red", lwd = 2)

# BC polygon

# calculate monthly mean
obs.n.pac.sst <- rowMeans(SST, na.rm = T)

# and annual observed means
n.pac.ann.obs.sst <- data.frame(year = 1854:2021,
                                ersst = tapply(obs.n.pac.sst, as.numeric(as.character(yr)), mean),
                                observations = "ERSST")