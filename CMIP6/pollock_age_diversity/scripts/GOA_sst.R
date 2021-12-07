# extract and summarize ERSST sst data

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

# load ERSST - 1960 through 2020

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1900-01-01):1:(2020-12-01T00:00:00Z)][(0.0):1:(0.0)][(52):1:(62)][(198):1:(226)]", "~temp")

# load and process SST data

nc <- nc_open("./data/nceiErsstv5_bf33_3a9d_7ad9.nc")

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


# need to drop Bristol Bay cells and trim southern extent
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

# and annual means
ann.sst <- tapply(obs.sst, as.numeric(as.character(yr)), mean)


# combine annual, 2-yr running mean, and 3-yr running mean SST
goa.sst <- data.frame(year = 1900:2020,
                      annual.sst = ann.sst,
                      two.yr.sst = rollmean(ann.sst, 2, fill = NA, align = "right"),
                      three.yr.sst = rollmean(ann.sst, 3, fill = NA, align = "right"))

plot.check <- goa.sst %>% 
  pivot_longer(cols = - year)

ggplot(plot.check, aes(year, value, color = name)) +
  geom_line()


# save
write.csv(goa.sst, "./data/goa.sst.csv")

# and annual anomalies wrt 1900-1999
goa.sst.anom <- goa.sst 

for(j in 2:4){
mean.20.century <- mean(goa.sst[goa.sst$year %in% 1900:1999,j], na.rm = T)
sd.20.century <- sd(goa.sst[goa.sst$year %in% 1900:1999,j], na.rm = T)
goa.sst.anom[,j] <- (goa.sst.anom[,j] - mean.20.century) / sd.20.century
}

# save
write.csv(goa.sst.anom, "./data/goa.sst.anom.csv")
