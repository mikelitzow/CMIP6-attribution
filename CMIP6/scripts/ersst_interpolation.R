# interpolate ERSST onto CMIP6 grid

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)
library(data.table)
library(metR)

# set palettes
new.col <- oceColorsPalette(64)
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set theme
theme_set(theme_bw())

## interpolate ERSST to CMIP6 grid ------------------------

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(2022-12-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(66)][(110):1:(250)]", "~temp")

# full North Pacific, 1854-2022
nc <- nc_open("./CMIP6/data/nceiErsstv5_89b2_444a_5b23.nc")

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

SST <- aperm(SST, 3:1)

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))
lon <- rep(x, each = length(y))
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check!
temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "")
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

ersst <- data.table(lat = lat,
                    lon = lon,
                    date = rep(lubridate::as_date(d), each = ncol(SST)),
                    sst = c(t(SST))
                    )

x.out <- seq(110.5, 249.5, by = 1)
y.out <- seq(20.5, 65.5, by = 1)

# Interpolate values to a new grid
interpolated <- ersst[, Interpolate(sst ~ lat + lon, y.out, x.out), by = date]

# compare raw and interpolated side by side

# get overall mean at each cell!
plot_interp <- interpolated %>%
  mutate(index = paste(lon, lat, sep = "-")) %>%
  group_by(index) %>%
  summarise(sst = mean(sst))
  
png("./CMIP6/figs/raw and interpolated ERSST.png", width = 8, height = 4, units = "in", res = 300)
par(mfrow = c(1,2))

temp.mean <- colMeans(SST)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", main = "Raw")
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")


# temp.mean <- interpolated$sst[interpolated$date == interpolated$date[1]]
temp.mean <- plot_interp$sst
z <- t(matrix(temp.mean,length(y.out)))
image.plot(x.out,y.out,z, col=oceColorsPalette(64), xlab = "", ylab = "", main = "Interpolated")
contour(x.out, y.out, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

dev.off()

### plot regions of interest-------------------

# define polygons

# EBS polygon 
ebs.x <- c(183, 183, 203, 203, 191) 
ebs.y <- c(53, 65, 65, 57.5, 53)

# GOA polygon
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)

# BC polygon
bc.x <- c(231, 238, 228, 225, 225)
bc.y <- c(55, 49, 49, 53, 55)

# northern CCE polygon
ncc.x <- c(238, 238, 233, 233, 238)
ncc.y <- c(49, 41, 41, 49, 49)

# southern CCE polygon
scc.x <- c(238, 243, 237, 233, 233, 238)
scc.y <- c(41, 33, 33, 39, 41, 41)


png("./CMIP6/figs/raw and interpolated ERSST NE Pacific regions.png", width = 8, height = 4, units = "in", res = 300)
par(mfrow = c(1,2))

temp.mean <- colMeans(SST)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", main = "Raw", xlim = c(170,245), ylim = c(30,66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")
polygon(ebs.x, ebs.y, border = "red", lwd = 2)
polygon(goa.x, goa.y, border = "red", lwd = 2)
polygon(bc.x, bc.y, border = "red", lwd = 2)
polygon(ncc.x, ncc.y, border = "red", lwd = 2)
polygon(scc.x, scc.y, border = "red", lwd = 2)

# temp.mean <- interpolated$sst[interpolated$date == interpolated$date[1]]
temp.mean <- plot_interp$sst
z <- t(matrix(temp.mean,length(y.out)))
image.plot(x.out,y.out,z, col=oceColorsPalette(64), xlab = "", ylab = "", main = "Interpolated", xlim = c(170,245), ylim = c(30,66))
contour(x.out, y.out, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")
polygon(ebs.x, ebs.y, border = "red", lwd = 2)
polygon(goa.x, goa.y, border = "red", lwd = 2)
polygon(bc.x, bc.y, border = "red", lwd = 2)
polygon(ncc.x, ncc.y, border = "red", lwd = 2)
polygon(scc.x, scc.y, border = "red", lwd = 2)

dev.off()


## save interpolated ERSST data set

write.csv(interpolated, "./CMIP6/data/interpolated_ERSST.csv", row.names = F)
