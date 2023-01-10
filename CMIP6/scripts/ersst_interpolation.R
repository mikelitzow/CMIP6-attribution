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


## load and process ERSST ------------------------

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(2022-12-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(66)][(110):1:(250)]", "~temp")

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(2021-12-01T00:00:00Z)][(0.0):1:(0.0)][(54):1:(66)][(180):1:(230)]", "~temp")
# nc <- nc_open("./CMIP6/data/nceiErsstv5_c5fc_6a40_5e5b.nc")

# full North Pacific, 1854-2022
nc <- nc_open("./CMIP6/data/nceiErsstv5_89b2_444a_5b23.nc")

# # using AK toy example to figure out plotting
# 
# nc <- nc_open("./CMIP6/data/nceiErsstv5_0341_0869_534b.nc")

# process

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


# put SST data into data.table

# ersst <- data.table(lat = lat,
#                     lon = lon,
#                     date = rep(chron::dates(d), each = ncol(SST)),
#                     sst = c(t(SST)),
#                     index = paste("N", lat, "E", lon, sep=""))

ersst <- data.table(lat = lat,
                    lon = lon,
                    date = rep(lubridate::as_date(d), each = ncol(SST)),
                    sst = c(t(SST))
                    )


# ersst <- data.table(lat = lat,
#                     lon = lon,
#                     sst = colMeans(SST))
# 
# ersst <- data.frame(lat = lat,
#                     lon = lon,
#                     sst = colMeans(SST))


# ersst.list <- list(x = x, y = y, z = SST)

# # get CMIP6 grid
# files.new <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585")
# path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.new[i], sep="")

# # load file
# nc <- nc_open(path)
# 
# x.cmip <- ncvar_get(nc, "lon")
# y.cmip <- ncvar_get(nc, "lat")
# 
# lat.cmip <- rep(y.cmip, length(x.cmip))
# lon.cmip <- rep(x.cmip, each = length(y.cmip))
# 
# new.list <- list(x = x.cmip, y = y.cmip)
# 
# new <- interp.surface.grid(ersst.list, new.list)
# 
# image.plot(  new)
# image.plot(ersst.list)
  
# new grid
# x.out <- c(min(x), seq(110.5, 249.5, by = 1), max(x))
# y.out <- c(min(y), seq(20.5, 65.5, by = 1), max(y))

# x.out <- seq(179.5, 230.5, by = 1)
# y.out <- seq(53.5, 66.5, by = 1)

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

### plot 


# Interpolate values to a new grid
# interpolated <- ersst[, Interpolate(sst ~ lon + lat, x.out, y.out)]
library(RColorBrewer)

interpolated <- Interpolate(data = ersst, sst ~ lon + lat, x.out, y.out)


plot <- data.frame(lon = lon,
                   lat = lat,
                   sst = colMeans(SST), type = "original")

plot <- rbind(plot,
               data.frame(lon = interpolated$lon,
                   lat = interpolated$lat,
                   sst = interpolated$sst, type = "interpolated"))

ggplot(plot, aes(lon, lat), fill = sst) +
  xlim(179,231) +
  ylim(53,67) +
  geom_tile(aes(fill = sst)) +
  scale_fill_distiller(palette = "Spectral") +
  facet_wrap(~type)


plot <- data.frame(lon = interpolated$lon,
                   lat = interpolated$lat,
                   sst = interpolated$sst)

ggplot(plot, aes(lon, lat), fill = sst) +
  geom_raster(aes(fill = sst))

# land <- is.na(colMeans(SST))

temp.mean <- interpolated$sst
# temp.mean <- ersst$sst
# temp.mean <- colMeans(SST)

z <- t(matrix(temp.mean,length(y.out)))
# z <- interpolated$sst
# x <- interpolated$lon  
# y <- interpolated$lat  

image.plot(x, y, z)

image.plot(x.out,y.out,z, col=oceColorsPalette(64), xlab = "", ylab = "")
contour(x.out, y.out, z, add=T)
map

interpolated <- ersst[, ersst.new := Interpolate(sst ~ lon + lat, lon, lat,
                                     data = interpolated, grid = FALSE)$ersst]


x.out <- seq(0, 360, by = 10)
y.out <- seq(-90, 0, by = 10)

interpolated <- geopotential[, Interpolate(gh | gh.t.w ~ lon + lat, x.out, y.out), 
                             by = date]

plot(ersst$lat, lat)
plot(ersst$lon, lon)





# Interpolate multiple values
geopotential[, c("u", "v") := GeostrophicWind(gh, lon, lat)]
interpolated <- geopotential[, Interpolate(u | v ~ lon + lat, x.out, y.out)]

# Interpolate values following a path
lats <- c(-34, -54, -30)   # start and end latitudes
lons <- c(302, 290, 180)   # start and end longituded
path <- geopotential[, Interpolate(gh ~ lon + lat, as.path(lons, lats))]
