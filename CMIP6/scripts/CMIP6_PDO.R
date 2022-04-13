# subset CMIP6 by region and evaluate models

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)
library(FactoMineR)

# set palette
new.col <- oceColorsPalette(64)

# set theme
theme_set(theme_bw())

## set up processing  ---------------------

# list of files
files.new <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585")

# make function to calculate cell weights
cell.weight <- function(x)  sqrt(cos(x*pi/180))


## begin with full North Pacific grid --------------------------

# create blank df to hold EOF loadings and PC scores
CMIP6.EOF1 <- CMIP6.PC1 <- data.frame()

for(i in 1:length(files.new)){ # start i loop (each CMIP6 model)

  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.new[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # extract dates

  d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
  m <- months(d)
  yr <- years(d)
  
  # extract lat / long
  
  x <- ncvar_get(nc, "lon")
  y <- ncvar_get(nc, "lat")
  
  # extract experiments and choose historical
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  experiment.keep <- grep("hist_ssp585", experiments)
  
  # extract SST

  SST <- ncvar_get(nc, "tos", verbose = F, start = c(experiment.keep,1,1,1), count = c(1,-1,-1,-1))
  
  SST <- aperm(SST, 3:1)
  
  SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))
  
  # check for values in degrees K and convert if needed
  if(max(colMeans(SST), na.rm = T) > 200) {
    
    SST <- SST - 273.15
    
  }
  
  
  # Keep track of corresponding latitudes and longitudes of each column:
  lat <- rep(y, length(x))
  lon <- rep(x, each = length(y))
  dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
  
  # 1950 - 2014
  SST <- SST[yr %in% 1950:2014, ]
  m <- m[yr %in% 1950:2014]
  yr <- yr[yr %in% 1950:2014]

  
  # get anomalies (remove seasonal signal)
  f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
  mu <- apply(SST, 2, f)	# compute monthly means for each time series (cell)
  mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location
  
  sst.anom <- SST - mu
  
  # detrend each cell
  detrend_steps <- 1:nrow(SST) # time steps for detrending
  detrend <- function(x) (x - lm(x ~ detrend_steps)$fitted.values)
  
  # identify land
  land <- is.na(colMeans(SST))
  
  # and detrend non-land cells
  detrended.sst <- apply(sst.anom[,!land], 2, detrend)
  
  # and cell weight for this model 
  weights <- cell.weight(lat[!land])
  
  # calculate EOF  
  EOF <- svd.triplet(cov(detrended.sst), col.w=weights) #weighting the columns
  pc1 <- as.matrix(detrended.sst) %*% EOF$U[,1]
  
  # add to data frames
  CMIP6.EOF1 <- rbind(CMIP6.EOF1, 
                      data.frame(model = str_remove(files.new[i], ".nc"),
                                 lat = lat[!land],
                                 lon = lon[!land],
                                 loading = scale(EOF$U[,1])))
    
  CMIP6.PC1 <- rbind(CMIP6.PC1, 
                        data.frame(model = str_remove(files.new[i], ".nc"),
                                   year = yr,
                                   month = m,
                                   pc1 = scale(pc1)))
  
} # close i loop (models)
  


## save ---------------------------------
write.csv(CMIP6.EOF1, "./CMIP6/summaries/CMIP6_EOF1_eigenvalues.csv", row.names = F)
write.csv(CMIP6.PC1, "./CMIP6/summaries/CMIP6_PC1.csv", row.names = F)


CMIP6.EOF1 <- read.csv("./CMIP6/summaries/CMIP6_EOF1_eigenvalues.csv")
CMIP6.PC1 <- read.csv("./CMIP6/summaries/CMIP6_PC1.csv")

## plot CMIP6 loadings ---------------------------------

# get SST dimensions again - for plotting
files.new <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585")
path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.new[1], sep="")

# load file
nc <- nc_open(path)

# extract dates
d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
m <- months(d)
yr <- years(d)

# extract lat / long
x <- ncvar_get(nc, "lon")
y <- ncvar_get(nc, "lat")

# extract experiments and choose historical
experiments <-  ncvar_get(nc, "experiment", verbose = F)
experiment.keep <- grep("hist_ssp585", experiments)

# extract SST
SST <- ncvar_get(nc, "tos", verbose = F, start = c(experiment.keep,1,1,1), count = c(1,-1,-1,-1))

SST <- aperm(SST, 3:1)

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))
lon <- rep(x, each = length(y))
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# identify land
land <- is.na(colMeans(SST))

# get range for plotting
limits <- CMIP6.EOF1 %>%
  group_by(model) %>%
  summarise(min = abs(min(loading)),
            max = abs(max(loading)))

limits$lim.1 <- limits$lim.2 <- NA

for(i in 1:nrow(limits)){
  
limits$lim.1[i] <- -(max(limits[i,2:3]))

limits$lim.2[i] <- max(limits[i,2:3])

}

plot.limits <- limits %>%
  select(-min, -max)

plot.dat <- CMIP6.EOF1 %>%
  left_join(., plot.limits)




models <- unique(CMIP6.EOF1$model)

png("./CMIP6/figs/CMIP6_and_ERSST_EOF1_loadings.png", width = 7, height = 11, units = 'in', res = 300)
par(mfrow = c(6, 4), mar = c(0.5, 0.5, 1, 0.5))

for(i in 1:length(models)){

temp.plot <- plot.dat %>%
  filter(model == models[i])
  
  z <- rep(NA, ncol(SST))
  z[!land] <- temp.plot$loading
  z <- t(matrix(z, length(y))) 
  image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", 
        zlim=c(mean(temp.plot$lim.1), mean(temp.plot$lim.2)))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
  contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
  map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia', 'Philippines'), 
      fill=T,add=T, lwd=1, col="lightyellow3")
  mtext(temp.plot$model[1], cex=0.8)

  
}
  # # plot time series
  # plot.pc <- data.frame(dec.yr = as.numeric(as.character(yr)) + 
  #                         (as.numeric(m) - 0.5)/12,
  #                       pc1 = pc1)
  # 
  # ggplot(plot.pc, aes(dec.yr, pc1)) +
  #   geom_line()
  
 
  
  ## add ERSST EOF1 -------------------
  
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
  
  # temp.mean <- colMeans(SST, na.rm=T)
  # z <- t(matrix(temp.mean,length(y)))
  # image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "")
  # contour(x, y, z, add=T)
  # map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")
  # 
  # 1950 - 2014
  SST <- SST[yr %in% 1950:2014, ]
  m <- m[yr %in% 1950:2014]
  yr <- yr[yr %in% 1950:2014]
  
  
  # get anomalies (remove seasonal signal)
  f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
  mu <- apply(SST, 2, f)	# compute monthly means for each time series (cell)
  mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location
  
  sst.anom <- SST - mu
  
  # detrend each cell
  detrend_steps <- 1:nrow(SST) # time steps for detrending
  detrend <- function(x) (x - lm(x ~ detrend_steps)$fitted.values)
  
  # identify land
  land <- is.na(colMeans(SST))
  
  # and detrend non-land cells
  detrended.sst <- apply(sst.anom[,!land], 2, detrend)
  
  # and cell weight for this model 
  weights <- cell.weight(lat[!land])
  
  # calculate EOF  
  EOF <- svd.triplet(cov(detrended.sst), col.w=weights) #weighting the columns
  pc1 <- as.matrix(detrended.sst) %*% EOF$U[,1]
  
  # plot loadings
  lim <- range(EOF$U[,1])
  
  # png("./CMIP6/figs/ERSST_EOF1_loading.png", width = 6, height = 4, units = 'in', res = 300)
  # par(mar = c(1,1,1,1))
  z <- rep(NA, ncol(SST))
  z[!land] <- EOF$U[,1]
  z <- t(matrix(z, length(y))) 
  image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
  contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
  map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia', 'Philippines'), 
      fill=T,add=T, lwd=1, col="lightyellow3") 
  mtext("ERSST", cex=0.8)
  
  dev.off()