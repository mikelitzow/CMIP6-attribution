# subset CMIP6 by region and evaluate models

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

## set up processing  ---------------------

# load polygons for subsetting each model
regional.polygons <- read.csv("./CMIP6/summaries/regional_polygons.csv")

# get vector of regions for indexing
regions <- unique(regional.polygons$region)

# loop through each file, save the file name and experiment,
# capture the time series of raw temps for area of interest for each file-experiment comparison

# objects for saving dates and temps and file-experiment ID from each realization
dates <- temps  <- matrix()
experiment.file <- data.frame()
files.new <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585")

length(files.new)

# make function to calculate cell weights
cell.weight <- function(x)  sqrt(cos(x*pi/180))

# and function to calculate weighted means with these weights 
weighted.cell.mean <- function(x) weighted.mean(x, weights, na.rm = T)

# load vector of clean region names to match with ERSST summaries
region.names <- read.csv("./CMIP6/summaries/clean_region_names.csv") 

## begin with full North Pacific grid --------------------------

# create blank df to hold time series of sst and anomalies wrt 1950-1999
CMIP6.sst.time.series <- CMIP6.anomaly.time.series <- data.frame()

for(i in 1:length(files.new)){ # start i loop (each CMIP6 model)

  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.new[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # extract one experiment at a time
  
  # get list of experiments
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  
  for(j in 1:length(experiments)){ # start j loop (each experiment)

  # extract dates

  d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
  m <- months(d)
  yr <- years(d)
  
  # extract lat / long
  
  x <- ncvar_get(nc, "lon")
  y <- ncvar_get(nc, "lat")

  # extract SST

  SST <- ncvar_get(nc, "tos", verbose = F, start = c(j,1,1,1), count = c(1,-1,-1,-1))
  
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
  
  # and cell weight for this model 
  weights <- cell.weight(lat)
  
  # extract the following time series for each regional subset:
  # annual unsmoothed, annual.two.yr.running.mean, annual.three.yr.running.mean
  # winter.unsmoothed, winter.two.yr.running.mean, winter.three.yr.running.mean
  
  # begin with full model grid 1850-2099
  
  # calculate monthly mean temp weighted by area  
  temp.monthly.sst <- apply(SST, 1, weighted.cell.mean) 
  
  # calculate annual means
  temp.annual <- tapply(temp.monthly.sst, yr, mean) # again, use full time series for entire North Pacific field
  
  temp.2yr <- rollmean(temp.annual, 2, fill = NA, align = "left") # for salmon - year of and year after ocean entry
  
  temp.3yr <- rollmean(temp.annual, 3, fill = NA, align = "center") # for salmon - year before, year of, and year after ocean entry
  
  # calculate winter means
  
  # first, define winter year
  winter.year <- if_else(m %in% c("Nov", "Dec"), as.numeric(as.character(yr))+1, as.numeric(as.character(yr)))
  winter.year <- winter.year[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
  
  # now select only winter sst data
  temp.winter.monthly.sst <- temp.monthly.sst[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
  
  # and calculate annual winter means
  temp.winter <- tapply(temp.winter.monthly.sst, winter.year, mean)
  
  # change leading and trailing years to NA b/c these are incomplete
  leading <- min(names(temp.winter))
  trailing <- max(names(temp.winter))
  
  temp.winter[names(temp.winter) %in% c(leading, trailing)] <- NA
  
  temp.winter.2yr <- rollmean(temp.winter, 2, fill = NA, align = "left")
  
  temp.winter.3yr <- rollmean(temp.winter, 3, fill = NA, align = "center")
  
  # combine into data frame of time series for each model
  CMIP6.sst.time.series <- rbind(CMIP6.sst.time.series,
                            data.frame(region = region.names[1,1], # fixed at "North_Pacific" for this first summary
                                       model = str_remove(files.new[i], ".nc"),
                                       experiment = experiments[j],
                                       year = 1850:2099,
                                       annual.unsmoothed = temp.annual,
                                       annual.two.yr.running.mean = temp.2yr,
                                       annual.three.yr.running.mean = temp.3yr,
                                       winter.unsmoothed = temp.winter[names(temp.winter) %in% 1850:2099], # need to drop winter 2100
                                       winter.two.yr.running.mean = temp.winter.2yr[names(temp.winter) %in% 1850:2099],
                                       winter.three.yr.running.mean = temp.winter.3yr[names(temp.winter) %in% 1850:2099]))
  
  ## now calculate the data as anomalies wrt 1850:1949
  # calculate annual anomalies
  annual.climatology.mean <- mean(temp.annual[names(temp.annual) %in% 1850:1949])
  
  annual.climatology.sd <- sd(temp.annual[names(temp.annual) %in% 1850:1949])
  
  temp.annual.anom <- (temp.annual - annual.climatology.mean) / annual.climatology.sd
  
  temp.anom.2yr <- rollmean(temp.annual.anom, 2, fill = NA, align = "left") # for salmon - year of and year after ocean entry
  
  temp.anom.3yr <- rollmean(temp.annual.anom, 3, fill = NA, align = "center") # for salmon - year before, year of, and year after ocean entry
  
  # calculate winter anomalies
  winter.climatology.mean <- mean(temp.winter[names(temp.winter) %in% 1850:1949], na.rm = T)
  
  winter.climatology.sd <- sd(temp.winter[names(temp.winter) %in% 1850:1949], na.rm = T)
  
  temp.winter.anom <- (temp.winter - winter.climatology.mean) / winter.climatology.sd
  
  temp.winter.anom.2yr <- rollmean(temp.winter.anom, 2, fill = NA, align = "left")
  
  temp.winter.anom.3yr <- rollmean(temp.winter.anom, 3, fill = NA, align = "center")
  
  # combine into data frame of anomalies
  CMIP6.anomaly.time.series <- rbind(CMIP6.anomaly.time.series,
                                    data.frame(region = region.names[1,1], # fixed at "North_Pacific" for this first summary
                                               model = str_remove(files.new[i], ".nc"),
                                               experiment = experiments[j],
                                               year = 1850:2099,
                                               annual.unsmoothed = temp.annual.anom,
                                               annual.two.yr.running.mean = temp.anom.2yr,
                                               annual.three.yr.running.mean = temp.anom.3yr,
                                               winter.unsmoothed = temp.winter.anom[names(temp.winter.anom) %in% 1850:2099], # need to drop winter 2100
                                               winter.two.yr.running.mean = temp.winter.anom.2yr[names(temp.winter.anom) %in% 1850:2099],
                                               winter.three.yr.running.mean = temp.winter.anom.3yr[names(temp.winter.anom) %in% 1850:2099]))


  ## query each of the regional time series
      
  for(k in 1:length(unique(regional.polygons$region))){ # open k loop (smaller regions)

    # keep track of progress
    print(paste(i, j, k, sep = "-"))
    
    # define regional polygon
    temp.polygon <- regional.polygons %>%
      filter(region == unique(regional.polygons$region)[k])
    
    # select sst cells in this polygon; all others become NA
    xp <- cbind(temp.polygon$x, temp.polygon$y)
    loc <- cbind(lon, lat)
    check <- in.poly(loc, xp=xp)
    
    temp.sst <- SST
    temp.sst[,!check] <- NA
    
    # # plot to check
    # temp.mean <- colMeans(temp.sst, na.rm=T)
    # z <- t(matrix(temp.mean,length(y)))
    # image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(210, 245), ylim = c(30, 55))
    # contour(x, y, z, add=T)
    # map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")
    
    # process time series for this polygon
    
    # note that we will include full time series for 1850-2099 
    # b/c we don't have the same data concerns as for some pre-1950 ERSST data

    ## now create time series of annual means
    ## unsmoothed, and two- and three-year running means
      
    # calculate monthly mean temp weighted by area 
    temp.monthly.sst <- apply(temp.sst, 1, weighted.cell.mean) 
      
    # calculate annual means
    temp.annual <- tapply(temp.monthly.sst, yr, mean) 
      
    temp.2yr <- rollmean(temp.annual, 2, fill = NA, align = "left") # for salmon - year of and year after ocean entry
      
    temp.3yr <- rollmean(temp.annual, 3, fill = NA, align = "center") # for salmon - year before, year of, and year after ocean entry
      
    # calculate winter means
    temp.winter.monthly.sst <- temp.monthly.sst[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
      
    temp.winter <- tapply(temp.winter.monthly.sst, winter.year, mean)
      
    # change leading and trailing years to NA b/c these are incomplete
    leading <- min(names(temp.winter))
    trailing <- max(names(temp.winter))
      
    temp.winter[names(temp.winter) %in% c(leading, trailing)] <- NA
      
    temp.winter.2yr <- rollmean(temp.winter, 2, fill = NA, align = "left")
      
    temp.winter.3yr <- rollmean(temp.winter, 3, fill = NA, align = "center")
      
    # combine into data frame of time series by region
    CMIP6.sst.time.series <- rbind(CMIP6.sst.time.series,
                                   data.frame(region = region.names[(k+1),1], # advance k by one to account for full North Pac region (done separately above)
                                              model = str_remove(files.new[i], ".nc"),
                                              experiment = experiments[j],
                                              year = 1850:2099,
                                              annual.unsmoothed = temp.annual,
                                              annual.two.yr.running.mean = temp.2yr,
                                              annual.three.yr.running.mean = temp.3yr,
                                              winter.unsmoothed = temp.winter[names(temp.winter) %in% 1850:2099], # need to drop winter 2100
                                              winter.two.yr.running.mean = temp.winter.2yr[names(temp.winter) %in% 1850:2099],
                                              winter.three.yr.running.mean = temp.winter.3yr[names(temp.winter) %in% 1850:2099]))
      
    ## now calculate the data as anomalies wrt 1950-1999
    # calculate annual anomalies
    annual.climatology.mean <- mean(temp.annual[names(temp.annual) %in% 1850:1949])
      
    annual.climatology.sd <- sd(temp.annual[names(temp.annual) %in% 1850:1949])
      
    temp.annual.anom <- (temp.annual - annual.climatology.mean) / annual.climatology.sd
      
    temp.anom.2yr <- rollmean(temp.annual.anom, 2, fill = NA, align = "left") # for salmon - year of and year after ocean entry
      
    temp.anom.3yr <- rollmean(temp.annual.anom, 3, fill = NA, align = "center") # for salmon - year before, year of, and year after ocean entry
      
    # calculate winter anomalies
    winter.climatology.mean <- mean(temp.winter[names(temp.winter) %in% 1850:1949], na.rm = T)
      
    winter.climatology.sd <- sd(temp.winter[names(temp.winter) %in% 1850:1949], na.rm = T)
      
    temp.winter.anom <- (temp.winter - winter.climatology.mean) / winter.climatology.sd
      
    temp.winter.anom.2yr <- rollmean(temp.winter.anom, 2, fill = NA, align = "left")
      
    temp.winter.anom.3yr <- rollmean(temp.winter.anom, 3, fill = NA, align = "center")
      
    # combine into data frame of time series by region
    CMIP6.anomaly.time.series <- rbind(CMIP6.anomaly.time.series,
                                       data.frame(region = region.names[(k+1),1],
                                                  model = str_remove(files.new[i], ".nc"),
                                                  experiment = experiments[j],
                                                  year = 1850:2099,
                                                  annual.unsmoothed = temp.annual.anom,
                                                  annual.two.yr.running.mean = temp.anom.2yr,
                                                  annual.three.yr.running.mean = temp.anom.3yr,
                                                  winter.unsmoothed = temp.winter.anom[names(temp.winter.anom) %in% 1850:2099], # need to drop winter 2100
                                                  winter.two.yr.running.mean = temp.winter.anom.2yr[names(temp.winter.anom) %in% 1850:2099],
                                                  winter.three.yr.running.mean = temp.winter.anom.3yr[names(temp.winter.anom) %in% 1850:2099]))
    
    } # close k loop (smaller regions)
  
  } # close j loop (experiments)
  
} # close i loop (models)

# save summaries
write.csv(CMIP6.sst.time.series, "./CMIP6/summaries/CMIP6.sst.time.series.csv", row.names = F)
write.csv(CMIP6.anomaly.time.series, "./CMIP6/summaries/CMIP6.anomaly.time.series.csv", row.names = F)
