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
dates <- temps <- experiment.file <- data.frame()


for(i in 1:length(files)){

  i <- 1
  
  path <- paste("./model_outputs/", files[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # nc
  
  # extract one experiment at a time and save with file name
  
  for(j in 1:2){
  
  j <- 1
  
  this.one <- data.frame(file = files[i], 
                         experiment = ncvar_get(nc, "experiment", start = j, count = 1))
  
  experiment.file <- rbind(experiment.file,
                           this.one)
  
  # extract dates

  
  d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))

  these.dates <- data.frame(this.one = d)
  dates <- rbind(dates, these.dates)

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


  # and plot 
  SST.mean <- colMeans(SST)
  z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  image(x,y,z, col=new.col)
contour(x, y, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)
