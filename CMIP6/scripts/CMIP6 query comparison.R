# explore offset in preindustrial and historical SST runs - 
# compare initial 1900-2099 NE Pacific query with 1850-2099 full N Pacific query

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

## process and plot 1850-2099 model outputs ---------------------

# loop through each file, save the file name and experiment,
# capture the time series of raw temps for area of interest for each file-experiment comparison

# objects for saving sst
preind.sst.2 <- hist.ssp585.sst.2  <- matrix()

files.2 <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585")

for(i in 1:length(files.2)){

  # i <- 1
  
  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.2[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # nc
  
  # extract one experiment at a time
  
  # get list of experiments
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  
  for(j in 1:length(experiments)){
    
    # j <- 2
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "hist_ssp585"){
      
      # extract spatial data
      x <- ncvar_get(nc, "lon")
      y <- ncvar_get(nc, "lat")
      
      # Keep track of corresponding latitudes and longitudes of each column:
      lat <- rep(y, length(x))   
      lon <- rep(x, each = length(y)) 
      
      SST <- ncvar_get(nc, "tos", verbose = F, start = c(j,1,1,1), count = c(1,-1,-1,-1))
      
      # Change data from a 3-D array to a matrix of monthly data by grid point:
      # First, reverse order of dimensions ("transpose" array)
      SST <- aperm(SST, 3:1)
      
      # Change to matrix with column for each grid point, rows for monthly means
      SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3])) 
      
      # remove  values south of 20N
      drop <- lat < 20
      
      SST[,drop] <- NA

      these.temps <- data.frame(xx = rowMeans(SST, na.rm = T))
      names(these.temps) <-  str_remove(files.2[i], ".nc")
      
      hist.ssp585.sst.2 <- cbind(hist.ssp585.sst.2, these.temps)
      
    } # close if 
    
    
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "piControl"){
  
      # extract spatial data
      x <- ncvar_get(nc, "lon")
      y <- ncvar_get(nc, "lat")
      
      # Keep track of corresponding latitudes and longitudes of each column:
      lat <- rep(y, length(x))   
      lon <- rep(x, each = length(y)) 
      
      SST <- ncvar_get(nc, "tos", verbose = F, start = c(j,1,1,1), count = c(1,-1,-1,-1))
      
      # Change data from a 3-D array to a matrix of monthly data by grid point:
      # First, reverse order of dimensions ("transpose" array)
      SST <- aperm(SST, 3:1)
      
      # Change to matrix with column for each grid point, rows for monthly means
      SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3])) 
      
      # remove  values south of 20N
      drop <- lat < 20
      
      SST[,drop] <- NA
      
      these.temps <- data.frame(xx = rowMeans(SST, na.rm = T))
      names(these.temps) <-  str_remove(files.2[i], ".nc")
      
      preind.sst.2 <- cbind(preind.sst.2, these.temps)
      
    } # close if 
    
  } # close j
} # close i

# remove leading column of NAs
piControl.temps.2 <- preind.sst.2[,-1]
hist.585.temps.2 <- hist.ssp585.sst.2[,-1]


## now the same for original 1900-2099 query----------------------

preind.sst.1 <- hist.ssp585.sst.1  <- matrix()

files.1.all <- list.files("./CMIP6/CMIP6_outputs/1900-2099_runs/")

# drop ssp245
drop <- grep("245", files.1.all)

files.1 <- files.1.all[-drop]

for(i in 1:length(files.1)){
  
  # i <- 1
  
  path <- paste("./CMIP6/CMIP6_outputs/1900-2099_runs/", files.1[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # nc
  
  # extract one experiment at a time
  
  # get list of experiments
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  
  for(j in 1:length(experiments)){
    
    # j <- 2
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "hist_ssp585"){
      
      # extract spatial data
      x <- ncvar_get(nc, "lon")
      y <- ncvar_get(nc, "lat")
      
      # Keep track of corresponding latitudes and longitudes of each column:
      lat <- rep(y, length(x))   
      lon <- rep(x, each = length(y)) 
      
      SST <- ncvar_get(nc, "tos", verbose = F, start = c(j,1,1,1), count = c(1,-1,-1,-1))
      
      # Change data from a 3-D array to a matrix of monthly data by grid point:
      # First, reverse order of dimensions ("transpose" array)
      SST <- aperm(SST, 3:1)
      
      # Change to matrix with column for each grid point, rows for monthly means
      SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3])) 
      
      # remove  values south of 20N
      drop <- lat < 20
      
      SST[,drop] <- NA
      
      these.temps <- data.frame(xx = rowMeans(SST, na.rm = T))
      names(these.temps) <-  str_remove(files.1[i], ".nc")
      
      hist.ssp585.sst.1 <- cbind(hist.ssp585.sst.1, these.temps)
      
    } # close if 
    
    
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "piControl"){
      
      # extract spatial data
      x <- ncvar_get(nc, "lon")
      y <- ncvar_get(nc, "lat")
      
      # Keep track of corresponding latitudes and longitudes of each column:
      lat <- rep(y, length(x))   
      lon <- rep(x, each = length(y)) 
      
      SST <- ncvar_get(nc, "tos", verbose = F, start = c(j,1,1,1), count = c(1,-1,-1,-1))
      
      # Change data from a 3-D array to a matrix of monthly data by grid point:
      # First, reverse order of dimensions ("transpose" array)
      SST <- aperm(SST, 3:1)
      
      # Change to matrix with column for each grid point, rows for monthly means
      SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3])) 
      
      # remove  values south of 20N
      drop <- lat < 20
      
      SST[,drop] <- NA
      
      these.temps <- data.frame(xx = rowMeans(SST, na.rm = T))
      names(these.temps) <-  str_remove(files.1[i], ".nc")
      
      preind.sst.1 <- cbind(preind.sst.1, these.temps)
      
    } # close if 
    
  } # close j
} # close i

# remove leading column of NAs
piControl.temps.1 <- preind.sst.1[,-1]
hist.585.temps.1 <- hist.ssp585.sst.1[,-1]


## process and combine the two query sets -------------------------------

# get a year vector for each time domain

# first for 1850-2099
path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.2[1], sep="")
nc <- nc_open(path)
d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
years.2 <- as.numeric(as.character(years(d)))

# now 1900-2099
path <- paste("./CMIP6/CMIP6_outputs/1900-2099_runs/", files.1[i], sep="")
nc <- nc_open(path)
d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
years.1 <- as.numeric(as.character(years(d)))

# get annual means
ff <- function(x) tapply(x, years.1, mean)

hist.585.annual.1 <- as.data.frame(apply(hist.585.temps.1, 2, ff))
piControl.annual.1 <- as.data.frame(apply(piControl.temps.1, 2, ff))

ff <- function(x) tapply(x, years.2, mean)

hist.585.annual.2 <- as.data.frame(apply(hist.585.temps.2, 2, ff))
piControl.annual.2 <- as.data.frame(apply(piControl.temps.2, 2, ff))

# change Kelvin to Celsius
change <- colMeans(hist.585.annual.1) > 200
hist.585.annual.1[,change] <- hist.585.annual.1[,change] - 273.15

change <- colMeans(hist.585.annual.2) > 200
hist.585.annual.2[,change] <- hist.585.annual.2[,change] - 273.15

change <- colMeans(piControl.annual.1) > 200
piControl.annual.1[,change] <- piControl.annual.1[,change] - 273.15

change <- colMeans(piControl.annual.2) > 200
piControl.annual.2[,change] <- piControl.annual.2[,change] - 273.15

## combine each and plot

# create year column
hist.585.annual.1$year <- piControl.annual.1$year <- 1900:2099
hist.585.annual.2$year <- piControl.annual.2$year <- 1850:2099

# and experiment column
piControl.annual.1$experiment <- piControl.annual.2$experiment <- "Preindustrial"
hist.585.annual.1$experiment <- hist.585.annual.2$experiment <- "Historical/ssp585"

piControl.annual.1 <- piControl.annual.1 %>%
  pivot_longer(cols = c(-year, -experiment))

piControl.annual.2 <- piControl.annual.2 %>%
  pivot_longer(cols = c(-year, -experiment))

hist.585.annual.1 <- hist.585.annual.1 %>%
  pivot_longer(cols = c(-year, -experiment))

hist.585.annual.2 <- hist.585.annual.2 %>%
  pivot_longer(cols = c(-year, -experiment))

plot.query.1 <- rbind(piControl.annual.1, hist.585.annual.1)

ggplot(plot.query.1, aes(year, value, color = experiment)) +
  geom_line() +
  facet_wrap(~name) +
  ggtitle("1900-2099 NE Pacific query")

ggsave("./CMIP6/figs/1900-2099_NE_Pacific_query_preindustrial_hist-ssp585_comparison.png", width = 10, height = 8, units = 'in')

plot.query.2 <- rbind(piControl.annual.2, hist.585.annual.2)

ggplot(plot.query.2, aes(year, value, color = experiment)) +
  geom_line() +
  facet_wrap(~name) +
  ggtitle("1850-2099 N Pacific query")

ggsave("./CMIP6/figs/1850-2099_NE_Pacific_query_preindustrial_hist-ssp585_comparison.png", width = 10, height = 8, units = 'in')
