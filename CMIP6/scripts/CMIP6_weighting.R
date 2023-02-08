# weight CMIP6 by climatology, change in temperature, annual SD, and annual autocorrelation

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
# make function to calculate cell weights
cell.weight <- function(x)  sqrt(cos(x*pi/180))

# load interpolated ERSST data
ersst <- read.csv("./CMIP6/data/interpolated_ERSST.csv")

ersst <- ersst %>%
  mutate(year = lubridate::year(date)) %>%
  filter(year %in% 1950:2014) %>%
  mutate(lat_long = paste("N", lat, "_E", lon, sep = "")) %>%
  group_by(lat_long) %>%
  summarize(ersst_mean_sst = mean(sst, na.rm = T)) %>%
  mutate(lat = as.numeric(str_sub(lat_long, start = 2, end = 5)),
         long = as.numeric(str_sub(lat_long, start = 8, end = 12)),
         weight = cell.weight(lat))


# loop through each file, save the file name and experiment,
# subset by region, and calculate RMSE 1900-1949 climatology

# objects for saving dates and temps and file-experiment ID from each realization
climatology_rmse <- data.frame()

files.new <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585")

length(files.new)

# load polygons for subsetting each model
regional.polygons <- read.csv("./CMIP6/summaries/regional_polygons.csv")

# clean up names
regional.polygons <- regional.polygons %>%
  mutate(region = case_when(
    region == "ebs" ~ "Eastern Bering Sea",
    region == "goa" ~ "Gulf of Alaska",
    region == "bc" ~ "British Columbia Coast",
    region == "ncc" ~ "Northern California Current",
    region == "scc" ~ "Southern California Current"
  ))

# and function to calculate weighted means with these weights 
weighted.cell.mean <- function(x) weighted.mean(x, weights, na.rm = T)

## loop through model-region combinations for climatology---------------

for(i in 1:length(files.new)){ # start i loop (each CMIP6 model)
# i <- 1
  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.new[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # extract one experiment at a time
  
  # get list of experiments and identify which is historical / ssp585
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  experiment_use <- str_which(experiments, "hist_ssp585")
  
  # extract dates
  d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
  m <- months(d)
  yr <- years(d)
  
  # extract lat / long
  x <- ncvar_get(nc, "lon")
  y <- ncvar_get(nc, "lat")

  # extract SST
  SST <- ncvar_get(nc, "tos", verbose = F, start = c(experiment_use,1,1,1), count = c(1,-1,-1,-1))
  
  SST <- aperm(SST, 3:1)
  
  SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))
  
  # check for values in degrees K and convert if needed
  if(max(colMeans(SST), na.rm = T) > 200) {
    
    SST <- SST - 273.15
    
  }
  
  # Keep track of corresponding latitudes and longitudes of each column:
  lat <- rep(y, length(x))
  lon <- rep(x, each = length(y))
  dimnames(SST) <- list(as.character(d), paste("N", lat, "_E", lon, sep=""))
  
  # and cell weight for this model 
  weights <- cell.weight(lat)
  
  # restrict to 1950:2014
  use <- yr %in% 1950:2014
  
  SST <- SST[use,]
  
  cmip.temp <- data.frame(lat_long = colnames(SST),
                          cmip_mean_sst = colMeans(SST))
  
  # cmip.temp <- cmip.temp %>%
  #   mutate(lat = as.numeric(str_sub(lat_long, start = 2, end = 5)),
  #                                long = as.numeric(str_sub(lat_long, start = 8, end = 12))) 
  
  compare <- left_join(ersst, cmip.temp)
  
  # get sample size
  sample <- compare %>%
    filter(!is.na(ersst_mean_sst),
           !is.na(cmip_mean_sst))

  # get RMSE for North Pacific climatology
  mod <- lm(ersst_mean_sst ~ cmip_mean_sst, weights = weight, data = compare)
  
  # record results
  climatology_rmse <- rbind(climatology_rmse,
                            data.frame(region = "North Pacific",
                                       climatology = "1950-2014",
                                       cmip_model = str_remove(files.new[i], ".nc"),
                                       experiment = "hist_ssp585",
                                       metric = "bias",
                                       RMSE = sqrt(mean(mod$residuals^2)),
                                       n = nrow(sample)))
  
  ## query each of the regional time series
      
  for(k in 1:length(unique(regional.polygons$region))){ # open k loop (smaller regions)
    # k <- 1
    # keep track of progress
    print(paste(i, k, sep = "-"))
    
    # define regional polygon
    temp.polygon <- regional.polygons %>%
      filter(region == unique(regional.polygons$region)[k])
    
    # select ersst cells in this polygon; all others become NA
    xp <- cbind(temp.polygon$x, temp.polygon$y)
    loc <- cbind(compare$long, compare$lat)
    check <- in.poly(loc, xp=xp)
    
    temp.compare <- compare
    temp.compare$ersst_mean_sst[!check] <- NA
    temp.compare$cmip_mean_sst[!check] <- NA
    
    # get sample size
    sample <- temp.compare %>%
      filter(!is.na(ersst_mean_sst),
             !is.na(cmip_mean_sst))
    
    # get RMSE for climatology
    mod <- lm(ersst_mean_sst ~ cmip_mean_sst, weights = weight, data = temp.compare)
    
    # record results
    climatology_rmse <- rbind(climatology_rmse,
                              data.frame(region = unique(regional.polygons$region)[k],
                                         climatology = "1950-2014",
                                         cmip_model = str_remove(files.new[i], ".nc"),
                                         experiment = "hist_ssp585",
                                         metric = "bias",
                                         RMSE = sqrt(mean(mod$residuals^2)),
                                         n = nrow(sample)))   
} # close k loop (regions)
  } # close i loop (models)
    
## loop through model-region combinations for interannual SD and AR(1)------------------------------------------

sd_rmse <- ar_rmse <- data.frame()

# get annual means for ERSST cells
# reload interpolated ERSST data
ersst <- read.csv("./CMIP6/data/interpolated_ERSST.csv")

ersst.ann <- ersst %>%
  mutate(year = lubridate::year(date)) %>%
  filter(year %in% 1950:2014) %>%
  mutate(lat_long = paste("N", lat, "_E", lon, sep = "")) %>%
  group_by(lat_long, year) %>%
  summarize(ersst_ann_sst = tapply(sst, year, mean, na.rm = T)) %>%
  mutate(lat = as.numeric(str_sub(lat_long, start = 2, end = 5)),
         long = as.numeric(str_sub(lat_long, start = 8, end = 12)),
         weight = cell.weight(lat))

for(i in 1:length(files.new)){ # start i loop (each CMIP6 model)
  # i <- 1
  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.new[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # extract one experiment at a time
  
  # get list of experiments and identify which is historical / ssp585
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  experiment_use <- str_which(experiments, "hist_ssp585")
  
  # extract dates
  d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
  m <- months(d)
  yr <- years(d)
  
  # extract lat / long
  x <- ncvar_get(nc, "lon")
  y <- ncvar_get(nc, "lat")
  
  # extract SST
  SST <- ncvar_get(nc, "tos", verbose = F, start = c(experiment_use,1,1,1), count = c(1,-1,-1,-1))
  
  SST <- aperm(SST, 3:1)
  
  SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))
  
  # check for values in degrees K and convert if needed
  if(max(colMeans(SST), na.rm = T) > 200) {
    
    SST <- SST - 273.15
    
  }
  
  # Keep track of corresponding latitudes and longitudes of each column:
  lat <- rep(y, length(x))
  lon <- rep(x, each = length(y))
  dimnames(SST) <- list(as.character(d), paste("N", lat, "_E", lon, sep=""))
  
  # restrict to 1950:2014
  use <- yr %in% 1950:2014
  
  SST <- SST[use,]
  dates <- rep(chron::dates(row.names(SST)), each = ncol(SST))
  
  cmip.ann <- as.data.frame(SST) %>%
    pivot_longer(cols = 1:ncol(SST)) %>%
    mutate(date = dates,
           year = lubridate::year(date),
           year = if_else(year > 2020, year - 100, year)) %>%
    rename(lat_long = name,
           sst = value) %>%
    group_by(lat_long, year) %>%
    summarize(cmip_ann_sst = tapply(sst, year, mean, na.rm = T))
    
  compare <- left_join(ersst.ann, cmip.ann)
  
  
  # calculate SD
  compare_sd <- compare %>%
    group_by(lat_long) %>%
    summarize(ersst_sd = sd(ersst_ann_sst),
              cmip_sd = sd(cmip_ann_sst)) %>%
    mutate(lat = as.numeric(str_sub(lat_long, start = 2, end = 5)),
           long = as.numeric(str_sub(lat_long, start = 8, end = 12)),
           weight = cell.weight(lat))
  
  # get sample size
  sample_sd <- compare_sd %>%
    filter(!is.na(ersst_sd),
           !is.na(cmip_sd))
  
  # get RMSE for North Pacific sd
  mod <- lm(ersst_sd ~ cmip_sd, weights = weight, data = compare_sd)
  
  # record results
  sd_rmse <- rbind(sd_rmse,
                            data.frame(region = "North Pacific",
                                       climatology = "1950-2014",
                                       cmip_model = str_remove(files.new[i], ".nc"),
                                       experiment = "hist_ssp585",
                                       metric = "SD",
                                       RMSE = sqrt(mean(mod$residuals^2)),
                                       n = nrow(sample_sd)))
  
  # calculate AR(1)
  compare_ar <- na.omit(compare) %>% # need to remove NA for ar function
    group_by(lat_long) %>%
    summarize(ersst_ar = ar(as.vector(ersst_ann_sst), aic = F, order.max = 1)$ar,
              cmip_ar = ar(as.vector(cmip_ann_sst), aic = F, order.max = 1)$ar) %>%
    mutate(lat = as.numeric(str_sub(lat_long, start = 2, end = 5)),
           long = as.numeric(str_sub(lat_long, start = 8, end = 12)),
           weight = cell.weight(lat))
  
  # get sample size
  sample_ar <- compare_ar %>%
    filter(!is.na(ersst_ar),
           !is.na(cmip_ar))
  
  # get RMSE for North Pacific sd
  mod <- lm(ersst_ar ~ cmip_ar, weights = weight, data = compare_ar)
  
  # record results
  ar_rmse <- rbind(ar_rmse,
                   data.frame(region = "North Pacific",
                              climatology = "1950-2014",
                              cmip_model = str_remove(files.new[i], ".nc"),
                              experiment = "hist_ssp585",
                              metric = "AR(1)",
                              RMSE = sqrt(mean(mod$residuals^2)),
                              n = nrow(sample_ar)))
  
  ## query each of the regional time series
  
  for(k in 1:length(unique(regional.polygons$region))){ # open k loop (smaller regions)
    # k <- 1
    # keep track of progress
    print(paste(i, k, "sd", sep = "-"))
    
    # define regional polygon
    temp.polygon <- regional.polygons %>%
      filter(region == unique(regional.polygons$region)[k])
    
    # select ersst cells in this polygon; all others become NA
    xp <- cbind(temp.polygon$x, temp.polygon$y)
    loc <- cbind(compare_sd$long, compare_sd$lat)
    check <- in.poly(loc, xp=xp)
    
    temp_compare_sd <- compare_sd
    temp_compare_sd$ersst_sd[!check] <- NA
    temp_compare_sd$cmip_sd[!check] <- NA
    
    # get sample size
    sample_temp_sd <- temp_compare_sd %>%
      filter(!is.na(ersst_sd),
             !is.na(cmip_sd))
    
    # get RMSE for climatology
    mod <- lm(ersst_sd ~ cmip_sd, weights = weight, data = temp_compare_sd)
    
    # record results
    sd_rmse <- rbind(sd_rmse,
                              data.frame(region = unique(regional.polygons$region)[k],
                                         climatology = "1950-2014",
                                         cmip_model = str_remove(files.new[i], ".nc"),
                                         experiment = "hist_ssp585",
                                         metric = "SD",
                                         RMSE = sqrt(mean(mod$residuals^2)),
                                         n = nrow(sample_temp_sd)))   
    
    ## regional ar(1)
    # keep track of progress
    print(paste(i, k, "ar", sep = "-"))
    
    # select ersst cells in this polygon; all others become NA
    xp <- cbind(temp.polygon$x, temp.polygon$y)
    loc <- cbind(compare_ar$long, compare_ar$lat)
    check <- in.poly(loc, xp=xp)
    
    temp_compare_ar <- compare_ar
    temp_compare_ar$ersst_ar[!check] <- NA
    temp_compare_ar$cmip_ar[!check] <- NA
    
    # get sample size
    sample_temp_ar <- temp_compare_ar %>%
      filter(!is.na(ersst_ar),
             !is.na(cmip_ar))
    
    # get RMSE for climatology
    mod <- lm(ersst_ar ~ cmip_ar, weights = weight, data = temp_compare_ar)
    
    # record results
    ar_rmse <- rbind(ar_rmse,
                     data.frame(region = unique(regional.polygons$region)[k],
                                climatology = "1950-2014",
                                cmip_model = str_remove(files.new[i], ".nc"),
                                experiment = "hist_ssp585",
                                metric = "AR(1)",
                                RMSE = sqrt(mean(mod$residuals^2)),
                                n = nrow(sample_temp_ar)))  
    
  } # close k loop (regions)
} # close i loop (models)



## loop through model-region combinations for climatology change (1854-1949 vs. 1993:2022) ----------------
climatology_change_rmse <- data.frame()

# load interpolated ERSST data
ersst <- read.csv("./CMIP6/data/interpolated_ERSST.csv")

ersst_change <- ersst %>%
  mutate(year = lubridate::year(date)) %>%
  filter(year %in% c(1854:1949, 1993:2022)) %>%
  mutate(lat_long = paste("N", lat, "_E", lon, sep = ""),
         era = if_else(year %in% 1854:1949, "Era 1854-1949", "Era 1993-2022")) %>%
  group_by(lat_long, era) %>%
  summarize(ersst_mean_sst = mean(sst, na.rm = T)) %>%
  pivot_wider(names_from = era, values_from = ersst_mean_sst) %>%
  mutate(ersst_era_diff = `Era 1993-2022` - `Era 1854-1949`,
          lat = as.numeric(str_sub(lat_long, start = 2, end = 5)),
         long = as.numeric(str_sub(lat_long, start = 8, end = 12)),
         weight = cell.weight(lat))

# ggplot(ersst_change, aes(ersst_era_diff)) +
#   geom_histogram(bins = 50, fill = "dark grey", color = "black") +
#   theme_bw()

for(i in 1:length(files.new)){ # start i loop (each CMIP6 model)
  # i <- 1
  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files.new[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # extract one experiment at a time
  
  # get list of experiments and identify which is historical / ssp585
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  experiment_use <- str_which(experiments, "hist_ssp585")
  
  # extract dates
  d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
  m <- months(d)
  yr <- years(d)
  
  # extract lat / long
  x <- ncvar_get(nc, "lon")
  y <- ncvar_get(nc, "lat")
  
  # extract SST
  SST <- ncvar_get(nc, "tos", verbose = F, start = c(experiment_use,1,1,1), count = c(1,-1,-1,-1))
  
  SST <- aperm(SST, 3:1)
  
  SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))
  
  # check for values in degrees K and convert if needed
  if(max(colMeans(SST), na.rm = T) > 200) {
    
    SST <- SST - 273.15
    
  }
  
  # Keep track of corresponding latitudes and longitudes of each column:
  lat <- rep(y, length(x))
  lon <- rep(x, each = length(y))
  dimnames(SST) <- list(as.character(d), paste("N", lat, "_E", lon, sep=""))

  dates <- rep(d, each = ncol(SST))
  
  # # check!
  # check_year <- years(d)
  # 
  # table(check_year)
  
  cmip_era1 <- as.data.frame(SST) %>%
    pivot_longer(cols = 1:ncol(SST)) %>%
    mutate(date = dates,
           year = years(date)) %>%
    filter(year %in% 1854:1949) %>%
    rename(lat_long = name,
           sst = value) %>%
    group_by(lat_long) %>%
    summarize(cmip_era1_sst = mean(sst)) 
  
  cmip_era2 <- as.data.frame(SST) %>%
    pivot_longer(cols = 1:ncol(SST)) %>%
    mutate(date = dates,
           year = years(date)) %>%
    filter(year %in% 1993:2022) %>%
    rename(lat_long = name,
           sst = value) %>%
    group_by(lat_long) %>%
    summarize(cmip_era2_sst = mean(sst)) 
  
  cmip_change <- left_join(cmip_era1, cmip_era2) %>%
    mutate(cmip_era_diff = cmip_era2_sst - cmip_era1_sst)
  
  # ggplot(cmip_change, aes(cmip_era_diff)) +
  #   geom_histogram(bins = 50, fill = "dark grey", color = "black") +
  #   theme_bw()

  compare <- left_join(ersst_change, cmip_change) %>%
    select(lat_long, lat, long, ersst_era_diff, weight, cmip_era_diff)
  
  ggplot(compare, aes(cmip_era_diff, ersst_era_diff)) +
    geom_point() +
    theme_bw() + 
    geom_smooth(method = "lm", se = F)
  
  # get sample size
  sample <- compare %>%
    filter(!is.na(ersst_era_diff),
           !is.na(cmip_era_diff))
  
  # get RMSE for North Pacific climatology
  mod <- lm(ersst_era_diff ~ cmip_era_diff, weights = weight, data = compare)
  
  # record results
  climatology_change_rmse <- rbind(climatology_change_rmse,
                            data.frame(region = "North Pacific",
                                       climatology = "1854-1949 vs 1993-2022",
                                       cmip_model = str_remove(files.new[i], ".nc"),
                                       experiment = "hist_ssp585",
                                       metric = "SST change",
                                       RMSE = sqrt(mean(mod$residuals^2)),
                                       n = nrow(sample)))
  
  ## query each of the regional time series
  
  for(k in 1:length(unique(regional.polygons$region))){ # open k loop (smaller regions)
    # k <- 1
    # keep track of progress
    print(paste(i, k, sep = "-"))
    
    # define regional polygon
    temp.polygon <- regional.polygons %>%
      filter(region == unique(regional.polygons$region)[k])
    
    # select ersst cells in this polygon; all others become NA
    xp <- cbind(temp.polygon$x, temp.polygon$y)
    loc <- cbind(compare$long, compare$lat)
    check <- in.poly(loc, xp=xp)
    
    temp.compare <- compare
    temp.compare$ersst_era_diff[!check] <- NA
    temp.compare$cmip_era_diff[!check] <- NA
    
    # get sample size
    sample <- temp.compare %>%
      filter(!is.na(ersst_era_diff),
             !is.na(cmip_era_diff))
    
    # get RMSE for climatology
    mod <- lm(ersst_era_diff ~ cmip_era_diff, weights = weight, data = temp.compare)
    
    # record results
    climatology_change_rmse <- rbind(climatology_change_rmse,
                              data.frame(region = unique(regional.polygons$region)[k],
                                         climatology = "1854-1949 vs 1993-2022",
                                         cmip_model = str_remove(files.new[i], ".nc"),
                                         experiment = "hist_ssp585",
                                         metric = "SST change",
                                         RMSE = sqrt(mean(mod$residuals^2)),
                                         n = nrow(sample)))   
    
  } # close k loop (regions)
} # close i loop (models)