## evaluate North Pacific-wide warming rate for CMIP6 models

## compare historical / ssp585 with pre-industrial 

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


# start from the beginning - get list of file names
files <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/")

# loop through each file, save the file name and experiment,
# capture the time series of raw temps for area of interest for each file-experiment comparison


# objects for saving dates and temps and file-experiment ID from each realization
piControl.temps <- hist.585.temps  <- matrix()

for(i in 1:length(files)){
  
  # i <- 1
  
  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585/", files[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # nc
  
  # extract one experiment at a time and save with file name
  
  # get list of experiments
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  
  for(j in 1:length(experiments)){
    
    # j <- 1
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "hist_ssp585"){
      
      temp <- data.frame(model =  files[i],
                         experiment = ncvar_get(nc, "experiment", start = j, count = 1))
      
      
      
      # extract dates
      
      d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
      yr <- as.numeric(as.character(years(d)))
      
      # extract spatial area
      x <- ncvar_get(nc, "lon")
      y <- ncvar_get(nc, "lat")
      
      SST <- ncvar_get(nc, "tos", verbose = F, start = c(j,1,1,1), count = c(1,-1,-1,-1))
      
      # Change data from a 3-D array to a matrix of monthly data by grid point:
      # First, reverse order of dimensions ("transpose" array)
      SST <- aperm(SST, 3:1)
      
      # Change to matrix with column for each grid point, rows for monthly means
      SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3])) 
      dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
      
      # remove  values south of 20N
      drop <- lat < 20
      
      SST[,drop] <- NA
      
      # plot to check 
      png(paste("./CMIP6/figs/", files[i], "_hist_spatial_means_1850-1949.png", sep = ""), 6, 4.5, units = "in", res = 300)
      SST.mean <- colMeans(SST[yr < 1950,])
      z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
      image(x,y,z, col=new.col)
      contour(x, y, z, add=T, col="white")
      map('world2Hires',fill=F,add=T, lwd=2)
      
      dev.off()
      
      these.temps <- data.frame(xx = rowMeans(SST, na.rm = T))
      names(these.temps)<-  str_remove(files[i], ".nc")
      
      hist.585.temps <- cbind(hist.585.temps, these.temps)
      
    } # close if 
    
    
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "piControl"){
      
      # extract dates
      
      d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
      yr <- as.numeric(as.character(years(d)))
      
      # extract spatial area
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
      dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
      
      # remove  values south of 20N
      drop <- lat < 20
      
      SST[,drop] <- NA
      
      # plot to check 
      png(paste("./CMIP6/figs/", files[i], "_piControl_spatial_means.png", sep = ""), 6, 4.5, units = "in", res = 300)
      SST.mean <- colMeans(SST)
      z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
      image(x,y,z, col=new.col)
      contour(x, y, z, add=T, col="white")
      map('world2Hires',fill=F,add=T, lwd=2)
      
      dev.off()
      
      
      these.temps <- data.frame(xx = rowMeans(SST, na.rm = T))
      names(these.temps)<-  str_remove(files[i], ".nc")
      
      piControl.temps <- cbind(piControl.temps, these.temps)
      
    } # close if 
    
  } # close j
} # close i

# remove leading column of NAs
piControl.temps <- piControl.temps[,-1]
hist.585.temps <- hist.585.temps[,-1]

# get annual means
ff <- function(x) tapply(x, as.numeric(as.character(years((d)))), mean)

hist.585.annual <- apply(hist.585.temps, 2, ff)
piControl.annual <- apply(piControl.temps, 2, ff)



# and ersst for the same area

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
image.plot(x,y,z, col=oceColorsPalette(64))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

# calculate monthly mean
obs.n.pac.sst <- rowMeans(SST, na.rm = T)

# and annual observed means
n.pac.ann.obs.sst <- data.frame(year = 1854:2021,
                                ersst = tapply(obs.n.pac.sst, as.numeric(as.character(yr)), mean),
                                observations = "ERSST")

# and observed warming rate
n.pac.obs.warming <- data.frame(year = 1854:2021,
                                ersst.warming = n.pac.ann.obs.sst$ersst - 
                                  mean(n.pac.ann.obs.sst$ersst[n.pac.ann.obs.sst$year %in% 1854:1949]))
# switch to C from K!
K_C <- colMeans(hist.585.annual) > 200

hist.585.annual[,K_C] <- hist.585.annual[,K_C] - 273.15
piControl.annual[,K_C] <- piControl.annual[,K_C] - 273.15

# compare time series!

check.hist <- as.data.frame(hist.585.annual) %>%
  mutate(experiment = "hist.ssp585",
         year = 1850:2099) %>%
  pivot_longer(cols = c(-year, -experiment))

check.control <- as.data.frame(piControl.annual) %>%
  mutate(experiment = "piControl",
         year = 1850:2099) %>%
  pivot_longer(cols = c(-year, -experiment))

check <- rbind(check.control, check.hist)

ggplot(check, aes(year, value, color = experiment)) +
  geom_line() +
  geom_line(data = n.pac.ann.obs.sst, aes(year, ersst, color = observations), lwd = 1) +
  facet_wrap(~name)

ggsave("./CMIP6/figs/control_hist585_comparison.png", width = 8, height = 6, units = 'in')

# now calculate warming rate
warming.rate <- hist.585.annual

for(j in 1:ncol(hist.585.annual)){
  
  warming.rate[,j] <- hist.585.annual[,j] - colMeans(hist.585.annual[1:100,])[j]
  
}

colMeans(hist.585.annual[1:100,])

warming.rate <- as.data.frame(warming.rate) %>%
  mutate(year = 1850:2099) %>%
  pivot_longer(cols = -year)


ggplot(warming.rate, aes(year, value, color = name)) +
  geom_line() +
  geom_line(data = n.pac.obs.warming, aes(year, ersst.warming), color = "black", lwd = 1) +
  ggtitle("Annual SST 1850-2099 (ERSST observations in black)") +
  ylab("SST (°C)") +
  theme(axis.title.x = element_blank())

ggsave("./CMIP6/figs/n_pac_model_estimated_warming_rate.png", width = 7, height = 5, units = 'in')

# save these values for future analysis
write.csv(warming.rate, "./CMIP6/summaries/north_pacific_annual_modeled_sst.csv")
write.csv(n.pac.obs.warming, "./CMIP6/summaries/north_pacific_annual_observed_sst.csv")

# need to:

# 1) fit loess to full time series for each model and get estimated time
# for 0.5, 1, 1.5, 2 degrees of warming wrt 1850-1949

# 2) fit loess to ersst to calculate warming wrt 1854-1949 through 2021

# 3) regress annual observed on annual modeled to weight each model

# 4) fit Bayesian regression model to estimate overall warming event


warming.rate <- warming.rate %>%
  pivot_wider(names_from = name, values_from = value)
warming.rate <- as.matrix(warming.rate)

# create an object to catch
n.pac.warming.timing <- warming.rate

for(j in 2:ncol(warming.rate)){
  
  mod <- loess(warming.rate[,j] ~ warming.rate[,1])
  n.pac.warming.timing[,j] <- predict(mod)
  
} 

# plot to check
check.plot <- as.data.frame(n.pac.warming.timing) %>%
  pivot_longer(cols = -year)

ggplot(check.plot, aes(year, value)) +
  geom_line() + 
  facet_wrap(~name) +
  geom_hline(yintercept = c(0.5, 1, 1.5, 2), lty = 2, color = "red") +
  coord_cartesian(xlim = c(2000,2050), ylim = c(0,4))

# looks good!

# get the year that each warming level is reached for each model
levels <- c(0.5, 1, 1.5, 2)

timing <- data.frame()
n.pac.warming.timing <- as.data.frame(n.pac.warming.timing)

for(j in 2:ncol(n.pac.warming.timing)){
  # j <- 2
  temp <- NA
  
  for(i in 1:length(levels)){
    # i <- 1
    
    temp[i] <- min(n.pac.warming.timing$year[n.pac.warming.timing[,j] >= levels[i]])
    
  }
  
  timing <- rbind(timing,
                  data.frame(
                    model = colnames(n.pac.warming.timing[j]),
                    level = levels,
                    year = temp))
}

# and plot
ggplot(timing, aes(year)) +
  geom_histogram(bins = 6) +
  facet_wrap(~level, scales = "free_x")

ggplot(timing, aes(level, year, color = model)) +
  geom_point() +
  labs(x = "N. Pacific warming wrt 1850-1949 (°C)",
       y = "Year first reached") +
  scale_y_continuous(breaks = seq(1960, 2090, by = 10))

ggsave("./CMIP6/figs/N_Pac_warming_rate_by_model.png", width = 6, height = 4, units = 'in')

timing$level <- as.factor(timing$level)

obs.timing <- data.frame()


# fit loess to ersst to get observed warming rate
mod <- loess(n.pac.obs.warming[,2] ~ n.pac.obs.warming[,1])
n.pac.obs.warming$trend <- predict(mod)

# plot to check
ggplot(n.pac.obs.warming, aes(year, trend)) +
  geom_line()

# looks right


# now get observed timing of different warming levels - 
# we've only reached 1, which is 0.5 

i <- 1
temp <- min(n.pac.obs.warming$year[n.pac.obs.warming$trend >= levels[i]]) 

obs.timing <- data.frame(level = levels,
                         timing = c(2003, NA, NA, NA))

# 2003


obs.timing$level <- as.factor(obs.timing$level)

ggplot(timing, aes(y = year, level)) +
  geom_boxplot() +
  labs(x = "N Pacific warming wrt 1850-1949 (°C)",
       y = "Year first reached",
       title = "Red dot = ERSSTv5 warming") +
  geom_point(data = obs.timing, aes(y = timing, as.factor(level)), color = "red") 

ggsave("./CMIP6/figs/ne_pacific_warming_rate_models_obs.png", width = 5, height = 4, units = 'in')



# save timing
write.csv(timing, "./CMIP6/summaries/model.ne.pacific.warming.timing.csv")


# evaluate ability of different models to predict warming for 1950-2021

# reload unsmoothed warming data
model.warming.rate <- read.csv("./CMIP6/summaries/ne_pacific_annual_modeled_sst.csv", row.names = 1)

obs.warming.rate <- read.csv("./CMIP6/summaries/ne_pacific_annual_observed_sst.csv", row.names = 1)

predict.timing <- model.warming.rate %>%
  filter(year %in% 1972:2021) 

response.timing <- obs.warming.rate %>%
  filter(year %in% 1972:2021) %>%
  rename(ersst = ersst.warming)

predict.timing <- left_join(predict.timing, response.timing)

# plot to check
ggplot(predict.timing, aes(value, ersst)) +
  geom_point() +
  facet_wrap(~name) # avoids a 'fishook' around declining rate of 1950s


model.warming.evaluation <- data.frame()

models <- unique(predict.timing$name)

for(i in 1:length(models)){
  
  # i <- 1
  
  temp <- predict.timing %>%
    filter(name == models[i])
  
  linear.fit <- lm(ersst ~ value, data = temp)
  
  RSS <- c(crossprod(linear.fit$residuals))
  Pearson.resid <- RSS / linear.fit$df.residual
  
  
  model.warming.evaluation <- rbind(model.warming.evaluation,
                                    data.frame(model = models[i],
                                               coeff = coefficients(linear.fit)[2],
                                               Pearson.resid = Pearson.resid))
  
}


model.warming.evaluation # should use coefficients! (inverse of difference from 1)

model.warming.evaluation$coeff.from.one <- 1-model.warming.evaluation$coeff
model.warming.evaluation$weight <- abs(1/model.warming.evaluation$coeff.from.one)

ggplot(model.warming.evaluation, aes(weight)) +
  geom_histogram(bins = 8, fill = "grey", color = "black") # seems right!


ggplot(model.warming.evaluation, aes(abs(coeff.from.one), weight)) +
  geom_point()

# save 
write.csv(model.warming.evaluation, "./CMIP6/summaries/N_Pac_warming_model_weights.csv", row.names = F)


# next, brms estimates of warming timing
