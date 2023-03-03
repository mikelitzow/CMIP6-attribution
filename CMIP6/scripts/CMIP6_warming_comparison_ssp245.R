## evaluate North Pacific-wide warming rate for CMIP6 models

## compare historical / ssp245 with pre-industrial 

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./CMIP6/scripts/stan_utils.R")

# set palette
new.col <- oceColorsPalette(64)

# set theme
theme_set(theme_bw())

# make function to calculate cell weights
cell.weight <- function(x)  sqrt(cos(x*pi/180))

# and function to calculate weighted means with these weights 
weighted.cell.mean <- function(x) weighted.mean(x, weights, na.rm = T)

# start from the beginning - get list of file names
files <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp245/")

# loop through each file, save the file name and experiment,
# capture the time series of raw temps for area of interest for each file-experiment comparison


# objects for saving dates and temps and file-experiment ID from each realization
piControl.temps <- hist.245.temps  <- matrix()

for(i in 1:length(files)){
  
  # i <- 1
  
  path <- paste("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp245/", files[i], sep="")
  
  # load file
  nc <- nc_open(path)
  
  # nc
  
  # extract one experiment at a time and save with file name
  
  # get list of experiments
  experiments <-  ncvar_get(nc, "experiment", verbose = F)
  
  for(j in 1:length(experiments)){
    
    # j <- 1
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "hist_ssp245"){
      
      temp <- data.frame(model =  files[i],
                         experiment = ncvar_get(nc, "experiment", start = j, count = 1))
      
      
      
      # extract dates
      
      d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
      yr <- as.numeric(as.character(years(d)))
      
      # extract spatial area
      x <- ncvar_get(nc, "lon")
      y <- ncvar_get(nc, "lat")
      
      # Keep track of corresponding latitudes and longitudes of each column:
      lat <- rep(y, length(x))   
      lon <- rep(x, each = length(y))   
      
      # calculate cell weights
      weights <- cell.weight(lat)
      
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
      png(paste("./CMIP6/figs/", files[i], "_hist_spatial_means_1850-1949_245_run.png", sep = ""), 6, 4.5, units = "in", res = 300)
      SST.mean <- colMeans(SST[yr < 1950,])
      z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
      image(x,y,z, col=new.col)
      contour(x, y, z, add=T, col="white")
      map('world2Hires',fill=F,add=T, lwd=2)
      
      dev.off()
      
      these.temps <- data.frame(xx = apply(SST, 1, weighted.cell.mean))
      names(these.temps)<-  str_remove(files[i], ".nc")
      
      hist.245.temps <- cbind(hist.245.temps, these.temps)
      
    } # close if 
    
    
    if(ncvar_get(nc, "experiment", start = j, count = 1) == "piControl"){ # no control for these files, returns empty
      
      # extract dates
      
      d <- dates(ncvar_get(nc, "time"), origin = c(1,1,1970))
      yr <- as.numeric(as.character(years(d)))
      
      # extract spatial area
      x <- ncvar_get(nc, "lon")
      y <- ncvar_get(nc, "lat")
      
      # Keep track of corresponding latitudes and longitudes of each column:
      lat <- rep(y, length(x))   
      lon <- rep(x, each = length(y))   
      
      # calculate cell weights
      weights <- cell.weight(lat)
      
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
      
      
      these.temps <- data.frame(xx = apply(SST, 1, weighted.cell.mean))
      names(these.temps)<-  str_remove(files[i], ".nc")
      
      piControl.temps <- cbind(piControl.temps, these.temps)
      
    } # close if 
    
  } # close j
} # close i

# remove leading column of NAs
hist.245.temps <- hist.245.temps[,-1]

# get annual means
ff <- function(x) tapply(x, as.numeric(as.character(years((d)))), mean)

hist.245.annual <- apply(hist.245.temps, 2, ff)



# and ersst for the same area

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(2022-12-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(66)][(110):1:(250)]", "~temp")

# load and process SST data
# nc <- nc_open("~temp")

# full North Pacific, 1854-2022
nc <- nc_open("./CMIP6/data/nceiErsstv5_89b2_444a_5b23.nc")

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


# calculate cell weights
weights <- cell.weight(lat)

# plot to check
temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

# calculate monthly mean
obs.n.pac.sst <- apply(SST, 1, weighted.cell.mean)

# and annual observed means
n.pac.ann.obs.sst <- data.frame(year = 1854:2022,
                                ersst = tapply(obs.n.pac.sst, as.numeric(as.character(yr)), mean),
                                observations = "ERSST")

# and observed warming rate
n.pac.obs.warming <- data.frame(year = 1854:2022,
                                ersst.warming = n.pac.ann.obs.sst$ersst - 
                                  mean(n.pac.ann.obs.sst$ersst[n.pac.ann.obs.sst$year %in% 1854:1949]))
# check for °K!
K_C <- colMeans(hist.245.annual) > 200
K_C # none


# compare time series!

check.hist <- as.data.frame(hist.245.annual) %>%
  mutate(experiment = "hist.ssp245",
         year = 1850:2099) %>%
  pivot_longer(cols = c(-year, -experiment))

check.control <- as.data.frame(piControl.annual) %>%
  mutate(experiment = "piControl",
         year = 1850:2099) %>%
  pivot_longer(cols = c(-year, -experiment))

check <- rbind(check.control, check.hist)

ggplot(check.hist, aes(year, value)) +
  geom_line(color = "dark blue") +
  geom_line(data = n.pac.ann.obs.sst, aes(year, ersst), color = "black", lwd = 1) +
  facet_wrap(~name)

ggsave("./CMIP6/figs/ersst_hist245_comparison.png", width = 8, height = 6, units = 'in')

# now calculate warming rate
warming.rate <- hist.245.annual

for(j in 1:ncol(hist.245.annual)){
  
  warming.rate[,j] <- hist.245.annual[,j] - colMeans(hist.245.annual[1:100,])[j]
  
}

colMeans(hist.245.annual[1:100,])

warming.rate <- as.data.frame(warming.rate) %>%
  mutate(year = 1850:2099) %>%
  pivot_longer(cols = -year)


ggplot(warming.rate, aes(year, value, color = name)) +
  geom_line() +
  geom_line(data = n.pac.obs.warming, aes(year, ersst.warming), color = "black", lwd = 1) +
  ggtitle("Annual SST 1850-2099 (ERSST observations in black)") +
  ylab("SST (°C)") +
  theme(axis.title.x = element_blank())

ggsave("./CMIP6/figs/n_pac_model_estimated_warming_rate_ssp245.png", width = 7, height = 5, units = 'in')

# save these values for future analysis
write.csv(warming.rate, "./CMIP6/summaries/north_pacific_annual_modeled_sst_wrt_1850-1949_ssp245.csv")

# need to:

# 1) fit brms model to full time series for each model and get estimated time
# for 0.5, 1, 1.5, 2 degrees of warming wrt 1850-1949

# 2) fit brms to ersst to calculate warming wrt 1854-1949 through 2021

# 3) regress annual observed on annual modeled to weight each model

# 4) fit Bayesian regression model to estimate overall warming event

warming.rate <- read.csv("./CMIP6/summaries/north_pacific_annual_modeled_sst_wrt_1850-1949_ssp245.csv", row.names = 1) %>%
  mutate(name = str_remove_all(name, "_245")) %>% # removing trailing _245 from model names
  rename(model = name,
         warming = value)

# fitting one brms model with 23 separate smooths (one for each CMIP6 model) is not practicable
# instead, loop through the models and fit each separately

models <- unique(warming.rate$model)

## fit brms model to each CMIP6 time series to estimate warming for each year------

# warming_formula <-  bf(warming ~ s(year))
# 
# for(m in 22:length(models)){
# # m <- 12
# temp.dat <- warming.rate %>%
#   filter(model == models[m]) 
# 
# warming_brm <- brm(warming_formula,
#                    data = temp.dat,
#                    cores = 4, chains = 4, iter = 5000,
#                    save_pars = save_pars(all = TRUE),
#                    control = list(adapt_delta = 0.99, max_treedepth = 16))
# 
# saveRDS(warming_brm, file = paste("./CMIP6/brms_output/warming_brm_", models[m], ".rds", sep = ""))
# 
# }
# 
# # run diagnostics for each model
# for(m in 1:length(models)){
#   # m <- 23
#   mod_check <- readRDS(file = paste("./CMIP6/brms_output/warming_brm_", models[m],".rds", sep = ""))
#   
#   print(paste("Model #", m, ": ", models[m], sep = ""))
#   
#   check_hmc_diagnostics(mod_check$fit)
#   
#   print(neff_lowest(mod_check$fit))
#   
#   print(rhat_highest(mod_check$fit))
# }
# 
# # finally, loop through each model and record predicted warming by year
# 
# warming_out <- data.frame()
# 
# new.dat <- data.frame(year = 1850:2099)
# 
# for(m in 1:length(models)){
#   # m <- 23
#   mod <- readRDS(file = paste("./CMIP6/brms_output/warming_brm_", models[m],".rds", sep = ""))
#   
#   print(paste("Model #", m, ": ", models[m], sep = ""))
#  
#   pred <- posterior_epred(mod, newdata = new.dat)
#   
#   warming_out <- rbind(warming_out,
#                        data.frame(model = models[m],
#                                   year = 1850:2099,
#                                   warming = colMeans(pred)))
#   
# }
# 
# # save
# write.csv(warming_out, "./CMIP6/summaries/CMIP6_brms_warming_rate_ssp585.csv", row.names = F)

## fit inverse model - year as a function of warming -----------------------------


# set brms formula
inverse_formula <-  bf(year ~ s(warming))


for(m in 1:length(models)){
  # m <- 20
  temp.dat <- warming.rate %>%
    filter(model == models[m],
           year >= 1973) # limit to 1973-on to ease fitting 

inverse_warming_brm <- brm(inverse_formula,
                           data = temp.dat, 
                           cores = 4, chains = 4, iter = 4000,
                           save_pars = save_pars(all = TRUE),
                           control = list(adapt_delta = 0.999, max_treedepth = 16))

saveRDS(inverse_warming_brm, file = paste("./CMIP6/brms_output/inverse_warming_brm_", models[m],"_ssp245.rds", sep = ""))

}

## now loop through again and run model diagnostics

for(m in 1:length(models)){
  # m <- 1
  mod_check <- readRDS(file = paste("./CMIP6/brms_output/inverse_warming_brm_", models[m],".rds", sep = ""))
  
  print(paste("Model #", m, ": ", models[m], sep = ""))

  check_hmc_diagnostics(mod_check$fit)

  print(neff_lowest(mod_check$fit))

  print(rhat_highest(mod_check$fit))
}

# [1] "Model #20: MIROC6"
# 
# Divergences:
#   1 of 8000 iterations ended with a divergence (0.0125%).
# Try increasing 'adapt_delta' to remove the divergences.
# 
# Tree depth:
#   0 of 8000 iterations saturated the maximum tree depth of 16.
# 
# Energy:
#   E-BFMI indicated no pathological behavior.
# lp__ sds_swarming_1      zs_1_1[3]  bs_swarming_1 
# 1390.439       1639.767       1649.904       1850.701 
# bs_swarming_1       zs_1_1[3]  sds_swarming_1            lp__ s_swarming_1[3] 
# 1.002567        1.001854        1.001832        1.001711        1.001613


## brms estimates of warming timing for each model-----------------------------
warming.timing <- data.frame()

for(m in 1:length(models)){
  # m <- 1
  inverse_warming_brm <- readRDS(file = paste("./CMIP6/brms_output/inverse_warming_brm_", models[m],".rds", sep = ""))
  
  print(paste("Model #", m, ": ", models[m], sep = ""))
  
  ce1s_1 <- conditional_effects(inverse_warming_brm, effect = "warming", re_formula = NA,
                                probs = c(0.025, 0.975), resolution = 10000)
  
  index <- ce1s_1$warming$warming
  
  choose <- c(which.min(abs(index - 0.5)),
              which.min(abs(index - 1.0)),
              which.min(abs(index - 1.5)),
              which.min(abs(index - 2.0)))
  
  warming.timing <- rbind(warming.timing,
                          data.frame(model = models[m],
                            warming = c(0.5, 1.0, 1.5, 2.0),
                          year = ce1s_1$warming$estimate__[choose],
                          UCI = ce1s_1$warming$upper__[choose],
                          LCI = ce1s_1$warming$lower__[choose]))
}
## warming model weights ----------------------------------
# evaluate ability of different models to predict warming for 1973-2022

# reload unsmoothed warming data
model.warming.rate <- read.csv("./CMIP6/summaries/ne_pacific_annual_modeled_sst.csv", row.names = 1)

obs.warming.rate <- read.csv("./CMIP6/summaries/north_pacific_annual_observed_sst.csv", row.names = 1)

predict.timing <- model.warming.rate %>%
  filter(year %in% 1973:2022) 

response.timing <- obs.warming.rate %>%
  filter(year %in% 1973:2022) %>%
  rename(ersst = ersst.warming)

predict.timing <- left_join(predict.timing, response.timing)

# plot to check
ggplot(predict.timing, aes(value, ersst)) +
  geom_point() +
  facet_wrap(~name) +
  geom_abline(color = "red")# avoids a 'fishook' around declining rate of 1950s


model.warming.evaluation <- data.frame()

models <- unique(predict.timing$name)

for(i in 1:length(models)){
  
  # i <- 1
  
  temp <- predict.timing %>%
    filter(name == models[i])
  
  linear.fit <- lm(ersst ~ value, data = temp)

  model.warming.evaluation <- rbind(model.warming.evaluation,
                                    data.frame(model = models[i],
                                               RMSE = sqrt(mean(linear.fit$residuals^2))))
  
}


model.warming.evaluation 

ggplot(model.warming.evaluation, aes(RMSE)) +
  geom_histogram(bins = 8, fill = "grey", color = "black") # seems right!


# save 
write.csv(model.warming.evaluation, "./CMIP6/summaries/N_Pac_warming_model_weights.csv", row.names = F)


