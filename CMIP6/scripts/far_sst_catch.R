# Combine FAR and catch for GOA sockeye example

source("./scripts/load.R")
dir.create("./figures/sst_catch", showWarnings = FALSE)


## Clean catch data ----------------------------------------

## Read in salmon data
goa_catch <- read.csv("./data/goa.catch.csv")
goa_age   <- read.csv("./data/goa_age.csv")


## Exclude 2022 for now (3-yr runnning mean = NA!)
# goa_catch <- goa_catch[goa_catch$year != 2022, ]


## Add region
goa_catch$region <- "GOA"

## Combine
catch_wide <- goa_catch

## Read in FAR data
goa_far <- read.csv("./data/complete_FAR_RR_time_series_with_uncertainty.csv") %>%
    filter(region == "Gulf_of_Alaska",
           window == "3yr_running_mean")

## Subset years
catch_wide <- catch_wide[catch_wide$year >= 1965, ]


## Remove south peninsula: from Greg R. -- Very very very few sockeye originate
## from the southern Alaska Peninsula management area, but the catch in this
## area likely dominates the GOA trend line because of catches from Bristol Bay
## origin sockeye
catch_wide <- catch_wide[catch_wide$area != "south.peninsula", ]

## Sum across catch regions
sock <- plyr::ddply(catch_wide, .(region, year), summarize, catch = sum(sockeye))

## Add species
sock$species <- "Sockeye"

## Combine in long format
catch <- sock
catch$species_fac <- factor(catch$species, levels = unique(catch$species))
catch$region_fac  <- factor(catch$region, levels = unique(catch$region))

## Add era variable
catch$era <- "1977-1988"
catch$era <- ifelse(catch$year >= 1989, "1989-2022", catch$era)
catch$era <- ifelse(catch$year <= 1976, "1965-1976", catch$era)


## Combine catch + sst -------------------------------------
## "winter" = November-March (year corresponding to January);
##            annual = January - December both winter and annual data scaled wrt
##            20th century (1901-1999 for winter; 1900-1999 for annual)
## .2 denotes two-yr running mean (year of and year following year of interest)
## .3 denotes three-yr running mean (year before, year of, and year following year of interest)

## Get mean age props for sockeye
age_wgt_goa <- c(mean(goa_age$R_ocean_1),
                 mean(goa_age$R_ocean_2),
                 mean(goa_age$R_ocean_3),
                 mean(goa_age$R_ocean_4),
                 mean(goa_age$R_ocean_5))


## Add FAR
catch$annual_far_3 <- NA

for(i in 1:nrow(catch)) {
    # i <- 1
    rg <- catch$region[i]
    sp <- catch$species[i]
    yr <- catch$year[i]
    if(rg == "GOA") {
        far_dat <- goa_far
        age_wgt_i <- age_wgt_goa
    }

    if(sp == "Sockeye") {
        far_i <- far_dat[far_dat$year %in% (yr - 1):(yr - 5), ]
        annual_far_3 <- weighted.mean(far_i$FAR, w = rev(age_wgt_i))
    }

    catch$annual_far_3[i] <- annual_far_3

}


## Create different time periods ---------------------------
catch_1965 <- catch
catch_1965 <- plyr::ddply(catch_1965, .(region, species), transform,
                          log_catch = log(catch),
                          log_catch_stnd = scale(log(catch)),
                          catch_stnd = scale(catch))


catch <- catch[catch$year >= 1989, ]
catch <- plyr::ddply(catch, .(region, species), transform,
                     log_catch = log(catch),
                     log_catch_stnd = scale(log(catch)),
                     catch_stnd = scale(catch))


write.csv(catch, file = "./data/GOA_sockeye_catch_far.csv", row.names = FALSE)
