## Create Fig 1: combined map + FAR time series

library(tidyverse)
library(oce)
library(chron)
library(ncdf4)
library(sf)
library(rnaturalearthdata)
library(rnaturalearth)
library(reshape2)
library(patchwork)

# set palettes
new.col <- oceColorsPalette(64)
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set theme
theme_set(theme_bw())


## convert_coords ------------------------------------------
convert_coords <- function(data,
                           proj.orig,
                           proj.new,
                           xvar.orig = "lon",
                           yvar.orig = "lat",
                           xvar.new = "easting_m",
                           yvar.new = "northing_m",
                           plot = FALSE, ...) {
    ## Add transformed x,y coordinates to the input data.frame
    ##
    ## This function takes as input a data.frame with lon and lat columns
    ## specified by `xvar.orig` and `yvar.orig` and a projection string for the
    ## lon/lat (`proj.orig`) and (1) converts the coordinates to the projection
    ## given by `proj.new` and adds columns `xvar.new` and `yvar.new` to
    ## the input data.frame. The output of the function is the input `data` with
    ## two additional columns added whose names are specified by `xvar.new`
    ## and `yvar.new`.
    ##
    ## data = a data.frame with lon and lat columns
    ## proj.orig = string giving the projection of the lon/lat coords in the
    ##             input data
    ## proj.new = string giving the projection to convert to
    ## xvar.orig = longitude column name in `data`
    ## yvar.orig = latitude column name in `data`
    ## xvar.new = name of easting/lon converted column in output data.frame
    ## yvar.new = name of northing/lat converted column in output data.frame
    ## plot = logical, should the converted coords be plotted
    ## ... = passed to plot()

    if(xvar.new %in% names(data))
        stop("`xvar.new` already exists in data")
    if(yvar.new %in% names(data))
        stop("`yvar.new` already exists in data")

    if(!xvar.orig %in% names(data))
        stop("`xvar.orig` not found")
    if(!yvar.orig %in% names(data))
        stop("`yvar.orig` not found")

    dat <- data.table::copy(data)

    ## Add index column
    dat$index <- 1:nrow(dat)

    ## Subset lon and lat columns
    df.sub <- subset(dat, select = c(xvar.orig, yvar.orig, "index"))

    ## Convert to sf
    sf.sub <- sf::st_as_sf(df.sub, coords = c(1, 2), crs = proj.orig)

    ## Convert from proj.orig to proj.new
    sf.proj <- sf::st_transform(sf.sub, proj.new)
    ind     <- st_drop_geometry(sf.proj)
    coord   <- as.data.frame(st_coordinates(sf.proj))
    df.proj <- data.frame(index = ind$index,
                          x = coord$X,
                          y = coord$Y)
    names(df.proj) <- c("index", xvar.orig, yvar.orig)

    names(df.proj)[names(df.proj) == xvar.orig] <- xvar.new
    names(df.proj)[names(df.proj) == yvar.orig] <- yvar.new

    ## Add eastings and northings to input df
    df.out <- merge(x = dat, y = df.proj, by = "index")
    df.out <- df.out[order(df.out$index), ]
    df.out[["index"]] <- NULL

    if(plot) {
        plot(df.out[[xvar.new]], df.out[[yvar.new]], ...)
    }

    return(df.out)
}



## load and process ERSST ----------------------------------
# load and process SST data
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



# image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "")
# contour(x, y, z, add=T)
# map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

# extract study area
# 54-66 deg. N, 188-202 deg. E
ebs.x <- c(183, 183, 203, 203, 191) 
ebs.y <- c(53, 65, 65, 57.5, 53)

# polygon(ebs.x, ebs.y, border = "red", lwd = 2)

# GOA polygon
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)

# polygon(goa.x, goa.y, border = "red", lwd = 2)

# BC polygon
bc.x <- c(231, 238, 228, 225, 225)
bc.y <- c(55, 49, 49, 53, 55)

# polygon(bc.x, bc.y, border = "red", lwd = 2)

# northern CCE polygon
ncc.x <- c(238, 238, 233, 233, 238)
ncc.y <- c(49, 41, 41, 49, 49)

# polygon(ncc.x, ncc.y, border = "red", lwd = 2)

# southern CCE polygon
scc.x <- c(238, 243, 237, 233, 233, 238)
scc.y <- c(41, 33, 33, 39, 41, 41)

# polygon(scc.x, scc.y, border = "red", lwd = 2)


## ggplot map ----------------------------------------------
## See here for Pacific centered map
## https://stackoverflow.com/questions/69970029/st-crop-on-pacific-centred-naturalearth-world-map

proj.wgs <- sf::st_crs(4326)
worldMap <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()

target_crs <- st_crs("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=133")

# define a long & slim polygon that overlaps the meridian line & set its CRS to match
# that of world
# Centered in lon 133
offset <- 180 - 133


polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)


# modify world dataset to remove overlapping portions with world's polygons
world2 <- worldMap %>% st_difference(polygon)

# Transform
world3 <- world2 %>% st_transform(crs = target_crs)

world4 <- st_crop(
  x = world3, 
  y = st_as_sfc(
    st_bbox(c(xmin= 100, ymin = 10, xmax = 270, ymax = 80), crs = 4326)
  ) %>% st_transform(target_crs)
)
# g <- ggplot() +
#     geom_sf(data = world4,
#             fill = NA, color = "grey10", size = 0.2)
# print(g)


## Reshape SST data
row.names(z) <- x
colnames(z) <- y
df <- melt(z)
names(df) <- c("x", "y", "z")
dfc <- convert_coords(df, proj.wgs, target_crs, "x", "y")

## Convert polygon coordinates
poly.ebs = data.frame(x = ebs.x, y = ebs.y)
poly.goa = data.frame(x = goa.x, y = goa.y)
poly.bc = data.frame(x = bc.x, y = bc.y)
poly.ncc = data.frame(x = ncc.x, y = ncc.y)
poly.scc = data.frame(x = scc.x, y = scc.y)
poly.ebs <- convert_coords(poly.ebs, proj.wgs, target_crs, "x", "y")
poly.goa <- convert_coords(poly.goa, proj.wgs, target_crs, "x", "y")
poly.bc <- convert_coords(poly.bc, proj.wgs, target_crs, "x", "y")
poly.ncc <- convert_coords(poly.ncc, proj.wgs, target_crs, "x", "y")
poly.scc <- convert_coords(poly.scc, proj.wgs, target_crs, "x", "y")

xlim <- range(dfc$easting_m)
ylim <- range(dfc$northing_m)
ylim[1] <- ylim[1] + (0.02 * ylim[1])
ylim[2] <- ylim[2] + (0.02 * ylim[2])

gmap <- ggplot(dfc) +
    geom_raster(aes(x = easting_m, y = northing_m, fill = z)) +
    geom_contour(aes(x = easting_m, y = northing_m, z = z), color = "black", na.rm = TRUE, size = 0.1) +
    geom_sf(data = world4,
            fill = "lightyellow3", color = "grey30", size = 0.2) +
    geom_polygon(data = poly.ebs, aes(x = easting_m, y = northing_m), color = "red", fill = NA, size = 1.0) +
    geom_polygon(data = poly.goa, aes(x = easting_m, y = northing_m), color = "red", fill = NA, size = 1.0) +
    geom_polygon(data = poly.bc, aes(x = easting_m, y = northing_m), color = "red", fill = NA, size = 1.0) +
    geom_polygon(data = poly.ncc, aes(x = easting_m, y = northing_m), color = "red", fill = NA, size = 1.0) +
    geom_polygon(data = poly.scc, aes(x = easting_m, y = northing_m), color = "red", fill = NA, size = 1.0) +
    scale_x_continuous(limits = xlim, expand = c(-0.01, 0)) +
    scale_y_continuous(limits = ylim, expand = c(-0.01, 0)) +
    scale_fill_gradientn(colors = new.col) +
    coord_sf(expand = TRUE) +
    # labs(x = "", y = "", fill = expression(paste('', degree, C))) +
    labs(x = "", y = "", fill = "") +
    # theme(legend.position="top")
    theme(legend.justification = c(0, 0),
          legend.position = c(0.00, 0.35),
          legend.title.align = 0.5,
          legend.key.size = unit(0.35, "cm"),
          legend.background = element_blank(),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6))
print(gmap)


## Setup time series data
far_pred <- read.csv("./CMIP6/summaries/fig1_far_pred.csv")
far_pred$region_plot <- factor(far_pred$region_plot, levels = c("North Pacific",
                                                                "Eastern Bering Sea",
                                                                "Gulf of Alaska",
                                                                "British Columbia Coast",
                                                                "Northern California Current",
                                                                "Southern California Current"))

gts <- ggplot(far_pred) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_line(aes(x = year, y = prob, color = window_plot), size = 0.25) +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = window_plot), alpha = 0.15) +
    facet_wrap(~region_plot, ncol = 2) +
    labs(y = "Fraction of Attributable Risk", color = "", fill = "") +
    scale_color_manual(values = cb[c(2,6)]) +
    scale_fill_manual(values = cb[c(2,6)]) +
    theme(axis.title.x = element_blank()) +
    theme(legend.position="bottom")
print(gts)


# gc <- gmap / gts + plot_layout(heights = unit(c(3, 1), c('in', 'null')))
gc <- gmap / gts
print(gc)
ggsave("./CMIP6/figs/fig1_combined.png", width = 5, height = 7)
