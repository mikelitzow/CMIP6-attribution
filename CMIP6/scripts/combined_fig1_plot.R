# combine map and FAR plots for Fig 1

library(magick)

map <- image_read("./CMIP6/figs/study_site.png")

# map <- image_scale(map, "1300")

far <- image_read("./CMIP6/figs/FAR_rolling_window_time_series_annual_3yrlabelled.png")
img <- c(map, far)
image_info(img)

fig <- image_append(img)

image_write(fig, path = "./CMIP6/figs/combined_fig1.png", format = "png", density = 300)

