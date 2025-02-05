rm(list = ls())
library(sf)
library(fields)
library(sp)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(spBayes)
library(elevatr)
no_nbs <- c(57, 170, 236, 269, 343, 685, 946, 947, 989, 1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203)
all_y <- readRDS(here::here("snow_cleaned.Rda"))[-no_nbs,]
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]

# First year
setwd("D:/77/Research/temp/snow_trend")
weekly_ini_without_lat_alt <- readRDS("weekly_ini_without_lat+out.Rda")
weekly_final_without_lat_alt <- readRDS("weekly_final_without_lat+out.Rda")
weekly_ini_lat_alt <- readRDS("weekly_ini_with_lat+out.Rda")
weekly_final_lat_alt <- readRDS("weekly_final_with_lat+out.Rda")

p_snow_diff <- weekly_final_without_lat_alt[,,,2] - weekly_ini_without_lat_alt[,,,2]
diff_mean_without_lat_alt <- apply(p_snow_diff,c(2,3),mean)
diff_sd_without_lat_alt <- apply(p_snow_diff,c(2,3),sd)

p_snow_diff <- weekly_final_lat_alt[,,,2] - weekly_ini_lat_alt[,,,2]
diff_mean_lat_alt <- apply(p_snow_diff,c(2,3),mean)
diff_sd_lat_alt <- apply(p_snow_diff,c(2,3),sd)





# ----------------------------------------------------------------------------Plots
world <- ne_countries(scale = "medium", returnclass = "sf")
world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]

# Convert the data frame to an sf object
aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"

# Trend
trend_data_without_lat_alt <- st_as_sf(data.frame(cbind(coords, diff_mean_without_lat_alt, diff_sd_without_lat_alt)), coords = c("LON", "LAT"), crs = 4326)
trend_aeqd_without_lat_alt <- st_transform(trend_data_without_lat_alt, crs = aeqd_proj)

trend_data_lat_alt <- st_as_sf(data.frame(cbind(coords, diff_mean_lat_alt, diff_sd_lat_alt)), coords = c("LON", "LAT"), crs = 4326)
trend_aeqd_lat_alt <- st_transform(trend_data_lat_alt, crs = aeqd_proj)


world_aeqd <- st_transform(world_north, crs = aeqd_proj)


equator_points <- data.frame(
  lon = seq(-180, 180, length.out = 200),
  lat = rep(0, 100)  # All points at latitude = 0
)

# Transform equator points to the azimuthal equidistant projection
equator_sf <- st_as_sf(equator_points, coords = c("lon", "lat"), crs = 4326) 
equator_aeqd <- st_transform(equator_sf, crs = aeqd_proj)


library(magick)
# Trend Plot
for (i in 2:52) {
  print(i)
  # Create the plot
  plot1 <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = trend_aeqd_without_lat_alt, aes(color = trend_aeqd_without_lat_alt[[i]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_gradient2(
      low = "red",          # Color for negative values
      mid = "white",  # Color for zero
      high = "blue",         # Color for positive values
      midpoint = 0,          # Set midpoint at zero
      limits = range(c(diff_mean_lat_alt[,-1], diff_mean_without_lat_alt[,-1])), # Set the range of the legend
      guide = "colourbar"
    ) +
    theme_minimal() +
    labs(
      title = paste("Trend Mean for Week (No Covariates)",i),
      color = ""
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) + 
    guides(
      color = guide_colorbar(
        barwidth = 20,   # Adjust the width of the color bar
        barheight = 0.5  # Adjust the height of the color bar
      )
    )
  
  plot2 <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = trend_aeqd_without_lat_alt, aes(color = trend_aeqd_without_lat_alt[[i+52]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_viridis_c(option = "C", direction = -1, limits = c(0, 0.17)) + # Vibrant color palette
    theme_minimal() +
    labs(
      title = paste("Trend sd for Week (No Covariates)",i),
      color = ""
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) + 
    guides(
      color = guide_colorbar(
        barwidth = 20,   # Adjust the width of the color bar
        barheight = 0.5  # Adjust the height of the color bar
      )
    )
  
  plot3 <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = trend_aeqd_lat_alt, aes(color = trend_aeqd_lat_alt[[i]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_gradient2(
      low = "red",          # Color for negative values
      mid = "white",  # Color for zero
      high = "blue",         # Color for positive values
      midpoint = 0,          # Set midpoint at zero
      limits = range(c(diff_mean_lat_alt[,-1], diff_mean_without_lat_alt[,-1])), # Set the range of the legend
      guide = "colourbar"
    ) +
    theme_minimal() +
    labs(
      title = paste("Trend Mean (with Lat+Alt) for Week",i),
      color = ""
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) + 
    guides(
      color = guide_colorbar(
        barwidth = 20,   # Adjust the width of the color bar
        barheight = 0.5  # Adjust the height of the color bar
      )
    )
  
  plot4 <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = trend_aeqd_lat_alt, aes(color = trend_aeqd_lat_alt[[i+52]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_viridis_c(option = "C", direction = -1, limits = c(0,0.17)) + # Vibrant color palette
    theme_minimal() +
    labs(
      title = paste("Trend sd (with Lat+Alt) for Week",i),
      color = ""
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) + 
    guides(
      color = guide_colorbar(
        barwidth = 20,   # Adjust the width of the color bar
        barheight = 0.5  # Adjust the height of the color bar
      )
    )
  
  plot = cowplot::plot_grid(plot1,plot2, plot3, plot4, nrow = 2)
  # Save the plot
  ggsave(
    filename = paste0("plot_",i, ".png"), # Save as plot_1.png, plot_2.png, ...
    plot = plot,                          # Specify the plot object
    width = 15, height = 15,                # Set width and height
    dpi = 300                             # High resolution
  )
} 




png_files <- list.files(pattern = "plot_\\d+\\.png") 
png_files <- png_files[order(as.numeric(gsub("\\D", "", png_files)))]# Find all saved PNGs
gif <- image_read(png_files)                       # Read images
gif <- image_animate(gif, fps = 5)                 # Set frames per second (5 fps)
image_write(gif, "trend_animation.gif")            # Save the GIF

