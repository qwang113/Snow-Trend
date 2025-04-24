rm(list = ls())
library(sf)
library(fields)
library(sp)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(spBayes)
library(elevatr)

setwd(here::here())
no_nbs <- c(57, 170, 236, 269, 343, 685, 946, 947, 989, 1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203)
all_y <- readRDS("snow_cleaned.Rda")[-no_nbs,]
all_y_nnbs <- readRDS("snow_cleaned.Rda")[no_nbs,]
y <- rbind(all_y, all_y_nnbs)[,-c(1,2)]
coords <- rbind(all_y, all_y_nnbs)[,1:2]

# First year
setwd("D:/77/Research/temp/snow")
weekly_ini_without_lat_alt <- readRDS("ini_year_bym.Rda")
weekly_final_without_lat_alt <- readRDS("fin_year_bym.Rda")

weekly_ini_lat_alt <- readRDS("ini_year_bym+sc.Rda")
weekly_final_lat_alt <- readRDS("fin_year_bym+sc.Rda")


weekly_ini_INDEP <- readRDS("ini_year_ind.Rda")
weekly_final_INDEP <- readRDS("fin_year_ind.Rda")



p_snow_diff <- weekly_final_without_lat_alt[,,,2] - weekly_ini_without_lat_alt[,,,2]
diff_mean_without_lat_alt <- apply(p_snow_diff,c(2,3),mean)
diff_sd_without_lat_alt <- apply(p_snow_diff,c(2,3),sd)

p_snow_diff <- weekly_final_lat_alt[,,,2] - weekly_ini_lat_alt[,,,2]
diff_mean_lat_alt <- apply(p_snow_diff,c(2,3),mean)
diff_sd_lat_alt <- apply(p_snow_diff,c(2,3),sd)

p_snow_diff <- weekly_final_INDEP[,,,2] - weekly_ini_INDEP[,,,2]
diff_mean_INDEP <- apply(p_snow_diff,c(2,3),mean)
diff_sd_INDEP <- apply(p_snow_diff,c(2,3),sd)

# p_snow_diff <- weekly_final_lat_alt_prior_100[,,,2] - weekly_ini_lat_alt_prior_100[,,,2]
# diff_mean_lat_alt_prior_100 <- apply(p_snow_diff,c(2,3),mean)
# diff_sd_lat_alt_prior_100 <- apply(p_snow_diff,c(2,3),sd)




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

trend_data_INDEP <- st_as_sf(data.frame(cbind(coords, diff_mean_INDEP, diff_sd_INDEP)), coords = c("LON", "LAT"), crs = 4326)
trend_aeqd_INDEP <- st_transform(trend_data_INDEP, crs = aeqd_proj)


# trend_data_lat_alt_prior_100 <- st_as_sf(data.frame(cbind(coords, diff_mean_lat_alt_prior_100, diff_sd_lat_alt_prior_100)), coords = c("LON", "LAT"), crs = 4326)
# trend_aeqd_lat_alt_prior_100 <- st_transform(trend_data_lat_alt_prior_100, crs = aeqd_proj)


world_aeqd <- st_transform(world_north, crs = aeqd_proj)


equator_points <- data.frame(
  lon = seq(-180, 180, length.out = 200),
  lat = rep(0, 100)  # All points at latitude = 0
)

# Transform equator points to the azimuthal equidistant projection
equator_sf <- st_as_sf(equator_points, coords = c("lon", "lat"), crs = 4326) 
equator_aeqd <- st_transform(equator_sf, crs = aeqd_proj)

wk_names <- substr(names(all_y[3:54]),2,6) 
wk_names <- sub("^([0-9]+\\.[0-9]+).*", "\\1", wk_names)

month_raw <- substr(wk_names, 1, 2)
month_num <- as.integer(gsub("\\.", "", month_raw))
month_abbr <- month.abb[month_num]

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
      limits = range(c(diff_mean_lat_alt[,-1], diff_mean_without_lat_alt[,-1], diff_mean_INDEP[,-1])), # Set the range of the legend
      guide = "colourbar"
    ) +
    theme_minimal() +
    labs(
      title = paste("Trend Mean for Week (No Covariates)",wk_names[i]),
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
    scale_color_viridis_c(option = "C", direction = -1, limits = c(0, 0.18)) + # Vibrant color palette
    theme_minimal() +
    labs(
      title = paste("Trend sd for Week (No Covariates)",wk_names[i]),
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
  trend_aeqd_lat_alt$color_val <- pmax(pmin(trend_aeqd_lat_alt[[i]], 0.5), -0.5)
  plot3 <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
    geom_sf(data = trend_aeqd_lat_alt, aes(color = color_val), size = 4, shape = 18) +
    geom_sf(data = world_aeqd, fill = NA, color = "black") +
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) +
    scale_color_gradient2(
      low = "red",         
      mid = "white",       
      high = "blue",        
      midpoint = 0,
      limits = c(-0.5, 0.5),
      guide = "colourbar"
    ) +
    theme_minimal() +
    labs(
      title = paste("Trend Mean (with Lat+Alt) for Week", wk_names[i]),
      color = ""
    ) +
    annotate("text", x = -7000000, y = -1000000, label = month_abbr[i], size = 20) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 20)
    ) +
    guides(
      color = guide_colorbar(
        barwidth = 50,
        barheight = 0.5
      )
    )

  plot4 <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = trend_aeqd_lat_alt, aes(color = trend_aeqd_lat_alt[[i+52]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_viridis_c(option = "C", direction = -1, limits = c(0,0.18)) + # Vibrant color palette
    theme_minimal() +
    labs(
      title = paste("Trend sd (with Lat+Alt) for Week",wk_names[i]),
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
  plot5 <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = trend_aeqd_lat_alt, aes(color = trend_aeqd_INDEP[[i]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_gradient2(
      low = "red",          # Color for negative values
      mid = "white",  # Color for zero
      high = "blue",         # Color for positive values
      midpoint = 0,          # Set midpoint at zero
      limits = range(c(diff_mean_lat_alt[,-1], diff_mean_without_lat_alt[,-1], diff_mean_INDEP[,-1])), # Set the range of the legend
      guide = "colourbar"
    ) +
    theme_minimal() +
    labs(
      title = paste("Trend Mean (Independent) for Week",wk_names[i]),
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

  plot6 <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = trend_aeqd_lat_alt, aes(color = trend_aeqd_INDEP[[i+52]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_viridis_c(option = "C", direction = -1, limits = c(0,0.18)) + # Vibrant color palette
    theme_minimal() +
    labs(
      title = paste("Trend sd (Independent) for Week",wk_names[i]),
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
  
  
  plot = cowplot::plot_grid(plot5,plot1, nrow = 2)
  # plot = plot3
  # Save the plot
  ggsave(
    filename = paste0("plot_",i, ".png"), # Save as plot_1.png, plot_2.png, ...
    plot = plot,                          # Specify the plot object
    width = 15, height = 15,                # Set width and height
    dpi = 100                             # High resolution
  )
}

setwd("D:/77/Research/temp/snow")
png_files <- list.files(pattern = "plot_\\d+\\.png$") 
png_files <- png_files[order(as.numeric(gsub("\\D", "", png_files)))]# Find all saved PNGs
gif <- image_read(png_files)                       # Read images
gif <- image_animate(gif, fps = 5)                 # Set frames per second (5 fps)
image_write(gif, "trend_animation_bym+.gif")            # Save the GIF

special <- c( which(wk_names %in% c("11.27", "1.15","4.1","6.3")))

trend_aeqd_lat_alt$color_val <- pmax(pmin(trend_aeqd_lat_alt[[special[1]]], 0.3), -0.3)
p1 <- ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
  geom_sf(data = trend_aeqd_lat_alt, aes(color = color_val), size = 4, shape = 18) +
  geom_sf(data = world_aeqd, fill = NA, color = "black") +
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) +
  scale_color_gradient2(
    low = "red",         
    mid = "white",       
    high = "blue",        
    midpoint = 0,
    limits = c(-0.3, 0.3),
    guide = "colourbar"
  ) +
  theme_minimal() +
  labs(
    title = paste("Trend Mean (with Lat+Alt) for Week", wk_names[special[1]]),
    color = ""
  ) +
  annotate("text", x = -7000000, y = -1000000, label = month_abbr[special[1]], size = 20) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(
    color = guide_colorbar(
      barwidth = 20,
      barheight = 0.5
    )
  )


trend_aeqd_without_lat_alt$color_val <- pmax(pmin(trend_aeqd_without_lat_alt[[special[1]]], 0.3), -0.3)
p2 <- ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
  geom_sf(data = trend_aeqd_without_lat_alt, aes(color = color_val), size = 4, shape = 18) +
  geom_sf(data = world_aeqd, fill = NA, color = "black") +
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) +
  scale_color_gradient2(
    low = "red",         
    mid = "white",       
    high = "blue",        
    midpoint = 0,
    limits = c(-0.3, 0.3),
    guide = "colourbar"
  ) +
  theme_minimal() +
  labs(
    title = paste("Trend Mean (with Lat+Alt) for Week", wk_names[special[2]]),
    color = ""
  ) +
  annotate("text", x = -7000000, y = -1000000, label = month_abbr[special[2]], size = 20) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(
    color = guide_colorbar(
      barwidth = 20,
      barheight = 0.5
    )
  )

trend_aeqd_lat_alt$color_val <- pmax(pmin(trend_aeqd_lat_alt[[special[3]]], 0.3), -0.3)
p3 <- ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
  geom_sf(data = trend_aeqd_lat_alt, aes(color = color_val), size = 4, shape = 18) +
  geom_sf(data = world_aeqd, fill = NA, color = "black") +
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) +
  scale_color_gradient2(
    low = "red",         
    mid = "white",       
    high = "blue",        
    midpoint = 0,
    limits = c(-0.3, 0.3),
    guide = "colourbar"
  ) +
  theme_minimal() +
  labs(
    title = paste("Trend Mean (with Lat+Alt) for Week", wk_names[special[3]]),
    color = ""
  ) +
  annotate("text", x = -7000000, y = -1000000, label = month_abbr[special[3]], size = 20) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(
    color = guide_colorbar(
      barwidth = 20,
      barheight = 0.5
    )
  )


trend_aeqd_lat_alt$color_val <- pmax(pmin(trend_aeqd_lat_alt[[special[4]]], 0.3), -0.3)
p4 <- ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
  geom_sf(data = trend_aeqd_lat_alt, aes(color = color_val), size = 4, shape = 18) +
  geom_sf(data = world_aeqd, fill = NA, color = "black") +
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) +
  scale_color_gradient2(
    low = "red",         
    mid = "white",       
    high = "blue",        
    midpoint = 0,
    limits = c(-0.3, 0.3),
    guide = "colourbar"
  ) +
  theme_minimal() +
  labs(
    title = paste("Trend Mean (with Lat+Alt) for Week", wk_names[special[4]]),
    color = ""
  ) +
  annotate("text", x = -7000000, y = -1000000, label = month_abbr[special[4]], size = 20) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(
    color = guide_colorbar(
      barwidth = 20,
      barheight = 0.5
    )
  )

plot(density(diff_mean_INDEP[,-1]))
lines(density(diff_mean_without_lat_alt[,-1]), col = "blue")
plot(density(diff_mean_lat_alt[,-1]), col = "red")

cowplot::plot_grid(p1,p2, nrow = 1)
cowplot::plot_grid(p3,p4, nrow = 1)
