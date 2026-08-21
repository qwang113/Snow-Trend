rm(list = ls())
library(sf)
library(fields)
library(sp)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
# Set this path to the directory containing snow_all_groups.Rda.
setwd("path/to/snow/data")
snow_dat <- readRDS("snow_all_groups.Rda")
snow_matrix <- as.matrix(snow_dat[,-(1:2)])
snow_long <- snow_dat$LON
snow_lat <- snow_dat$LAT

# saveRDS(snow_dat,"snow_cleaned.Rda")
world <- ne_countries(scale = "medium", returnclass = "sf")
world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]


# Convert the data frame to an sf object
sf_data <- st_as_sf(snow_dat, coords = c("LON", "LAT"), crs = 4326)

aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"
world_aeqd <- st_transform(world_north, crs = aeqd_proj)
sf_data_aeqd <- st_transform(sf_data, crs = aeqd_proj)

equator_points <- data.frame(
  lon = seq(-180, 180, length.out = 200),
  lat = rep(0, 100)  # All points at latitude = 0
)

# Transform equator points to the azimuthal equidistant projection
equator_sf <- st_as_sf(equator_points, coords = c("lon", "lat"), crs = 4326) 
equator_aeqd <- st_transform(equator_sf, crs = aeqd_proj)


idx = which(colnames(snow_matrix)=="1990-03-05")

# Plot for a single time
varname <- names(sf_data_aeqd)[idx]

sf_data_aeqd$plot_snow <- ifelse(is.na(sf_data_aeqd[[varname]]), 0, sf_data_aeqd[[varname]])
sf_data_aeqd$plot_snow <- factor(sf_data_aeqd$plot_snow,
                                 levels = c(0, 1),
                                 labels = c("0 or NA", "Snow"))

library(ggplot2)

ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
  geom_sf(data = sf_data_aeqd, aes(color = plot_snow), size = 1.5, shape = 18,show.legend = FALSE) +
  scale_color_manual(
    values = c("0 or NA" = "lightgray", "Snow" = "#1f78b4"),
    breaks = c("0 or NA", "Snow"),
    name = "Snow Existence"
  ) +
  geom_sf(data = world_aeqd, fill = NA, color = "black") +
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) +
  theme_minimal() +
  labs(title = "Snow Presence On Mar.5, 1990") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
