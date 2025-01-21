library(terra)
library(dplyr)
library(sf)
library(fields)
library(sp)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(spBayes)
library(elevatr)
setwd("D:/77/Research/temp/snow_trend/elevations")
no_nbs <- c(57, 170, 236, 269, 343, 685, 946, 947, 989, 1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203)
all_y <- readRDS(here::here("snow_cleaned.Rda"))[-no_nbs,]
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]

curr_elev <- read.table("1.txt", sep = "")[,-1]
for (i in 2:6) {
  print(i)
  curr_elev <- rbind(curr_elev, read.table(paste(i,".txt", sep = ""), sep = "")[,-1])
}
colnames(curr_elev) <- c("lat","long","elev")
full_dat <- curr_elev[,c(2,1,3)]
world <- ne_countries(scale = "medium", returnclass = "sf")
world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]
world_aeqd <- st_transform(world_north, crs = aeqd_proj)
# Convert the data frame to an sf object
sf_data <- st_as_sf(full_dat, coords = c("long", "lat"), crs = 4326)
aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"
sf_data_aeqd <- st_transform(sf_data, crs = aeqd_proj)
equator_points <- data.frame(
  lon = seq(-180, 180, length.out = 200),
  lat = rep(0, 100)  # All points at latitude = 0
)

# Transform equator points to the azimuthal equidistant projection
equator_sf <- st_as_sf(equator_points, coords = c("lon", "lat"), crs = 4326) 
equator_aeqd <- st_transform(equator_sf, crs = aeqd_proj)



ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
  geom_sf(data = sf_data_aeqd, aes(color = elev ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_viridis_c(option = "C", direction = -1) + # Vibrant color palette
  theme_minimal() +
  labs(
    title = names(sf_data_aeqd)[1],
    color = names(sf_data_aeqd)[1]
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
