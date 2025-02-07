rm(list = ls())
library(sf)
library(fields)
library(sp)
library(spdep)
library(fields)
library(spatstat)
library(reshape2)
library(ggplot2)
library(rstan)
library(maps)
library(spBayes)
library(BayesLogit)
library(Matrix)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)


setwd(here::here())
all_y <- readRDS("snow_cleaned.Rda")
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]

ggplot() +
  geom_point(aes( x = coords[,1], y = coords[,2]))


S <- nrow(y)
TT <- ncol(y)
sf_coords <- st_as_sf(data.frame(coords), coords = c("LON", "LAT"), crs = 4326)
sf_data_aeqd <- st_transform(sf_coords, crs = "+proj=aeqd +lat_0=90 +lon_0=-100")

ggplot(data = sf_data_aeqd[1:3,]) +
  geom_sf() +
  theme_minimal() +
  labs(title = "Spatial Data in AEQD Projection",
       x = "Projected X",
       y = "Projected Y")

transformed_coordinates <- st_coordinates(sf_data_aeqd)

dif <- st_coordinates(sf_data_aeqd[2,]) - st_coordinates(sf_data_aeqd[3,])
theta <- atan2(dif[2], dif[1])

rotate_points <- function(coords, angle) {
  # Create a rotation matrix for 2D rotation
  rotation_matrix <- matrix(c(cos(angle), -sin(angle), 
                              sin(angle), cos(angle)), 
                            ncol = 2)
  
  # Apply the rotation to the coordinates
  rotated_coords <- t(rotation_matrix %*% t(coords))
  
  return(rotated_coords)
}

rotated_coordinates <- data.frame(rotate_points(transformed_coordinates, theta))/1e6


ggplot() +
  geom_point(aes( x = rotated_coordinates[,1], y = rotated_coordinates[,2]))+
  coord_fixed(ratio = 1) +
  xlim( c(-10,10)  ) +
  ylim( c(-10,10)  )
  

Distances <- pairdist(rotated_coordinates)
Omg <- Matrix(0, nrow = nrow(coords), ncol = nrow(coords),sparse = TRUE)
Omg[which(Distances <= 0.22)] = 1

num_nbs <- rowSums(Omg) - 1
# 
# ggplot() +
#   geom_point(aes( x = rotated_coordinates[,1], y = rotated_coordinates[,2], color = num_nbs))+
#   coord_fixed(ratio = 1) +
#   xlim( c(-10,10)  ) +
#   ylim( c(-10,10)  ) +
#   scale_color_viridis_c () +  # You can change the option to control color range
#   labs(color = "Number of Neighbors") +  # Optional: Add color legend title
#   theme_minimal()





# saveRDS(snow_dat,"snow_cleaned.Rda")
world <- ne_countries(scale = "medium", returnclass = "sf")
world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]



# Convert the data frame to an sf object
sf_data <- st_as_sf(data.frame(cbind(coords, num_nbs))
                    , coords = c("LON", "LAT"), crs = 4326)


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

ggplot() +
  geom_sf(data = sf_data_aeqd, aes(color = factor(num_nbs), shape = factor(num_nbs)), size = 2) + # Data points
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map borders
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) +
  theme_minimal() +
  labs(
    title = "number_of_neighbors",
    color = "number_of_neighbors"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggplot() +
  geom_sf(data = sf_data_aeqd, aes(color = factor(num_nbs), shape = factor(num_nbs)), size = 2) + # Data points
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map borders
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) +  # Color scale
  scale_shape_manual(values = c(15, 3, 17, 18, 19)) +  # Shape scale
  theme_minimal() +
  labs(
    title = "Number of Neighbors",
    color = "Number of Neighbors",
    shape = "Number of Neighbors"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) 



