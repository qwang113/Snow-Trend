rm(list = ls())
library(sf)
library(fields)
library(sp)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
idx <- which(readRDS("github_DF")$Group=="4")
snow_dat <- read.csv(here::here("snow_dat.csv"))[idx,-1]
snow_matrix <- as.matrix(snow_dat[,-(1:2)])
snow_long <- snow_dat$LON
snow_lat <- snow_dat$LAT


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


idx = 1000

# Plot for a single time
ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = "NA") + # World map
  geom_sf(data = sf_data_aeqd, aes(color = factor(sf_data_aeqd[[names(sf_data_aeqd)[idx]]])), size = 1.5, shape = 18) + # Data points with color mapped
  scale_color_manual(values = c("0" = "orange", "1" = "skyblue"),  # Custom color mapping
                     labels = c("No Snow" = "0", "Snow" = "1"),
                     name = "Snow Existance") +
  geom_sf(data = world_aeqd, fill = "NA", color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1)+
  theme_minimal() +
  labs(title = paste("Snow Existence Map from the North Pole","-",substr(names(sf_data_aeqd)[idx],start = 2, stop = 50) ),
       color = "Snow Existence") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))



