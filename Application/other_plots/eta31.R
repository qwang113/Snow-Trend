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

elev <- scale(c(read.csv(here::here("curr_elev.csv"))[,4]
                ,read.csv(here::here("nnbs_elev.csv"), sep = "\t", row.names = NULL)[,4]))
lats <- scale(coords[,2])

sample_idx <- seq(from = 204, to = 1000, by = 4)

setwd("D:/77/Research/temp/snow/")

theta <- readRDS("theta01_bym+notime.Rda")[,sample_idx]
theta_s <- readRDS("theta10_bym+notime.Rda")[,sample_idx]
S <- nrow(y)

theta_01_all <- theta[1:S,]
theta_02_all <- theta[(S+1):(2*S),]
theta_11_all <- theta[(2*S+1):(3*S),]
theta_12_all <- theta[(3*S+1):(4*S),]
theta_21_all <- theta[(4*S+1):(5*S),]
theta_22_all <- theta[(5*S+1):(6*S),]
theta_a1_all <- theta[(6*S+1):(7*S),]
theta_a2_all <- theta[(7*S+1):(8*S),]

theta_0L_all <- theta[8*S+1,]
theta_1L_all <- theta[8*S+2,]
theta_2L_all <- theta[8*S+3,]
theta_aL_all <- theta[8*S+4,]

theta_0A_all <- theta[8*S+5,]
theta_1A_all <- theta[8*S+6,]
# theta_2A_all <- theta[8*S+7,]
# theta_aA_all <- theta[8*S+8,]

theta_01s_all <- theta_s[1:S,]
theta_02s_all <- theta_s[(S+1):(2*S),]
theta_11s_all <- theta_s[(2*S+1):(3*S),]
theta_12s_all <- theta_s[(3*S+1):(4*S),]
theta_21s_all <- theta_s[(4*S+1):(5*S),]
theta_22s_all <- theta_s[(5*S+1):(6*S),]
theta_a1s_all <- theta_s[(6*S+1):(7*S),]
theta_a2s_all <- theta_s[(7*S+1):(8*S),]

theta_0Ls_all <- theta_s[8*S+1,]
theta_1Ls_all <- theta_s[8*S+2,]
theta_2Ls_all <- theta_s[8*S+3,]
theta_aLs_all <- theta_s[8*S+4,]

theta_0As_all <- theta_s[8*S+5,]
theta_1As_all <- theta_s[8*S+6,]
# theta_2As_all <- theta_s[8*S+7,]
# theta_aAs_all <- theta_s[8*S+8,]



aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"


# Trend
eta31 <- st_as_sf(data.frame(cbind(coords, apply(theta_a1_all,1,mean))), coords = c("LON", "LAT"), crs = 4326)
eta31 <- st_transform(eta31, crs = aeqd_proj)

eta31s <- st_as_sf(data.frame(cbind(coords, apply(theta_a1s_all,1,mean))), coords = c("LON", "LAT"), crs = 4326)
eta31s <- st_transform(eta31s, crs = aeqd_proj)


world <- ne_countries(scale = "medium", returnclass = "sf")
world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]
world_aeqd <- st_transform(world_north, crs = aeqd_proj)




equator_points <- data.frame(
  lon = seq(-180, 180, length.out = 200),
  lat = rep(0, 100)  # All points at latitude = 0
)

# Transform equator points to the azimuthal equidistant projection
equator_sf <- st_as_sf(equator_points, coords = c("lon", "lat"), crs = 4326) 
equator_aeqd <- st_transform(equator_sf, crs = aeqd_proj)

plot1 <- ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
  geom_sf(data = eta31, aes(color = eta31[[1]] ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_gradient2( low = "red",          # Color for negative values
                         mid = "white",  # Color for zero
                         high = "blue") +
  theme_minimal() +
  labs(
    title ="eta31 Map",
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
  geom_sf(data = eta31s, aes(color = eta31s[[1]] ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_gradient2(limit = c(-0.002,0.002), low = "red",          # Color for negative values
                        mid = "white",  # Color for zero
                        high = "blue") +
  theme_minimal() +
  labs(
    title ="eta31s Map",
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

cowplot::plot_grid(plot1, plot2)

