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

sample_idx <- seq(from = 5005, to = 10000, by = 5)

setwd("D:/77/Research/temp/snow_trend")

theta <- readRDS("self_theta.Rda")[,sample_idx]
all_tau <- readRDS("self_tau.Rda")[,sample_idx]

theta_s <- readRDS("self_theta_10.Rda")[,sample_idx]
all_taus <- readRDS("self_tau_10.Rda")[,sample_idx]

theta_mean_01 <- apply(theta, 1, mean)
theta_mean_10 <- apply(theta_s, 1, mean)

tau <- apply(all_tau, 1, mean)
taus <- apply(all_taus, 1,mean)

S <- nrow(y)

theta_01_all <- theta[1:S,]
theta_02_all <- theta[(S+1):(2*S),]
theta_11_all <- theta[(2*S+1):(3*S),]
theta_12_all <- theta[(3*S+1):(4*S),]
theta_21_all <- theta[(4*S+1):(5*S),]
theta_22_all <- theta[(5*S+1):(6*S),]
theta_a1_all <- theta[(6*S+1):(7*S),]
theta_a2_all <- theta[(7*S+1):(8*S),]

theta_01s_all <- theta_s[1:S,]
theta_02s_all <- theta_s[(S+1):(2*S),]
theta_11s_all <- theta_s[(2*S+1):(3*S),]
theta_12s_all <- theta_s[(3*S+1):(4*S),]
theta_21s_all <- theta_s[(4*S+1):(5*S),]
theta_22s_all <- theta_s[(5*S+1):(6*S),]
theta_a1s_all <- theta_s[(6*S+1):(7*S),]
theta_a2s_all <- theta_s[(7*S+1):(8*S),]


beta_0_hat <- theta_01_all + theta_02_all
beta_1_hat <- theta_11_all + theta_12_all
beta_2_hat <- theta_21_all + theta_22_all
alpha_hat <- theta_a1_all + theta_a2_all

beta_0s_hat <- theta_01s_all + theta_02s_all
beta_1s_hat <- theta_11s_all + theta_12s_all
beta_2s_hat <- theta_21s_all + theta_22s_all
alpha_s_hat <- theta_a1s_all + theta_a2s_all

# Transaction Matrix

SS <- dim(y)[1]
TT <- dim(y)[2]
period <- 52


P <- array(NA, dim = c(SS, TT-1, 2, 2))
inv_logit <- function(x){return(1/(1+exp(-x)))}

# Winter begins: Week 20 (December 21st).
# Winter midpoint: Week 27 (February 4th).
# Winter ends: Week 33 (March 20th).
winter_weeks <- c(20, 27, 33)
weekly_pred <- array(NA, dim = c(length(sample_idx),SS,length(winter_weeks),2))

for (week_idx in 1:length(winter_weeks)) {
  print(paste("Now predicting week",winter_weeks[week_idx]))
  for (idx in 1:length(sample_idx)) {
    print(idx)
    for (time in 1:(TT-1)) {
      P[, time, 1, 2] <- inv_logit(beta_0_hat[,idx] + beta_1_hat[,idx] * cos(2*pi*time/period) + beta_2_hat[,idx] * sin(2*pi*time/period) + alpha_hat[,idx] * time)
      P[, time, 1, 1] <- 1 - P[, time, 1, 2]
      P[, time, 2, 1] <- inv_logit(beta_0s_hat[,idx] + beta_1s_hat[,idx] * cos(2*pi*time/period) + beta_2s_hat[,idx] * sin(2*pi*time/period)+ alpha_s_hat[,idx] * time)
      P[, time, 2, 2] <- 1 - P[, time, 2, 1]
    }
    
    for (s in 1:SS) {
      curr_p0 <- matrix(c(y[s,1] == 0, y[s,1] == 1), ncol = 1)
      for(t in 1:(winter_weeks[week_idx]+53*52-1)){
        curr_p0 <- P[s,t,,]%*%(curr_p0)
      }
      weekly_pred[idx, s, week_idx, ] <- curr_p0
    }
  }
}



# Select a week and calculate the trend
# Select the second dim on the third dim of array to get P(x_t = 1), i.e., E(x_t)

diffs <- matrix(NA, nrow = SS, ncol = length(winter_weeks))
for (nu in 1:52) {
  week_idx <- seq(from = nu, to = TT, by = 52)
  yearly_pred = yearly_trend <- weekly_pred[,week_idx,2]
  for (yrs in 2:ncol(yearly_pred)) {
    yearly_trend[,yrs] = yearly_pred[,yrs] - yearly_pred[,1]
  }
  diffs[,nu] <- yearly_trend[,ncol(yearly_trend)]
}
colnames(diffs) <- paste0("week", 1:52)


world <- ne_countries(scale = "medium", returnclass = "sf")
world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]

# Convert the data frame to an sf object
sf_data <- st_as_sf(data.frame(cbind(coords, theta_01, theta_02, theta_11, theta_12, theta_21, theta_22, theta_a1, theta_a2,
                                     theta_01s, theta_02s, theta_11s, theta_12s, theta_21s, theta_22s, theta_a1s, theta_a2s, diffs))
                    , coords = c("LON", "LAT"), crs = 4326)

aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"
world_aeqd <- st_transform(world_north, crs = aeqd_proj)
sf_data_aeqd <- st_transform(sf_data, crs = aeqd_proj)
elev_coords <- as.data.frame(cbind(coords, 1:nrow(y)))
colnames(elev_coords) <- c("x", "y","Loc")
elevations <- get_elev_point(locations = elev_coords, prj = 4326)$elevation


equator_points <- data.frame(
  lon = seq(-180, 180, length.out = 200),
  lat = rep(0, 100)  # All points at latitude = 0
)

# Transform equator points to the azimuthal equidistant projection
equator_sf <- st_as_sf(equator_points, coords = c("lon", "lat"), crs = 4326) 
equator_aeqd <- st_transform(equator_sf, crs = aeqd_proj)

p_elev <-
ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
  geom_sf(data = sf_data_aeqd, aes(color = elevations ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_viridis_c(option = "C", direction = -1) + # Vibrant color palette
  theme_minimal() +
  labs(
    title = 'elevation',
    color = 'elevation'
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
i = 1
p_theta_1 <- ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
  geom_sf(data = sf_data_aeqd, aes(color = sf_data_aeqd[[i]] ), size = 3, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_viridis_c(option = "C", direction = -1) + # Vibrant color palette
  theme_minimal() +
  labs(
    title = names(sf_data_aeqd)[i],
    color = names(sf_data_aeqd)[i]
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
cowplot::plot_grid(p_elev, p_theta_1)


# Plot for a single time
for (i in 1:16) {
  # Create the plot
  plot <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = sf_data_aeqd, aes(color = sf_data_aeqd[[i]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_viridis_c(option = "C", direction = -1) + # Vibrant color palette
    theme_minimal() +
    labs(
      title = names(sf_data_aeqd)[i],
      color = names(sf_data_aeqd)[i]
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save the plot
  ggsave(
    filename = paste0(names(sf_data_aeqd)[i], ".png"), # Save as plot_1.png, plot_2.png, ...
    plot = plot,                          # Specify the plot object
    width = 8, height = 6,                # Set width and height
    dpi = 300                             # High resolution
  )
} 

# Trend plot:
# i = 19
# ggplot() +
#   geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
#   geom_sf(data = sf_data_aeqd, aes(color = sf_data_aeqd[[i]]/53), size = 2, shape = 18) + # Data points with color mapped
#   geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
#   geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
#   scale_color_viridis_c(option = "C", direction = -1) + # Vibrant color palette
#   theme_minimal() +
#   labs(
#     title = names(sf_data_aeqd)[i],
#     color = names(sf_data_aeqd)[i]
#   ) +
#   theme(
#     legend.position = "bottom",
#     plot.title = element_text(hjust = 0.5)
#   )


# Directory to store the individual frames
# dir.create("gif_frames")
# color_limits <- range(diffs) / 53
# # Loop to generate plots and save as PNG
# for (i in 1:52 + 16) {
#   # Create the plot
#   plot <- ggplot() +
#     geom_sf(data = sf_data_aeqd, aes(color = sf_data_aeqd[[i]] / 53), size = 2, shape = 18) + # Data points
#     geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map borders
#     geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
#     scale_color_gradient2(
#       low = "blue",         # Negative values
#       high = "red",         # Positive values
#       midpoint = 0,         # Center at 0
#       limits = color_limits, # Use consistent limits for color scale
#       guide = guide_colorbar(barwidth = 20, barheight = 0.5) # Adjust legend bar size
#     ) +
#     theme_minimal() +
#     labs(
#       title = names(sf_data_aeqd)[i],
#       color = names(sf_data_aeqd)[i]
#     ) +
#     theme(
#       legend.position = "bottom",
#       plot.title = element_text(hjust = 0.5)
#     )
#   
#   # Save the plot as a PNG file
#   ggsave(filename = sprintf("gif_frames/frame_%02d.png", i), plot = plot, width = 8, height = 6, dpi = 300)
# }
# 
# 
# # Combine PNGs into a GIF using magick
# frames <- image_read(sprintf("gif_frames/frame_%02d.png", 1:52+16))
# gif <- image_animate(frames, fps = 2)  # Adjust fps for speed of animation
# 
# # Save the GIF
# image_write(gif, path = "animated_map.gif")
# 
# # Cleanup: Remove temporary files (optional)
# unlink("gif_frames", recursive = TRUE)

