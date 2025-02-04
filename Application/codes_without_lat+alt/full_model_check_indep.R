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
elev <- read.csv(here::here("curr_elev.csv"))[,3]
lats <- coords[,2]
sample_idx <- seq(from = 1005, to = 2000, by = 5)

setwd("D:/77/Research/temp/snow_trend")
theta <- readRDS("self_theta_INDEP.Rda")[,sample_idx]
theta_s <- readRDS("self_theta_10_INDEP.Rda")[,sample_idx]
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



# Transaction Matrix ----------------------------------------------------------------------------Trend

SS <- dim(y)[1]
TT <- dim(y)[2]
period <- 52


P <- array(NA, dim = c(SS, TT-1, 2, 2))
inv_logit <- function(x){return(1/(1+exp(-x)))}

# Winter begins: Week 20 (December 21st).
# Winter midpoint: Week 27 (February 4th).
# Winter ends: Week 33 (March 20th).
winter_weeks <- c(20, 27, 33)

# Prediction for the first year


weekly_ini <- array(NA, dim = c(length(sample_idx),SS,52,2))
weekly_final <- array(NA, dim = c(length(sample_idx),SS,52,2))
for (week_idx in 1:52) {
  for (idx in 1:length(sample_idx)) {
    print(paste("Now doing week", week_idx, "Sample index",idx))
    for (time in 1:(TT-1)) {
      P[, time, 1, 2] <- inv_logit(
        beta_0_hat[,idx] + beta_1_hat[,idx] * cos(2*pi*time/period) + beta_2_hat[,idx] * sin(2*pi*time/period) + alpha_hat[,idx] * time)
      P[, time, 1, 1] <- 1 - P[, time, 1, 2]
      P[, time, 2, 1] <- inv_logit(beta_0s_hat[,idx] + beta_1s_hat[,idx] * cos(2*pi*time/period) + beta_2s_hat[,idx] * sin(2*pi*time/period)+ alpha_s_hat[,idx] * time)
      P[, time, 2, 2] <- 1 - P[, time, 2, 1]
    }
    for (s in 1:SS) {
      # For location s, calculate the first year probability for week week_idx
      curr_p0 <- t(c(y[s,1] == 0, y[s,1] == 1))
      if(week_idx == 1){
        weekly_ini[idx, s, week_idx, ] <- curr_p0
      }else{
        for(t in 1:(week_idx-1)){
          curr_p0 <- curr_p0 %*% P[s,t,,]
        }
        weekly_ini[idx, s, week_idx, ] <- curr_p0
      }
      # For location s, calculate the last year probability for week week_idx
      curr_p0 <- t(c(y[s,1] == 0, y[s,1] == 1))
      for(t in 1:(week_idx+53*52-1)){
        curr_p0 <- curr_p0 %*% P[s,t,,]
      }
      weekly_final[idx, s, week_idx, ] <- curr_p0
    }
  }
}

saveRDS(weekly_ini, "weekly_ini_without_lat+out.Rda")
saveRDS(weekly_final, "weekly_final_without_lat+out.Rda")


# First year
weekly_ini <- readRDS("weekly_ini_without_lat+out.Rda")
weekly_final <- readRDS("weekly_final_without_lat+out.Rda")
p_snow_diff <- weekly_final[,,,2] - weekly_ini[,,,2]
diff_mean <- apply(p_snow_diff,c(2,3),mean)
diff_sd <- apply(p_snow_diff,c(2,3),sd)


# ----------------------------------------------------------------------------Plots
world <- ne_countries(scale = "medium", returnclass = "sf")
world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]

# Convert the data frame to an sf object
sf_data <- st_as_sf(data.frame(cbind(coords, beta_0_hat, beta_1_hat, beta_2_hat, alpha_hat, beta_0s_hat, beta_1s_hat, beta_2s_hat, alpha_s_hat))
                    , coords = c("LON", "LAT"), crs = 4326)

aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"

# Trend
trend_data <- st_as_sf(data.frame(cbind(coords, diff_mean, diff_sd)), coords = c("LON", "LAT"), crs = 4326)
trend_aeqd <- st_transform(trend_data, crs = aeqd_proj)
world_aeqd <- st_transform(world_north, crs = aeqd_proj)
sf_data_aeqd <- st_transform(sf_data, crs = aeqd_proj)

equator_points <- data.frame(
  lon = seq(-180, 180, length.out = 200),
  lat = rep(0, 100)  # All points at latitude = 0
)

# Transform equator points to the azimuthal equidistant projection
equator_sf <- st_as_sf(equator_points, coords = c("lon", "lat"), crs = 4326) 
equator_aeqd <- st_transform(equator_sf, crs = aeqd_proj)


# Plot for a single time
for (i in 1:8) {
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
    ) + 
    guides(
      color = guide_colorbar(
        barwidth = 20,   # Adjust the width of the color bar
        barheight = 0.5  # Adjust the height of the color bar
      )
    )
  
  # Save the plot
  ggsave(
    filename = paste0(names(sf_data_aeqd)[i], ".png"), # Save as plot_1.png, plot_2.png, ...
    plot = plot,                          # Specify the plot object
    width = 8, height = 6,                # Set width and height
    dpi = 300                             # High resolution
  )
} 
library(magick)
# Trend Plot
for (i in 1:52) {
  print(i)
  # Create the plot
  plot1 <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = trend_aeqd, aes(color = trend_aeqd[[i]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_viridis_c(option = "C", direction = -1) + # Vibrant color palette
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
    geom_sf(data = trend_aeqd, aes(color = trend_aeqd[[i+52]] ), size = 2, shape = 18) + # Data points with color mapped
    geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
    geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
    scale_color_viridis_c(option = "C", direction = -1) + # Vibrant color palette
    theme_minimal() +
    labs(
      title = paste("Trend sd for Week0",i),
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
  plot = cowplot::plot_grid(plot1,plot2, nrow = 1)
  # Save the plot
  ggsave(
    filename = paste0("plot_",i, ".png"), # Save as plot_1.png, plot_2.png, ...
    plot = plot,                          # Specify the plot object
    width = 10, height = 6,                # Set width and height
    dpi = 300                             # High resolution
  )
} 
png_files <- png_files[order(as.numeric(gsub("\\D", "", png_files)))]# Find all saved PNGs
gif <- image_read(png_files)                       # Read images
gif <- image_animate(gif, fps = 5)                 # Set frames per second (5 fps)
image_write(gif, "trend_animation.gif")            # Save the GIF




