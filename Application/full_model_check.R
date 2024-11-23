rm(list = ls())

library(sf)
library(fields)
library(sp)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)



all_y <- readRDS(here::here("snow_cleaned.Rda"))
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]

setwd("D:/77/Research/temp/snow_trend")
theta <- readRDS("self_theta.Rda")
all_tau <- readRDS("self_tau.Rda")

theta_s <- readRDS("self_theta_10.Rda")
all_taus <- readRDS("self_tau_10.Rda")

theta_mean_01 <- apply(theta, 1, mean)
theta_mean_10 <- apply(theta_s, 1, mean)

tau <- apply(all_tau, 1, mean)
taus <- apply(all_taus, 1,mean)

S <- nrow(y)

theta_01 <- theta_mean_01[1:S]
theta_02 <- theta_mean_01[(S+1):(2*S)]
theta_11 <- theta_mean_01[(2*S+1):(3*S)]
theta_12 <- theta_mean_01[(3*S+1):(4*S)]
theta_21 <- theta_mean_01[(4*S+1):(5*S)]
theta_22 <- theta_mean_01[(5*S+1):(6*S)]
theta_a1 <- theta_mean_01[(6*S+1):(7*S)]
theta_a2 <- theta_mean_01[(7*S+1):(8*S)]

theta_01s <- theta_mean_10[1:S]
theta_02s <- theta_mean_10[(S+1):(2*S)]
theta_11s <- theta_mean_10[(2*S+1):(3*S)]
theta_12s <- theta_mean_10[(3*S+1):(4*S)]
theta_21s <- theta_mean_10[(4*S+1):(5*S)]
theta_22s <- theta_mean_10[(5*S+1):(6*S)]
theta_a1s <- theta_mean_10[(6*S+1):(7*S)]
theta_a2s <- theta_mean_10[(7*S+1):(8*S)]


beta_0_hat <- theta_01 + theta_02
beta_1_hat <- theta_11 + theta_12
beta_2_hat <- theta_21 + theta_22
alpha_hat <- theta_a1 + theta_a2

beta_0s_hat <- theta_01s + theta_02s
beta_1s_hat <- theta_11s + theta_12s
beta_2s_hat <- theta_21s + theta_22s
alpha_s_hat <- theta_a1s + theta_a2s



# Transaction Matrix

SS <- dim(y)[1]
TT <- dim(y)[2]
period <- 52

P <- array(NA, dim = c(SS, TT-1, 2, 2))
inv_logit <- function(x){return(1/(1+exp(-x)))}
for (time in 1:(TT-1)) {
  P[, time, 1, 2] <- inv_logit(beta_0_hat + beta_1_hat * cos(2*pi*time/period) + beta_2_hat * sin(2*pi*time/period) + alpha_hat * time)
  P[, time, 1, 1] <- 1 - P[, time, 1, 2]
  P[, time, 2, 1] <- inv_logit(beta_0s_hat + beta_1s_hat * cos(2*pi*time/period) + beta_2s_hat * sin(2*pi*time/period)+ alpha_s_hat * time)
  P[, time, 2, 2] <- 1 - P[, time, 2, 1]
}

# Calculate weekly prediction
weekly_pred <- array(NA, dim = c(dim(y),2))
for (s in 1:SS) {
  weekly_pred[s,1,] = c(y[s,1] == 0, y[s,1] == 1)
}

for (s in 1:SS) {
  for(t in 1:(TT-1)){
    weekly_pred[s,t+1,] <- t(weekly_pred[s,t,])%*%P[s,t,,]
  }
}

# Select a week and calculate the trend
# Select the second dim on the third dim of array to get P(x_t = 1), i.e., E(x_t)
nu = 1
week_idx <- seq(from = nu, to = TT, by = 52)
yearly_pred = yearly_trend <- weekly_pred[,week_idx,2]
for (yrs in 2:ncol(yearly_pred)) {
  yearly_trend[,yrs] = yearly_pred[,yrs] - yearly_pred[,1]
}
last_year <- yearly_trend[,ncol(yearly_trend)]

diffs <- last_year - yearly_pred[,1]
diffs[which(diffs <= -0.1)]


# saveRDS(snow_dat,"snow_cleaned.Rda")
world <- ne_countries(scale = "medium", returnclass = "sf")
world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]



# Convert the data frame to an sf object
sf_data <- st_as_sf(data.frame(cbind(coords, theta_01, theta_02, theta_11, theta_12, theta_21, theta_22, theta_a1, theta_a2,
                                     theta_01s, theta_02s, theta_11s, theta_12s, theta_21s, theta_22s, theta_a1s, theta_a2s))
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

# Plot for a single time
for (i in 1:16) {
  # Create the plot
  plot <- ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
    geom_sf(data = sf_data_aeqd, aes(color = sf_data_aeqd[[i]]), size = 2, shape = 18) + # Data points with color mapped
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




ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
  geom_sf(data = sf_data_aeqd, aes(color = sf_data_aeqd[[i]]), size = 2, shape = 18) + # Data points with color mapped
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





