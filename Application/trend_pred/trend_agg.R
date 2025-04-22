rm(list = ls())
library(sf)
library(fields)
library(sp)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(spBayes)
library(elevatr)
library(abind)
no_nbs <- c(57, 170, 236, 269, 343, 685, 946, 947, 989, 1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203)
all_y <- readRDS(here::here("snow_cleaned.Rda"))[-no_nbs,]
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]
coords_nnbs <- readRDS(here::here("snow_cleaned.Rda"))[no_nbs,1:2]
coords <- rbind(coords, coords_nnbs)
# First year
setwd("D:/77/Research/temp/snow")
weekly_ini_without_lat_alt <- readRDS("ini_year_bym.Rda")
weekly_final_without_lat_alt <- readRDS("fin_year_bym.Rda")

weekly_ini_lat_alt <- readRDS("ini_year_bym+.Rda")
weekly_final_lat_alt <- readRDS("fin_year_bym+.Rda")

weekly_ini_INDEP <- readRDS("ini_year_ind.Rda")
weekly_final_INDEP <- readRDS("fin_year_ind.Rda")

# weekly_ini_lat_alt_prior_100 <- readRDS("weekly_ini_with_lat+out_prior_100.Rda")
# weekly_final_lat_alt_prior_100 <- readRDS("weekly_final_with_lat+out_prior_100.Rda")

dpd_snow_diff <- (apply(weekly_final_without_lat_alt[,,,2],c(1,2),sum) - apply(weekly_ini_without_lat_alt[,,,2], c(1,2),sum))/53
dpd_diff_mean_without_lat_alt <- apply(dpd_snow_diff,2,mean)
dpd_diff_sd_without_lat_alt <- apply(dpd_snow_diff,2,sd)

dpd_snow_diff <- (apply(weekly_final_lat_alt[1:50,,,2],c(1,2),sum) - apply(weekly_ini_lat_alt[1:50,,,2], c(1,2),sum))/53
dpd_diff_mean_lat_alt <- apply(dpd_snow_diff,2,mean)
dpd_diff_sd_lat_alt <- apply(dpd_snow_diff,2,sd)

dpd_snow_diff <- (apply(weekly_final_INDEP[,,,2],c(1,2),sum) - apply(weekly_ini_INDEP[,,,2], c(1,2),sum))/53
dpd_diff_mean_INDEP <- apply(dpd_snow_diff,2,mean)
dpd_diff_sd_INDEP <- apply(dpd_snow_diff,2,sd)

# ----------------------------------------------------------------------------Plots
# Convert the data frame to an sf object
aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"


# Trend
dpd_without_lat_alt <- st_as_sf(data.frame(cbind(coords, dpd_diff_mean_without_lat_alt, dpd_diff_sd_without_lat_alt)), coords = c("LON", "LAT"), crs = 4326)
dpd_aeqd_without_lat_alt <- st_transform(dpd_without_lat_alt, crs = aeqd_proj)


dpd_lat_alt <- st_as_sf(data.frame(cbind(coords, dpd_diff_mean_lat_alt, dpd_diff_sd_lat_alt)), coords = c("LON", "LAT"), crs = 4326)
dpd_aeqd_lat_alt <- st_transform(dpd_lat_alt, crs = aeqd_proj)

dpd_INDEP <- st_as_sf(data.frame(cbind(coords, dpd_diff_mean_INDEP, dpd_diff_sd_INDEP)), coords = c("LON", "LAT"), crs = 4326)
dpd_aeqd_INDEP <- st_transform(dpd_INDEP, crs = aeqd_proj)

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

rg_mean <- range(dpd_diff_mean_INDEP, dpd_diff_mean_lat_alt,dpd_diff_mean_without_lat_alt)
rg_sd <- log(range(dpd_diff_sd_INDEP, dpd_diff_sd_lat_alt, dpd_diff_sd_without_lat_alt))

plot1 <- ggplot() +
  geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
  geom_sf(data = dpd_aeqd_without_lat_alt, aes(color = dpd_aeqd_without_lat_alt[[1]] ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_gradient2(
    low = "red",          # Color for negative values
    mid = "white",  # Color for zero
    high = "blue",         # Color for positive values
    midpoint = 0,          # Set midpoint at zero
    # Set the range of the legend
    limit = rg_mean,
    guide = "colourbar"
  ) +
  theme_minimal() +
  labs(
    title ="Aggregated Trend - SP",
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
  geom_sf(data = dpd_aeqd_lat_alt, aes(color = dpd_aeqd_lat_alt[[1]] ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_gradient2(
    low = "red",          # Color for negative values
    mid = "white",  # Color for zero
    high = "blue",         # Color for positive values
    midpoint = 0,          # Set midpoint at zero
    # Set the range of the legend
    limit = rg_mean,
    guide = "colourbar"
  ) +
  theme_minimal() +
  labs(
    title ="Aggregated Trend - SP + Covariates",
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
  geom_sf(data = dpd_aeqd_INDEP, aes(color = dpd_aeqd_INDEP[[1]] ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_gradient2(
    low = "red",          # Color for negative values
    mid = "white",  # Color for zero
    high = "blue",         # Color for positive values
    midpoint = 0,          # Set midpoint at zero
    # Set the range of the legend
    limit = rg_mean,
    guide = "colourbar"
  ) +
  theme_minimal() +
  labs(
    title ="Aggregated Trend - IND",
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
  geom_sf(data = dpd_aeqd_without_lat_alt, aes(color = log(dpd_aeqd_without_lat_alt[[2]]) ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_viridis_c(option = "C", direction = 1, limits = rg_sd)+

  theme_minimal() +
  labs(
    title ="Aggregated Trend log of SD - SP",
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
  geom_sf(data = dpd_aeqd_lat_alt, aes(color = log(dpd_aeqd_lat_alt[[2]]) ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_viridis_c(option = "C", direction = 1, limits = rg_sd)+

  theme_minimal() +
  labs(
    title ="Aggregated Trend log of SD - SP + Covariates",
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
  geom_sf(data = dpd_aeqd_INDEP, aes(color = log(dpd_aeqd_INDEP[[2]]) ), size = 2, shape = 18) + # Data points with color mapped
  geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
  geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
  scale_color_viridis_c(option = "C", direction = 1, limits = rg_sd) +
  theme_minimal() +
  labs(
    title ="Aggregated Trend log of SD - IND",
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

plot = cowplot::plot_grid(plot3,plot1, plot2, plot6, plot4, plot5, nrow = 2)

plot = cowplot::plot_grid(plot3, plot1, plot2, nrow = 1)
plot = cowplot::plot_grid(plot6, plot4, plot5, nrow = 1)
plot

# diff_mean_indep_bym <- dpd_diff_mean_without_lat_alt - dpd_diff_mean_INDEP
# diff_sd_indep_bym <- (dpd_diff_sd_INDEP - dpd_diff_sd_lat_alt)/dpd_diff_sd_INDEP
# 
# diff_mean_indep_bym_sf <- st_as_sf(data.frame(cbind(coords, diff_mean_indep_bym)), coords = c("LON", "LAT"), crs = 4326)
# diff_mean_indep_bym_sf  <- st_transform(diff_mean_indep_bym_sf , crs = aeqd_proj)
# 
# diff_sd_indep_bym_sf <- st_as_sf(data.frame(cbind(coords, diff_sd_indep_bym)), coords = c("LON", "LAT"), crs = 4326)
# diff_sd_indep_bym_sf  <- st_transform(diff_sd_indep_bym_sf , crs = aeqd_proj)
# 
# 
# 
# ggplot() +
#   geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
#   geom_sf(data = diff_mean_indep_bym_sf, aes(color = diff_mean_indep_bym_sf[[1]] ), size = 2, shape = 18) + # Data points with color mapped
#   geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
#   geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
#   scale_color_gradient2(
#     low = "red",          # Color for negative values
#     mid = "white",  # Color for zero
#     high = "blue",         # Color for positive values
#     midpoint = 0,    # Set midpoint at zero
#     guide = "colourbar"
#   ) +
#   theme_minimal() +
#   labs(
#     title ="Different in Predicted Mean Snowy Days per Decade (Spatial Model - Independent Model) ",
#     color = ""
#   ) +
#   theme(
#     legend.position = "bottom",
#     plot.title = element_text(hjust = 0.5)
#   ) + 
#   guides(
#     color = guide_colorbar(
#       barwidth = 20,   # Adjust the width of the color bar
#       barheight = 0.5  # Adjust the height of the color bar
#     )
#   )
# 
# 
# ggplot() +
#   geom_sf(data = world_aeqd, fill = "lightgray", color = NA) + # World map
#   geom_sf(data = diff_sd_indep_bym_sf, aes(color = diff_sd_indep_bym_sf[[1]] ), size = 2, shape = 18) + # Data points with color mapped
#   geom_sf(data = world_aeqd, fill = NA, color = "black") + # World map
#   geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1) + # Equator
#   scale_color_gradient2(
#     low = "red",          # Color for negative values
#     mid = "white",  # Color for zero
#     high = "blue",         # Color for positive values
#     midpoint = 0,    # Set midpoint at zero
#     guide = "colourbar"
#   ) +
#   theme_minimal() +
#   labs(
#     title ="Change in Relative SD",
#     color = ""
#   ) +
#   theme(
#     legend.position = "bottom",
#     plot.title = element_text(hjust = 0.5)
#   ) + 
#   guides(
#     color = guide_colorbar(
#       barwidth = 20,   # Adjust the width of the color bar
#       barheight = 0.5  # Adjust the height of the color bar
#     )
#   )
# 
# 
# data <- data.frame(
#   x = dpd_diff_sd_INDEP,
#   y = dpd_diff_sd_without_lat_alt
# )
# 
# # Plot using ggplot2
# ggplot(data, aes(x = x, y = y)) +
#   geom_point() +  # Scatter plot
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 2, alpha = 0.6) +  # y = x line
#   labs(x = "SD for Independent Model Trend Estimation", y = "SD for Spatial Model Trend Estimation",
#        title = "") +
#   theme_minimal()

# Tables
tb <-
t(
  cbind(
c(table(dpd_diff_mean_INDEP < 0),
  round(mean(dpd_diff_mean_INDEP),4),
  round(median(dpd_diff_mean_INDEP),4), 
  round(mean(dpd_diff_sd_INDEP),4), 
  round(mean(dpd_diff_mean_INDEP/dpd_diff_sd_INDEP),4),
  length(which(dpd_diff_mean_INDEP/dpd_diff_sd_INDEP>2)),
  length(which(dpd_diff_mean_INDEP/dpd_diff_sd_INDEP < -2))
),

c(table(dpd_diff_mean_without_lat_alt < 0),
  round(mean(dpd_diff_mean_without_lat_alt),4),
  round(median(dpd_diff_mean_without_lat_alt),4),
  round(mean(dpd_diff_sd_without_lat_alt),4),
  round(mean(dpd_diff_mean_without_lat_alt/dpd_diff_sd_without_lat_alt),4),
  length(which(dpd_diff_mean_without_lat_alt/dpd_diff_sd_without_lat_alt > 2)),
  length(which(dpd_diff_mean_without_lat_alt/dpd_diff_sd_without_lat_alt < -2))
),

c(table(dpd_diff_mean_lat_alt < 0),
  round(mean(dpd_diff_mean_lat_alt),4),
  round(median(dpd_diff_mean_lat_alt),4),
  round(mean(dpd_diff_sd_lat_alt),4),
  round(mean(dpd_diff_mean_lat_alt/dpd_diff_sd_lat_alt),4),
  length(which(dpd_diff_mean_lat_alt/dpd_diff_sd_lat_alt>2)),
  length(which(dpd_diff_mean_lat_alt/dpd_diff_sd_lat_alt < -2))
  )
))

colnames(tb) <- c("Increase","Decrease","Mean","Median","SD","Z-Score Mean","z>2","z<-2")
rownames(tb) <- c("IND","BYM","BYM+")
knitr::kable(tb, align = 'c', format = "latex")
