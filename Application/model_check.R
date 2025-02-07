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

setwd(here::here())
all_y <- readRDS("snow_cleaned.Rda")

states <- c("Idaho","Washington","Oregon","Montana","Wyoming")
test_idx <- NULL
for (i in 1:length(states)) {
  tem = spBayes::pointsInPoly(as.matrix(map_data("state", region = states[i])[,1:2]),cbind(all_y$LON, all_y$LAT))
  test_idx <- unique(c(test_idx, tem))
  
}

test_area <- data.frame(map_data("state", region = states))
points_inside <- data.frame(all_y[test_idx, 1:2])
test_period <- 1000


ggplot() +
  # Add Idaho polygon
  geom_polygon(data = test_area, aes(x = long, y = lat, group = group), fill = "lightblue", color = "black") +
  # Add points that are inside Idaho
  geom_point(data = points_inside, aes(x = LON, y = LAT), color = "red", size =15, shape = 18) +
  # Set plot title and labels
  ggtitle("Map of Test Data with Points Inside") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(ratio = 1) +
  theme_minimal()

y <- all_y[test_idx,3:(2+test_period)]
coords <- all_y[test_idx,1:2]
samples <- readRDS(here::here("test_app_samples.Rda"))


beta_0 <- apply(samples$beta_0, 2, mean)
beta_0s <- apply(samples$beta_0s, 2, mean)
beta_1 <- apply(samples$beta_1, 2, mean)
beta_1s <- apply(samples$beta_1s, 2, mean)
beta_2 <- apply(samples$beta_2, 2, mean)
beta_2s <- apply(samples$beta_2s, 2, mean)
alpha <- apply(samples$alpha, 2, mean)
alpha_s <- apply(samples$alpha_s, 2, mean)

p0 <-
ggplot() +
  geom_polygon(data = test_area, aes(x = long, y = lat, group = group), fill = "#A6CEE3", color = "black") +
  geom_point(data = points_inside, aes(x = LON, y = LAT, color = beta_0), size =15, shape = 18) +
  scale_color_gradient(low = "#FFEDA0", high = "#F03B20") +
  ggtitle("Posterior Mean for Beta_0") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p0s <-
  ggplot() +
  geom_polygon(data = test_area, aes(x = long, y = lat, group = group), fill = "#A6CEE3", color = "black") +
  geom_point(data = points_inside, aes(x = LON, y = LAT, color = beta_0s), size =15, shape = 18) +
  scale_color_gradient(low = "#FFEDA0", high = "#F03B20") +
  ggtitle("Posterior Mean for Beta_0s") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p1 <-
  ggplot() +
  geom_polygon(data = test_area, aes(x = long, y = lat, group = group), fill = "#A6CEE3", color = "black") +
  geom_point(data = points_inside, aes(x = LON, y = LAT, color = beta_1), size =15, shape = 18) +
  scale_color_gradient(low = "#FFEDA0", high = "#F03B20") +
  ggtitle("Posterior Mean for Beta_1") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p1s <-
  ggplot() +
  geom_polygon(data = test_area, aes(x = long, y = lat, group = group), fill = "#A6CEE3", color = "black") +
  geom_point(data = points_inside, aes(x = LON, y = LAT, color = beta_1s), size =15, shape = 18) +
  scale_color_gradient(low = "#FFEDA0", high = "#F03B20") +
  ggtitle("Posterior Mean for Beta_1s") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <-
  ggplot() +
  geom_polygon(data = test_area, aes(x = long, y = lat, group = group), fill = "#A6CEE3", color = "black") +
  geom_point(data = points_inside, aes(x = LON, y = LAT, color = beta_2), size =15, shape = 18) +
  scale_color_gradient(low = "#FFEDA0", high = "#F03B20") +
  ggtitle("Posterior Mean for Beta_2") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2s <-
  ggplot() +
  geom_polygon(data = test_area, aes(x = long, y = lat, group = group), fill = "#A6CEE3", color = "black") +
  geom_point(data = points_inside, aes(x = LON, y = LAT, color = beta_2s), size =15, shape = 18) +
  scale_color_gradient(low = "#FFEDA0", high = "#F03B20") +
  ggtitle("Posterior Mean for Beta_2s") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

pa <-
  ggplot() +
  geom_polygon(data = test_area, aes(x = long, y = lat, group = group), fill = "#A6CEE3", color = "black") +
  geom_point(data = points_inside, aes(x = LON, y = LAT, color = alpha), size =15, shape = 18) +
  scale_color_gradient(low = "#FFEDA0", high = "#F03B20") +
  ggtitle("Posterior Mean for Alpha") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

pas <-
  ggplot() +
  geom_polygon(data = test_area, aes(x = long, y = lat, group = group), fill = "#A6CEE3", color = "black") +
  geom_point(data = points_inside, aes(x = LON, y = LAT, color = alpha_s), size =15, shape = 18) +
  scale_color_gradient(low = "#FFEDA0", high = "#F03B20") +
  ggtitle("Posterior Mean for Alpha_s") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


cowplot::plot_grid(p0,p0s,p1,p1s,p2,p2s,pa,pas, nrow = 4)

