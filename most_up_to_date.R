rm(list = ls())
library(ncdf4)
library(terra)
library(dplyr)
setwd("D:/77/Research/temp/snow/")
snow_rast <- rast("snow_raw.nc")
snow_df <- as.data.frame(snow_rast, xy = TRUE, na.rm = FALSE)
xy_df <- data.frame(x = snow_df$x, y = snow_df$y)
xy_vect <- vect(xy_df, geom = c("x", "y"), crs = crs(snow_rast))
xy_lonlat <- project(xy_vect, "EPSG:4326")
coords <- as.data.frame(round(geom(xy_lonlat)[, c("x", "y")],4))
coords_yisu <- readRDS("snow_cleaned.Rda")[,1:2]
colnames(coords) <- c("LON", "LAT")
yisu_mat  <- as.matrix(coords_yisu[, c("LON", "LAT")])
coords_mat <- as.matrix(coords[, c("LON", "LAT")])

yisu_mat   <- as.matrix(coords_yisu[, c("LON", "LAT")])
coords_mat <- as.matrix(coords[, c("LON", "LAT")])

get_nearest <- function(pt) {
  dists <- rowSums((t(t(coords_mat) - pt))^2)
  min_idx <- which.min(dists)
  min_coords <- coords_mat[min_idx, ]
  c(min_coords[1], min_coords[2], sqrt(dists[min_idx]), min_idx)
}

yisu_mat   <- as.matrix(coords_yisu[, c("LON", "LAT")])
coords_mat <- as.matrix(coords[, c("LON", "LAT")])

get_nearest <- function(pt) {
  dists <- rowSums((t(t(coords_mat) - pt))^2)
  min_idx <- which.min(dists)
  min_coords <- coords_mat[min_idx, ]
  c(min_coords[1], min_coords[2], sqrt(dists[min_idx]), min_idx)
}
nearest_results <- t(apply(yisu_mat, 1, get_nearest))
result_df <- cbind(
  coords_yisu,
  nearest_lon = nearest_results[, 1],
  nearest_lat = nearest_results[, 2],
  distance    = nearest_results[, 3],
  coords_row  = nearest_results[, 4]
)

result_df <- as.data.frame(result_df)
result_df$coords_row <- as.integer(result_df$coords_row)
result_df$distance   <- as.numeric(result_df$distance)


snow_df[,1:2] <- coords
colnames(snow_df) <- c("LON", "LAT", as.character(seq(from = as.Date("1966-10-10"), by = "1 week", length.out = 3057)))

start_date <- as.Date("1967-08-07")
weekly_dates <- as.character(seq(from = start_date, by = "1 week", length.out = 57*52))
cleaned_df <- snow_df[,c("LON","LAT",weekly_dates)]



