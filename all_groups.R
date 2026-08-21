rm(list = ls())
library(ncdf4)
library(terra)
library(dplyr)
library(lubridate)
# Set this path to the directory containing snow_raw.nc.
setwd("path/to/raw/snow/data")
snow_rast <- rast("snow_raw.nc")
snow_df <- as.data.frame(snow_rast, xy = TRUE, na.rm = FALSE)
xy_df <- data.frame(x = snow_df$x, y = snow_df$y)
xy_vect <- vect(xy_df, geom = c("x", "y"), crs = crs(snow_rast))
xy_lonlat <- project(xy_vect, "EPSG:4326")
coords <- as.data.frame(round(geom(xy_lonlat)[, c("x", "y")],4))
# Set this path to the local Snow-Trend repository.
setwd("path/to/Snow-Trend")
coords_yisu <- readRDS(file.path("data", "snow_cleaned.Rda"))[,1:2]
colnames(coords) <- c("LON", "LAT")
yisu_mat  <- as.matrix(coords_yisu[, c("LON", "LAT")])
coords_mat <- as.matrix(coords[, c("LON", "LAT")])

yisu_mat   <- as.matrix(coords_yisu[, c("LON", "LAT")])
coords_mat <- as.matrix(coords[, c("LON", "LAT")])

snow_df[,1:2] <- coords
snow_df <- snow_df[which(snow_df$land==1),]
snow_df <- snow_df[,-c(3,4)]
colnames(snow_df) <- c("LON", "LAT", as.character(seq(from = as.Date("1966-10-10"), by = "1 week", length.out = 3070)))

start_date <- as.Date("1972-08-01")
end_date <- as.Date("2024-07-31")
weekly_dates <- as.character(seq(from = as.Date("1967-08-07"), by = "1 week", to = as.Date("2025-05-05")))
dates_truncated <- weekly_dates[weekly_dates >= start_date & weekly_dates <= end_date]


get_annual_label <- function(date) {
  y <- year(date)
  m <- month(date)
  if (m >= 8) {
    return(paste0(y, "-", y + 1))
  } else {
    return(paste0(y - 1, "-", y))
  }
}


df <- data.frame(
  date = dates_truncated,
  year = year(dates_truncated),
  month = month(dates_truncated),
  cycle = sapply(dates_truncated, get_annual_label)
)

cycles_with_53 <- df %>%
  group_by(cycle) %>%
  summarise(n = n()) %>%
  filter(n == 53) %>%
  pull(cycle)

last_july_weeks <- df %>%
  filter(cycle %in% cycles_with_53, month == 7) %>%
  group_by(cycle) %>%
  filter(date == max(date)) %>%
  pull(date)

dates_cleaned <- df %>%
  filter(!(date %in% last_july_weeks)) %>%
  pull(date)

snow_df <- snow_df[,c("LON", "LAT", as.character(dates_cleaned))]
# snow_df[,1:2] <- coords_yisu
saveRDS(snow_df, here::here("snow_all_groups.Rda"))
