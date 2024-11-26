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
library(future.apply)
library(pbapply)
setwd(here::here())
all_y <- readRDS("snow_cleaned.Rda")
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]

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
slopes <- (sf_data_aeqd[2,]- sf_data_aeqd[3,])





