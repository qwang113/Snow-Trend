rm(list = ls())
setwd("D:/77/Research/temp/snow_trend")
thetas_01_nnbs <- readRDS("self_theta_01_lat+alt_nnbs.Rda")
thetas_01 <- readRDS("self_theta_lat+alt_prior_100.Rda")
last_thing_nnbs <- thetas_01[(nrow(thetas_01_nnbs)-7):nrow(thetas_01_nnbs),]
last_thing <- thetas_01[(nrow(thetas_01)-7):nrow(thetas_01),]
