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

elev <- read.csv(here::here("curr_elev.csv"))
elev_nnbs <- read.csv(here::here("nnbs_elev.csv"), sep = "\t", row.names = NULL)[,-5]
colnames(elev_nnbs) <- colnames(elev)
elev <- rbind(elev, elev_nnbs)[,4]
lats <- coords[,2]

sample_idx <- 1:2000

setwd("D:/77/Research/temp/snow/")

theta_bymp <- readRDS("theta01_bym+.Rda")[,sample_idx]
thetas_bymp <- readRDS("theta10_bym+.Rda")[,sample_idx]

theta_bymp2 <- readRDS("D:/77/Research/temp/snow_trend/self_theta_lat+alt_prior_100.Rda")[,sample_idx]
thetas_bymp2 <- readRDS("D:/77/Research/temp/snow_trend/self_theta_10_lat+alt_prior_100.Rda")[,sample_idx]
S <- nrow(y)

theta_01_all <- theta_bymp[1:S,]
theta_02_all <- theta_bymp[(S+1):(2*S),]
theta_11_all <- theta_bymp[(2*S+1):(3*S),]
theta_12_all <- theta_bymp[(3*S+1):(4*S),]
theta_21_all <- theta_bymp[(4*S+1):(5*S),]
theta_22_all <- theta_bymp[(5*S+1):(6*S),]
theta_a1_all <- theta_bymp[(6*S+1):(7*S),]
theta_a2_all <- theta_bymp[(7*S+1):(8*S),]

theta_0L_all <- theta_bymp[8*S+1,]
theta_1L_all <- theta_bymp[8*S+2,]
theta_2L_all <- theta_bymp[8*S+3,]
theta_aL_all <- theta_bymp[8*S+4,]

theta_0A_all <- theta_bymp[8*S+5,]
theta_1A_all <- theta_bymp[8*S+6,]
theta_2A_all <- theta_bymp[8*S+7,]
theta_aA_all <- theta_bymp[8*S+8,]


theta_01s_all <- thetas_bymp[1:S,]
theta_02s_all <- thetas_bymp[(S+1):(2*S),]
theta_11s_all <- thetas_bymp[(2*S+1):(3*S),]
theta_12s_all <- thetas_bymp[(3*S+1):(4*S),]
theta_21s_all <- thetas_bymp[(4*S+1):(5*S),]
theta_22s_all <- thetas_bymp[(5*S+1):(6*S),]
theta_a1s_all <- thetas_bymp[(6*S+1):(7*S),]
theta_a2s_all <- thetas_bymp[(7*S+1):(8*S),]

theta_0Ls_all <- thetas_bymp[8*S+1,]
theta_1Ls_all <- thetas_bymp[8*S+2,]
theta_2Ls_all <- thetas_bymp[8*S+3,]
theta_aLs_all <- thetas_bymp[8*S+4,]

theta_0As_all <- thetas_bymp[8*S+5,]
theta_1As_all <- thetas_bymp[8*S+6,]
theta_2As_all <- thetas_bymp[8*S+7,]
theta_aAs_all <- thetas_bymp[8*S+8,]


S <- nrow(y) - nrow(all_y_nnbs)

theta_01dd_all <- theta_bymp2[1:S,]
theta_02dd_all <- theta_bymp2[(S+1):(2*S),]
theta_11dd_all <- theta_bymp2[(2*S+1):(3*S),]
theta_12dd_all <- theta_bymp2[(3*S+1):(4*S),]
theta_21dd_all <- theta_bymp2[(4*S+1):(5*S),]
theta_22dd_all <- theta_bymp2[(5*S+1):(6*S),]
theta_a1dd_all <- theta_bymp2[(6*S+1):(7*S),]
theta_a2dd_all <- theta_bymp2[(7*S+1):(8*S),]

theta_0Ldd_all <- theta_bymp2[8*S+1,]
theta_1Ldd_all <- theta_bymp2[8*S+2,]
theta_2Ldd_all <- theta_bymp2[8*S+3,]
theta_aLdd_all <- theta_bymp2[8*S+4,]

theta_0Add_all <- theta_bymp2[8*S+5,]
theta_1Add_all <- theta_bymp2[8*S+6,]
theta_2Add_all <- theta_bymp2[8*S+7,]
theta_aAdd_all <- theta_bymp2[8*S+8,]


theta_01sdd_all <- thetas_bymp2[1:S,]
theta_02sdd_all <- thetas_bymp2[(S+1):(2*S),]
theta_11sdd_all <- thetas_bymp2[(2*S+1):(3*S),]
theta_12sdd_all <- thetas_bymp2[(3*S+1):(4*S),]
theta_21sdd_all <- thetas_bymp2[(4*S+1):(5*S),]
theta_22sdd_all <- thetas_bymp2[(5*S+1):(6*S),]
theta_a1sdd_all <- thetas_bymp2[(6*S+1):(7*S),]
theta_a2sdd_all <- thetas_bymp2[(7*S+1):(8*S),]

theta_0Lsdd_all <- thetas_bymp2[8*S+1,]
theta_1Lsdd_all <- thetas_bymp2[8*S+2,]
theta_2Lsdd_all <- thetas_bymp2[8*S+3,]
theta_aLsdd_all <- thetas_bymp2[8*S+4,]

theta_0Asdd_all <- thetas_bymp2[8*S+5,]
theta_1Asdd_all <- thetas_bymp2[8*S+6,]
theta_2Asdd_all <- thetas_bymp2[8*S+7,]
theta_aAsdd_all <- thetas_bymp2[8*S+8,]


checks <- function(A,B){
  par(mfrow = c(1,1))
  hist(A-B)
  hist(A)
  return(list(mean(A-B),sd(A-B),mean(A),mean(B),sd(A),sd(B)))
}

checks(theta_aA_all, theta_aAdd_all)
plot(theta_aA_all, type = 'l', ylim = c(-0.0001,0.0001))
lines(theta_aAdd_all, col = "red")