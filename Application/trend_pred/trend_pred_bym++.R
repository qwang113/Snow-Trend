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
all_y <- readRDS("snow_cleaned_full.Rda")[-no_nbs,]
all_y_nnbs <- readRDS("snow_cleaned_full.Rda")[no_nbs,]
y <- rbind(all_y, all_y_nnbs)[,-c(1,2)]
coords <- rbind(all_y, all_y_nnbs)[,1:2]

elev <- scale(c(read.csv(here::here("curr_elev.csv"))[,4]
        ,read.csv(here::here("nnbs_elev.csv"), sep = "\t", row.names = NULL)[,4]))
lats <- scale(coords[,2])

sample_idx <- seq(from = 1204, to = 2000, by = 4)

setwd("D:/77/Research/temp/snow/")

theta <- readRDS("theta01_bym++.Rda")[,sample_idx]
theta_s <- readRDS("theta10_bym++.Rda")[,sample_idx]
S <- nrow(y)

theta_01_all <- theta[1:S,]
theta_02_all <- theta[(S+1):(2*S),]
theta_11_all <- theta[(2*S+1):(3*S),]
theta_12_all <- theta[(3*S+1):(4*S),]
theta_21_all <- theta[(4*S+1):(5*S),]
theta_22_all <- theta[(5*S+1):(6*S),]
theta_a1_all <- theta[(6*S+1):(7*S),]
theta_a2_all <- theta[(7*S+1):(8*S),]
theta_E_all <- theta[(8*S+1),]
theta_L_all <- theta[(8*S+2),]



theta_01s_all <- theta_s[1:S,]
theta_02s_all <- theta_s[(S+1):(2*S),]
theta_11s_all <- theta_s[(2*S+1):(3*S),]
theta_12s_all <- theta_s[(3*S+1):(4*S),]
theta_21s_all <- theta_s[(4*S+1):(5*S),]
theta_22s_all <- theta_s[(5*S+1):(6*S),]
theta_a1s_all <- theta_s[(6*S+1):(7*S),]
theta_a2s_all <- theta_s[(7*S+1):(8*S),]
theta_Es_all <- theta_s[(8*S+1),]
theta_Ls_all <- theta_s[(8*S+2),]



beta_0_hat <- theta_01_all + theta_02_all
beta_1_hat <- theta_11_all + theta_12_all
beta_2_hat <- theta_21_all + theta_22_all
alpha_hat <- theta_a1_all + theta_a2_all
beta_E_hat <- theta_E_all
beta_L_hat <- theta_L_all


beta_0s_hat <- theta_01s_all + theta_02s_all 
beta_1s_hat <- theta_11s_all + theta_12s_all
beta_2s_hat <- theta_21s_all + theta_22s_all
alpha_s_hat <- theta_a1s_all + theta_a2s_all
beta_Es_hat <- theta_Es_all
beta_Ls_hat <- theta_Ls_all


# Transaction Matrix ----------------------------------------------------------------------------Trend

SS <- dim(y)[1]
TT <- dim(y)[2]
period <- 52


P <- array(NA, dim = c(SS, TT-1, 2, 2))
inv_logit <- function(x){return(1/(1+exp(-x)))}

# Prediction for the first year
weekly_ini_lat_alt <- array(NA, dim = c(length(sample_idx),SS,52,2))
weekly_final_lat_alt <- array(NA, dim = c(length(sample_idx),SS,52,2))
curr_idx <- 1
for (idx in curr_idx:length(sample_idx)) {
  for (time in 1:(TT-1)) {
    P[, time, 1, 2] <- inv_logit(
      beta_0_hat[,idx] +
        beta_1_hat[,idx] * cos(2*pi*time/period) +
        beta_2_hat[,idx] * sin(2*pi*time/period) +
        alpha_hat[,idx] * time +
        lats*beta_L_hat[idx] +
        elev*beta_E_hat[idx])
    
    P[, time, 1, 1] <- 1 - P[, time, 1, 2]
    
    P[, time, 2, 1] <- inv_logit(
      beta_0s_hat[,idx] +
        beta_1s_hat[,idx] * cos(2*pi*time/period) +
        beta_2s_hat[,idx] * sin(2*pi*time/period) +
        alpha_s_hat[,idx] * time + 
        lats*beta_Ls_hat[idx] +
        elev*beta_Es_hat[idx])
    
    P[, time, 2, 2] <- 1 - P[, time, 2, 1]
    
  }
  
  for (week_idx in 1:52) {
    print(paste("Now doing week", week_idx, "Sample index",idx))
    for (s in 1:SS) {
      # For location s, calculate the first year probability for week week_idx
      curr_p0 <- t(c(y[s,1] == 0, y[s,1] == 1))
      if(week_idx == 1){
        weekly_ini_lat_alt[idx, s, week_idx, ] <- curr_p0
      }else{
        for(t in 1:(week_idx-1)){
          curr_p0 <- curr_p0 %*% P[s,t,,]
        }
        weekly_ini_lat_alt[idx, s, week_idx, ] <- curr_p0
      }
      # For location s, calculate the last year probability for week week_idx
      curr_p0 <- t(c(y[s,1] == 0, y[s,1] == 1))
      for(t in 1:(week_idx+51*52-1)){
        curr_p0 <- curr_p0 %*% P[s,t,,]
      }
      weekly_final_lat_alt[idx, s, week_idx, ] <- curr_p0
    }
  }
}
setwd("D:/77/Research/temp/snow/")
saveRDS(weekly_ini_lat_alt, "ini_year_bym++.Rda")
saveRDS(weekly_final_lat_alt, "fin_year_bym++.Rda")
