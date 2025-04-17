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
all_y <- readRDS("snow_cleaned.Rda")[-no_nbs,]
all_y_nnbs <- readRDS("snow_cleaned.Rda")[no_nbs,]
y <- rbind(all_y, all_y_nnbs)[,-c(1,2)]
coords <- rbind(all_y, all_y_nnbs)[,1:2]

sample_idx <- seq(from = 1005, to = 2000, by = 5)

setwd("D:/77/Research/temp/snow/")
theta <- readRDS("theta01_bym.Rda")[,sample_idx]
theta_s <- readRDS("theta10_bym.Rda")[,sample_idx]

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

saveRDS(weekly_ini_lat_alt, "ini_year_bym.Rda")
saveRDS(weekly_final_lat_alt, "fin_year_bym.Rda")