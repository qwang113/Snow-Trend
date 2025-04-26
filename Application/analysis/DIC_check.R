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

elev <- c(read.csv(here::here("curr_elev.csv"))[,4]
          ,read.csv(here::here("nnbs_elev.csv"), sep = "\t", row.names = NULL)[,4])
lats <- coords[,2]

sample_idx <- seq(from = 204, to = 1000, by = 4)

setwd("D:/77/Research/temp/snow/")

# No time likelihood
theta <- readRDS("theta01_bym+notime.Rda")[,sample_idx]
theta_s <- readRDS("theta10_bym+notime.Rda")[,sample_idx]
S <- nrow(y)

theta_01_all <- theta[1:S,]
theta_02_all <- theta[(S+1):(2*S),]
theta_11_all <- theta[(2*S+1):(3*S),]
theta_12_all <- theta[(3*S+1):(4*S),]
theta_21_all <- theta[(4*S+1):(5*S),]
theta_22_all <- theta[(5*S+1):(6*S),]
theta_a1_all <- theta[(6*S+1):(7*S),]
theta_a2_all <- theta[(7*S+1):(8*S),]

theta_0L_all <- theta[8*S+1,]
theta_1L_all <- theta[8*S+2,]
theta_2L_all <- theta[8*S+3,]
theta_0A_all <- theta[8*S+4,]
theta_1A_all <- theta[8*S+5,]
theta_2A_all <- theta[8*S+6,]


theta_01s_all <- theta_s[1:S,]
theta_02s_all <- theta_s[(S+1):(2*S),]
theta_11s_all <- theta_s[(2*S+1):(3*S),]
theta_12s_all <- theta_s[(3*S+1):(4*S),]
theta_21s_all <- theta_s[(4*S+1):(5*S),]
theta_22s_all <- theta_s[(5*S+1):(6*S),]
theta_a1s_all <- theta_s[(6*S+1):(7*S),]
theta_a2s_all <- theta_s[(7*S+1):(8*S),]

theta_0Ls_all <- theta_s[8*S+1,]
theta_1Ls_all <- theta_s[8*S+2,]
theta_2Ls_all <- theta_s[8*S+3,]
theta_0As_all <- theta_s[8*S+4,]
theta_1As_all <- theta_s[8*S+5,]
theta_2As_all <- theta_s[8*S+6,]



beta_0_hat <- theta_01_all + theta_02_all + lats%*%t(theta_0L_all) + elev%*%t(theta_0A_all)
beta_1_hat <- theta_11_all + theta_12_all + lats%*%t(theta_1L_all) + elev%*%t(theta_1A_all)
beta_2_hat <- theta_21_all + theta_22_all + lats%*%t(theta_2L_all) + elev%*%t(theta_2A_all)
alpha_hat <- theta_a1_all + theta_a2_all

beta_0s_hat <- theta_01s_all + theta_02s_all + lats%*%t(theta_0Ls_all) + elev%*%t(theta_0As_all)
beta_1s_hat <- theta_11s_all + theta_12s_all + lats%*%t(theta_1Ls_all) + elev%*%t(theta_1As_all)
beta_2s_hat <- theta_21s_all + theta_22s_all + lats%*%t(theta_2Ls_all) + elev%*%t(theta_2As_all)
alpha_s_hat <- theta_a1s_all + theta_a2s_all

inv_logit <- function(x){return(1/(1+exp(-x)))}
my_param_notime <- my_param <- cbind(rowMeans(beta_0_hat),
                  rowMeans(beta_1_hat),
                  rowMeans(beta_2_hat),
                  rowMeans(alpha_hat),
                  rowMeans(beta_0s_hat),
                  rowMeans(beta_1s_hat),
                  rowMeans(beta_2s_hat),
                  rowMeans(alpha_s_hat)
                  )

SS <- dim(y)[1]
TT <- dim(y)[2]
all_llh <- rep(NA, SS)
period <- 52
  P <- array(NA, dim = c(SS, TT-1, 2, 2))
  for (time in 1:(TT-1)) {
    P[, time, 1, 2] <- inv_logit(
      my_param[,1] + my_param[,2]* cos(2*pi*time/period) + my_param[,3] * sin(2*pi*time/period) + my_param[,4] * time)
    P[, time, 1, 1] <- 1 - P[, time, 1, 2]
    P[, time, 2, 1] <- inv_logit(
      my_param[,5] + my_param[,6] * cos(2*pi*time/period) + my_param[,7] * sin(2*pi*time/period)+ my_param[,8] * time)
    P[, time, 2, 2] <- 1 - P[, time, 2, 1]
  }
for (s in 1:SS) {
  print(s)
  P_curr <- P[s,,,]
  location_time_0 <- which(y[s,-ncol(y)]==0, arr.ind =  TRUE)
  next_y <- y[cbind(s, location_time_0[,2]+1)]
  
  llh <- 0
  for (i in 1:length(next_y)) {
    llh <- llh + ifelse(next_y[i]==0, log(P_curr[location_time_0[i,2],1,1]), log(P_curr[location_time_0[i,2],1,2]))
  }
  
  
  location_time_1<- which(y[s,-ncol(y)]==1, arr.ind =  TRUE)
  next_y <- y[cbind(s, location_time_1[,2]+1)]
  for (i in 1:length(next_y)) {
    llh <- llh + ifelse(next_y[i]==0, log(P_curr[location_time_1[i,2],2,1]), log(P_curr[location_time_1[i,2],2,2]))
  }
  all_llh[s] <- llh
}

notime_llh <- all_llh

# Time included llh

theta <- readRDS("theta01_bym+withtime.Rda")[,sample_idx]
theta_s <- readRDS("theta10_bym+withtime.Rda")[,sample_idx]
S <- nrow(y)

theta_01_all <- theta[1:S,]
theta_02_all <- theta[(S+1):(2*S),]
theta_11_all <- theta[(2*S+1):(3*S),]
theta_12_all <- theta[(3*S+1):(4*S),]
theta_21_all <- theta[(4*S+1):(5*S),]
theta_22_all <- theta[(5*S+1):(6*S),]
theta_a1_all <- theta[(6*S+1):(7*S),]
theta_a2_all <- theta[(7*S+1):(8*S),]

theta_0L_all <- theta[8*S+1,]
theta_1L_all <- theta[8*S+2,]
theta_2L_all <- theta[8*S+3,]
theta_aL_all <- theta[8*S+4,]

theta_0A_all <- theta[8*S+5,]
theta_1A_all <- theta[8*S+6,]
theta_2A_all <- theta[8*S+7,]
theta_aA_all <- theta[8*S+8,]

theta_01s_all <- theta_s[1:S,]
theta_02s_all <- theta_s[(S+1):(2*S),]
theta_11s_all <- theta_s[(2*S+1):(3*S),]
theta_12s_all <- theta_s[(3*S+1):(4*S),]
theta_21s_all <- theta_s[(4*S+1):(5*S),]
theta_22s_all <- theta_s[(5*S+1):(6*S),]
theta_a1s_all <- theta_s[(6*S+1):(7*S),]
theta_a2s_all <- theta_s[(7*S+1):(8*S),]

theta_0Ls_all <- theta_s[8*S+1,]
theta_1Ls_all <- theta_s[8*S+2,]
theta_2Ls_all <- theta_s[8*S+3,]
theta_aLs_all <- theta_s[8*S+4,]

theta_0As_all <- theta_s[8*S+5,]
theta_1As_all <- theta_s[8*S+6,]
theta_2As_all <- theta_s[8*S+7,]
theta_aAs_all <- theta_s[8*S+8,]


beta_0_hat <- theta_01_all + theta_02_all + lats%*%t(theta_0L_all) + elev%*%t(theta_0A_all)
beta_1_hat <- theta_11_all + theta_12_all + lats%*%t(theta_1L_all) + elev%*%t(theta_1A_all)
beta_2_hat <- theta_21_all + theta_22_all + lats%*%t(theta_2L_all) + elev%*%t(theta_2A_all)
alpha_hat <- theta_a1_all + theta_a2_all + lats%*%t(theta_aL_all) + elev%*%t(theta_aA_all)

beta_0s_hat <- theta_01s_all + theta_02s_all + lats%*%t(theta_0Ls_all) + elev%*%t(theta_0As_all)
beta_1s_hat <- theta_11s_all + theta_12s_all + lats%*%t(theta_1Ls_all) + elev%*%t(theta_1As_all)
beta_2s_hat <- theta_21s_all + theta_22s_all + lats%*%t(theta_2Ls_all) + elev%*%t(theta_2As_all)
alpha_s_hat <- theta_a1s_all + theta_a2s_all + lats%*%t(theta_aLs_all) + elev%*%t(theta_aAs_all)

inv_logit <- function(x){return(1/(1+exp(-x)))}
my_param_withtime <- my_param <- cbind(rowMeans(beta_0_hat),
                  rowMeans(beta_1_hat),
                  rowMeans(beta_2_hat),
                  rowMeans(alpha_hat),
                  rowMeans(beta_0s_hat),
                  rowMeans(beta_1s_hat),
                  rowMeans(beta_2s_hat),
                  rowMeans(alpha_s_hat)
)

SS <- dim(y)[1]
TT <- dim(y)[2]
all_llh <- rep(NA, SS)
period <- 52
P <- array(NA, dim = c(SS, TT-1, 2, 2))
for (time in 1:(TT-1)) {
  P[, time, 1, 2] <- inv_logit(
    my_param[,1] + my_param[,2]* cos(2*pi*time/period) + my_param[,3] * sin(2*pi*time/period) + my_param[,4] * time)
  P[, time, 1, 1] <- 1 - P[, time, 1, 2]
  P[, time, 2, 1] <- inv_logit(
    my_param[,5] + my_param[,6] * cos(2*pi*time/period) + my_param[,7] * sin(2*pi*time/period)+ my_param[,8] * time)
  P[, time, 2, 2] <- 1 - P[, time, 2, 1]
}
for (s in 1:SS) {
  print(s)
  P_curr <- P[s,,,]
  location_time_0 <- which(y[s,-ncol(y)]==0, arr.ind =  TRUE)
  next_y <- y[cbind(s, location_time_0[,2]+1)]
  
  llh <- 0
  for (i in 1:length(next_y)) {
    llh <- llh + ifelse(next_y[i]==0, log(P_curr[location_time_0[i,2],1,1]), log(P_curr[location_time_0[i,2],1,2]))
  }
  
  
  location_time_1<- which(y[s,-ncol(y)]==1, arr.ind =  TRUE)
  next_y <- y[cbind(s, location_time_1[,2]+1)]
  for (i in 1:length(next_y)) {
    llh <- llh + ifelse(next_y[i]==0, log(P_curr[location_time_1[i,2],2,1]), log(P_curr[location_time_1[i,2],2,2]))
  }
  all_llh[s] <- llh
}
withtime_llh <- all_llh

