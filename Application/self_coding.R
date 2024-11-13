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
setwd(here::here())
all_y <- readRDS("snow_cleaned.Rda")
# states <- c("Idaho","Washington","Oregon","Montana","Wyoming")
# test_idx <- NULL
# library(mvtnorm)
# for (i in 1:length(states)) {
#   tem = spBayes::pointsInPoly(as.matrix(map_data("state", region = states[i])[,1:2]),cbind(all_y$LON, all_y$LAT))
#   test_idx <- unique(c(test_idx, tem))
#   
# }
# 
# test_area <- data.frame(map_data("state", region = states))
# points_inside <- data.frame(all_y[test_idx, 1:2])
# test_period = 1000
# y <- all_y[test_idx,3:(2+test_period)]
# coords <- all_y[test_idx,1:2]

y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]

S <- nrow(y)
TT <- ncol(y)

Distances <- pairdist(coords)
Omg <- Matrix(0, nrow = nrow(coords), ncol = nrow(coords),sparse = TRUE)
Omg[which(Distances <= 2)] = 1
Omg <- Omg - diag(nrow = nrow(coords))
# D <- diag(rowSums(Omg))
D <- Matrix(diag(rowSums(Omg)) ,sparse = TRUE)
sigma = 1
eps = 0.0001
period = 52

location_time_0 <- which(y[,-ncol(y)]==0, arr.ind =  TRUE)
next_y <- y[cbind(location_time_0[,1], location_time_0[,2]+1)]
design_mat <- Matrix(0, nrow = length(next_y), ncol = 8*S, sparse = TRUE)
location_idx <- sparse.model.matrix(~ factor(row) - 1, data = data.frame(location_time_0))

# pb <- txtProgressBar(min = 0, max = nrow(design_mat), style = 3)
covariates <- 
cbind(1,1,
      cos(2*pi*location_time_0[,2]/period),
      cos(2*pi*location_time_0[,2]/period),
      sin(2*pi*location_time_0[,2]/period), 
      sin(2*pi*location_time_0[,2]/period), 
      location_time_0[,2], location_time_0[,2])


for (i in 1:nrow(design_mat)) {
  print(i)
  design_mat[i,] <- as.vector(outer( location_idx[i,],
                                          covariates[i,], "*"))
  # setTxtProgressBar(pb, i)
}

tot_samples <- 1000

all_theta <- matrix(NA, nrow = 8*S, ncol = tot_samples)
all_tau <- matrix(NA, nrow = 8, ncol = tot_samples)
curr_theta_vec <- rep(0, 8*S)
curr_tau_vec <- rep(0.1,8)
a_tau <- 0.001
b_tau <- 0.001
curr_idx <- 0
save_idx <- 0
burn = 1000
thin = 2
while(save_idx < tot_samples) {
  curr_idx = curr_idx + 1
  print(curr_idx)
  # Sample current Omega
  kappas <- next_y - 1/2
  curr_phi <- Matrix(design_mat, sparse = TRUE) %*% curr_theta_vec
  curr_omega <- mapply(function(b, z) rpg(1, h = b, z = z), rep(1, length(next_y)), as.numeric(curr_phi))
  
  # Sample current theta
  curr_B <- bdiag(
    curr_tau_vec[1]*solve(D-Omg+diag(eps, S)),
    curr_tau_vec[2]*diag(1,S),
    curr_tau_vec[3]*solve(D-Omg+diag(eps, S)),
    curr_tau_vec[4]*diag(1,S),
    curr_tau_vec[5]*solve(D-Omg+diag(eps, S)),
    curr_tau_vec[6]*diag(1,S),
    curr_tau_vec[7]*solve(D-Omg+diag(eps, S)),
    curr_tau_vec[8]*diag(1,S)
  )                  
  
  pos_sigma <- solve(t(design_mat)%*%(design_mat * curr_omega) + solve(curr_B) )
  pos_mu <- pos_sigma%*%(t(design_mat)%*%kappas)
  L <- chol(pos_sigma)
  sth <- rnorm(length(pos_mu))
  curr_theta_vec <- t(L) %*% sth + pos_mu
  
 # Sample taus
  curr_tau_vec[1] <- 1/rgamma(1, shape = a_tau+S/2, rate = b_tau + t(curr_theta_vec[1:S])%*%(D-Omg+diag(eps, S))%*%curr_theta_vec[1:S]/2 )
  curr_tau_vec[2] <- 1/rgamma(1, shape = a_tau+S/2, rate = b_tau + t(curr_theta_vec[(S+1):(2*S)])%*%curr_theta_vec[(S+1):(2*S)]/2 )
  
  curr_tau_vec[3] <- 1/rgamma(1, shape = a_tau+S/2, rate = b_tau + t(curr_theta_vec[(2*S+1):(3*S)])%*%(D-Omg+diag(eps, S))%*%curr_theta_vec[(2*S+1):(3*S)]/2 )
  curr_tau_vec[4] <- 1/rgamma(1, shape = a_tau+S/2, rate = b_tau + t(curr_theta_vec[(3*S+1):(4*S)])%*%curr_theta_vec[(3*S+1):(4*S)]/2 )
  
  curr_tau_vec[5] <- 1/rgamma(1, shape = a_tau+S/2, rate = b_tau + t(curr_theta_vec[(4*S+1):(5*S)])%*%(D-Omg+diag(eps, S))%*%curr_theta_vec[(4*S+1):(5*S)] /2)
  curr_tau_vec[6] <- 1/rgamma(1, shape = a_tau+S/2, rate = b_tau + t(curr_theta_vec[(5*S+1):(6*S)])%*%curr_theta_vec[(5*S+1):(6*S)]/2 )
 
  curr_tau_vec[7] <- 1/rgamma(1, shape = a_tau+S/2, rate = b_tau + t(curr_theta_vec[(6*S+1):(7*S)])%*%(D-Omg+diag(eps, S))%*%curr_theta_vec[(6*S+1):(7*S)]/2 )
  curr_tau_vec[8] <- 1/rgamma(1, shape = a_tau+S/2, rate = b_tau + t(curr_theta_vec[(7*S+1):(8*S)])%*%curr_theta_vec[(7*S+1):(8*S)]/2 )

  if((curr_idx > burn) & (curr_idx %% thin == 0) ){
    save_idx <- save_idx + 1
    all_theta[,save_idx] <- as.vector(curr_theta_vec)
    all_tau[,save_idx] <-  as.vector(curr_tau_vec)
  }
  
}
theta_mean <- apply(all_theta, 1, mean)
theta_01 <- theta_mean[1:S]
theta_02 <- theta_mean[(S+1):(2*S)]
theta_11 <- theta_mean[(2*S+1):(3*S)]
theta_12 <- theta_mean[(3*S+1):(4*S)]
theta_21 <- theta_mean[(4*S+1):(5*S)]
theta_22 <- theta_mean[(5*S+1):(6*S)]
theta_a1 <- theta_mean[(6*S+1):(7*S)]
theta_a2 <- theta_mean[(7*S+1):(8*S)]

saveRDS(all_theta, "self_theta.Rda")
saveRDS(all_tau, "self_tau.Rda")
# Compare with stan
# samples <- readRDS(here::here("test_app_samples.Rda"))
# theta_01_stan <- apply(samples$theta_01, 2, mean)
# theta_02_stan <- apply(samples$theta_02, 2, mean)
# theta_11_stan <- apply(samples$theta_11, 2, mean)
# theta_12_stan <- apply(samples$theta_12, 2, mean)
# theta_21_stan <- apply(samples$theta_21, 2, mean)
# theta_22_stan <- apply(samples$theta_22, 2, mean)
# theta_a1_stan <- apply(samples$theta_a1, 2, mean)
# theta_a2_stan <- apply(samples$theta_a2, 2, mean)
# 
# par(mfrow = c(4,2))
# boxplot(cbind(theta_01, theta_01_stan))
# boxplot(cbind(theta_02, theta_02_stan))
# boxplot(cbind(theta_11, theta_11_stan))
# boxplot(cbind(theta_12, theta_12_stan))
# boxplot(cbind(theta_21, theta_21_stan))
# boxplot(cbind(theta_22, theta_22_stan))
# boxplot(cbind(theta_a1, theta_a1_stan))
# boxplot(cbind(theta_a2, theta_a2_stan))

