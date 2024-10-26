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
  set.seed(0)
  
  inv_logit <- function(x){return(1/(1+exp(-x)))}
  
  eps = 0.0001
  long <- runif(500)
  lat <- runif(500)
  period <- 50
  coords <- cbind(long, lat)
  Distances <- pairdist(coords)
  Omg <- matrix(0, nrow = nrow(coords), ncol = nrow(coords))
  Omg[which(Distances <= 0.05)] = 1
  Omg <- Omg - diag(nrow = nrow(coords))
  # D <- diag(rowSums(Omg))
  D <- diag(rowSums(Omg)) + diag(eps, nrow = length(rowSums(Omg)))
  
  
  BYM_mat <- (D-Omg)
  tau_01 = .001
  tau_02 = .01
  
  tau_01s = .001
  tau_02s = .01
  
  tau_11 = .001
  tau_12 = .01
  
  tau_11s = .001
  tau_12s = .01
  
  tau_21 = .001
  tau_22 = .01
  
  tau_21s = .001
  tau_22s = .01
  
  tau_a1 = 1e-7
  tau_a2 = 1e-7
  
  tau_a1s = 1e-7
  tau_a2s = 1e-7
  
  
  theta_01 <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_01)
  theta_01s <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_01s)
  
  theta_02 <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_02)
  theta_02s <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_02s)
  # // Define BYM model parameters for cosine
  theta_11 <- rmvnorm(1, rep(2, nrow(coords)), sigma = solve(BYM_mat)*tau_11)
  theta_11s <- rmvnorm(1, rep(-2, nrow(coords)), sigma = solve(BYM_mat)*tau_11s)
  
  theta_12 <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_12)
  theta_12s <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_12s)
  
  # // Define BYM model parameters for sine
  theta_21 <- rmvnorm(1, rep(5, nrow(coords)), sigma = solve(BYM_mat)*tau_21)
  theta_21s <- rmvnorm(1, rep(-5, nrow(coords)), sigma = solve(BYM_mat)*tau_21s)
  
  theta_22 <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_22)
  theta_22s <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_22s)
  
  # // Fixed time trend coefficient
  theta_a1 <- rmvnorm(1, rep(0.01, nrow(coords)), sigma = solve(BYM_mat)*tau_a1)
  theta_a1s <- rmvnorm(1, rep(-0.01, nrow(coords)), sigma = solve(BYM_mat)*tau_a1s)
  
  theta_a2 <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_a2)
  theta_a2s <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_a2s)
  
  beta_0 <-  theta_01 + theta_02
  beta_0s <-  theta_01s + theta_02s
  beta_1 <-  theta_11 + theta_12
  beta_1s <-  theta_11s + theta_12s
  beta_2 <-  theta_21 + theta_22
  beta_2s <-  theta_21s + theta_22s
  alpha_0 <-  theta_a1 + theta_a2
  alpha_0s <-  theta_a1s + theta_a2s
  
  SS <- nrow(coords)
  time_step <- 1:500
  TT <- length(time_step)
  
  P <- array(NA, dim = c(SS, TT-1, 2, 2))
  
  for (time in 1:(TT-1)) {
    P[, time, 1, 2] <- inv_logit(beta_0 + beta_1 * cos(2*pi*time/period) + beta_2 * sin(2*pi*time/period) + alpha_0 * time)
    P[, time, 1, 1] <- 1 - P[, time, 1, 2]
    P[, time, 2, 1] <- inv_logit(beta_0s + beta_1s * cos(2*pi*time/period) + beta_2s * sin(2*pi*time/period)+ alpha_0s * time)
    P[, time, 2, 2] <- 1 - P[, time, 2, 1]
  }
  
  y <- matrix(NA, nrow = SS, ncol = TT)
  y[,1] <- rbinom(SS, 1, 0.5)
  for (time in 2:TT) {
    y[,time] <- rbinom(SS, 1, (y[,time-1]==0) * P[, time-1, 1, 2] +  (y[,time-1]==1) * P[, time-1, 2, 2])
  }
  par(mfrow = c(4,1))
  plot(y[1,], type = 'l')
  plot(y[2,], type = 'l')
  plot(y[3,], type = 'l')
  plot(y[4,], type = 'l')
  
  
  
  
  # # Test whether theta is spatially correlated
  # test_thetas <- data.frame("long" = coords[,1], "lat" = coords[,2], "value" = as.vector(alpha_0))
  # ggplot(test_thetas, aes(x = long, y = lat, color = value)) +
  #   geom_point() +
  #   scale_color_gradient(low = "blue", high = "yellow") +
  #   labs(title = "Heatmap of Values by Location", x = "Longitude", y = "Latitude", color = "Value") +
  #   theme_minimal()
  
