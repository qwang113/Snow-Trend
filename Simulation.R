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
eps = 0.0001
long <- runif(20)
lat <- runif(20)
coords <- expand.grid(long, lat)
Distances <- pairdist(coords)
Omg <- matrix(0, nrow = nrow(coords), ncol = nrow(coords))
Omg[which(Distances <= 0.05)] = 1
Omg <- Omg - diag(nrow = nrow(coords))
D <- diag(rowSums(Omg)) + diag(eps, nrow = length(rowSums(Omg)))

BYM_mat <- D-Omg
tau_0 = 1
tau_0s = 1
tau_1 = 2
tau_1s = 2
tau_2 = 3
tau_2s = 3
tau_a = 4
tau_as = 4

theta_01 <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_0)
theta_01s <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_0s)
theta_02 <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_0)
theta_02s <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_0s)
# // Define BYM model parameters for cosine
theta_11 <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_1)
theta_11s <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_1s)
theta_12 <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_1)
theta_12s <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_1s)

# // Define BYM model parameters for sine
theta_21 <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_2)
theta_21s <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_2s)
theta_22 <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_2)
theta_22s <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_2s)

# // Fixed time trend coefficient
theta_a1 <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_a)
theta_a1s <- rmvnorm(1, rep(0, nrow(coords)), sigma = solve(BYM_mat)*tau_as)
theta_a2 <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_a)
theta_a2s <- rmvnorm(1, rep(0, nrow(coords)), sigma = diag(nrow(coords))*tau_as)

# // Hyperparameter for hierarchical prior
