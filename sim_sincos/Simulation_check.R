library(ggplot2)
rm(list = ls())
all_y <- readRDS(here::here("sim_y.Rda"))
y <- all_y[,-c(1,2)]
beta_0 <- readRDS("beta_0.Rda")
beta_1 <- readRDS("beta_1.Rda")
beta_0s <- readRDS("beta_0s.Rda")
beta_1s <- readRDS("beta_1s.Rda")
beta_2 <- readRDS("beta_2.Rda")
beta_2s <- readRDS("beta_2s.Rda")
alpha_0 <- readRDS("alpha_0.Rda")
alpha_0s <- readRDS("alpha_0s.Rda")
samples <- readRDS("samples.Rda")
beta_0_hat <- apply(samples$beta_0, 2, mean)
beta_0s_hat <- apply(samples$beta_0s, 2, mean)
beta_1_hat <- apply(samples$beta_1, 2, mean)
beta_1s_hat <- apply(samples$beta_1s, 2, mean)
beta_2_hat <- apply(samples$beta_2, 2, mean)
beta_2s_hat <- apply(samples$beta_2s, 2, mean)
alpha_0_hat <- apply(samples$alpha, 2, mean)
alpha_0s_hat <- apply(samples$alpha_s, 2, mean)
SS <- dim(y)[1]
TT <- dim(y)[2]
period <- 20

P <- array(NA, dim = c(SS, TT-1, 2, 2))
inv_logit <- function(x){return(1/(1+exp(-x)))}
for (time in 1:(TT-1)) {
  P[, time, 1, 2] <- inv_logit(beta_0_hat + beta_1_hat * cos(2*pi*time/period) + beta_2_hat * sin(2*pi*time/period) + alpha_0_hat * time)
  P[, time, 1, 1] <- 1 - P[, time, 1, 2]
  P[, time, 2, 1] <- inv_logit(beta_0s_hat + beta_1s_hat * cos(2*pi*time/period) + beta_2s_hat * sin(2*pi*time/period)+ alpha_0s_hat * time)
  P[, time, 2, 2] <- 1 - P[, time, 2, 1]
}
y_hat <- matrix(NA, nrow = SS, ncol = TT)
y_hat[,1] <- y[,1]
for (time in 2:TT) {
  y_hat[,time] <- rbinom(SS, 1, (y_hat[,time-1]==0) * P[, time-1, 1, 2] + (y_hat[,time-1]==1) * P[, time-1, 2, 2])
}

ggplot() +
  geom_path(aes(x = 1:length(beta_0), y = t(beta_2))) +
  geom_path(aes(x = 1:length(beta_0), y = beta_2_hat), color = "blue")

