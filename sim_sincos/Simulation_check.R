rm(list = ls())
library(ggplot2)
setwd(here::here("./sim_sincos/"))
all_y <- readRDS("sim_y.Rda")
y <- all_y[,-c(1,2)]
beta_0 <- readRDS("beta_0.Rda")
beta_1 <- readRDS("beta_1.Rda")
beta_0s <- readRDS("beta_0s.Rda")
beta_1s <- readRDS("beta_1s.Rda")
beta_2 <- readRDS("beta_2.Rda")
beta_2s <- readRDS("beta_2s.Rda")
alpha <- readRDS("alpha.Rda")
alpha_s <- readRDS("alpha_s.Rda")
samples <- readRDS("samples.Rda")
beta_0_hat <- apply(samples$beta_0, 2, mean)
beta_0s_hat <- apply(samples$beta_0s, 2, mean)
beta_1_hat <- apply(samples$beta_1, 2, mean)
beta_1s_hat <- apply(samples$beta_1s, 2, mean)
beta_2_hat <- apply(samples$beta_2, 2, mean)
beta_2s_hat <- apply(samples$beta_2s, 2, mean)
alpha_hat <- apply(samples$alpha, 2, mean)
alpha_s_hat <- apply(samples$alpha_s, 2, mean)
SS <- dim(y)[1]
TT <- dim(y)[2]
period <- 20

# P <- array(NA, dim = c(SS, TT-1, 2, 2))
# inv_logit <- function(x){return(1/(1+exp(-x)))}
# for (time in 1:(TT-1)) {
#   P[, time, 1, 2] <- inv_logit(beta_0_hat + beta_1_hat * cos(2*pi*time/period) + beta_2_hat * sin(2*pi*time/period) + alpha_hat * time)
#   P[, time, 1, 1] <- 1 - P[, time, 1, 2]
#   P[, time, 2, 1] <- inv_logit(beta_0s_hat + beta_1s_hat * cos(2*pi*time/period) + beta_2s_hat * sin(2*pi*time/period)+ alpha_s_hat * time)
#   P[, time, 2, 2] <- 1 - P[, time, 2, 1]
# }
# y_hat <- matrix(NA, nrow = SS, ncol = TT)
# y_hat[,1] <- y[,1]
# for (time in 2:TT) {
#   y_hat[,time] <- rbinom(SS, 1, (y_hat[,time-1]==0) * P[, time-1, 1, 2] + (y_hat[,time-1]==1) * P[, time-1, 2, 2])
# }

# Beta_0 is got by beta_0_hat
# Beta_0s is got by beta_2_hat
# Beta_1 is got by beta_0s_hat
# Beta_1s is got by beta_2s_hat
# Beta_2 is got by beta_1_hat
# Beta_2s is got by alpha_hat
# Alpha is got by beta_1s_hat
# Alpha_s is got by alpha_s_hat

all_dat <- data.frame("beta_0"= as.vector(beta_0), "beta_0s" = as.vector(beta_0s),
                      "beta_1"= as.vector(beta_1), "beta_1s" = as.vector(beta_1s),
                      "beta_2"= as.vector(beta_2), "beta_2s" = as.vector(beta_2s),
                      "alpha"= as.vector(alpha), "alpha_s" = as.vector(alpha_s),
                      
                      "beta_0_hat"= as.vector(beta_0_hat), "beta_0s_hat" = as.vector(beta_0s_hat),
                      "beta_1_hat"= as.vector(beta_1_hat), "beta_1s_hat" = as.vector(beta_1s_hat),
                      "beta_2_hat"= as.vector(beta_2_hat), "beta_2s_hat" = as.vector(beta_2s_hat),
                      "alpha_hat"= as.vector(alpha_hat), "alpha_s_hat" = as.vector(alpha_s_hat))


p1 <- ggplot(data = all_dat) +
  geom_path(aes(x = 1:length(beta_0), y = beta_0, color = "beta_0")) +
  geom_path(aes(x = 1:length(beta_0), y = beta_0_hat, color = "beta_0_hat")) +
  geom_ribbon(aes(x = 1:length(beta_0), ymin = apply(samples$beta_0,2, quantile, 0.05), 
                  ymax = apply(samples$beta_0,2, quantile, 0.95), fill = "beta_0_hat"), alpha = 0.3)+
  labs(color = "Parameter", y = "Value", x = "") +
  theme(legend.position = "bottom")

p2 <- ggplot(data = all_dat) +
  geom_path(aes(x = 1:length(beta_0), y = beta_0s, color = "beta_0s")) +
  geom_path(aes(x = 1:length(beta_0), y = beta_0s_hat, color = "beta_0s_hat")) +
  geom_ribbon(aes(x = 1:length(beta_0), ymin = apply(samples$beta_0s,2, quantile, 0.05), 
                  ymax = apply(samples$beta_0s,2, quantile, 0.95), fill = "beta_0s_hat"), alpha = 0.3)+
  labs(color = "Parameter", y = "Value", x = "") +
  theme(legend.position = "bottom")

p3 <- ggplot(data = all_dat) +
  geom_path(aes(x = 1:length(beta_0), y = beta_1, color = "beta_1")) +
  geom_path(aes(x = 1:length(beta_0), y = beta_1_hat, color = "beta_1_hat")) +
  geom_ribbon(aes(x = 1:length(beta_0), ymin = apply(samples$beta_1,2, quantile, 0.05), 
                  ymax = apply(samples$beta_1,2, quantile, 0.95), fill = "beta_1_hat"), alpha = 0.3)+
  labs(color = "Parameter", y = "Value", x = "") +
  theme(legend.position = "bottom")

p4 <- ggplot(data = all_dat) +
  geom_path(aes(x = 1:length(beta_0), y = beta_1s, color = "beta_1s")) +
  geom_path(aes(x = 1:length(beta_0), y = beta_1s_hat, color = "beta_1s_hat")) +
  geom_ribbon(aes(x = 1:length(beta_0), ymin = apply(samples$beta_1s,2, quantile, 0.05), 
                  ymax = apply(samples$beta_1s,2, quantile, 0.95), fill = "beta_1s_hat"), alpha = 0.3)+
  labs(color = "Parameter", y = "Value", x = "") +
  theme(legend.position = "bottom")

p5 <- ggplot(data = all_dat) +
  geom_path(aes(x = 1:length(beta_0), y = beta_2, color = "beta_2")) +
  geom_path(aes(x = 1:length(beta_0), y = beta_2_hat, color = "beta_2_hat")) +
  geom_ribbon(aes(x = 1:length(beta_0), ymin = apply(samples$beta_2,2, quantile, 0.05), 
                  ymax = apply(samples$beta_2,2, quantile, 0.95), fill = "beta_2_hat"), alpha = 0.3)+
  labs(color = "Parameter", y = "Value", x = "") +
  theme(legend.position = "bottom")

p6 <- ggplot(data = all_dat) +
  geom_path(aes(x = 1:length(beta_0), y = beta_2s, color = "beta_2s")) +
  geom_path(aes(x = 1:length(beta_0), y = beta_2s_hat, color = "beta_2s_hat")) +
  geom_ribbon(aes(x = 1:length(beta_0), ymin = apply(samples$beta_2s,2, quantile, 0.05), 
                  ymax = apply(samples$beta_2s,2, quantile, 0.95), fill = "beta_2s_hat"), alpha = 0.3)+
  labs(color = "Parameter", y = "Value", x = "") +
  theme(legend.position = "bottom")

p7 <- ggplot(data = all_dat) +
  geom_path(aes(x = 1:length(beta_0), y = alpha, color = "alpha")) +
  geom_path(aes(x = 1:length(beta_0), y = alpha_hat, color = "alpha_hat")) +
  geom_ribbon(aes(x = 1:length(beta_0), ymin = apply(samples$alpha,2, quantile, 0.05), 
                  ymax = apply(samples$alpha,2, quantile, 0.95), fill = "alpha_hat"), alpha = 0.3)+
  labs(color = "Parameter", y = "Value", x = "") +
  theme(legend.position = "bottom")


p8 <- ggplot(data = all_dat) +
  geom_path(aes(x = 1:length(beta_0), y = alpha_s, color = "alpha_s")) +
  geom_path(aes(x = 1:length(beta_0), y = alpha_s_hat, color = "alpha_s_hat")) +
  geom_ribbon(aes(x = 1:length(beta_0), ymin = apply(samples$alpha_s,2, quantile, 0.05), 
                  ymax = apply(samples$alpha_s,2, quantile, 0.95), fill = "alpha_s_hat"), alpha = 0.3)+
  labs(color = "Parameter", y = "Value", x = "") +
  theme(legend.position = "bottom")

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
ggplot()+
  geom_point(aes(x = all_y[,1], y = all_y[,2], color = t(beta_0)))
