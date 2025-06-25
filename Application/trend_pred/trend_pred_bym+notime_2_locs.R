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

target_points <- list(
  # Red trend (negative values)
  russia_novosibirsk = c(82.9204, 55.0302),
  china_urumqi = c(87.6168, 43.8256),
  greenland_nuuk = c(-51.7216, 64.1835),
  norway_tromso = c(18.9553, 69.6496),
  usa_fairbanks = c(-147.7164, 64.8378),
  
  # Blue trend (positive values)
  russia_yakutsk = c(129.7331, 62.0355),
  canada_fort_mcmurray = c(-111.3800, 56.7260),
  japan_sapporo = c(141.3545, 43.0621),
  china_dandong = c(124.3547, 40.0005),
  usa_boston = c(-71.0589, 42.3601)
)


nearest_indices <- sapply(target_points, function(loc) {
  dists <- apply(coords, 1, function(row) sqrt(sum((row - loc)^2)))
  which.min(dists)
})



ids <- nearest_indices


sample_idx <- seq(from = 204, to = 1000, by = 4)

setwd("D:/77/Research/temp/snow/")

theta <- readRDS("theta01_bym+notime.Rda")[,sample_idx]
theta_s <- readRDS("theta10_bym+notime.Rda")[,sample_idx]
S <- nrow(y)

theta_01_all <- theta[ids,]
theta_02_all <- theta[S+ids,]
theta_11_all <- theta[2*S+ids,]
theta_12_all <- theta[3*S+ids,]
theta_21_all <- theta[4*S+ids,]
theta_22_all <- theta[5*S+ids,]
theta_a1_all <- theta[6*S+ids,]
theta_a2_all <- theta[7*S+ids,]

theta_0L_all <- theta[8*S+1,]
theta_1L_all <- theta[8*S+2,]
theta_2L_all <- theta[8*S+3,]
theta_0A_all <- theta[8*S+4,]
theta_1A_all <- theta[8*S+5,]
theta_2A_all <- theta[8*S+6,]


theta_01s_all <- theta_s[ids,]
theta_02s_all <- theta_s[S+ids,]
theta_11s_all <- theta_s[2*S+ids,]
theta_12s_all <- theta_s[3*S+ids,]
theta_21s_all <- theta_s[4*S+ids,]
theta_22s_all <- theta_s[5*S+ids,]
theta_a1s_all <- theta_s[6*S+ids,]
theta_a2s_all <- theta_s[7*S+ids,]

theta_0Ls_all <- theta_s[8*S+1,]
theta_1Ls_all <- theta_s[8*S+2,]
theta_2Ls_all <- theta_s[8*S+3,]
theta_0As_all <- theta_s[8*S+4,]
theta_1As_all <- theta_s[8*S+5,]
theta_2As_all <- theta_s[8*S+6,]



beta_0_hat <- theta_01_all + theta_02_all + lats[ids]%*%t(theta_0L_all) + elev[ids]%*%t(theta_0A_all)
beta_1_hat <- theta_11_all + theta_12_all + lats[ids]%*%t(theta_1L_all) + elev[ids]%*%t(theta_1A_all)
beta_2_hat <- theta_21_all + theta_22_all + lats[ids]%*%t(theta_2L_all) + elev[ids]%*%t(theta_2A_all)
alpha_hat <- theta_a1_all + theta_a2_all

beta_0s_hat <- theta_01s_all + theta_02s_all + lats[ids]%*%t(theta_0Ls_all) + elev[ids]%*%t(theta_0As_all)
beta_1s_hat <- theta_11s_all + theta_12s_all + lats[ids]%*%t(theta_1Ls_all) + elev[ids]%*%t(theta_1As_all)
beta_2s_hat <- theta_21s_all + theta_22s_all + lats[ids]%*%t(theta_2Ls_all) + elev[ids]%*%t(theta_2As_all)
alpha_s_hat <- theta_a1s_all + theta_a2s_all



# Transaction Matrix ----------------------------------------------------------------------------Trend

SS <- length(ids)
TT <- dim(y)[2]
period <- 52


P <- array(NA, dim = c(SS, TT-1, 2, 2))
inv_logit <- function(x){return(1/(1+exp(-x)))}

# Prediction for the first year
weekly_pred_lat_alt <- array(NA, dim = c(length(sample_idx),SS,TT))
weekly_pred_lat_alt[,,1] <- y[ids,1]

curr_idx <- 1
for (idx in curr_idx:length(sample_idx)) {
  for (time in 1:(TT-1)) {
    P[, time, 1, 2] <- inv_logit(
      beta_0_hat[,idx] + beta_1_hat[,idx] * cos(2*pi*time/period) + beta_2_hat[,idx] * sin(2*pi*time/period) + alpha_hat[,idx] * time)
    P[, time, 1, 1] <- 1 - P[, time, 1, 2]
    P[, time, 2, 1] <- inv_logit(beta_0s_hat[,idx] + beta_1s_hat[,idx] * cos(2*pi*time/period) + beta_2s_hat[,idx] * sin(2*pi*time/period)+ alpha_s_hat[,idx] * time)
    P[, time, 2, 2] <- 1 - P[, time, 2, 1]
  }
  for (s in 1:SS) {
    curr_idx <- ids[s]
    curr_p0 <- t(c(y[curr_idx,1] == 0, y[curr_idx,1] == 1))
    for(t in 1:(TT-1)){
      curr_p0 <- curr_p0 %*% P[s,t,,]
      weekly_pred_lat_alt[idx,s,t+1] <- curr_p0[,2]
    }
  }
}

num_groups <- dim(weekly_pred_lat_alt)[3] / 52

summed_array <- array(0, dim = c(200, length(ids), num_groups))

for (i in 1:num_groups) {
  idx_start <- (i - 1) * 52 + 1
  idx_end <- i * 52
  summed_array[,,i] <- apply(weekly_pred_lat_alt[,,idx_start:idx_end], c(1, 2), sum)
}


y_sub <- y[ids, ]
summed_y <- matrix(0, nrow = nrow(y_sub), ncol = num_groups)

for (i in 1:num_groups) {
  idx_start <- (i - 1) * 52 + 1
  idx_end <- i * 52
  summed_y[, i] <- rowSums(y_sub[, idx_start:idx_end])
}

library(ggplot2)
library(dplyr)
library(tidyr)

# Combine both locations
locations <- names(target_points)
plots <- list()

for (i in 1:SS) {
  print(i)
  pred_samples <- summed_array[, i, ]
  pred_df <- as.data.frame(t(pred_samples))
  colnames(pred_df) <- paste0("draw_", 1:200)
  pred_df$year <- 1:52
  
  summary_df <- pred_df %>%
    pivot_longer(cols = starts_with("draw_"), names_to = "draw", values_to = "value") %>%
    group_by(year) %>%
    summarise(
      mean = mean(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    )
  
  summary_df$truth <- summed_y[i, ]
  
  p <- ggplot(summary_df, aes(x = year)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
    geom_line(aes(y = mean), color = "blue", linewidth = 1) +
    geom_line(aes(y = truth), color = "red", size = 0.8) +
    labs(
      title = paste(locations[i]),
      x = "Year",
      y = "Number of Snowy Weeks"
    ) +
    # lims(y = c(0,52)) +
    theme_minimal()
  
  plots[[i]] <- p
}

# Display both plots
cowplot::plot_grid(plotlist = plots, nrow = 2)






setwd("D:/77/Research/temp/snow/")
saveRDS(weekly_pred_lat_alt, "fb_hk_pred.Rda")

