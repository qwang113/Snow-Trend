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
  # china_urumqi = c(87.6168, 43.8256),
  # greenland_nuuk = c(-51.7216, 64.1835),
  # norway_tromso = c(18.9553, 69.6496),
  # usa_fairbanks = c(-147.7164, 64.8378),
  
  # Blue trend (positive values)
  # russia_yakutsk = c(129.7331, 62.0355),
  # canada_fort_mcmurray = c(-111.3800, 56.7260),
  # japan_sapporo = c(141.3545, 43.0621),
  china_dandong = c(124.3547, 40.0005)
  # usa_boston = c(-71.0589, 42.3601)
)


nearest_indices <- sapply(target_points, function(loc) {
  dists <- apply(coords, 1, function(row) sqrt(sum((row - loc)^2)))
  which.min(dists)
})



ids <- nearest_indices


sample_idx <- seq(from = 204, to = 1000, by = 4)

setwd("D:/77/Research/temp/snow/")


theta <- readRDS("theta01_bym.Rda")[,sample_idx]
theta_s <- readRDS("theta10_bym.Rda")[,sample_idx]
S <- nrow(y)

theta_01_all <- theta[ids,]
theta_02_all <- theta[S+ids,]
theta_11_all <- theta[2*S+ids,]
theta_12_all <- theta[3*S+ids,]
theta_21_all <- theta[4*S+ids,]
theta_22_all <- theta[5*S+ids,]
theta_a1_all <- theta[6*S+ids,]
theta_a2_all <- theta[7*S+ids,]



theta_01s_all <- theta_s[ids,]
theta_02s_all <- theta_s[S+ids,]
theta_11s_all <- theta_s[2*S+ids,]
theta_12s_all <- theta_s[3*S+ids,]
theta_21s_all <- theta_s[4*S+ids,]
theta_22s_all <- theta_s[5*S+ids,]
theta_a1s_all <- theta_s[6*S+ids,]
theta_a2s_all <- theta_s[7*S+ids,]



beta_0_hat <- theta_01_all + theta_02_all
beta_1_hat <- theta_11_all + theta_12_all
beta_2_hat <- theta_21_all + theta_22_all
alpha_hat <- theta_a1_all + theta_a2_all




beta_0s_hat <- theta_01s_all + theta_02s_all 
beta_1s_hat <- theta_11s_all + theta_12s_all
beta_2s_hat <- theta_21s_all + theta_22s_all
alpha_s_hat <- theta_a1s_all + theta_a2s_all


coefs <- round(
  rbind(
    cbind(
      apply(beta_0_hat, 1, median),
      apply(beta_1_hat, 1, median),
      apply(beta_2_hat, 1, median),
      apply(alpha_hat, 1, median),
      apply(beta_0s_hat, 1, median),
      apply(beta_1s_hat, 1, median),
      apply(beta_2s_hat, 1, median),
      apply(alpha_s_hat, 1, median)
    ),
    cbind(
      apply(beta_0_hat, 1, quantile,0.025 ),
      apply(beta_1_hat, 1, quantile,0.025 ),
      apply(beta_2_hat, 1, quantile,0.025 ),
      apply(alpha_hat, 1, quantile,0.025 ),
      apply(beta_0s_hat, 1, quantile,0.025 ),
      apply(beta_1s_hat, 1, quantile,0.025 ),
      apply(beta_2s_hat, 1, quantile,0.025 ),
      apply(alpha_s_hat, 1, quantile,0.025 )
    ),
    cbind(
      apply(beta_0_hat, 1, quantile,0.975 ),
      apply(beta_1_hat, 1, quantile,0.975 ),
      apply(beta_2_hat, 1, quantile,0.975 ),
      apply(alpha_hat, 1, quantile,0.975 ),
      apply(beta_0s_hat, 1, quantile,0.975 ),
      apply(beta_1s_hat, 1, quantile,0.975 ),
      apply(beta_2s_hat, 1, quantile,0.975 ),
      apply(alpha_s_hat, 1, quantile,0.975 )
    )
  )
,6)

colnames(coefs) <- c(1:8)
knitr::kable(coefs, format = "latex")

# Transition Matrix ----------------------------------------------------------------------------Trend

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
      beta_0_hat[,idx] + beta_1_hat[,idx] * cos(2*pi*time/period) + beta_2_hat[,idx] * sin(2*pi*time/period) + alpha_hat[,idx] * time )
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
locations <- c("Russia - Novosibirsk","China - Dandong")
plots <- list()

for (i in 1:SS) {
  pred_samples <- summed_array[, i, ]
  pred_df <- as.data.frame(t(pred_samples))
  colnames(pred_df) <- paste0("draw_", 1:200)
  pred_df$year <- 1972:(1972+52-1)
  
  summary_df <- pred_df %>%
    pivot_longer(cols = starts_with("draw_"), names_to = "draw", values_to = "value") %>%
    group_by(year) %>%
    summarise(
      mean = mean(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    )
  
  summary_df$truth <- summed_y[i, ]
  
  summary_df_long <- summary_df %>%
    pivot_longer(cols = c(mean, truth), names_to = "Type", values_to = "Value")
  
  p <- ggplot(summary_df_long, aes(x = year)) +
    geom_ribbon(
      data = summary_df,
      aes(x = year, ymin = lower, ymax = upper),
      inherit.aes = FALSE, fill = "blue", alpha = 0.2
    ) +
    geom_line(aes(y = Value, color = Type), linewidth = 1) +
    scale_color_manual(
      values = c("truth" = "red", "mean" = "blue"),
      labels = c("truth" = "Truth", "mean" = "Prediction")
    ) +
    labs(
      title = paste(locations[i]),
      x = "Year",
      y = "Number of Snowy Weeks",
      color = NULL
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (i != SS) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "bottom")
  }
  
  plots[[i]] <- p
}

patchwork::wrap_plots(plots, nrow = 2)






setwd("D:/77/Research/temp/snow/")
saveRDS(weekly_pred_lat_alt, "fb_hk_pred.Rda")

