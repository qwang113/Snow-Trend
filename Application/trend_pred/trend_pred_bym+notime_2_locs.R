# rm(list = ls())
# 
# library(sf)
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(reticulate)
# library(patchwork)
# 
# # =====================================================
# # CONFIG
# # =====================================================
# BASE_DIR <- "D:/77/Research/temp/snow"
# setwd(BASE_DIR)
# 
# period <- 52
# thin <- 15
# chains <- 0:9
# 
# no_nbs <- c(
#   57,170,236,269,343,685,946,947,989,
#   1037,1084,1090,1109,1118,1127,1176,1203
# )
# 
# # =====================================================
# # LOAD DATA
# # =====================================================
# library(sf)
# library(igraph)
# 
# snow <- readRDS("snow_cleaned_full.Rda")
# 
# coords_full <- as.matrix(snow[,1:2])
# y_full <- as.matrix(snow[,-c(1,2)])
# 
# # -------------------------
# # 1️⃣ remove isolated
# # -------------------------
# coords_tmp <- coords_full[-no_nbs,]
# y_tmp <- y_full[-no_nbs,]
# 
# # -------------------------
# # 2️⃣ AEQD projection
# # -------------------------
# coords_sf <- st_as_sf(
#   data.frame(lon=coords_tmp[,1], lat=coords_tmp[,2]),
#   coords=c("lon","lat"),
#   crs=4326
# )
# 
# coords_proj <- st_transform(coords_sf, "+proj=aeqd +lat_0=90 +lon_0=-100")
# 
# xy <- st_coordinates(coords_proj) / 1e6
# 
# # -------------------------
# # 3️⃣ adjacency
# # -------------------------
# D <- as.matrix(dist(xy))
# W <- (D <= 0.22) * 1
# diag(W) <- 0
# 
# # -------------------------
# # 4️⃣ connected components
# # -------------------------
# g <- graph_from_adjacency_matrix(W, mode="undirected")
# 
# comp <- components(g)
# sizes <- comp$csize
# order <- order(sizes, decreasing=TRUE)
# 
# keep <- which(comp$membership %in% order[1:2])
# 
# # -------------------------
# # 5️⃣ final data
# # -------------------------
# coords <- coords_tmp[keep,]
# y <- y_tmp[keep,]
# 
# S <- nrow(y)
# TT <- ncol(y)
# 
# cat("Using S =", S, "\n")
# 
# 
# # covariates
# lat <- scale(coords[,2])[,1]
# 
# elev_all <- read.csv("curr_elev.csv")[,4]
# elev_tmp <- elev_all[-no_nbs]
# elev <- scale(elev_tmp[keep])[,1]
# 
# load("snow_temp_full.Rda")
# snow_temp <- sce_temp
# temp_all <- as.matrix(snow_temp[-no_nbs, -c(1,2)])
# temp <- scale(temp_all[keep,])
# 
# t_full <- 1:TT
# t_scaled <- scale(t_full)[,1]
# 
# # =====================================================
# # TARGET LOCATIONS
# # =====================================================
# target_points <- list(
#   Novosibirsk = c(82.9204, 55.0302),
#   Dandong     = c(124.3547, 40.0005),
#   Winnipeg    = c(-97.1384, 49.8951),
#   Minneapolis = c(-93.2650, 44.9778)
# )
# nearest_indices <- sapply(target_points, function(loc){
#   dists <- apply(coords, 1, function(row) sqrt(sum((row - loc)^2)))
#   which.min(dists)
# })
# 
# ids <- nearest_indices
# 
# SS <- length(ids)
# 
# # =====================================================
# # PYTHON LOADER
# # =====================================================
# use_python("D:/anaconda3/envs/CPD/python.exe", required = TRUE)
# np <- import("numpy")
# pickle <- import("pickle")
# builtins <- import_builtins()
# 
# load_bym <- function(prefix, use_cov=FALSE){
#   
#   eta_list <- list()
#   tau_list <- list()
#   
#   for(c in chains){
#     
#     fname <- if(use_cov){
#       sprintf("%s_chain%d+cov.pkl", prefix, c)
#     } else {
#       sprintf("%s_chain%d.pkl", prefix, c)
#     }
#     
#     path <- file.path(BASE_DIR, fname)
#     
#     f <- builtins$open(path, "rb")
#     d <- pickle$load(f)
#     f$close()
#     
#     eta_raw <- d$eta
#     tau_raw <- d$tau
#     
#     if (inherits(eta_raw, "python.builtin.object")) {
#       eta <- py_to_r(eta_raw)
#       tau <- py_to_r(tau_raw)
#     } else {
#       eta <- eta_raw
#       tau <- tau_raw
#     }
#     
#     # thinning
#     eta <- eta[, seq(1, ncol(eta), by=thin)]
#     tau <- tau[, seq(1, ncol(tau), by=thin)]
#     
#     eta_list[[length(eta_list)+1]] <- eta
#     tau_list[[length(tau_list)+1]] <- tau
#   }
#   
#   list(eta=eta_list, tau=tau_list)
# }
# 
# load_ind <- function(prefix){
#   
#   eta_list <- list()
#   
#   for(c in chains){
#     
#     d <- np$load(sprintf("%s_chain%d.npz", prefix, c))
#     key <- d$files[[1]]
#     
#     eta_raw <- d[[key]]
#     
#     if (inherits(eta_raw, "python.builtin.object")) {
#       eta <- py_to_r(eta_raw)
#     } else {
#       eta <- eta_raw
#     }
#     
#     eta <- eta[, seq(1, ncol(eta), by=thin)]
#     
#     eta_list[[length(eta_list)+1]] <- eta
#   }
#   
#   eta_list
# }
# 
# # =====================================================
# # LOAD POSTERIOR
# # =====================================================
# bym01_cov <- load_bym("bym01", TRUE)
# bym10_cov <- load_bym("bym10", TRUE)
# 
# # bym01 <- load_bym("bym01", FALSE)
# # bym10 <- load_bym("bym10", FALSE)
# # 
# # ind01 <- load_ind("ind01")
# # ind10 <- load_ind("ind10")
# 
# # =====================================================
# # HELPER
# # =====================================================
# inv_logit <- function(x) 1/(1+exp(-x))
# 
# one_step <- function(phi01, phi10, y_prev){
#   
#   p <- ifelse(y_prev==0, inv_logit(phi01), 1-inv_logit(phi10))
#   
#   # rbinom(1, 1, p)
#   return(p)
# }
# 
# # =====================================================
# # IID
# # =====================================================
# compute_iid <- function(eta01_list, eta10_list){
#   
#   eta01_all <- do.call(cbind, eta01_list)
#   eta10_all <- do.call(cbind, eta10_list)
#   
#   M <- ncol(eta01_all)
#   out <- array(NA, c(M, SS, TT-1))
#   
#   for(m in 1:M){
#     for(i in 1:SS){
#       s <- ids[i]
#       
#       for(t in 1:(TT-1)){
#         
#         cov4 <- c(1,
#                   cos(2*pi*t/period),
#                   sin(2*pi*t/period),
#                   t_scaled[t])
#         
#         phi01 <- sum(eta01_all[((0:3)*S + s), m] * cov4)
#         phi10 <- sum(eta10_all[((0:3)*S + s), m] * cov4)
#         
#         out[m,i,t] <- one_step(phi01, phi10, y[s,t])
#       }
#     }
#   }
#   
#   out
# }
# 
# # =====================================================
# # BYM (no cov)
# # =====================================================
# compute_bym <- function(eta01, tau01, eta10, tau10){
#   
#   eta01_all <- do.call(cbind, eta01$eta)
#   tau01_all <- do.call(cbind, eta01$tau)
#   
#   eta10_all <- do.call(cbind, eta10$eta)
#   tau10_all <- do.call(cbind, eta10$tau)
#   
#   M <- ncol(eta01_all)
#   out <- array(NA, c(M, SS, TT-1))
#   
#   for(m in 1:M){
#     for(i in 1:SS){
#       
#       s <- ids[i]
#       
#       for(t in 1:(TT-1)){
#         
#         w <- (t-1) %% 52 + 1
#         cov4 <- c(1, cos(2*pi*t/period), sin(2*pi*t/period), t_scaled[t])
#         
#         phi01 <- 0
#         phi10 <- 0
#         
#         for(k in 1:8){
#           
#           phi01 <- phi01 +
#             eta01_all[(k-1)*S+s, m] *
#             tau01_all[(k-1)*52+w, m] *
#             cov4[(k-1)%/%2+1]
#           
#           phi10 <- phi10 +
#             eta10_all[(k-1)*S+s, m] *
#             tau10_all[(k-1)*52+w, m] *
#             cov4[(k-1)%/%2+1]
#         }
#         
#         out[m,i,t] <- one_step(phi01, phi10, y[s,t])
#       }
#     }
#   }
#   
#   out
# }
# 
# # =====================================================
# # BYM + COV
# # =====================================================
# compute_bym_cov <- function(eta01, tau01, eta10, tau10){
#   
#   eta01_all <- do.call(cbind, eta01$eta)
#   tau01_all <- do.call(cbind, eta01$tau)
#   
#   eta10_all <- do.call(cbind, eta10$eta)
#   tau10_all <- do.call(cbind, eta10$tau)
#   
#   M <- ncol(eta01_all)
#   out <- array(NA, c(M, SS, TT-1))
#   
#   for(m in 1:M){
#     for(i in 1:SS){
#       
#       s <- ids[i]
#       
#       for(t in 1:(TT-1)){
#         
#         w <- (t-1) %% 52 + 1
#         cov4 <- c(1, cos(2*pi*t/period), sin(2*pi*t/period), t_scaled[t])
#         
#         phi01 <- 0
#         phi10 <- 0
#         
#         for(k in 1:8){
#           
#           phi01 <- phi01 +
#             eta01_all[(k-1)*S+s, m] *
#             tau01_all[(k-1)*52+w, m] *
#             cov4[(k-1)%/%2+1]
#           
#           phi10 <- phi10 +
#             eta10_all[(k-1)*S+s, m] *
#             tau10_all[(k-1)*52+w, m] *
#             cov4[(k-1)%/%2+1]
#         }
#         
#         gamma01 <- eta01_all[(8*S+1):(8*S+3), m]
#         gamma10 <- eta10_all[(8*S+1):(8*S+3), m]
#         
#         fac <- c(
#           t_scaled[t]*lat[s],
#           t_scaled[t]*elev[s],
#           t_scaled[t]*temp[s,t]
#         )
#         
#         phi01 <- phi01 + sum(gamma01 * fac)
#         phi10 <- phi10 + sum(gamma10 * fac)
#         
#         out[m,i,t] <- one_step(phi01, phi10, y[s,t])
#       }
#     }
#   }
#   
#   out
# }
# 
# 
# compute_iid_sim <- function(eta01_list, eta10_list){
#   
#   eta01_all <- do.call(cbind, eta01_list)
#   eta10_all <- do.call(cbind, eta10_list)
#   
#   M <- ncol(eta01_all)
#   out <- array(NA, c(M, SS, TT-1))
#   
#   for(m in 1:M){
#     for(i in 1:SS){
#       
#       s <- ids[i]
#       y_curr <- y[s,1]
#       
#       for(t in 1:(TT-1)){
#         
#         cov4 <- c(1,
#                   cos(2*pi*t/period),
#                   sin(2*pi*t/period),
#                   t_scaled[t])
#         
#         phi01 <- sum(eta01_all[((0:3)*S + s), m] * cov4)
#         phi10 <- sum(eta10_all[((0:3)*S + s), m] * cov4)
#         
#         p <- ifelse(y_curr==0,
#                     inv_logit(phi01),
#                     1-inv_logit(phi10))
#         
#         y_next <- rbinom(1,1,p)
#         
#         out[m,i,t] <- y_next
#         y_curr <- y_next
#       }
#     }
#   }
#   
#   out
# }
# 
# compute_bym_sim <- function(eta01, tau01, eta10, tau10){
#   
#   eta01_all <- do.call(cbind, eta01$eta)
#   tau01_all <- do.call(cbind, eta01$tau)
#   
#   eta10_all <- do.call(cbind, eta10$eta)
#   tau10_all <- do.call(cbind, eta10$tau)
#   
#   M <- ncol(eta01_all)
#   out <- array(NA, c(M, SS, TT-1))
#   
#   for(m in 1:M){
#     for(i in 1:SS){
#       
#       s <- ids[i]
#       y_curr <- y[s,1]
#       
#       for(t in 1:(TT-1)){
#         
#         w <- (t-1) %% 52 + 1
#         cov4 <- c(1, cos(2*pi*t/period), sin(2*pi*t/period), t_scaled[t])
#         
#         phi01 <- 0
#         phi10 <- 0
#         
#         for(k in 1:8){
#           
#           phi01 <- phi01 +
#             eta01_all[(k-1)*S+s, m] *
#             tau01_all[(k-1)*52+w, m] *
#             cov4[(k-1)%/%2+1]
#           
#           phi10 <- phi10 +
#             eta10_all[(k-1)*S+s, m] *
#             tau10_all[(k-1)*52+w, m] *
#             cov4[(k-1)%/%2+1]
#         }
#         
#         p <- ifelse(y_curr==0,
#                     inv_logit(phi01),
#                     1-inv_logit(phi10))
#         
#         y_next <- rbinom(1,1,p)
#         
#         out[m,i,t] <- y_next
#         y_curr <- y_next
#       }
#     }
#   }
#   
#   out
# }
# 
# 
# compute_bym_cov_sim <- function(eta01, tau01, eta10, tau10){
#   
#   eta01_all <- do.call(cbind, eta01$eta)
#   tau01_all <- do.call(cbind, eta01$tau)
#   
#   eta10_all <- do.call(cbind, eta10$eta)
#   tau10_all <- do.call(cbind, eta10$tau)
#   
#   M <- ncol(eta01_all)
#   out <- array(NA, c(M, SS, TT-1))
#   
#   for(m in 1:M){
#     for(i in 1:SS){
#       
#       s <- ids[i]
#       y_curr <- y[s,1]
#       
#       for(t in 1:(TT-1)){
#         
#         w <- (t-1) %% 52 + 1
#         cov4 <- c(1, cos(2*pi*t/period), sin(2*pi*t/period), t_scaled[t])
#         
#         phi01 <- 0
#         phi10 <- 0
#         
#         for(k in 1:8){
#           
#           phi01 <- phi01 +
#             eta01_all[(k-1)*S+s, m] *
#             tau01_all[(k-1)*52+w, m] *
#             cov4[(k-1)%/%2+1]
#           
#           phi10 <- phi10 +
#             eta10_all[(k-1)*S+s, m] *
#             tau10_all[(k-1)*52+w, m] *
#             cov4[(k-1)%/%2+1]
#         }
#         
#         gamma01 <- eta01_all[(8*S+1):(8*S+3), m]
#         gamma10 <- eta10_all[(8*S+1):(8*S+3), m]
#         
#         fac <- c(
#           t_scaled[t]*lat[s],
#           t_scaled[t]*elev[s],
#           t_scaled[t]*temp[s,t]
#         )
#         
#         phi01 <- phi01 + sum(gamma01 * fac)
#         phi10 <- phi10 + sum(gamma10 * fac)
#         
#         p <- ifelse(y_curr==0,
#                     inv_logit(phi01),
#                     1-inv_logit(phi10))
#         
#         y_next <- rbinom(1,1,p)
#         
#         out[m,i,t] <- y_next
#         y_curr <- y_next
#       }
#     }
#   }
#   
#   out
# }
# 
# 
# 
# 
# # =====================================================
# # RUN
# # =====================================================
# # iid_pred <- compute_iid_sim(ind01, ind10)
# # bym_pred <- compute_bym_sim(bym01, bym01, bym10, bym10)
# bym_cov_pred <- compute_bym_cov(bym01_cov, bym01_cov, bym10_cov, bym10_cov)
# 
# # =====================================================
# # YEAR AGGREGATION
# # =====================================================
# agg_year <- function(arr){
#   
#   num_year <- dim(arr)[3] / 52
#   out <- array(NA, c(dim(arr)[1], dim(arr)[2], num_year))
#   
#   for(i in 1:num_year){
#     idx <- ((i-1)*52+1):(i*52)
#     out[,,i] <- apply(arr[,,idx], c(1,2), sum)
#   }
#   
#   out
# }
# 
# # iid_year <- agg_year(iid_pred)
# # bym_year <- agg_year(bym_pred)
# bym_cov_year <- agg_year(bym_cov_pred)
# 
# # =====================================================
# # PLOT
# # =====================================================
# names(ids) <- names(target_points)
# plot_model <- function(arr, name){
#   
#   plots <- list()
#   
#   for(i in 1:SS){
#     
#     pred <- arr[,i,]
#     df <- as.data.frame(t(pred))
#     df$year <- 1972:(1972+ncol(pred)-1)
#     
#     summary <- df |>
#       pivot_longer(-year) |>
#       group_by(year) |>
#       summarise(mean=mean(value),
#                 low=quantile(value,0.025),
#                 up=quantile(value,0.975),
#                 .groups="drop")
#     
#     truth <- sapply(1:nrow(summary), function(j){
#       idx <- ((j-1)*52+1):(j*52)
#       sum(y[ids[i], idx])
#     })
#     
#     summary$truth <- truth
#     
#     p <- ggplot(summary, aes(year)) +
#       geom_ribbon(aes(ymin=low,ymax=up), fill="blue", alpha=0.2) +
#       geom_line(aes(y=mean), color="blue") +
#       geom_line(aes(y=truth), color="red") +
#       ggtitle(paste(names(ids)[i])) +
#       theme_minimal() +
#       theme(plot.title = element_text(hjust = 0.5))
#     
#     plots[[i]] <- p
#   }
#   
#   wrap_plots(plots)
# }
# 
# # p1 <- plot_model(iid_year, "IID")
# # p2 <- plot_model(bym_year, "BYM")
# p3 <- plot_model(bym_cov_year, "BYM+Cov")
# 
# # p1 / p2 / p3
# p3
# # 
# # compute_metrics <- function(arr){
# #   
# #   # arr: M × SS × T_year
# #   
# #   M <- dim(arr)[1]
# #   SS <- dim(arr)[2]
# #   TT <- dim(arr)[3]
# #   
# #   results <- list()
# #   
# #   for(i in 1:SS){
# #     
# #     pred <- arr[,i,]   # M × TT
# #     
# #     # posterior summary
# #     mean_pred <- apply(pred, 2, mean)
# #     low <- apply(pred, 2, quantile, 0.025)
# #     up  <- apply(pred, 2, quantile, 0.975)
# #     
# #     # truth
# #     truth <- sapply(1:TT, function(j){
# #       idx <- ((j-1)*52+1):(j*52)
# #       sum(y[ids[i], idx])
# #     })
# #     
# #     # ---------------------------
# #     # MSE
# #     # ---------------------------
# #     mse <- mean((mean_pred - truth)^2)
# #     
# #     # ---------------------------
# #     # Coverage
# #     # ---------------------------
# #     covered <- (truth >= low) & (truth <= up)
# #     coverage <- mean(covered)
# #     
# #     results[[i]] <- data.frame(
# #       location = names(ids)[i],
# #       MSE = mse,
# #       Coverage = coverage
# #     )
# #   }
# #   
# #   bind_rows(results)
# # }
# # 
# # iid_metrics     <- compute_metrics(iid_year)
# # bym_metrics     <- compute_metrics(bym_year)
# # bym_cov_metrics <- compute_metrics(bym_cov_year)
# # 
# # iid_metrics$model <- "IID"
# # bym_metrics$model <- "BYM"
# # bym_cov_metrics$model <- "BYM+Cov"
# # 
# # all_metrics <- bind_rows(
# #   iid_metrics,
# #   bym_metrics,
# #   bym_cov_metrics
# # )
# # 
# # print(all_metrics)
# # ids



rm(list = ls())

library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reticulate)
library(patchwork)

# =====================================================
# CONFIG
# =====================================================
BASE_DIR <- "D:/77/Research/temp/snow"
setwd(BASE_DIR)

period <- 52
thin <- 15
chains <- 0:9

# =====================================================
# LOAD DATA (FULL DATA, NO DELETION)
# =====================================================
snow <- readRDS("snow_cleaned_full.Rda")

coords <- as.matrix(snow[, 1:2])
y <- as.matrix(snow[, -c(1, 2)])

S <- nrow(y)
TT <- ncol(y)

cat("Using FULL S =", S, "\n")

t_full <- 1:TT
t_scaled <- scale(t_full)[, 1]

# =====================================================
# TARGET LOCATIONS
# =====================================================
target_points <- list(
  Novosibirsk = c(82.9204, 55.0302),
  Dandong     = c(124.3547, 40.0005),
  Winnipeg    = c(-97.1384, 49.8951),
  Minneapolis = c(-93.2650, 44.9778)
)

nearest_indices <- sapply(target_points, function(loc) {
  dists <- apply(coords, 1, function(row) sqrt(sum((row - loc)^2)))
  which.min(dists)
})

ids <- nearest_indices
SS <- length(ids)

# =====================================================
# PYTHON LOADER
# =====================================================
use_condaenv("CPD", required = TRUE)
np <- import("numpy")
pickle <- import("pickle")
builtins <- import_builtins()

load_bym <- function(prefix) {
  
  eta_list <- list()
  tau_list <- list()
  
  for (c in chains) {
    
    fname <- sprintf("%s_chain%d.pkl", prefix, c)
    path <- file.path(BASE_DIR, fname)
    
    f <- builtins$open(path, "rb")
    d <- pickle$load(f)
    f$close()
    
    eta_raw <- d$eta
    tau_raw <- d$tau
    
    eta <- py_to_r(eta_raw)
    tau <- py_to_r(tau_raw)
    
    if (is.null(dim(eta))) {
      eta_shape <- as.integer(unlist(py_to_r(eta_raw$shape)))
      eta <- matrix(eta, nrow = eta_shape[1], ncol = eta_shape[2])
    } else {
      eta <- as.matrix(eta)
    }
    
    if (is.null(dim(tau))) {
      tau_shape <- as.integer(unlist(py_to_r(tau_raw$shape)))
      tau <- matrix(tau, nrow = tau_shape[1], ncol = tau_shape[2])
    } else {
      tau <- as.matrix(tau)
    }
    
    eta <- eta[, seq(1, ncol(eta), by = thin), drop = FALSE]
    tau <- tau[, seq(1, ncol(tau), by = thin), drop = FALSE]
    
    eta_list[[length(eta_list) + 1]] <- eta
    tau_list[[length(tau_list) + 1]] <- tau
  }
  
  list(eta = eta_list, tau = tau_list)
}

load_ind <- function(prefix) {
  
  eta_list <- list()
  
  for (c in chains) {
    
    d <- np$load(file.path(BASE_DIR, sprintf("%s_chain%d.npz", prefix, c)), allow_pickle = TRUE)
    key <- d$files[[1]]
    eta_raw <- d[[key]]
    
    eta <- py_to_r(eta_raw)
    
    if (is.null(dim(eta))) {
      shape <- as.integer(unlist(py_to_r(eta_raw$shape)))
      eta <- matrix(eta, nrow = shape[1], ncol = shape[2])
    } else {
      eta <- as.matrix(eta)
    }
    
    idx <- seq(1, ncol(eta), by = thin)
    eta <- eta[, idx, drop = FALSE]
    
    eta_list[[length(eta_list) + 1]] <- eta
  }
  
  eta_list
}

# =====================================================
# LOAD POSTERIOR (ONLY IID + BYM)
# =====================================================
ind01 <- load_ind("p01_ind_all")
ind10 <- load_ind("p10_ind_all")

bym01 <- load_bym("p01_weekly")
bym10 <- load_bym("p10_weekly")

# =====================================================
# HELPER
# =====================================================
inv_logit <- function(x) 1 / (1 + exp(-x))

one_step_prob <- function(phi01, phi10, y_prev) {
  ifelse(y_prev == 0, inv_logit(phi01), 1 - inv_logit(phi10))
}

# =====================================================
# IID
# =====================================================
compute_iid_sim <- function(eta01_list, eta10_list) {
  
  eta01_all <- do.call(cbind, eta01_list)
  eta10_all <- do.call(cbind, eta10_list)
  
  M <- ncol(eta01_all)
  out <- array(NA_real_, c(M, SS, TT - 1))
  
  for (m in 1:M) {
    for (i in 1:SS) {
      
      s <- ids[i]
      y_curr <- y[s, 1]
      
      for (t in 1:(TT - 1)) {
        
        cov4 <- c(
          1,
          cos(2 * pi * t / period),
          sin(2 * pi * t / period),
          t_scaled[t]
        )
        
        phi01 <- sum(eta01_all[((0:3) * S + s), m] * cov4)
        phi10 <- sum(eta10_all[((0:3) * S + s), m] * cov4)
        
        p <- one_step_prob(phi01, phi10, y_curr)
        y_next <- rbinom(1, 1, p)
        
        out[m, i, t] <- y_next
        y_curr <- y_next
      }
    }
  }
  
  out
}

# =====================================================
# BYM (NO COV)
# =====================================================
compute_bym_sim <- function(eta01_obj, eta10_obj) {
  
  eta01_all <- do.call(cbind, eta01_obj$eta)
  tau01_all <- do.call(cbind, eta01_obj$tau)
  
  eta10_all <- do.call(cbind, eta10_obj$eta)
  tau10_all <- do.call(cbind, eta10_obj$tau)
  
  M <- ncol(eta01_all)
  out <- array(NA_real_, c(M, SS, TT - 1))
  
  for (m in 1:M) {
    for (i in 1:SS) {
      
      s <- ids[i]
      y_curr <- y[s, 1]
      
      for (t in 1:(TT - 1)) {
        
        w <- (t - 1) %% 52 + 1
        cov4 <- c(
          1,
          cos(2 * pi * t / period),
          sin(2 * pi * t / period),
          t_scaled[t]
        )
        
        phi01 <- 0
        phi10 <- 0
        
        for (k in 1:8) {
          
          phi01 <- phi01 +
            eta01_all[(k - 1) * S + s, m] *
            tau01_all[(k - 1) * 52 + w, m] *
            cov4[(k - 1) %/% 2 + 1]
          
          phi10 <- phi10 +
            eta10_all[(k - 1) * S + s, m] *
            tau10_all[(k - 1) * 52 + w, m] *
            cov4[(k - 1) %/% 2 + 1]
        }
        
        p <- one_step_prob(phi01, phi10, y_curr)
        y_next <- rbinom(1, 1, p)
        
        out[m, i, t] <- y_next
        y_curr <- y_next
      }
    }
  }
  
  out
}

# =====================================================
# RUN
# =====================================================
iid_pred <- compute_iid_sim(ind01, ind10)
bym_pred <- compute_bym_sim(bym01, bym10)

# =====================================================
# YEAR AGGREGATION
# =====================================================
agg_year <- function(arr) {
  
  num_year <- dim(arr)[3] / 52
  out <- array(NA_real_, c(dim(arr)[1], dim(arr)[2], num_year))
  
  for (i in 1:num_year) {
    idx <- ((i - 1) * 52 + 1):(i * 52)
    out[, , i] <- apply(arr[, , idx, drop = FALSE], c(1, 2), sum)
  }
  
  out
}

iid_year <- agg_year(iid_pred)
bym_year <- agg_year(bym_pred)

# =====================================================
# PLOT
# =====================================================
names(ids) <- names(target_points)

plot_model <- function(arr, name) {
  
  plots <- list()
  
  for (i in 1:SS) {
    
    pred <- arr[, i, ]
    df <- as.data.frame(t(pred))
    df$year <- 1972:(1972 + ncol(pred) - 1)
    
    summary <- df |>
      pivot_longer(-year) |>
      group_by(year) |>
      summarise(
        mean = mean(value),
        low  = quantile(value, 0.025),
        up   = quantile(value, 0.975),
        .groups = "drop"
      )
    
    truth <- sapply(1:nrow(summary), function(j) {
      idx <- ((j - 1) * 52 + 1):(j * 52)
      sum(y[ids[i], idx])
    })
    
    summary$truth <- truth
    
    p <- ggplot(summary, aes(year)) +
      geom_ribbon(aes(ymin = low, ymax = up), fill = "blue", alpha = 0.2) +
      geom_line(aes(y = mean), color = "blue") +
      geom_line(aes(y = truth), color = "red") +
      ggtitle(paste(name, "-", names(ids)[i])) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    plots[[i]] <- p
  }
  
  wrap_plots(plots)
}

p1 <- plot_model(iid_year, "IID")
p2 <- plot_model(bym_year, "BYM")

p1 / p2

# =====================================================
# OPTIONAL METRICS
# =====================================================
compute_metrics <- function(arr) {
  
  M <- dim(arr)[1]
  SS <- dim(arr)[2]
  TT_year <- dim(arr)[3]
  
  results <- list()
  
  for (i in 1:SS) {
    
    pred <- arr[, i, ]
    
    mean_pred <- apply(pred, 2, mean)
    low <- apply(pred, 2, quantile, 0.025)
    up  <- apply(pred, 2, quantile, 0.975)
    
    truth <- sapply(1:TT_year, function(j) {
      idx <- ((j - 1) * 52 + 1):(j * 52)
      sum(y[ids[i], idx])
    })
    
    mse <- mean((mean_pred - truth)^2)
    coverage <- mean((truth >= low) & (truth <= up))
    
    results[[i]] <- data.frame(
      location = names(ids)[i],
      MSE = mse,
      Coverage = coverage
    )
  }
  
  bind_rows(results)
}

iid_metrics <- compute_metrics(iid_year)
bym_metrics <- compute_metrics(bym_year)

iid_metrics$model <- "IID"
bym_metrics$model <- "BYM"

all_metrics <- bind_rows(iid_metrics, bym_metrics)
print(all_metrics)

ids
