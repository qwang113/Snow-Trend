rm(list = ls())

library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(reticulate)

# =====================================================
# CONFIG
# =====================================================
BASE_DIR <- "D:/77/Research/temp/snow"
setwd(BASE_DIR)

period <- 52
thin <- 15
chains <- 0:9

# =====================================================
# LOAD DATA
# =====================================================
snow <- readRDS("snow_cleaned_full.Rda")

coords <- as.matrix(snow[,1:2])
y <- as.matrix(snow[,-c(1,2)])

S <- nrow(y)
TT <- ncol(y)

cat("Using S =", S, "\n")

# =====================================================
# COVARIATES
# =====================================================

lat <- scale(coords[,2])[,1]

# ===== NEW: longitude split =====
lon_raw <- coords[,1]

region <- rep(0, S)
region[lon_raw >= -30] <- 1

lon_na <- rep(0, S)
lon_euas <- rep(0, S)

mask_na <- (region == 0)
mask_euas <- (region == 1)

lon_na[mask_na] <- scale(lon_raw[mask_na])[,1]
lon_euas[mask_euas] <- scale(lon_raw[mask_euas])[,1]

# elevation
elev_raw <- read.csv("curr_elev.csv")[,4]
nnbs_df <- read.csv("nnbs_elev.csv", sep="\t", row.names=NULL)
nnbs_elev <- nnbs_df[,3]

no_nbs <- c(
  57,170,236,269,343,685,946,947,989,
  1037,1084,1090,1109,1118,1127,1176,1203
)

mask <- rep(TRUE, S)
mask[no_nbs] <- FALSE

elev_all <- numeric(S)
elev_all[mask] <- elev_raw
elev_all[no_nbs] <- nnbs_elev

elev <- scale(elev_all)[,1]

# temperature
load("snow_temp_full.Rda")
snow_temp <- sce_temp
temp <- as.matrix(snow_temp[, -c(1,2)])
temp <- scale(temp)

# time
t_scaled <- scale(1:TT)[,1]

inv_logit <- function(x) 1/(1+exp(-x))

# =====================================================
# TARGET LOCATIONS
# =====================================================
target_points <- list(
  Novosibirsk = c(82.9204, 55.0302),
  Dandong     = c(124.3547, 40.0005),
  Winnipeg    = c(-97.1384, 49.8951),
  Minneapolis = c(-93.2650, 44.9778),
  Hokkaido    = c(141.3545, 43.0618),
  Quebec      = c(-71.2080, 46.8139),
  Chamonix    = c(6.8694, 45.9237),
  Alaska      = c(-147.7164, 64.8378)
)

nearest_indices <- sapply(target_points, function(loc){
  dists <- apply(coords, 1, function(row) sqrt(sum((row - loc)^2)))
  which.min(dists)
})

ids <- nearest_indices
SS <- length(ids)
names(ids) <- names(target_points)

# =====================================================
# NUMPY（npz）
# =====================================================
np <- import("numpy")

# =====================================================
# LOAD NPZ
# =====================================================

load_bym_cov <- function(prefix){
  
  eta_list <- list()
  tau_list <- list()
  
  for(c in chains){
    
    cat("Loading chain", c, "\n")
    
    eta <- as.matrix(readRDS(file.path(BASE_DIR,
                                       sprintf("%s_eta_chain%d.rds", prefix, c))))
    
    tau <- as.matrix(readRDS(file.path(BASE_DIR,
                                       sprintf("%s_tau_chain%d.rds", prefix, c))))
    
    storage.mode(eta) <- "numeric"
    storage.mode(tau) <- "numeric"

    eta <- eta[, seq(1, ncol(eta), by = thin), drop = FALSE]
    tau <- tau[, seq(1, ncol(tau), by = thin), drop = FALSE]
    
    cat("dim eta:", dim(eta), "\n")
    cat("dim tau:", dim(tau), "\n")
    
    eta_list[[length(eta_list)+1]] <- eta
    tau_list[[length(tau_list)+1]] <- tau
  }
  
  list(eta=eta_list, tau=tau_list)
}

# =====================================================
# LOAD POSTERIOR（已是longitude模型）
# =====================================================
bym01 <- load_bym_cov("p01_weekly_cov+lon")
bym10 <- load_bym_cov("p10_weekly_cov+lon")

# =====================================================
# SIMULATION（核心）
# =====================================================
compute_bym_cov <- function(eta01_list, tau01_list,
                            eta10_list, tau10_list){
  
  eta01_all <- do.call(cbind, eta01_list)
  tau01_all <- do.call(cbind, tau01_list)
  
  eta10_all <- do.call(cbind, eta10_list)
  tau10_all <- do.call(cbind, tau10_list)
  
  M <- ncol(eta01_all)
  
  out <- array(NA, c(M, SS, TT-1))
  
  pb <- txtProgressBar(min = 0, max = M, style = 3)
  
  for(m in 1:M){
    
    setTxtProgressBar(pb, m)
    
    for(i in 1:SS){
      
      s <- ids[i]
      y_curr <- y[s,1]
      
      for(t in 1:(TT-1)){
        
        w <- ((t-1) %% period) + 1
        
        cov4 <- c(
          1,
          cos(2*pi*t/period),
          sin(2*pi*t/period),
          t_scaled[t]
        )
        
        phi01 <- 0
        phi10 <- 0
        
        # BYM
        for(k in 1:8){
          phi01 <- phi01 +
            eta01_all[(k-1)*S + s, m] *
            tau01_all[(k-1)*period + w, m] *
            cov4[(k-1)%/%2 + 1]
          
          phi10 <- phi10 +
            eta10_all[(k-1)*S + s, m] *
            tau10_all[(k-1)*period + w, m] *
            cov4[(k-1)%/%2 + 1]
        }
        
        # ===== NEW γ=5 =====
        phi01 <- phi01 +
          t_scaled[t] * (
            eta01_all[8*S + 1, m] * lon_na[s] +
              eta01_all[8*S + 2, m] * lon_euas[s] +
              eta01_all[8*S + 3, m] * lat[s] +
              eta01_all[8*S + 4, m] * elev[s] +
              eta01_all[8*S + 5, m] * temp[s,t]
          )
        
        phi10 <- phi10 +
          t_scaled[t] * (
            eta10_all[8*S + 1, m] * lon_na[s] +
              eta10_all[8*S + 2, m] * lon_euas[s] +
              eta10_all[8*S + 3, m] * lat[s] +
              eta10_all[8*S + 4, m] * elev[s] +
              eta10_all[8*S + 5, m] * temp[s,t]
          )
        
        p01 <- inv_logit(phi01)
        p10 <- inv_logit(phi10)
        
        p <- (1 - y_curr) * p01 + y_curr * (1 - p10)
        
        out[m,i,t] <- p
        y_curr <- p
      }
    }
  }
  
  close(pb)
  out
}

# =====================================================
# RUN
# =====================================================
bym_cov_pred <- compute_bym_cov(
  bym01$eta, bym01$tau,
  bym10$eta, bym10$tau
)

# =====================================================
# YEAR AGG
# =====================================================
agg_year <- function(arr){
  
  num_year <- dim(arr)[3] / 52
  out <- array(NA, c(dim(arr)[1], dim(arr)[2], num_year))
  
  for(i in 1:num_year){
    idx <- ((i-1)*52+1):(i*52)
    out[,,i] <- apply(arr[,,idx], c(1,2), sum)
  }
  
  out
}

bym_cov_year <- agg_year(bym_cov_pred)

# =====================================================
# PLOT
# =====================================================
plot_model <- function(arr){
  
  plots <- list()
  
  for(i in 1:SS){
    
    pred <- arr[,i,]
    
    df <- as.data.frame(t(pred))
    df$year <- 1972:(1972+ncol(pred)-1)
    
    summary <- df |>
      pivot_longer(-year) |>
      group_by(year) |>
      summarise(mean=mean(value),
                low=quantile(value,0.025),
                up=quantile(value,0.975),
                .groups="drop")
    
    truth <- sapply(1:nrow(summary), function(j){
      idx <- ((j-1)*52+1):(j*52)
      sum(y[ids[i], idx])
    })
    
    summary$truth <- truth
    
    p <- ggplot(summary, aes(year)) +
      geom_ribbon(aes(ymin=low,ymax=up), fill="blue", alpha=0.2) +
      geom_line(aes(y=mean), color="blue") +
      geom_line(aes(y=truth), color="red") +
      ggtitle(names(ids)[i]) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    plots[[i]] <- p
  }
  
  wrap_plots(plots)
}

p_cov <- plot_model(bym_cov_year)
p_cov

