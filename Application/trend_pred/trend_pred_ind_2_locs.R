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
# LOAD DATA（不删点）
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

# lat
lat <- scale(coords[,2])[,1]

# =========================
# ELEVATION（关键修复）
# =========================
elev_raw <- read.csv("curr_elev.csv")[,4]

nnbs_df <- read.csv("nnbs_elev.csv", sep="\t", row.names = NULL)
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

# sanity
cat("elev NA count:", sum(is.na(elev_all)), "\n")

elev <- scale(elev_all)[,1]

# =========================
# TEMP（full）
# =========================
load("snow_temp_full.Rda")
snow_temp <- sce_temp

temp <- as.matrix(snow_temp[, -c(1,2)])
temp <- scale(temp)

# =========================
# TIME
# =========================
t_full <- 1:TT
t_scaled <- scale(t_full)[,1]

# =====================================================
# TARGET LOCATIONS
# =====================================================
target_points <- list(
  Novosibirsk = c(82.9204, 55.0302),
  Dandong     = c(124.3547, 40.0005),
  Winnipeg    = c(-97.1384, 49.8951),
  Minneapolis = c(-93.2650, 44.9778)
)

nearest_indices <- sapply(target_points, function(loc){
  dists <- apply(coords, 1, function(row) sqrt(sum((row - loc)^2)))
  which.min(dists)
})

ids <- nearest_indices
SS <- length(ids)

# =====================================================
# PYTHON
# =====================================================
np <- import("numpy")

# =====================================================
# LOAD IND（新文件名）
# =====================================================
load_ind <- function(prefix){
  
  eta_list <- list()
  
  for(c in chains){
    
    cat("Loading chain", c, "\n")
    
    d <- np$load(file.path(BASE_DIR, sprintf("%s_chain%d.npz", prefix, c)))
    key <- d$files[[1]]
    
    eta_raw <- d[[key]]
    
    eta <- py_to_r(eta_raw)
    eta <- as.matrix(eta)
    
    eta <- eta[, seq(1, ncol(eta), by=thin), drop=FALSE]
    
    cat("  dim:", dim(eta), "\n")
    
    eta_list[[length(eta_list)+1]] <- eta
  }
  
  eta_list
}

# =====================================================
# LOAD POSTERIOR
# =====================================================
ind01 <- load_ind("p01_ind_all")
ind10 <- load_ind("p10_ind_all")

# =====================================================
# HELPER
# =====================================================
inv_logit <- function(x) 1/(1+exp(-x))

# =====================================================
# IID SIM（原始结构）
compute_iid_sim <- function(eta01_list, eta10_list){
  
  eta01_all <- do.call(cbind, eta01_list)
  eta10_all <- do.call(cbind, eta10_list)
  
  M <- ncol(eta01_all)
  out <- array(NA, c(M, SS, TT-1))
  cat("\n===== START IID SIM =====\n")
  cat("Total samples:", M, "\n")
  
  for(m in 1:M){
    if(m %% 10 == 0){
      cat("\nsample", m, "/", ncol(eta01_all), "\n")
    }
    for(i in 1:SS){
      s <- ids[i]
      y_curr <- y[s,1]
      
      for(t in 1:(TT-1)){
        
        cov4 <- c(1,
                  cos(2*pi*t/period),
                  sin(2*pi*t/period),
                  t_scaled[t])
        
        phi01 <- sum(eta01_all[((0:3)*S + s), m] * cov4)
        phi10 <- sum(eta10_all[((0:3)*S + s), m] * cov4)
        
        p01 <- inv_logit(phi01)
        p10 <- inv_logit(phi10)
        p <- (1 - y_curr) * p01 + y_curr * (1 - p10)
        
        y_next <- p
        
        out[m,i,t] <- y_next
        y_curr <- y_next
      }
    }
  }
  
  cat("===== DONE IID SIM =====\n")
  
  out
}

# =====================================================
# RUN
# =====================================================
iid_pred <- compute_iid_sim(ind01, ind10)

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

iid_year <- agg_year(iid_pred)

# =====================================================
# PLOT
# =====================================================
names(ids) <- names(target_points)

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
      theme(plot.title = element_text(hjust = 0.5))
    
    plots[[i]] <- p
  }
  
  wrap_plots(plots)
}

p1 <- plot_model(iid_year)
p1
