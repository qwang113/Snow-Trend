rm(list = ls())

library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
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
# LOAD DATA（FULL S）
# =====================================================
snow <- readRDS("snow_cleaned_full.Rda")

coords <- as.matrix(snow[,1:2])
y <- as.matrix(snow[,-c(1,2)])

S <- nrow(y)
TT <- ncol(y)

cat("Using S =", S, "\n")

# =====================================================
# TIME BASIS（和 IID 一致）
# =====================================================
t_scaled <- scale(1:TT)[,1]

inv_logit <- function(x) 1/(1+exp(-x))

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
names(ids) <- names(target_points)

# =====================================================
# LOAD BYM（.rds version）
# =====================================================
library(reticulate)
py_config()
pickle <- import("pickle")
py <- import_builtins()

load_bym <- function(prefix){
  
  eta_list <- list()
  tau_list <- list()
  
  for(c in chains){
    
    cat("Loading chain", c, "\n")
    
    f <- py$open(file.path(BASE_DIR, sprintf("%s_chain%d.pkl", prefix, c)), "rb")
    d <- pickle$load(f)
    f$close()
    
    eta <- py_to_r(d$eta)
    tau <- py_to_r(d$tau)
    
    eta <- as.matrix(eta)
    tau <- as.matrix(tau)
    
    eta <- eta[, seq(1, ncol(eta), by=thin), drop=FALSE]
    tau <- tau[, seq(1, ncol(tau), by=thin), drop=FALSE]
    
    cat("  eta dim:", dim(eta), "\n")
    cat("  tau dim:", dim(tau), "\n")
    
    eta_list[[length(eta_list)+1]] <- eta
    tau_list[[length(tau_list)+1]] <- tau
  }
  
  list(eta=eta_list, tau=tau_list)
}

# =====================================================
# LOAD POSTERIOR
# =====================================================
bym01 <- load_bym("p01_weekly")
bym10 <- load_bym("p10_weekly")

# =====================================================
# BYM SIM（posterior mean version）
# =====================================================
compute_bym_sim <- function(eta01_list, tau01_list,
                            eta10_list, tau10_list){
  
  eta01_all <- do.call(cbind, eta01_list)
  tau01_all <- do.call(cbind, tau01_list)
  
  eta10_all <- do.call(cbind, eta10_list)
  tau10_all <- do.call(cbind, tau10_list)
  
  M <- ncol(eta01_all)
  
  out <- array(NA, c(M, SS, TT-1))
  
  cat("\n===== START BYM SIM =====\n")
  cat("Total samples:", M, "\n")
  
  for(m in 1:M){
    
    if(m %% 10 == 0){
      cat("\nsample", m, "/", M, "\n")
    }
    
    for(i in 1:SS){
      
      s <- ids[i]
      y_curr <- y[s,1]
      
      for(t in 1:(TT-1)){
        
        w <- ((t-1) %% period) + 1
        
        cov4 <- c(1,
                  cos(2*pi*t/period),
                  sin(2*pi*t/period),
                  t_scaled[t])
        
        phi01 <- 0
        phi10 <- 0
        
        # 8 components
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
        
        p01 <- inv_logit(phi01)
        p10 <- inv_logit(phi10)
        
        # 🔥 expectation recursion（关键）
        p <- (1 - y_curr) * p01 + y_curr * (1 - p10)
        
        y_next <- p
        
        out[m,i,t] <- y_next
        y_curr <- y_next
      }
    }
  }
  
  cat("===== DONE BYM SIM =====\n")
  
  out
}

# =====================================================
# RUN
# =====================================================
bym_pred <- compute_bym_sim(
  bym01$eta, bym01$tau,
  bym10$eta, bym10$tau
)

# =====================================================
# YEAR AGG（snow weeks）
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

bym_year <- agg_year(bym_pred)

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
      theme(plot.title = element_text(hjust = 0.5))
    
    plots[[i]] <- p
  }
  
  wrap_plots(plots)
}

p_bym <- plot_model(bym_year)
p_bym
