rm(list = ls())

library(reticulate)
library(ggplot2)
library(patchwork)
library(tidyr)
# Paths and settings
# Set this path to the local data and results directory.
BASE_DIR <- "path/to/snow/data-and-results"

period <- 52
thin2 <- 15

locations <- c(70, 1400)
weeks <- c(20, 35)
year <- 20

times <- (year * 52 + weeks - 1)
chains <- 0:9
# Data
snow <- readRDS("snow_cleaned_full.Rda")

coords <- as.matrix(snow[,1:2])
y <- as.matrix(snow[,-c(1,2)])

S <- nrow(y)
TT <- ncol(y)

cat("Using FULL S =", S, "TT =", TT, "\n")
# COVARIATES
lat <- scale(coords[,2])[,1]

# longitude split
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
no_nbs <- c(
  57,170,236,269,343,685,946,947,989,
  1037,1084,1090,1109,1118,1127,1176,1203
)

elev_raw <- read.csv(file.path(BASE_DIR,"curr_elev.csv"))[,4]
nnbs_elev <- read.csv(file.path(BASE_DIR,"nnbs_elev.csv"), sep="\t", header = TRUE, row.names = NULL)[,4]

elev_all <- rep(0, S)
mask <- rep(TRUE, S)
mask[no_nbs] <- FALSE

elev_all[mask] <- elev_raw
elev_all[no_nbs] <- nnbs_elev
elev <- scale(elev_all)[,1]

# temperature
load("snow_temp_full.Rda")
snow_temp <- sce_temp
temp_full <- as.matrix(snow_temp[,-c(1,2)])
temp_scaled <- scale(temp_full)

# time
t_full <- 1:TT
t_scaled_full <- scale(t_full)[,1]
# LOAD CHAINS
np <- import("numpy")

load_bym <- function(prefix){
  eta_list <- list()
  tau_list <- list()
  
  for(c in chains){
    d <- np$load(file.path(BASE_DIR,
                           sprintf("%s_chain%d.npz", prefix, c)))
    
    eta <- d[["eta"]]
    tau <- d[["tau"]]
    
    idx <- seq(1, dim(eta)[2], by = thin2)
    
    eta_list[[length(eta_list)+1]] <- eta[, idx, drop=FALSE]
    tau_list[[length(tau_list)+1]] <- tau[, idx, drop=FALSE]
  }
  
  list(eta=eta_list, tau=tau_list)
}

load_ind <- function(prefix){
  eta_list <- list()
  
  for(c in chains){
    d <- np$load(file.path(BASE_DIR,
                           sprintf("%s_chain%d.npz", prefix, c)))
    
    key <- d$files[[1]]
    eta <- d[[key]]
    
    idx <- seq(1, dim(eta)[2], by = thin2)
    eta_list[[length(eta_list)+1]] <- eta[, idx, drop=FALSE]
  }
  
  eta_list
}
# COMPUTE FUNCTIONS
compute_p_ind <- function(eta01_list, eta10_list){
  
  K <- 4
  results <- list()
  
  for(s in locations){
    for(t in times){
      
      t_s <- t_scaled_full[t]
      
      cov4 <- c(
        1,
        cos(2*pi*(t+1)/period),
        sin(2*pi*(t+1)/period),
        t_s
      )
      
      chains_out <- list()
      
      for(c in 1:length(eta01_list)){
        
        eta01 <- eta01_list[[c]]
        eta10 <- eta10_list[[c]]
        
        phi01 <- rep(0, ncol(eta01))
        phi10 <- rep(0, ncol(eta10))
        
        for(k in 1:K){
          phi01 <- phi01 + eta01[(k-1)*S + s,] * cov4[k]
          phi10 <- phi10 + eta10[(k-1)*S + s,] * cov4[k]
        }
        
        p01 <- 1/(1+exp(-phi01))
        p10 <- 1/(1+exp(-phi10))
        
        if(y[s,t]==0){
          p_final <- p01
        } else {
          p_final <- 1 - p10
        }
        
        chains_out[[c]] <- p_final
      }
      
      results[[paste(s,t,sep="_")]] <- chains_out
    }
  }
  
  results
}

compute_p_bym <- function(eta01_list, tau01_list,
                          eta10_list, tau10_list,
                          use_cov=TRUE){
  
  K_total <- 8
  results <- list()
  
  for(s in locations){
    for(t in times){
      
      w <- t %% 52
      t_s <- t_scaled_full[t]
      
      cov4 <- c(
        1,
        cos(2*pi*(t+1)/period),
        sin(2*pi*(t+1)/period),
        t_s
      )
      
      chains_out <- list()
      
      for(c in 1:length(eta01_list)){
        
        eta01 <- eta01_list[[c]]
        tau01 <- tau01_list[[c]]
        
        eta10 <- eta10_list[[c]]
        tau10 <- tau10_list[[c]]
        
        phi01 <- rep(0, ncol(eta01))
        phi10 <- rep(0, ncol(eta10))
        
        for(k in 1:K_total){
          phi01 <- phi01 +
            eta01[(k-1)*S + s,] *
            tau01[(k-1)*52 + w + 1,] *
            cov4[(k-1)%/%2 + 1]
          
          phi10 <- phi10 +
            eta10[(k-1)*S + s,] *
            tau10[(k-1)*52 + w + 1,] *
            cov4[(k-1)%/%2 + 1]
        }
        
        if(use_cov){
          
          gamma01 <- eta01[(K_total*S+1):(K_total*S+5),]
          gamma10 <- eta10[(K_total*S+1):(K_total*S+5),]
          
          phi01 <- phi01 +
            gamma01[1,]*t_s*lon_na[s] +
            gamma01[2,]*t_s*lon_euas[s] +
            gamma01[3,]*t_s*lat[s] +
            gamma01[4,]*t_s*elev[s] +
            gamma01[5,]*t_s*temp_scaled[s,t]
          
          phi10 <- phi10 +
            gamma10[1,]*t_s*lon_na[s] +
            gamma10[2,]*t_s*lon_euas[s] +
            gamma10[3,]*t_s*lat[s] +
            gamma10[4,]*t_s*elev[s] +
            gamma10[5,]*t_s*temp_scaled[s,t]
        }
        
        p01 <- 1/(1+exp(-phi01))
        p10 <- 1/(1+exp(-phi10))
        
        if(y[s,t]==0){
          p_final <- p01
        } else {
          p_final <- 1 - p10
        }
        
        chains_out[[c]] <- p_final
      }
      
      results[[paste(s,t,sep="_")]] <- chains_out
    }
  }
  
  results
}
# LOAD
bym01_cov <- load_bym("p01_weekly_cov+lon")
bym10_cov <- load_bym("p10_weekly_cov+lon")

bym01 <- load_bym("p01_weekly")
bym10 <- load_bym("p10_weekly")

ind01 <- load_ind("p01_ind_all")
ind10 <- load_ind("p10_ind_all")
# RUN
res_cov <- compute_p_bym(bym01_cov$eta, bym01_cov$tau,
                         bym10_cov$eta, bym10_cov$tau, TRUE)

res_bym <- compute_p_bym(bym01$eta, bym01$tau,
                         bym10$eta, bym10$tau, FALSE)

res_ind <- compute_p_ind(ind01, ind10)
# PLOT 1: MCMC TRACE
plot_all_models <- function(res_cov, res_bym, res_ind){
  
  build_df <- function(results, model_name){
    out <- list()
    
    for(name in names(results)){
      chains <- results[[name]]
      mat <- do.call(cbind, chains)
      
      df <- as.data.frame(mat)
      df$iter <- 1:nrow(df)
      
      df_long <- reshape(df,
                         varying=1:(ncol(df)-1),
                         v.names="value",
                         timevar="chain",
                         times=1:(ncol(df)-1),
                         direction="long")
      
      parts <- strsplit(name, "_")[[1]]
      s <- as.numeric(parts[1])
      t <- as.numeric(parts[2])
      
      df_long$location <- paste0("loc", s)
      df_long$week <- paste0("week", t %% 52 + 1)
      df_long$model <- model_name
      
      out[[length(out)+1]] <- df_long
    }
    
    do.call(rbind, out)
  }
  
  df_all <- rbind(
    build_df(res_bym, "BYM"),
    build_df(res_cov, "BYM+Lon"),
    build_df(res_ind, "Independent")
  )
  
  ggplot(df_all, aes(iter, value, color=factor(chain))) +
    geom_line(alpha=0.6) +
    facet_grid(model ~ location + week) +
    ylim(0,1) +
    scale_color_discrete(name = "Chain") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    )
}

plot_all_models(res_cov, res_bym, res_ind)
# PLOT 2: GAMMA TRACE
extract_gamma <- function(eta_list, label){
  out <- list()
  
  for(c in 1:length(eta_list)){
    eta <- eta_list[[c]]
    
    gamma <- eta[(8*S+1):(8*S+5), , drop=FALSE]
    
    df <- as.data.frame(t(gamma))
    colnames(df) <- c("lon_na","lon_euas","lat","elev","temp")
    
    df$iter <- 1:nrow(df)
    df$chain <- paste0("chain", c)
    df$transition <- label
    
    out[[c]] <- df
  }
  
  do.call(rbind, out)
}

gamma01 <- extract_gamma(bym01_cov$eta, "0→1")
gamma10 <- extract_gamma(bym10_cov$eta, "1→0")

gamma_all <- rbind(gamma01, gamma10)

gamma_long <- pivot_longer(
  gamma_all,
  cols=c(lon_na,lon_euas,lat,elev,temp),
  names_to="covariate",
  values_to="value"
)

ggplot(gamma_long,
       aes(iter, value, color = factor(chain))) +
  geom_line(alpha = 0.6) +
  facet_grid(covariate ~ transition, scales = "free_y") +
  scale_color_discrete(
    name = "Chain",
    labels = 1:length(unique(gamma_long$chain))   # Set chain labels here
  ) +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(
    legend.position = "bottom"
  )
