rm(list = ls())

library(reticulate)
library(ggplot2)
library(patchwork)

# =====================================================
# CONFIG
# =====================================================
BASE_DIR <- "D:/77/Research/temp/snow"

period <- 52
thin2 <- 15

locations <- c(70, 1400)
weeks <- c(20, 35)
year <- 20

times <- (year * 52 + weeks - 1)
chains <- 0:9

# =====================================================
# LOAD DATA
# =====================================================
snow <- readRDS("snow_cleaned_full.Rda")

coords <- as.matrix(snow[,1:2])
y <- as.matrix(snow[,-c(1,2)])

S <- nrow(y)
TT <- ncol(y)

cat("Using FULL S =", S, "TT =", TT, "\n")

# =====================================================
# COVARIATES（完全对齐Python）
# =====================================================

lat <- scale(coords[,2])[,1]

no_nbs <- c(
  57,170,236,269,343,685,946,947,989,
  1037,1084,1090,1109,1118,1127,1176,1203
)

elev_raw <- read.csv(file.path(BASE_DIR,"curr_elev.csv"))[,4]

nnbs_elev <- read.table(
  file.path(BASE_DIR,"nnbs_elev.csv"),
  sep="\t",
  header=TRUE,
  row.names=NULL
)[,3]

elev_all <- rep(0, S)
mask <- rep(TRUE, S)
mask[no_nbs] <- FALSE

elev_all[mask] <- elev_raw
elev_all[no_nbs] <- nnbs_elev

elev <- scale(elev_all)[,1]

load("snow_temp_full.Rda")
snow_temp <- sce_temp

temp_full <- as.matrix(snow_temp[,-c(1,2)])
temp_scaled <- scale(temp_full)

t_full <- 1:TT
t_scaled_full <- scale(t_full)[,1]

# =====================================================
# LOAD CHAINS（稳定R版）
# =====================================================

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
    
    eta <- eta[, idx, drop=FALSE]
    tau <- tau[, idx, drop=FALSE]
    
    eta_list[[length(eta_list)+1]] <- eta
    tau_list[[length(tau_list)+1]] <- tau
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
    eta <- eta[, idx, drop=FALSE]
    
    eta_list[[length(eta_list)+1]] <- eta
  }
  
  eta_list
}

# =====================================================
# POSTERIOR（完全复刻Python）
# =====================================================

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
          
          gamma01 <- eta01[(K_total*S+1):(K_total*S+3),]
          gamma10 <- eta10[(K_total*S+1):(K_total*S+3),]
          
          phi01 <- phi01 +
            gamma01[1,]*t_s*lat[s] +
            gamma01[2,]*t_s*elev[s] +
            gamma01[3,]*t_s*temp_scaled[s,t]
          
          phi10 <- phi10 +
            gamma10[1,]*t_s*lat[s] +
            gamma10[2,]*t_s*elev[s] +
            gamma10[3,]*t_s*temp_scaled[s,t]
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
        
        # ----------------------------
        # linear predictor
        # ----------------------------
        for(k in 1:K){
          
          phi01 <- phi01 + eta01[(k-1)*S + s,] * cov4[k]
          phi10 <- phi10 + eta10[(k-1)*S + s,] * cov4[k]
        }
        
        # ----------------------------
        # probability
        # ----------------------------
        p01 <- 1/(1+exp(-phi01))
        p10 <- 1/(1+exp(-phi10))
        
        # ----------------------------
        # Markov transition（关键）
        # ----------------------------
        if(y[s,t] == 0){
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

# =====================================================
# PLOT
# =====================================================

plot_model <- function(results, name){
  
  plots <- list()
  
  for(s in locations){
    for(t in times){
      
      chains <- results[[paste(s,t,sep="_")]]
      df <- as.data.frame(do.call(cbind, chains))
      
      df$iter <- 1:nrow(df)
      
      df_long <- reshape(df,
                         varying=1:(ncol(df)-1),
                         v.names="value",
                         timevar="chain",
                         times=1:(ncol(df)-1),
                         direction="long")
      
      p <- ggplot(df_long, aes(iter, value, color=factor(chain))) +
        geom_line(alpha=0.6) +
        ggtitle(paste(name,
                      "loc",s,
                      "week", t %% 52 + 1)) +
        ylim(0,1) +
        theme_minimal()
      
      plots[[length(plots)+1]] <- p
    }
  }
  
  wrap_plots(plots)
}

# =====================================================
# RUN
# =====================================================

cat("Loading chains...\n")

bym01_cov <- load_bym("p01_weekly_cov")
bym10_cov <- load_bym("p10_weekly_cov")

bym01 <- load_bym("p01_weekly")
bym10 <- load_bym("p10_weekly")

ind01 <- load_ind("p01_ind_all")
ind10 <- load_ind("p10_ind_all")

cat("Computing...\n")

res_cov <- compute_p_bym(
  bym01_cov$eta, bym01_cov$tau,
  bym10_cov$eta, bym10_cov$tau,
  TRUE
)

res_bym <- compute_p_bym(
  bym01$eta, bym01$tau,
  bym10$eta, bym10$tau,
  FALSE
)

res_ind <- compute_p_ind(
  ind01, ind10
)

cat("Plotting...\n")

cat("Plotting...\n")

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
    build_df(res_cov, "BYM+Cov"),
    build_df(res_ind, "Independent")
  )
  
  ggplot(df_all, aes(iter, value, color=factor(chain))) +
    geom_line(alpha=0.6) +
    facet_grid(model ~ location + week) +
    ylim(0,1) +
    scale_color_discrete(
      name = "Chain",
      labels = 0:9
    ) +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.background = element_rect(fill = "grey85", color = "black"),
      strip.text = element_text(size = 10, face = "bold"),
      panel.spacing = unit(0.8, "lines"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_line(color = "grey90", size = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      strip.placement = "outside"
    )
}

plot_all_models(res_cov, res_bym, res_ind)

extract_gamma <- function(eta_list, label){
  
  # gamma 在最后 3 行
  # (8*S+1):(8*S+3)
  
  out <- list()
  
  for(c in 1:length(eta_list)){
    
    eta <- eta_list[[c]]
    
    gamma <- eta[(8*S+1):(8*S+3), , drop=FALSE]   # 3 × M
    
    df <- as.data.frame(t(gamma))   # M × 3
    colnames(df) <- c("lat_trend", "elev_trend", "temp_trend")
    
    df$iter <- 1:nrow(df)
    df$chain <- paste0("chain", c)
    df$transition <- label
    
    out[[c]] <- df
  }
  
  do.call(rbind, out)
}

gamma01_df <- extract_gamma(bym01_cov$eta, "0→1")
gamma10_df <- extract_gamma(bym10_cov$eta, "1→0")

gamma_all <- rbind(gamma01_df, gamma10_df)

library(tidyr)

gamma_long <- pivot_longer(
  gamma_all,
  cols = c(lat_trend, elev_trend, temp_trend),
  names_to = "covariate",
  values_to = "value"
)

gamma_long$covariate <- factor(
  gamma_long$covariate,
  levels = c("lat_trend", "elev_trend", "temp_trend"),
  labels = c("lat × trend", "elev × trend", "temp × trend")
)
ggplot(gamma_long,
       aes(iter, value, color = factor(chain))) +
  
  geom_line(alpha = 0.6) +
  
  facet_grid(covariate ~ transition, switch = "y") +
  
  labs(
    title = "",
    x = "Iteration",
    y = "Coefficient value"
  ) +
  
  ylim(-0.35, 0.5) +
  
  scale_color_discrete(
    name = "Chain",
    labels = 0:9
  ) +
  
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    
    strip.background = element_rect(fill = "grey85", color = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    
    panel.spacing = unit(0.8, "lines"),
    
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_line(color = "grey90", size = 0.3),
    
    legend.position = "bottom",
    strip.placement = "outside",
    legend.key = element_blank()
  )

