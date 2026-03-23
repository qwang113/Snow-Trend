rm(list = ls())
library(sf)
library(dplyr)
library(ggplot2)
library(reticulate)
library(Matrix)
library(igraph)
BASE_DIR <- "D:/77/Research/temp/snow"
period <- 52
thin <- 15

locations <- c(70, 1400)
weeks <- c(20, 35)
year <- 20

times <- (year * 52 + weeks - 1)
chains <- 0:9

no_nbs <- c(
  57,170,236,269,343,685,946,947,989,
  1037,1084,1090,1109,1118,1127,1176,1203
)
snow_cleaned_full <- readRDS("snow_cleaned_full.Rda")

coords_full <- as.matrix(snow_cleaned_full[,1:2])
y_full <- as.matrix(snow_cleaned_full[,-c(1,2)])

coords_tmp <- coords_full[-no_nbs,]
y_tmp <- y_full[-no_nbs,]


coords_sf <- st_as_sf(
  data.frame(lon=coords_tmp[,1], lat=coords_tmp[,2]),
  coords=c("lon","lat"),
  crs=4326
)

coords_proj <- st_transform(coords_sf, "+proj=aeqd +lat_0=90 +lon_0=-100")

xy <- st_coordinates(coords_proj) / 1e6

D <- as.matrix(dist(xy))
W <- (D <= 0.22) * 1
diag(W) <- 0

g <- graph_from_adjacency_matrix(W, mode="undirected")

comp <- components(g)
sizes <- comp$csize
order <- order(sizes, decreasing=TRUE)

keep <- which(comp$membership %in% order[1:2])

coords <- coords_tmp[keep,]
y <- y_tmp[keep,]

S <- nrow(y)
TT <- ncol(y)

cat("Using S =", S, "\n")

# ---- lat
lat <- scale(coords[,2])[,1]

# ---- elevation
elev_raw <- read.csv(file.path(BASE_DIR,"curr_elev.csv"))[,4]
elev <- scale(elev_raw[keep])[,1]

# ---- temperature
load(file.path(BASE_DIR,"snow_temp_full.Rda"))
temp_full <- as.matrix(sce_temp[-no_nbs, -c(1,2)])
temp <- scale(temp_full[keep,])  # same global scaling

# ---- time
t_full <- 1:TT
t_scaled <- scale(t_full)[,1]

use_python("D:/anaconda3/envs/CPD/python.exe", required = TRUE)
np <- import("numpy")
pickle <- import("pickle")

load_bym <- function(prefix, use_cov=FALSE){
  
  eta_list <- list()
  tau_list <- list()
  
  for(c in chains){
    
    fname <- if(use_cov){
      sprintf("%s_chain%d+cov.pkl", prefix, c)
    } else {
      sprintf("%s_chain%d.pkl", prefix, c)
    }
    
    path <- file.path(BASE_DIR, fname)
    py_file <- import_builtins()$open(path, "rb")
    d <- pickle$load(py_file)
    py_file$close()
    
    eta <- d$eta[, seq(1, dim(d$eta)[2], by=thin)]
    tau <- d$tau[, seq(1, dim(d$tau)[2], by=thin)]
    
    eta_list[[length(eta_list)+1]] <- eta
    tau_list[[length(tau_list)+1]] <- tau
  }
  
  list(eta=eta_list, tau=tau_list)
}

load_ind <- function(prefix){
  
  eta_list <- list()
  
  for(c in chains){
    
    path <- file.path(BASE_DIR, sprintf("%s_chain%d.npz", prefix, c))
    
    d <- np$load(path)
    
    key <- d$files[[1]]
    eta <- d[[key]][ , seq(1, dim(d[[key]])[2], by=thin)]
    
    eta_list[[length(eta_list)+1]] <- eta
  }
  
  eta_list
}

compute_p_bym <- function(eta01_list, tau01_list,
                          eta10_list, tau10_list){
  
  results <- list()
  
  for(s in locations){
    for(t in times){
      
      w <- t %% 52
      t_s <- t_scaled[t]
      
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
        
        # ----------------------------
        # spatial weekly part
        # ----------------------------
        for(k in 1:8){
          
          beta01 <- eta01[(k-1)*S + s,] * tau01[(k-1)*52 + w,]
          beta10 <- eta10[(k-1)*S + s,] * tau10[(k-1)*52 + w,]
          
          phi01 <- phi01 + beta01 * cov4[(k-1)%/%2 + 1]
          phi10 <- phi10 + beta10 * cov4[(k-1)%/%2 + 1]
        }
        
        # ----------------------------
        # covariates
        # ----------------------------
        gamma01 <- eta01[(8*S+1):(8*S+3),]
        gamma10 <- eta10[(8*S+1):(8*S+3),]
        
        phi01 <- phi01 + gamma01[1,]*t_s*lat[s]
        phi01 <- phi01 + gamma01[2,]*t_s*elev[s]
        phi01 <- phi01 + gamma01[3,]*t_s*temp[s,t]
        
        phi10 <- phi10 + gamma10[1,]*t_s*lat[s]
        phi10 <- phi10 + gamma10[2,]*t_s*elev[s]
        phi10 <- phi10 + gamma10[3,]*t_s*temp[s,t]
        
        # ----------------------------
        # probability (CRITICAL)
        # ----------------------------
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

compute_p_bym_nocov <- function(eta01_list, tau01_list,
                                eta10_list, tau10_list){
  
  results <- list()
  
  for(s in locations){
    for(t in times){
      
      w <- t %% 52
      t_s <- t_scaled[t]
      
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
        
        for(k in 1:8){
          
          beta01 <- eta01[(k-1)*S + s,] * tau01[(k-1)*52 + w,]
          beta10 <- eta10[(k-1)*S + s,] * tau10[(k-1)*52 + w,]
          
          phi01 <- phi01 + beta01 * cov4[(k-1)%/%2 + 1]
          phi10 <- phi10 + beta10 * cov4[(k-1)%/%2 + 1]
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
  
  K <- 4  # cov4
  results <- list()
  
  for(s in locations){
    for(t in times){
      
      t_s <- t_scaled[t]
      
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
        
        # -------- IND linear predictor --------
        for(k in 1:K){
          
          beta01 <- eta01[(k-1)*S + s,]
          beta10 <- eta10[(k-1)*S + s,]
          
          phi01 <- phi01 + beta01 * cov4[k]
          phi10 <- phi10 + beta10 * cov4[k]
        }
        
        # -------- probability --------
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

build_trace_df <- function(results, model){
  
  df <- data.frame()
  
  for(name in names(results)){
    
    parts <- strsplit(name, "_")[[1]]
    s <- parts[1]
    t <- parts[2]
    
    chains_list <- results[[name]]
    
    for(c in 1:length(chains_list)){
      
      vals <- chains_list[[c]]
      
      df <- rbind(df, data.frame(
        iter = 1:length(vals),
        value = vals,
        chain = factor(paste0("chain", c)),
        location = paste0("loc", s),
        week = paste0("week", (as.numeric(t) %% 52)+1),
        model = model
      ))
    }
  }
  
  df
}



bym01_cov <- load_bym("bym01", TRUE)
bym10_cov <- load_bym("bym10", TRUE)

bym01 <- load_bym("bym01", FALSE)
bym10 <- load_bym("bym10", FALSE)

ind01 <- load_ind("ind01")
ind10 <- load_ind("ind10")

res_cov <- compute_p_bym(
  bym01_cov$eta, bym01_cov$tau,
  bym10_cov$eta, bym10_cov$tau
)

res_bym <- compute_p_bym_nocov(
  bym01$eta, bym01$tau,
  bym10$eta, bym10$tau
)
res_ind <- compute_p_ind(ind01, ind10)
df_ind <- build_trace_df(res_ind, "Independent")
df_cov <- build_trace_df(res_cov, "BYM+Cov")
df_bym <- build_trace_df(res_bym, "BYM")
df_all <- rbind(df_ind, df_cov, df_bym)

ggplot(df_all, aes(x=iter, y=value, color=chain)) +
  geom_line(alpha=0.6, size=0.3) +
  facet_grid(model ~ location + week) +
  theme_bw() +
  labs(
    x="Iteration",
    y="Predicted probability of locations covered by snow"
  ) +
  theme(
    legend.position="bottom",
    plot.title = element_text(hjust=0.5)
  )



extract_gamma_trace <- function(eta_list, S, label_prefix){
  
  cov_names <- c("lat × trend", "elev × trend", "temp × trend")
  
  df <- data.frame()
  
  for(c in 1:length(eta_list)){
    
    eta <- eta_list[[c]]
    
    gamma <- eta[(8*S+1):(8*S+3), ]   # 3 × M
    
    for(k in 1:3){
      
      vals <- gamma[k, ]
      
      df <- rbind(df, data.frame(
        iter = 1:length(vals),
        value = vals,
        chain = factor(paste0("chain", c)),
        covariate = cov_names[k],
        type = label_prefix   # "01" or "10"
      ))
    }
  }
  
  df
}
df_gamma01 <- extract_gamma_trace(bym01_cov$eta, S, "0→1")
df_gamma10 <- extract_gamma_trace(bym10_cov$eta, S, "1→0")

df_gamma_all <- rbind(df_gamma01, df_gamma10)

df_gamma_all$type <- factor(df_gamma_all$type, levels = c("0→1", "1→0"))
df_gamma_all$covariate <- factor(
  df_gamma_all$covariate,
  levels = c("lat × trend", "elev × trend", "temp × trend")
)

ggplot(df_gamma_all, aes(x=iter, y=value, color=chain)) +
  geom_line(alpha=0.6, size=0.3) +
  facet_grid(covariate ~ type) +
  theme_bw() +
  labs(
    x="Iteration",
    y="Coefficient value",
    title="Traceplots of Covariate Effects"
  ) +
  theme(
    legend.position="bottom",
    
    strip.text.y = element_text(size=13, face="bold"),
    
    strip.text.x = element_text(size=13, face="bold"),
    
    plot.title = element_text(size=14, hjust=0.5, face="bold")
  )
