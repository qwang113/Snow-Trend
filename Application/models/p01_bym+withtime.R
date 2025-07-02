rm(list = ls())
library(sf)
library(fields)
library(sp)
library(spdep)
library(fields)
library(spatstat)
library(reshape2)
library(ggplot2)
library(maps)
library(spBayes)
library(BayesLogit)
library(Matrix)
library(future.apply)
library(pbapply)
library(sparseMVN) 
setwd(here::here())
no_nbs <- c(57, 170, 236, 269, 343, 685, 946, 947, 989, 1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203)
all_y <- rbind(readRDS("snow_cleaned_full.Rda")[-no_nbs,],readRDS("snow_cleaned_full.Rda")[no_nbs,])
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]
elev <- scale(c(read.csv(here::here("curr_elev.csv"))[,4],
          read.csv(here::here("nnbs_elev.csv"),sep = "\t", row.names = NULL)[,4]))
lats <- scale(coords[,2])
S <- nrow(y)
TT <- ncol(y)
sf_coords <- st_as_sf(data.frame(coords), coords = c("LON", "LAT"), crs = 4326)
sf_data_aeqd <- st_transform(sf_coords, crs = "+proj=aeqd +lat_0=90 +lon_0=-100")
transformed_coordinates <- st_coordinates(sf_data_aeqd)

dif <- st_coordinates(sf_data_aeqd[2,]) - st_coordinates(sf_data_aeqd[3,])
theta <- atan2(dif[2], dif[1])

rotate_points <- function(coords, angle) {
  # Create a rotation matrix for 2D rotation
  rotation_matrix <- matrix(c(cos(angle), -sin(angle), 
                              sin(angle), cos(angle)), 
                            ncol = 2)
  
  # Apply the rotation to the coordinates
  rotated_coords <- t(rotation_matrix %*% t(coords))
  
  return(rotated_coords)
}

rotated_coordinates <- data.frame(rotate_points(transformed_coordinates, theta))/1e6

Distances <- pairdist(rotated_coordinates)
Omg <- Matrix(0, nrow = nrow(coords), ncol = nrow(coords),sparse = TRUE)
Omg[which(Distances <= 0.22)] = 1
Omg <- Omg - diag(nrow = nrow(coords))
# D <- diag(rowSums(Omg))
D <- Matrix(diag(rowSums(Omg)) ,sparse = TRUE)
sigma = 1
period = 52

location_time_0 <- which(y[,-ncol(y)]==0, arr.ind =  TRUE)
row_idx <- location_time_0[,1]
next_y <- y[cbind(location_time_0[,1], location_time_0[,2]+1)]
design_mat <- Matrix(0, nrow = length(next_y), ncol = 8*S, sparse = TRUE)
location_idx <- sparse.model.matrix(~ factor(row) - 1, data = data.frame(location_time_0))

# pb <- txtProgressBar(min = 0, max = nrow(design_mat), style = 3)
covariates <- Matrix(
  cbind(1,1,
        cos(2*pi*location_time_0[,2]/period),
        cos(2*pi*location_time_0[,2]/period),
        sin(2*pi*location_time_0[,2]/period), 
        sin(2*pi*location_time_0[,2]/period), 
        location_time_0[,2], location_time_0[,2]), sparse = TRUE )

# # Loop through chunks of rows
# chunk_size <- 10
# curr_row <- 1
# for (start_row in seq(curr_row , nrow(location_idx), by = chunk_size)) {
#   print(start_row)
#   # Define the end row for the current chunk
#   end_row <- min(start_row + chunk_size - 1, nrow(location_idx))
# 
#   # Extract the relevant rows from location_idx and covariates
#   loc_chunk <- location_idx[start_row:end_row, , drop = FALSE]
#   cov_chunk <- covariates[start_row:end_row, , drop = FALSE]
# 
#   # Apply the outer product to each pair of rows in the chunk
#   # mapply to compute the outer product for each pair of rows
#   result_chunk <- mapply(function(loc, cov) as.vector(outer(loc, cov, "*")),
#                          split(loc_chunk, row(loc_chunk)),
#                          split(cov_chunk, row(cov_chunk)))
# 
#   # Reshape result to match design matrix row structure
#   result_matrix <- matrix(result_chunk, nrow = end_row - start_row + 1, byrow = TRUE)
# 
#   # Convert result matrix to a sparse Matrix format and store in design_mat
#   design_mat[start_row:end_row, ] <- Matrix(result_matrix, sparse = TRUE)
# }
# 
# saveRDS(design_mat,"design_mat_01.Rda")
design_mat <- readRDS("D:/77/Research/temp/snow/design01.Rda")
lats_design <- lats[row_idx]
elev_design <- elev[row_idx]

other_covariates <- Matrix(
  cbind(lats_design, 
        lats_design*cos(2*pi*location_time_0[,2]/period),
        lats_design*sin(2*pi*location_time_0[,2]/period),
        lats_design*location_time_0[,2],
        elev_design,
        elev_design*cos(2*pi*location_time_0[,2]/period),
        elev_design*sin(2*pi*location_time_0[,2]/period),
        elev_design*location_time_0[,2]
        ), sparse = TRUE )

design_mat <- cbind(design_mat, other_covariates)
tot_samples <- 2000

all_theta <- matrix(NA, nrow = 8*S + 8, ncol = tot_samples)
all_tau <- matrix(NA, nrow = 8, ncol = tot_samples)
curr_theta_vec <- matrix(0, nrow = 1, ncol = 8*S + 8)
curr_tau_vec <- rep(1,8)
a_tau <- 0.001
b_tau <- 0.001
curr_idx <- 0
save_idx <- 0
burn = 0
thin = 1
prec <- D - Omg

while(save_idx < tot_samples) {
  curr_idx = curr_idx + 1
  print(curr_idx)
  # Sample current Omega
  kappas <- next_y - 1/2
  curr_phi <- design_mat %*% t(curr_theta_vec)
  curr_omega <-  rpg(length(next_y), h = 1, z = as.numeric(curr_phi))
  # Sample current theta
  curr_prec <- bdiag(
    1/curr_tau_vec[1]*prec,
    1/curr_tau_vec[2]*diag(1,S),
    1/curr_tau_vec[3]*prec,
    1/curr_tau_vec[4]*diag(1,S),
    1/curr_tau_vec[5]*prec,
    1/curr_tau_vec[6]*diag(1,S),
    1/curr_tau_vec[7]*prec,
    1/curr_tau_vec[8]*diag(1,S),
    1/10000*diag(1,8)
  )                  
  xtxomg <- t(design_mat)%*% Diagonal(length(curr_omega), curr_omega)%*%(design_mat)
  pos_prec <- xtxomg + curr_prec
  CH <- Cholesky(pos_prec, LDL = FALSE)
  b <- t(design_mat) %*% kappas
  # Solve pos_prec %*% x = b for x
  pos_mu <- solve(CH, b)
  # pos_mu_2 <- solve(pos_prec, b)
  curr_theta_vec <- rmvn.sparse(1, mu = pos_mu, CH = CH, prec = TRUE)
  
  # Sample taus
  curr_tau_vec[1] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[1:S])%*%prec%*%curr_theta_vec[1:S]/2) )
  curr_tau_vec[2] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(S+1):(2*S)])%*%curr_theta_vec[(S+1):(2*S)]/2 ) )
  
  curr_tau_vec[3] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(2*S+1):(3*S)])%*%prec%*%curr_theta_vec[(2*S+1):(3*S)]/2 ) )
  curr_tau_vec[4] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(3*S+1):(4*S)])%*%curr_theta_vec[(3*S+1):(4*S)]/2 ) )
  
  curr_tau_vec[5] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(4*S+1):(5*S)])%*%prec%*%curr_theta_vec[(4*S+1):(5*S)] /2) )
  curr_tau_vec[6] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(5*S+1):(6*S)])%*%curr_theta_vec[(5*S+1):(6*S)]/2 ) )
  
  curr_tau_vec[7] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(6*S+1):(7*S)])%*%prec%*%curr_theta_vec[(6*S+1):(7*S)]/2 ) )
  curr_tau_vec[8] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(7*S+1):(8*S)])%*%curr_theta_vec[(7*S+1):(8*S)]/2 ) )
  
  if((curr_idx > burn) & (curr_idx %% thin == 0) ){
    save_idx <- save_idx + 1
    all_theta[,save_idx] <- as.vector(curr_theta_vec)
    all_tau[,save_idx] <-  as.vector(curr_tau_vec)
  }
}

setwd("D:/77/Research/temp/snow/")
saveRDS(all_theta, "theta01_bym+withtime.Rda")
saveRDS(all_tau, "tau01_bym+withtime.Rda")

