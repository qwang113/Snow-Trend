rm(list = ls())
library(sf)
library(fields)
library(sp)
library(spdep)
library(fields)
library(spatstat)
library(reshape2)
library(ggplot2)
library(rstan)
library(maps)
library(spBayes)
library(BayesLogit)
library(Matrix)
library(future.apply)
library(pbapply)
library(sparseMVN) 
setwd(here::here())
no_nbs <- c(57, 170, 236, 269, 343, 685, 946, 947, 989, 1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203)
all_y <- readRDS("snow_cleaned.Rda")[-no_nbs,]
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]

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
eps = 1e-7
period = 52

location_time_0 <- which(y[,-ncol(y)]==0, arr.ind =  TRUE)
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
design_mat <- readRDS("D:/77/Research/temp/snow_trend/design_mat_01.Rda")
tot_samples <- 5000

all_theta <- matrix(NA, nrow = 8*S, ncol = tot_samples)
all_tau <- matrix(NA, nrow = 8, ncol = tot_samples)
curr_theta_vec <- matrix(0, nrow = 1, ncol = 8*S)
curr_tau_vec <- rep(0.1,8)
a_tau <- 0.001
b_tau <- 0.001
curr_idx <- 0
save_idx <- 0
burn = 0
thin = 1
 
cov_inv <- solve(D - Omg + diag(eps, S))
while(save_idx < tot_samples) {
    curr_idx = curr_idx + 1
    print(curr_idx)
    # Sample current Omega
    kappas <- next_y - 1/2
    curr_phi <- design_mat %*% t(curr_theta_vec)
    curr_omega <-  rpg(length(next_y), h = 1, z = as.numeric(curr_phi))
    # Sample current theta
    curr_B <- bdiag(
      curr_tau_vec[1]*cov_inv,
      curr_tau_vec[2]*diag(1,S),
      curr_tau_vec[3]*cov_inv,
      curr_tau_vec[4]*diag(1,S),
      curr_tau_vec[5]*cov_inv,
      curr_tau_vec[6]*diag(1,S),
      curr_tau_vec[7]*cov_inv,
      curr_tau_vec[8]*diag(1,S)
    )                  
    B_inv <- solve(curr_B)
    xtxomg <- t(design_mat)%*% Diagonal(length(curr_omega), curr_omega)%*%(design_mat)
    
    pos_prec <- xtxomg + B_inv
    CH <- Cholesky(pos_prec, LDL = FALSE)
    b <- t(design_mat) %*% kappas
    # Solve pos_prec %*% x = b for x
    pos_mu <- solve(CH, b)
    # pos_mu_2 <- solve(pos_prec, b)
    curr_theta_vec <- rmvn.sparse(1, mu = pos_mu, CH = CH, prec = TRUE)
    
    # curr_theta_vec[1:S] <- polya_aug_sampling(design_mat[1:S,1:S], xtxomg[1:S,1:S], B_inv[1:S,1:S], kappas[1:S])
    # 
    # curr_theta_vec[(S+1):(2*S)] <- polya_aug_sampling(design_mat[(S+1):(2*S),(S+1):(2*S)], xtxomg[(S+1):(2*S),(S+1):(2*S)],
    #                                                   B_inv[(S+1):(2*S),(S+1):(2*S)], kappas[(S+1):(2*S)])
    # 
    # curr_theta_vec[(2*S+1):(3*S)] <- polya_aug_sampling(design_mat[(2*S+1):(3*S),(2*S+1):(3*S)], xtxomg[(2*S+1):(3*S),(2*S+1):(3*S)],
    #                                                   B_inv[(2*S+1):(3*S),(2*S+1):(3*S)], kappas[(2*S+1):(3*S)])
    # 
    # curr_theta_vec[(3*S+1):(4*S)] <- polya_aug_sampling(design_mat[(3*S+1):(4*S),(3*S+1):(4*S)], xtxomg[(3*S+1):(4*S),(3*S+1):(4*S)],
    #                                                     B_inv[(3*S+1):(4*S),(3*S+1):(4*S)], kappas[(3*S+1):(4*S)])
    # 
    # curr_theta_vec[(4*S+1):(5*S)] <- polya_aug_sampling(design_mat[(4*S+1):(5*S),(4*S+1):(5*S)], xtxomg[(4*S+1):(5*S),(4*S+1):(5*S)],
    #                                                     B_inv[(4*S+1):(5*S),(4*S+1):(5*S)], kappas[(4*S+1):(5*S)])
    # 
    # curr_theta_vec[(5*S+1):(6*S)] <- polya_aug_sampling(design_mat[(5*S+1):(6*S),(5*S+1):(6*S)], xtxomg[(5*S+1):(6*S),(5*S+1):(6*S)],
    #                                                     B_inv[(5*S+1):(6*S),(5*S+1):(6*S)], kappas[(5*S+1):(6*S)])
    # 
    # curr_theta_vec[(6*S+1):(7*S)] <- polya_aug_sampling(design_mat[(6*S+1):(7*S),(6*S+1):(7*S)], xtxomg[(6*S+1):(7*S),(6*S+1):(7*S)],
    #                                                     B_inv[(6*S+1):(7*S),(6*S+1):(7*S)], kappas[(6*S+1):(7*S)])
    # 
    # curr_theta_vec[(7*S+1):(8*S)] <- polya_aug_sampling(design_mat[(7*S+1):(8*S),(7*S+1):(8*S)], xtxomg[(7*S+1):(8*S),(7*S+1):(8*S)],
    #                                                     B_inv[(7*S+1):(8*S),(7*S+1):(8*S)], kappas[(7*S+1):(8*S)])
    
    
    
    
    # Sample taus
    curr_tau_vec[1] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[1:S])%*%(D-Omg + diag(eps, S))%*%curr_theta_vec[1:S]/2) )
    curr_tau_vec[2] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(S+1):(2*S)])%*%curr_theta_vec[(S+1):(2*S)]/2 ) )
    
    curr_tau_vec[3] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(2*S+1):(3*S)])%*%(D-Omg + diag(eps, S))%*%curr_theta_vec[(2*S+1):(3*S)]/2 ) )
    curr_tau_vec[4] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(3*S+1):(4*S)])%*%curr_theta_vec[(3*S+1):(4*S)]/2 ) )
    
    curr_tau_vec[5] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(4*S+1):(5*S)])%*%(D-Omg + diag(eps, S))%*%curr_theta_vec[(4*S+1):(5*S)] /2) )
    curr_tau_vec[6] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(5*S+1):(6*S)])%*%curr_theta_vec[(5*S+1):(6*S)]/2 ) )
    
    curr_tau_vec[7] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(6*S+1):(7*S)])%*%(D-Omg + diag(eps, S))%*%curr_theta_vec[(6*S+1):(7*S)]/2 ) )
    curr_tau_vec[8] <- 1/rgamma(1, shape = a_tau+S/2, rate = as.numeric(b_tau + t(curr_theta_vec[(7*S+1):(8*S)])%*%curr_theta_vec[(7*S+1):(8*S)]/2 ) )
    
    if((curr_idx > burn) & (curr_idx %% thin == 0) ){
      save_idx <- save_idx + 1
      all_theta[,save_idx] <- as.vector(curr_theta_vec)
      all_tau[,save_idx] <-  as.vector(curr_tau_vec)
    }
}
theta_mean <- apply(all_theta, 1, mean)
theta_01 <- theta_mean[1:S]
theta_02 <- theta_mean[(S+1):(2*S)]
theta_11 <- theta_mean[(2*S+1):(3*S)]
theta_12 <- theta_mean[(3*S+1):(4*S)]
theta_21 <- theta_mean[(4*S+1):(5*S)]
theta_22 <- theta_mean[(5*S+1):(6*S)]
theta_a1 <- theta_mean[(6*S+1):(7*S)]
theta_a2 <- theta_mean[(7*S+1):(8*S)]

saveRDS(all_theta, "self_theta.Rda")
saveRDS(all_tau, "self_tau.Rda")
# Compare with stan
# samples <- readRDS(here::here("test_app_samples.Rda"))
# theta_01_stan <- apply(samples$theta_01, 2, mean)
# theta_02_stan <- apply(samples$theta_02, 2, mean)
# theta_11_stan <- apply(samples$theta_11, 2, mean)
# theta_12_stan <- apply(samples$theta_12, 2, mean)
# theta_21_stan <- apply(samples$theta_21, 2, mean)
# theta_22_stan <- apply(samples$theta_22, 2, mean)
# theta_a1_stan <- apply(samples$theta_a1, 2, mean)
# theta_a2_stan <- apply(samples$theta_a2, 2, mean)
# 
# par(mfrow = c(4,2))
# boxplot(cbind(theta_01, theta_01_stan))
# boxplot(cbind(theta_02, theta_02_stan))
# boxplot(cbind(theta_11, theta_11_stan))
# boxplot(cbind(theta_12, theta_12_stan))
# boxplot(cbind(theta_21, theta_21_stan))
# boxplot(cbind(theta_22, theta_22_stan))
# boxplot(cbind(theta_a1, theta_a1_stan))
# boxplot(cbind(theta_a2, theta_a2_stan))