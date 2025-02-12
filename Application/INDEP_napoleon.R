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
napoleon <- c(1433,2)
all_y <- readRDS("snow_cleaned.Rda")
all_y <- all_y[napoleon,]
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]
S <- nrow(y)
TT <- ncol(y)


period = 52

location_time_0 <- which(y[,-ncol(y)]==0, arr.ind =  TRUE)
row_idx <- location_time_0[,1]
next_y <- y[cbind(location_time_0[,1], location_time_0[,2]+1)]

design_mat <- Matrix(0, nrow = length(next_y), ncol = 4*S, sparse = TRUE)
location_idx <- sparse.model.matrix(~ factor(row) - 1, data = data.frame(location_time_0))

# pb <- txtProgressBar(min = 0, max = nrow(design_mat), style = 3)
covariates <- Matrix(
  cbind(1,
        cos(2*pi*location_time_0[,2]/period),
        sin(2*pi*location_time_0[,2]/period), 
        location_time_0[,2]), sparse = TRUE )

# # Loop through chunks of rows
chunk_size <- 10
curr_row <- 1
for (start_row in seq(curr_row , nrow(location_idx), by = chunk_size)) {
  print(start_row)
  # Define the end row for the current chunk
  end_row <- min(start_row + chunk_size - 1, nrow(location_idx))

  # Extract the relevant rows from location_idx and covariates
  loc_chunk <- location_idx[start_row:end_row, , drop = FALSE]
  cov_chunk <- covariates[start_row:end_row, , drop = FALSE]

  # Apply the outer product to each pair of rows in the chunk
  # mapply to compute the outer product for each pair of rows
  result_chunk <- mapply(function(loc, cov) as.vector(outer(loc, cov, "*")),
                         split(loc_chunk, row(loc_chunk)),
                         split(cov_chunk, row(cov_chunk)))

  # Reshape result to match design matrix row structure
  result_matrix <- matrix(result_chunk, nrow = end_row - start_row + 1, byrow = TRUE)

  # Convert result matrix to a sparse Matrix format and store in design_mat
  design_mat[start_row:end_row, ] <- Matrix(result_matrix, sparse = TRUE)
}
# 
# saveRDS(design_mat,"design_mat_01_INDEP.Rda")
# design_mat <- readRDS("D:/77/Research/temp/snow_trend/design_mat_01_INDEP.Rda")
tot_samples <- 1000

all_theta <- matrix(NA, nrow = 4*S, ncol = tot_samples)
curr_theta_vec <- matrix(0, nrow = 1, ncol = 4*S)
curr_idx <- 0
save_idx <- 0
burn = 200
thin = 2

while(save_idx < tot_samples) {
  curr_idx = curr_idx + 1
  print(curr_idx)
  # Sample current Omega
  kappas <- next_y - 1/2
  curr_phi <- design_mat %*% t(curr_theta_vec)
  curr_omega <-  rpg(length(next_y), h = rep(1, length(curr_phi)), z = as.numeric(curr_phi))
  # Sample current theta
  curr_prec <- bdiag(
    1/100*diag(1,S),
    1/100*diag(1,S),
    1/100*diag(1,S),
    1/100*diag(1,S)
  )                  
  xtxomg <- t(design_mat)%*%Diagonal(length(curr_omega), curr_omega)%*%(design_mat)
  pos_prec <- xtxomg + curr_prec
  CH <- Cholesky(pos_prec, LDL = FALSE)
  b <- t(design_mat) %*% kappas
  # Solve pos_prec %*% x = b for x
  pos_mu <- solve(CH, b)
  # pos_mu_2 <- solve(pos_prec, b)
  curr_theta_vec <- rmvn.sparse(1, mu = pos_mu, CH = CH, prec = TRUE)
  
  if((curr_idx > burn) & (curr_idx %% thin == 0) ){
    save_idx <- save_idx + 1
    all_theta[,save_idx] <- as.vector(curr_theta_vec)
  }
}

p01_params <- all_theta


# p10 -----------------------------------------


location_time_0 <- which(y[,-ncol(y)]==1, arr.ind =  TRUE)
row_idx <- location_time_0[,1]
next_y <- abs(y[cbind(location_time_0[,1], location_time_0[,2]+1)] - 1)


design_mat <- Matrix(0, nrow = length(next_y), ncol = 4*S, sparse = TRUE)
location_idx <- sparse.model.matrix(~ factor(row) - 1, data = data.frame(location_time_0))

# pb <- txtProgressBar(min = 0, max = nrow(design_mat), style = 3)
covariates <- Matrix(
  cbind(1,
        cos(2*pi*location_time_0[,2]/period),
        sin(2*pi*location_time_0[,2]/period), 
        location_time_0[,2]), sparse = TRUE )

# # Loop through chunks of rows

chunk_size <- 10
curr_row <- 1
for (start_row in seq(curr_row , nrow(location_idx), by = chunk_size)) {
  print(start_row)
  # Define the end row for the current chunk
  end_row <- min(start_row + chunk_size - 1, nrow(location_idx))
  
  # Extract the relevant rows from location_idx and covariates
  loc_chunk <- location_idx[start_row:end_row, , drop = FALSE]
  cov_chunk <- covariates[start_row:end_row, , drop = FALSE]
  
  # Apply the outer product to each pair of rows in the chunk
  # mapply to compute the outer product for each pair of rows
  result_chunk <- mapply(function(loc, cov) as.vector(outer(loc, cov, "*")),
                         split(loc_chunk, row(loc_chunk)),
                         split(cov_chunk, row(cov_chunk)))
  
  # Reshape result to match design matrix row structure
  result_matrix <- matrix(result_chunk, nrow = end_row - start_row + 1, byrow = TRUE)
  
  # Convert result matrix to a sparse Matrix format and store in design_mat
  design_mat[start_row:end_row, ] <- Matrix(result_matrix, sparse = TRUE)
}
 

tot_samples <- 1000

all_theta <- matrix(NA, nrow = 4*S, ncol = tot_samples)
curr_theta_vec <- matrix(0, nrow = 1, ncol = 4*S)
curr_idx <- 0
save_idx <- 0
burn = 200
thin = 2

while(save_idx < tot_samples) {
  curr_idx = curr_idx + 1
  print(curr_idx)
  # Sample current Omega
  kappas <- next_y - 1/2
  curr_phi <- design_mat %*% t(curr_theta_vec)
  curr_omega <-  rpg(length(next_y), h = rep(1, length(curr_phi)), z = as.numeric(curr_phi))
  # Sample current theta
  curr_prec <- bdiag(
    1/100*diag(1,S),
    1/100*diag(1,S),
    1/100*diag(1,S),
    1/100*diag(1,S)
  )                  
  xtxomg <- t(design_mat)%*% Diagonal(length(curr_omega), curr_omega)%*%(design_mat)
  pos_prec <- (xtxomg+t(xtxomg))/2 + curr_prec
  CH <- Cholesky(pos_prec, LDL = FALSE)
  b <- t(design_mat) %*% kappas
  # Solve pos_prec %*% x = b for x
  pos_mu <- solve(CH, b)
  curr_theta_vec <- rmvn.sparse(1, mu = pos_mu, CH = CH, prec = TRUE)
  
  if((curr_idx > burn) & (curr_idx %% thin == 0) ){
    save_idx <- save_idx + 1
    all_theta[,save_idx] <- as.vector(curr_theta_vec)
  }
}

p10_params <- all_theta

# -----------------------------------------compare the likelihood
theta <- p01_params
theta_s <- p10_params

beta_0_hat <- mean(theta[1,])
beta_1_hat <- mean(theta[3,])
beta_2_hat <- mean(theta[5,])
alpha_hat <- mean(theta[7,])
beta_0s_hat <- mean(theta_s[1,])
beta_1s_hat <- mean(theta_s[3,])
beta_2s_hat <- mean(theta_s[5,])
alpha_s_hat <- mean(theta_s[7,])

SS <- dim(y)[1]
TT <- dim(y)[2]
period <- 52

P <- array(NA, dim = c(SS, TT-1, 2, 2))
inv_logit <- function(x){return(1/(1+exp(-x)))}
for (time in 1:(TT-1)) {
  P[, time, 1, 2] <- inv_logit(
    beta_0_hat + beta_1_hat* cos(2*pi*time/period) + beta_2_hat * sin(2*pi*time/period) + alpha_hat * time)
  P[, time, 1, 1] <- 1 - P[, time, 1, 2]
  P[, time, 2, 1] <- inv_logit(beta_0s_hat + beta_1s_hat * cos(2*pi*time/period) + beta_2s_hat * sin(2*pi*time/period)+ alpha_s_hat * time)
  P[, time, 2, 2] <- 1 - P[, time, 2, 1]
}
P_nap <- P[1,,,]
my_P <- P_nap
location_time_0 <- which(y[1,-ncol(y)]==0, arr.ind =  TRUE)
next_y <- y[cbind(location_time_0[,1], location_time_0[,2]+1)]

llh <- 0
for (i in 1:length(next_y)) {
  llh <- llh + ifelse(next_y[i]==0, log(P_nap[location_time_0[i,2],1,1]), log(P_nap[location_time_0[i,2],1,2]))
}

location_time_1<- which(y[1,-ncol(y)]==1, arr.ind =  TRUE)
next_y <- y[cbind(location_time_1[,1], location_time_1[,2]+1)]
for (i in 1:length(next_y)) {
  llh <- llh + ifelse(next_y[i]==0, log(P_nap[location_time_1[i,2],2,1]), log(P_nap[location_time_1[i,2],2,2]))
}
my_llh <- llh



beta_0_hat <- -3.323236
beta_1_hat <- -4.204593
beta_2_hat <- 0.886424 
alpha_hat <- -3.464589e-05
beta_0s_hat <- 1.532151 
beta_1s_hat <- 3.76802
beta_2s_hat <- 0.33008
alpha_s_hat <- -0.0004243435

SS <- dim(y)[1]
TT <- dim(y)[2]
period <- 52

P <- array(NA, dim = c(SS, TT-1, 2, 2))
inv_logit <- function(x){return(1/(1+exp(-x)))}
for (time in 1:(TT-1)) {
  P[, time, 1, 2] <- inv_logit(
    beta_0_hat + beta_1_hat* cos(2*pi*time/period) + beta_2_hat * sin(2*pi*time/period) + alpha_hat * time)
  P[, time, 1, 1] <- 1 - P[, time, 1, 2]
  P[, time, 2, 1] <- inv_logit(beta_0s_hat + beta_1s_hat * cos(2*pi*time/period) + beta_2s_hat * sin(2*pi*time/period)+ alpha_s_hat * time)
  P[, time, 2, 2] <- 1 - P[, time, 2, 1]
}
P_nap <- P[1,,,]
location_time_0 <- which(y[1,-ncol(y)]==0, arr.ind =  TRUE)
next_y <- y[cbind(location_time_0[,1], location_time_0[,2]+1)]
yisu_P <- P_nap
llh <- 0
for (i in 1:length(next_y)) {
  llh <- llh + ifelse(next_y[i]==0, log(P_nap[location_time_0[i,2],1,1]), log(P_nap[location_time_0[i,2],1,2]))
}
location_time_1<- which(y[1,-ncol(y)]==1, arr.ind =  TRUE)
next_y <- y[cbind(location_time_1[,1], location_time_1[,2]+1)]
for (i in 1:length(next_y)) {
  llh <- llh + ifelse(next_y[i]==0, log(P_nap[location_time_1[i,2],2,1]), log(P_nap[location_time_1[i,2],2,2]))
}
yisu_llh <- llh

pred <- matrix(NA, nrow = 2, ncol = TT)

curr_p0 <- t(c(y[1,1] == 0, y[1,1] == 1))
pred[,1] <- curr_p0

for(t in 1:(TT-1)){
  curr_p0 <- curr_p0 %*% my_P[t,,]
  pred[,t+1] <- curr_p0
}

# Assume your vector is named vec
vec <- pred[2, ]  # Example data
mat <- matrix(vec, nrow = 52, ncol = 54)
my_esn <- colSums(mat)

curr_p0 <- t(c(y[1,1] == 0, y[1,1] == 1))
pred[,1] <- curr_p0
for(t in 1:(TT-1)){
  curr_p0 <- curr_p0 %*% yisu_P[t,,]
  pred[,t+1] <- curr_p0
}

# Assume your vector is named vec
vec <- pred[2, ]  # Example data
mat <- matrix(vec, nrow = 52, ncol = 54)
my_esn <- colSums(mat)


vec <- as.numeric(y[1,])
mat <- matrix(vec, nrow = 52, ncol = 54)
true_sn <- colSums(mat)



# apply(p01_params[c(1,3,5,7),], 1, mean)
# apply(p10_params[c(1,3,5,7),], 1, mean)
# my_param <- c(apply(p01_params[c(1,3,5,7),], 1, mean), apply(p10_params[c(1,3,5,7),], 1, mean))
# yisu_param <- c(-3.2016, -4.073, 0.791, 0.0000382, 1.7258, 3.669, -0.944, -0.0004935)
# round(cbind(my_param,yisu_param),7)
