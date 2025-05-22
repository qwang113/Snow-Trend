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

# all_y <- readRDS("snow_cleaned.Rda")[-no_nbs,]
# all_y_nnbs <- readRDS("snow_cleaned.Rda")[no_nbs,]

all_y <- readRDS("snow_cleaned_full.Rda")[-no_nbs,]
all_y_nnbs <- readRDS("snow_cleaned_full.Rda")[no_nbs,]

y <- rbind(all_y, all_y_nnbs)[,-c(1,2)]
coords <- rbind(all_y, all_y_nnbs)[,1:2]

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
chunk_size <- 100
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

setwd("D:/77/Research/temp/snow/")
saveRDS(design_mat,"design01_IND.Rda")
design_mat <- readRDS("design01_IND.Rda")

tot_samples <- 2000

all_theta <- matrix(NA, nrow = 4*S, ncol = tot_samples)
all_tau <- matrix(NA, nrow = 4, ncol = tot_samples)
curr_theta_vec <- matrix(0, nrow = 1, ncol = 4*S)
curr_tau_vec <- rep(1,4)
a_tau <- 0.001
b_tau <- 0.001
curr_idx <- 0
save_idx <- 0
burn = 0
thin = 1

while(save_idx < tot_samples) {
  curr_idx = curr_idx + 1
  print(curr_idx)
  # Sample current Omega
  kappas <- next_y - 1/2
  curr_phi <- design_mat %*% t(curr_theta_vec)
  curr_omega <-  rpg(length(next_y), h = 1, z = as.numeric(curr_phi))
  # Sample current theta
  curr_prec <- bdiag(
    1/25*diag(1,S),
    1/25*diag(1,S),
    1/25*diag(1,S),
    1/1*diag(1,S)
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
  if((curr_idx > burn) & (curr_idx %% thin == 0) ){
    save_idx <- save_idx + 1
    all_theta[,save_idx] <- as.vector(curr_theta_vec)
  }
}

saveRDS(all_theta, "theta01_IND.Rda")
