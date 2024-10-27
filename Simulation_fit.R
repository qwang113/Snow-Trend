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

all_y <- readRDS(here::here("sim_y.Rda"))
y <- all_y[,-c(1,2)]

coords <- all_y[,1:2]
Distances <- pairdist(coords)
Omg <- matrix(0, nrow = nrow(coords), ncol = nrow(coords))
Omg[which(Distances <= 0.1)] = 1
Omg <- Omg - diag(nrow = nrow(coords))
# D <- diag(rowSums(Omg))
D <- diag(rowSums(Omg)) 
sigma = 1
eps = 0.0001
period = 20

stan_code <- {"
  functions {
    // Define the matrix scale function
    
    matrix scale(matrix A) {
        int S = rows(A); // Get the number of rows
        vector[S] row_sums; // Vector to store the row sums
        matrix[S, cols(A)] scaled_A; // Initialize the scaled matrix

        for (s in 1:S) {
         scaled_A[s, ] = A[s, ]/A[s,s];
        };

        return scaled_A; // Return the scaled matrix
    }
    
    
    // Define Likelihood
    real llh(matrix y,int S, int TT, real period,
    vector beta_0, vector beta_1, vector beta_2, vector alpha,
    vector beta_0s, vector beta_1s, vector beta_2s, vector alpha_s) {
    
            real ll = 0; // Initialize log-likelihood
        
        for (s in 1:S) {
            for (t in 1:(TT-1)) {
                if(y[s,t]==0) {
                real p01 = inv_logit(beta_0[s] + beta_1[s]*cos(2*pi()*t/period) + beta_2[s]*sin(2*pi()*t/period) + alpha[s]*t);
                 ll += (1-y[s,t+1])*log(1-p01) + y[s,t+1]*log(p01);
                } else{
                real p10 = inv_logit(beta_0s[s] + beta_1s[s]*cos(2*pi()*t/period) + beta_2s[s]*sin(2*pi()*t/period) + alpha_s[s]*t);
                 ll += (1-y[s,t+1])*log(p10) + y[s,t+1]*log(1-p10);
                }
            }
        }
        return ll;
    }
  }

  data {
    int<lower=0> TT;          // Number of observations
    int<lower=0> S;          // Number of Locations
    matrix[S, TT] y;     // Data: observations by locations
    matrix<lower = 0>[S,S] D; // Diagonal matrix with number of neighbors
    matrix<lower = 0>[S,S] Omg; // Adjacency matrix Omega
    real<lower = 0> sigma; // Hyperparameters for tau's prior
    real<lower = 0> eps; // A small positive number to make sure BYM prior to be positive definite
    real<lower = 0> period;
}




parameters {
  // Define BYM model parameters for intercept
    vector[S] theta_01;
    vector[S] theta_01s;
    vector[S] theta_02;
    vector[S] theta_02s;
  // Define BYM model parameters for cosine
    vector[S] theta_11;
    vector[S] theta_11s;
    vector[S] theta_12;
    vector[S] theta_12s;
  
  // Define BYM model parameters for sine
    vector[S] theta_21;
    vector[S] theta_21s;
    vector[S] theta_22;
    vector[S] theta_22s;
    
  // Fixed time trend coefficient
    vector[S] theta_a1;
    vector[S] theta_a1s;
    vector[S] theta_a2;
    vector[S] theta_a2s;
    
  // Hyperparameter for hierarchical prior
    real<lower=0> tau_01;
    real<lower=0> tau_02;
    real<lower=0> tau_01s;
    real<lower=0> tau_02s;
    real<lower=0> tau_11;
    real<lower=0> tau_12;
    real<lower=0> tau_11s;
    real<lower=0> tau_12s;
    real<lower=0> tau_21;
    real<lower=0> tau_22;
    real<lower=0> tau_21s;
    real<lower=0> tau_22s;
    real<lower=0> tau_a1;
    real<lower=0> tau_a2;
    real<lower=0> tau_a1s;
    real<lower=0> tau_a2s;
    
}


transformed parameters {
    // Define regression coefficients 
    vector[S] beta_0 = theta_01 + theta_02;
    vector[S] beta_0s = theta_01s + theta_02s;
    vector[S] beta_1 = theta_11 + theta_12;
    vector[S] beta_1s = theta_11s + theta_12s;
    vector[S] beta_2 = theta_21 + theta_22;
    vector[S] beta_2s = theta_21s + theta_22s;
    vector[S] alpha = theta_a1 + theta_a2;
    vector[S] alpha_s = theta_a1s + theta_a2s;
}


model {
    
    // Define BYM normal priors
    theta_01 ~ multi_normal_prec( rep_vector(0,S), 1/tau_01 * (D-Omg+diag_matrix(rep_vector(eps, S))));
    theta_02 ~ multi_normal_prec( rep_vector(0,S), 1/tau_02 * diag_matrix(rep_vector(1, S)));
    theta_01s ~ multi_normal_prec( rep_vector(0,S), 1/tau_01s * (D-Omg+diag_matrix(rep_vector(eps, S))));
    theta_02s ~ multi_normal_prec( rep_vector(0,S), 1/tau_02s * diag_matrix(rep_vector(1, S)));  
    
    theta_11 ~ multi_normal_prec( rep_vector(0,S), 1/tau_11 * (D-Omg+diag_matrix(rep_vector(eps, S))));
    theta_12 ~ multi_normal_prec( rep_vector(0,S), 1/tau_12 * diag_matrix(rep_vector(1, S)));
    theta_11s ~ multi_normal_prec( rep_vector(0,S), 1/tau_11s * (D-Omg+diag_matrix(rep_vector(eps, S))));
    theta_12s ~ multi_normal_prec( rep_vector(0,S), 1/tau_12s * diag_matrix(rep_vector(1, S))); 
    
    theta_21 ~ multi_normal_prec( rep_vector(0,S), 1/tau_21 * (D-Omg+diag_matrix(rep_vector(eps, S))));
    theta_22 ~ multi_normal_prec( rep_vector(0,S), 1/tau_22 * diag_matrix(rep_vector(1, S)));
    theta_21s ~ multi_normal_prec( rep_vector(0,S), 1/tau_21s * (D-Omg+diag_matrix(rep_vector(eps, S))));
    theta_22s ~ multi_normal_prec( rep_vector(0,S), 1/tau_22s * diag_matrix(rep_vector(1, S))); 
    
    theta_a1 ~ multi_normal_prec( rep_vector(0,S), 1/tau_a1 * (D-Omg+diag_matrix(rep_vector(eps, S))));
    theta_a2 ~ multi_normal_prec( rep_vector(0,S), 1/tau_a2 * diag_matrix(rep_vector(1, S)));
    theta_a1s ~ multi_normal_prec( rep_vector(0,S), 1/tau_a1s * (D-Omg+diag_matrix(rep_vector(eps, S))));
    theta_a2s ~ multi_normal_prec( rep_vector(0,S), 1/tau_a2s * diag_matrix(rep_vector(1, S))); 
    
    // Define tau hierarchical priors
    tau_01 ~ cauchy(0,sigma);
    tau_01s ~ cauchy(0,sigma);
    tau_11 ~ cauchy(0,sigma);
    tau_11s ~ cauchy(0,sigma);
    tau_21 ~ cauchy(0,sigma);
    tau_21s ~ cauchy(0,sigma);
    tau_01 ~ cauchy(0,sigma);
    tau_01s ~ cauchy(0,sigma);
    tau_02 ~ cauchy(0,sigma);
    tau_02s ~ cauchy(0,sigma);
    tau_12 ~ cauchy(0,sigma);
    tau_12s ~ cauchy(0,sigma);
    tau_22 ~ cauchy(0,sigma);
    tau_22s ~ cauchy(0,sigma);
    tau_02 ~ cauchy(0,sigma);
    tau_02s ~ cauchy(0,sigma);
    tau_a1 ~ cauchy(0,sigma);
    tau_a2 ~ cauchy(0,sigma);
    tau_a1s ~ cauchy(0,sigma);
    tau_a2s ~ cauchy(0,sigma);
    
    
    // Target function
    
    target += llh(y, S, TT, period, beta_0, beta_0s, beta_1, beta_1s, beta_2, beta_2s, alpha, alpha_s);
}
  "}
data_list <- list(TT = dim(y)[2], S = dim(y)[1], y = y, D = D, Omg = Omg, sigma = sigma, eps = eps, period = period)
# stan_model_code <- stan_model(model_code = stan_code)

fit <- stan(model_code = stan_code, 
            data = data_list,
            chains = 6,             # Number of chains
            iter = 3000,            # Total iterations per chain
            warmup = 1000,          # Number of warmup iterations
            thin = 2,               # Thinning interval
            cores = 6,
            init = "0")             # Number of cores to use
samples <- extract(fit)
saveRDS(samples, "samples.Rda")
