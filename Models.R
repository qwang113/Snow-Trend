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

snow_dat <- readRDS("snow_cleaned.Rda")

snow_matrix <- as.matrix(snow_dat[,-(1:2)])

snow_long <- snow_dat$LON
snow_lat <- snow_dat$LAT
snow_long_miles <- snow_dat$LON*54.6
snow_lat_miles <- snow_dat$LAT*64
Distances <- pairdist(cbind(snow_long_miles, snow_lat_miles))
Omg <- ifelse(Distances > 200,0,1) - diag(nrow = nrow(snow_dat))
D <- diag(rowSums(Omg)+1)

# Hyperparameter for tau's prior
sigma <- 1.0  # Set this to your desired hyperparameter value

# Create a list for Stan


stan_code <- {"
  functions {
    // Define the matrix scale function
    
    matrix scale(matrix A) {
        int S = rows(A); // Get the number of rows
        vector[S] row_sums; // Vector to store the row sums
        matrix[S, cols(A)] scaled_A; // Initialize the scaled matrix

        for (s in 1:S) {
         scaled_A[s, ] = A[s, ]/A[s,s]
        }

        return scaled_A; // Return the scaled matrix
    }
    
    
    // Define Likelihood
    real llh(matrix y,int S, int TT,
    vector beta_0, vector beta_1, vector beta_2, vector alpha,
    vector beta_0s, vector beta_1s, vector beta_2s, vector alpha_s) {
    
            real ll = 0; // Initialize log-likelihood
        
        for (s in 1:S) {
            for (t in 1:(TT-1)) {
                if(y[s,t]==0) {
                real p01 = inv_logit(beta_0[s] + beta_1[s]*cos(2*pi()*t/52) + beta_2[s]*sin(2*pi()*t/52) + alpha[s]*t);
                 ll += (1-y[s,t+1])*log(1-p01) + y[s,t+1]*log(p01);
                } else{
                real p10 = inv_logit(beta_0s[s] + beta_1s[s]*cos(2*pi()*t/52) + beta_2s[s]*sin(2*pi()*t/52) + alpha_s[s]*t);
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
    real<lower=0> tau_0;
    real<lower=0> tau_0s;
    real<lower=0> tau_1;
    real<lower=0> tau_1s;
    real<lower=0> tau_2;
    real<lower=0> tau_2s;
    real<lower=0> tau_a;
    real<lower=0> tau_as;
    
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
    theta_01 ~ multi_normal_prec( rep_vector(0,S), 1/tau_0 * scale(D-Omg));
    theta_02 ~ multi_normal_prec( rep_vector(0,S), 1/tau_0 * diag_matrix(rep_vector(1, S)));
    theta_01s ~ multi_normal_prec( rep_vector(0,S), 1/tau_0s * scale(D-Omg));
    theta_02s ~ multi_normal_prec( rep_vector(0,S), 1/tau_0s * diag_matrix(rep_vector(1, S)));  
    
    theta_11 ~ multi_normal_prec( rep_vector(0,S), 1/tau_1 * scale(D-Omg));
    theta_12 ~ multi_normal_prec( rep_vector(0,S), 1/tau_1 * diag_matrix(rep_vector(1, S)));
    theta_11s ~ multi_normal_prec( rep_vector(0,S), 1/tau_1s * scale(D-Omg));
    theta_12s ~ multi_normal_prec( rep_vector(0,S), 1/tau_1s * diag_matrix(rep_vector(1, S))); 
    
    theta_21 ~ multi_normal_prec( rep_vector(0,S), 1/tau_2 * scale(D-Omg));
    theta_22 ~ multi_normal_prec( rep_vector(0,S), 1/tau_2 * diag_matrix(rep_vector(1, S)));
    theta_21s ~ multi_normal_prec( rep_vector(0,S), 1/tau_2s * scale(D-Omg));
    theta_22s ~ multi_normal_prec( rep_vector(0,S), 1/tau_2s * diag_matrix(rep_vector(1, S))); 
    
    theta_a1 ~ multi_normal_prec( rep_vector(0,S), 1/tau_a * scale(D-Omg));
    theta_a2 ~ multi_normal_prec( rep_vector(0,S), 1/tau_a * diag_matrix(rep_vector(1, S)));
    theta_a1s ~ multi_normal_prec( rep_vector(0,S), 1/tau_as * scale(D-Omg) );
    theta_a2s ~ multi_normal_prec( rep_vector(0,S), 1/tau_as * diag_matrix(rep_vector(1, S))); 
    
    // Define tau hierarchical priors
    tau_0 ~ cauchy(0,sigma);
    tau_0s ~ cauchy(0,sigma);
    tau_1 ~ cauchy(0,sigma);
    tau_1s ~ cauchy(0,sigma);
    tau_2 ~ cauchy(0,sigma);
    tau_2s ~ cauchy(0,sigma);
    tau_0 ~ cauchy(0,sigma);
    tau_0s ~ cauchy(0,sigma);
    
    
    // Target function
    
    target += llh(y, S, TT, beta_0, beta_0s, beta_1, beta_1s, beta_2, beta_2s, alpha, alpha_s);
}
  "}
data_list <- list(TT = dim(snow_matrix)[2], S = dim(snow_matrix)[1], y = snow_matrix, D = D, Omg = Omg, sigma = sigma)
# stan_model_code <- stan_model(model_code = stan_code)

fit <- stan(model_code = stan_code, 
                   data = data_list,
                   chains = 1,             # Number of chains
                   iter = 2000,            # Total iterations per chain
                   warmup = 1000,          # Number of warmup iterations
                   thin = 2,               # Thinning interval
                   cores = 2,
                   init = "0")             # Number of cores to use



# # Check the quality of neighborhood choice
# no_nbs <- which(rowSums(D)==1)
# nb_check <- data.frame("long" = snow_long, "lat" = snow_lat, "nonbs" = 0)
# nb_check$nonbs[no_nbs] <- 1
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]
# 
# 
# # Convert the data frame to an sf object
# sf_data <- st_as_sf(nb_check, coords = c("long", "lat"), crs = 4326)
# 
# aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"
# world_aeqd <- st_transform(world_north, crs = aeqd_proj)
# sf_data_aeqd <- st_transform(sf_data, crs = aeqd_proj)
# equator_points <- data.frame(
#   lon = seq(-180, 180, length.out = 200),
#   lat = rep(0, 100)  # All points at latitude = 0
# )
# 
# # Transform equator points to the azimuthal equidistant projection
# equator_sf <- st_as_sf(equator_points, coords = c("lon", "lat"), crs = 4326)
# equator_aeqd <- st_transform(equator_sf, crs = aeqd_proj)
# 
# 
# idx = 1
# # Plot for a single time
# ggplot() +
#   geom_sf(data = world_aeqd, fill = "lightgray", color = "NA") + # World map
#   geom_sf(data = sf_data_aeqd, aes(color = factor(nonbs)), size = 1.5, shape = 18) + # Data points with color mapped
#   geom_sf(data = world_aeqd, fill = "NA", color = "black") + # World map
#   geom_sf(data = equator_aeqd, color = "red", linetype = "dashed", size = 0.1)+
#   theme_minimal() +
#   labs(title = paste("Snow Existence Map from the North Pole","-",substr(names(sf_data_aeqd)[idx],start = 2, stop = 50) ),
#        color = "Snow Existence") +
#   theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))