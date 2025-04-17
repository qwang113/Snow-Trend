setwd("D:/77/Research/temp/snow_trend")
no_nbs <- c(57, 170, 236, 269, 343, 685, 946, 947, 989, 1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203)
all_y <- readRDS(here::here("snow_cleaned.Rda"))
all_y <- all_y[-no_nbs,]
y <- all_y[,-c(1,2)]
coords <- all_y[,1:2]
S <- nrow(y)
TT <- ncol(y)

theta_01_all <- apply(theta[1:S,],1,mean)
theta_11_all <- apply(theta[(S+1):(2*S),],1,mean)
theta_21_all <- apply(theta[(2*S+1):(3*S),],1,mean)
theta_a1_all <- apply(theta[(3*S+1):(4*S),],1,mean)
theta_01s_all <- apply(theta_s[1:S,],1,mean)
theta_11s_all <- apply(theta_s[(S+1):(2*S),],1,mean)
theta_21s_all <- apply(theta_s[(2*S+1):(3*S),],1,mean)
theta_a1s_all <- apply(theta_s[(3*S+1):(4*S),],1,mean)

my_param <- cbind(theta_01_all,theta_11_all, theta_21_all, theta_a1_all, theta_01s_all, theta_11s_all, theta_21s_all, theta_a1s_all)
yisu_param <- readRDS("llh.rds")[-no_nbs,3:10]

convert_cosine_form <- function(A, B) {
  C <- A * cos(2 * pi * B / 52)
  D <- A * sin(2 * pi * B / 52)
  return(list(C = C, D = D))
}

yisu_param_trans <- convert_cosine_form(A = yisu_param[,2], B = yisu_param[,3])
yisu_param_trans_s <- convert_cosine_form(A = yisu_param[,6], B = yisu_param[,7])
yisu_param[,2] <- yisu_param_trans$C
yisu_param[,3] <- yisu_param_trans$D
yisu_param[,6] <- yisu_param_trans_s$C
yisu_param[,7] <- yisu_param_trans_s$D

inv_logit <- function(x){return(1/(1+exp(-x)))}

SS <- dim(y)[1]
TT <- dim(y)[2]
all_llh <- rep(NA, SS)
period <- 52
  P <- array(NA, dim = c(SS, TT-1, 2, 2))
  for (time in 1:(TT-1)) {
    P[, time, 1, 2] <- inv_logit(
      my_param[s,1] + my_param[,2]* cos(2*pi*time/period) + my_param[,3] * sin(2*pi*time/period) + my_param[,4] * time)
    P[, time, 1, 1] <- 1 - P[, time, 1, 2]
    P[, time, 2, 1] <- inv_logit(my_param[,5] + my_param[,6] * cos(2*pi*time/period) + my_param[,7] * sin(2*pi*time/period)+ my_param[,8] * time)
    P[, time, 2, 2] <- 1 - P[, time, 2, 1]
  }
for (s in 1:SS) {
  print(s)
  P_curr <- P[s,,,]
  location_time_0 <- which(y[s,-ncol(y)]==0, arr.ind =  TRUE)
  next_y <- y[cbind(s, location_time_0[,2]+1)]
  
  llh <- 0
  for (i in 1:length(next_y)) {
    llh <- llh + ifelse(next_y[i]==0, log(P_curr[location_time_0[i,2],1,1]), log(P_curr[location_time_0[i,2],1,2]))
  }
  
  
  location_time_1<- which(y[s,-ncol(y)]==1, arr.ind =  TRUE)
  next_y <- y[cbind(s, location_time_1[,2]+1)]
  for (i in 1:length(next_y)) {
    llh <- llh + ifelse(next_y[i]==0, log(P_curr[location_time_1[i,2],2,1]), log(P_curr[location_time_1[i,2],2,2]))
  }
  all_llh[s] <- llh
  
}

my_llh <- all_llh



my_param <- yisu_param
P <- array(NA, dim = c(SS, TT-1, 2, 2))
for (time in 1:(TT-1)) {
  P[, time, 1, 2] <- inv_logit(
    my_param[s,1] + my_param[,2]* cos(2*pi*time/period) + my_param[,3] * sin(2*pi*time/period) + my_param[,4] * time)
  P[, time, 1, 1] <- 1 - P[, time, 1, 2]
  P[, time, 2, 1] <- inv_logit(my_param[,5] + my_param[,6] * cos(2*pi*time/period) + my_param[,7] * sin(2*pi*time/period)+ my_param[,8] * time)
  P[, time, 2, 2] <- 1 - P[, time, 2, 1]
}
for (s in 1:SS) {
  print(s)
  P_curr <- P[s,,,]
  location_time_0 <- which(y[s,-ncol(y)]==0, arr.ind =  TRUE)
  next_y <- y[cbind(s, location_time_0[,2]+1)]
  
  llh <- 0
  for (i in 1:length(next_y)) {
    llh <- llh + ifelse(next_y[i]==0, log(P_curr[location_time_0[i,2],1,1]), log(P_curr[location_time_0[i,2],1,2]))
  }
  
  
  location_time_1<- which(y[s,-ncol(y)]==1, arr.ind =  TRUE)
  next_y <- y[cbind(s, location_time_1[,2]+1)]
  for (i in 1:length(next_y)) {
    llh <- llh + ifelse(next_y[i]==0, log(P_curr[location_time_1[i,2],2,1]), log(P_curr[location_time_1[i,2],2,2]))
  }
  all_llh[s] <- llh
  
}

yisu_llh <- all_llh


library(ggplot2)

# Create a data frame for ggplot2
df <- data.frame(my_llh = my_llh, yisu_llh = yisu_llh)

# Plot using ggplot2
ggplot(df, aes(x = my_llh, y = yisu_llh)) +
  geom_point() +  # Scatter plot
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1, lty = 2) +  # Reference line y = x
  theme_minimal() +  # Use a clean theme
  labs(x = "my_llh", y = "yisu_llh", title = "Comparison Plot")  # Labels and title



