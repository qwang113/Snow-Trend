rm(list = ls())

library(sf)
library(ggplot2)
library(rnaturalearth)
library(cowplot)
library(reticulate)
library(knitr)
# Paths and settings
period <- 52
# 1. Load FULL snow
snow_full <- readRDS("snow_cleaned.Rda")

coords <- as.matrix(snow_full[, 1:2])
y <- as.matrix(snow_full[, -c(1, 2)])

S <- nrow(coords)
cat("FULL S =", S, "\n")
# 2. Covariates (FULL, aligned with full-data fitting)
no_nbs <- c(
  57, 170, 236, 269, 343, 685, 946, 947, 989,
  1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203
)

# latitude
lats <- as.numeric(scale(coords[, 2]))

# elevation: curr_elev + nnbs_elev
elev_raw <- read.csv("curr_elev.csv")[, 4]
nnbs_elev <- read.delim("nnbs_elev.csv", sep = "\t", row.names = NULL)[, 4]

mask <- rep(TRUE, S)
mask[no_nbs] <- FALSE

elev_all <- numeric(S)
elev_all[mask] <- elev_raw
elev_all[no_nbs] <- nnbs_elev

elev <- as.numeric(scale(elev_all))
# 3. NPZ summaries
use_condaenv("CPD", required = TRUE)
py_config()

# Set this path to the local data and results directory.
setwd("path/to/snow/data-and-results")
np <- import("numpy", convert = FALSE)

ind  <- np$load("trend_ind_summary.npz", allow_pickle = TRUE)
bym  <- np$load("trend_bym_summary.npz", allow_pickle = TRUE)
bymp <- np$load("trend_bymp+lon_summary.npz", allow_pickle = TRUE)

to_r_vec <- function(x) {
  as.numeric(py_to_r(x$tolist()))
}

trend_INDEP <- list(
  mean  = to_r_vec(ind$f[["mean"]]),
  sd    = to_r_vec(ind$f[["sd"]]),
  lower = to_r_vec(ind$f[["lower"]]),
  upper = to_r_vec(ind$f[["upper"]])
)

trend_BM <- list(
  mean  = to_r_vec(bym$f[["mean"]]),
  sd    = to_r_vec(bym$f[["sd"]]),
  lower = to_r_vec(bym$f[["lower"]]),
  upper = to_r_vec(bym$f[["upper"]])
)

trend_BMP <- list(
  mean  = to_r_vec(bymp$f[["mean"]]),
  sd    = to_r_vec(bymp$f[["sd"]]),
  lower = to_r_vec(bymp$f[["lower"]]),
  upper = to_r_vec(bymp$f[["upper"]])
)

stopifnot(length(trend_INDEP$mean) == S)
stopifnot(length(trend_BM$mean) == S)
stopifnot(length(trend_BMP$mean) == S)
# 4. Map prep
aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"

make_sf <- function(mean, sd) {
  st_transform(
    st_as_sf(
      data.frame(
        LON = coords[, 1],
        LAT = coords[, 2],
        mean = mean,
        sd = sd
      ),
      coords = c("LON", "LAT"),
      crs = 4326
    ),
    crs = aeqd_proj
  )
}

sf_INDEP <- make_sf(trend_INDEP$mean, trend_INDEP$sd)
sf_BM    <- make_sf(trend_BM$mean, trend_BM$sd)
sf_BMP   <- make_sf(trend_BMP$mean, trend_BMP$sd)

world <- ne_countries(scale = "medium", returnclass = "sf")
world_north <- world[st_coordinates(st_centroid(world))[, 2] > 0, ]
world_aeqd <- st_transform(world_north, crs = aeqd_proj)
# 5. Mean map
rg_mean <- range(
  trend_INDEP$mean,
  trend_BM$mean,
  trend_BMP$mean
)

plot_mean <- function(sf_obj, title) {
  ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
    geom_sf(data = sf_obj, aes(color = mean), size = 2, shape = 18) +
    geom_sf(data = world_aeqd, fill = NA, color = "black") +
    scale_color_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      limits = rg_mean,
      guide = guide_colorbar(
        barwidth = 20,
        barheight = 0.5
      )
    ) +
    theme_minimal() +
    labs(title = title, color = "Mean") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
}

p1 <- plot_mean(sf_INDEP, "Annual Trend - IND")
p2 <- plot_mean(sf_BM, "Annual Trend - Weekly BYM")
p3 <- plot_mean(sf_BMP, "Annual Trend - BYM+")

cowplot::plot_grid(p1, p2, p3, nrow = 1)
# 6. SD map
rg_sd <- range(
  log(trend_INDEP$sd),
  log(trend_BM$sd),
  log(trend_BMP$sd)
)

plot_sd <- function(sf_obj, title) {
  ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
    geom_sf(data = sf_obj, aes(color = log(sd)), size = 2, shape = 18) +
    geom_sf(data = world_aeqd, fill = NA, color = "black") +
    scale_color_viridis_c(
      option = "C",
      direction = 1,
      limits = rg_sd,
      guide = guide_colorbar(
        barwidth = 20,
        barheight = 0.5
      )
    ) +
    theme_minimal() +
    labs(title = title, color = "log(SD)") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
}

p4 <- plot_sd(sf_INDEP, "Annual Trend log(SD) - IND")
p5 <- plot_sd(sf_BM, "Annual Trend log(SD) - Weekly BYM")
p6 <- plot_sd(sf_BMP, "Annual Trend log(SD) - BYM+")

cowplot::plot_grid(p4, p5, p6, nrow = 1)
cowplot::plot_grid(p3, p6, nrow = 1)
# 7. Table
make_table_row <- function(trend) {
  inc_total <- sum(trend$mean > 0)
  dec_total <- sum(trend$mean < 0)
  
  inc_sig <- sum(trend$lower > 0)
  dec_sig <- sum(trend$upper < 0)
  
  c(
    paste0(inc_total, "(*", inc_sig, ")"),
    paste0(dec_total, "(*", dec_sig, ")"),
    sprintf("%.4f", mean(trend$mean)),
    sprintf("%.4f", median(trend$mean)),
    sprintf("%.4f", mean(trend$sd))
  )
}

tb <- rbind(
  make_table_row(trend_INDEP),
  make_table_row(trend_BM),
  make_table_row(trend_BMP)
)

colnames(tb) <- c("Increase", "Decrease", "Mean", "Median", "SD")
rownames(tb) <- c("IND", "BYM", "BYM+")

kable(
  tb,
  format = "latex",
  booktabs = FALSE,
  align = "c",
  col.names = colnames(tb)
)
# 8. Regression
M1 <- lm(trend_BM$mean ~ lats + elev)
knitr::kable(summary(M1)$coefficients, format = "latex", digits = 5)


df_scatter <- data.frame(
  lat = lats,
  trend_ind  = trend_INDEP$mean,
  trend_bym  = trend_BM$mean,
  trend_bymp = trend_BMP$mean
)

p_lat_bymp <- ggplot(df_scatter, aes(x = lat, y = trend_bymp)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Trend vs Latitude (BYM+)",
    x = "Latitude (scaled)",
    y = "Trend Estimate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_lat_bymp
# Scatter: Trend vs Elevation
df_scatter$elev <- elev
df_scatter$elev_raw <- scale(elev_all)

p_elev <- ggplot(df_scatter, aes(x = elev_raw, y = trend_bymp)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Trend vs Elevation (BYM+)",
    x = "Elevation (scaled)",
    y = "Trend Estimate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_elev

load("snow_temp_full.Rda")
snow_temp <- sce_temp
temp_full <- as.matrix(snow_temp[,-c(1,2)])
temp_scaled <- scale(temp_full)

period <- 52
TT <- ncol(temp_full)

w_t <- ((1:TT - 1) %% period) + 1
t_scaled <- as.numeric(scale(1:TT))

S <- nrow(temp_full)

alpha <- rep(NA, S)

for (s in 1:S) {
  
  y <- temp_full[s, ]
  
  if (any(is.na(y))) next
  
  df <- data.frame(
    y = y,
    t = t_scaled,
    w = factor(w_t)
  )
  
  fit <- lm(y ~ t + w, data = df)
  
  alpha[s] <- coef(fit)["t"]
}


df_temp <- data.frame(
  trend = trend_BMP$mean,
  alpha = alpha
)

df_temp <- df_temp[complete.cases(df_temp), ]

p_temp <- ggplot(df_temp, aes(x = alpha, y = trend)) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Trend vs Temperature Trend (BYM+)",
    x = "Temperature Trend",
    y = "Trend Estimate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_temp
# Scatter: Trend vs Longitude
lons <- as.numeric(scale(coords[, 1]))

df_scatter$lon <- lons

p_lon <- ggplot(df_scatter, aes(x = lon, y = trend_bymp)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Trend vs Longitude (BYM+)",
    x = "Longitude (scaled)",
    y = "Trend Estimate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_lon


df_scatter$group <- ifelse(df_scatter$lon < -0.6, "North America", "Asia and Europe")

p_lon2 <- ggplot(df_scatter, aes(x = lon, y = trend_bymp, color = group)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Trend vs Longitude (BYM+)",
    x = "Longitude (scaled)",
    y = "Trend Estimate",
    color = "Area"
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p_lon2
