rm(list = ls())

library(sf)
library(ggplot2)
library(rnaturalearth)
library(cowplot)
library(reticulate)
library(knitr)

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
period <- 52

# ------------------------------------------------------------
# 1️⃣ Load FULL snow
# ------------------------------------------------------------
snow_full <- readRDS("snow_cleaned.Rda")

coords <- as.matrix(snow_full[, 1:2])
y <- as.matrix(snow_full[, -c(1, 2)])

S <- nrow(coords)
cat("FULL S =", S, "\n")

# ------------------------------------------------------------
# 2️⃣ Covariates (FULL, aligned with full-data fitting)
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 3️⃣ Load NPZ summaries (FULL DATA)
# ------------------------------------------------------------
use_condaenv("CPD", required = TRUE)
py_config()

setwd("D:/77/Research/temp/snow")
np <- import("numpy", convert = FALSE)

ind  <- np$load("trend_ind_summary.npz", allow_pickle = TRUE)
bym  <- np$load("trend_bym_summary.npz", allow_pickle = TRUE)
bymp <- np$load("trend_bymp_summary.npz", allow_pickle = TRUE)

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

# ------------------------------------------------------------
# 4️⃣ Map prep
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 5️⃣ Mean map
# ------------------------------------------------------------
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

p1 <- plot_mean(sf_INDEP, "Aggregated Trend - IND")
p2 <- plot_mean(sf_BM, "Aggregated Trend - Weekly BYM")
p3 <- plot_mean(sf_BMP, "Aggregated Trend - Weekly BYM+")

cowplot::plot_grid(p1, p2, p3, nrow = 1)

# ------------------------------------------------------------
# 6️⃣ SD map
# ------------------------------------------------------------
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

p4 <- plot_sd(sf_INDEP, "Aggregated Trend log(SD) - IND")
p5 <- plot_sd(sf_BM, "Aggregated Trend log(SD) - Weekly BYM")
p6 <- plot_sd(sf_BMP, "Aggregated Trend log(SD) - Weekly BYM+")

cowplot::plot_grid(p4, p5, p6, nrow = 1)
cowplot::plot_grid(p3, p6, nrow = 1)

# ------------------------------------------------------------
# 7️⃣ Table
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 8️⃣ Regression
# ------------------------------------------------------------
M1 <- lm(trend_BM$mean ~ lats + elev)
knitr::kable(summary(M1)$coefficients, format = "latex", digits = 5)

