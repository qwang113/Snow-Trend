rm(list = ls())

library(sf)
library(igraph)
library(ggplot2)
library(rnaturalearth)
library(cowplot)
library(reticulate)
library(knitr)

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
DIST_TH <- 0.22
period <- 52

# ------------------------------------------------------------
# 1️⃣ Load snow + remove isolated
# ------------------------------------------------------------
no_nbs <- c(57, 170, 236, 269, 343, 685, 946, 947, 989,
            1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203)

snow_full <- readRDS("snow_cleaned.Rda")
snow_full <- snow_full[-no_nbs,]

coords_full <- as.matrix(snow_full[,1:2])
y_full <- as.matrix(snow_full[,-c(1,2)])

# ------------------------------------------------------------
# 2️⃣ Keep two largest connected components
# ------------------------------------------------------------
coords_sf <- st_as_sf(
  data.frame(LON=coords_full[,1], LAT=coords_full[,2]),
  coords=c("LON","LAT"),
  crs=4326
)

coords_aeqd <- st_transform(coords_sf,
                            "+proj=aeqd +lat_0=90 +lon_0=-100")

xy <- st_coordinates(coords_aeqd) / 1e6

Dmat <- as.matrix(dist(xy))
W <- (Dmat <= DIST_TH) * 1
diag(W) <- 0

g <- graph_from_adjacency_matrix(W, mode="undirected")
comp <- components(g)

sizes <- comp$csize
order_idx <- order(sizes, decreasing=TRUE)
keep_idx <- sort(which(comp$membership %in% order_idx[1:2]))

coords <- coords_full[keep_idx,]
y <- y_full[keep_idx,]

S <- nrow(coords)
cat("Final S =", S, "\n")

# ------------------------------------------------------------
# 3️⃣ Covariates (aligned)
# ------------------------------------------------------------
elev_all <- read.csv("curr_elev.csv")[,4]
elev <- scale(elev_all[-no_nbs][keep_idx])

lats <- scale(coords[,2])

# ------------------------------------------------------------
# 4️⃣ Load NPZ (no Rda)
# ------------------------------------------------------------
library(reticulate)

use_python("e:/anaconda3/envs/CPD/python.exe", required = TRUE)

py_config()

np <- import("numpy")
setwd("D:/77/Research/temp/snow")
ind_data  <- np$load("trend_ind.npz", allow_pickle=TRUE)
bym_data  <- np$load("trend_bym_weekly.npz", allow_pickle=TRUE)
bymp_data <- np$load("trend_weekly_bym+cov.npz", allow_pickle=TRUE)

weekly_ini_INDEP   <- ind_data$f[["weekly_ini"]]
weekly_final_INDEP <- ind_data$f[["weekly_final"]]

weekly_ini_BM   <- bym_data$f[["weekly_ini"]]
weekly_final_BM <- bym_data$f[["weekly_final"]]

weekly_ini_BMP   <- bymp_data$f[["weekly_ini"]]
weekly_final_BMP <- bymp_data$f[["weekly_final"]]

# ------------------------------------------------------------
# 5️⃣ Trend computation
# ------------------------------------------------------------
compute_trend <- function(ini, final) {
  
  dpd_snow_diff <-
    (apply(final[,,,2], c(1,2), sum) -
       apply(ini[,,,2], c(1,2), sum)) / 51
  
  list(
    mean  = apply(dpd_snow_diff, 2, mean),
    sd    = apply(dpd_snow_diff, 2, sd),
    lower = apply(dpd_snow_diff, 2, quantile, 0.025),
    upper = apply(dpd_snow_diff, 2, quantile, 0.975)
  )
}

trend_INDEP <- compute_trend(weekly_ini_INDEP, weekly_final_INDEP)
trend_BM    <- compute_trend(weekly_ini_BM, weekly_final_BM)
trend_BMP   <- compute_trend(weekly_ini_BMP, weekly_final_BMP)

# ------------------------------------------------------------
# 6️⃣ Map prep
# ------------------------------------------------------------
aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"

make_sf <- function(mean, sd) {
  st_transform(
    st_as_sf(
      data.frame(LON=coords[,1], LAT=coords[,2],
                 mean=mean, sd=sd),
      coords=c("LON","LAT"),
      crs=4326
    ),
    crs=aeqd_proj
  )
}

sf_INDEP <- make_sf(trend_INDEP$mean, trend_INDEP$sd)
sf_BM    <- make_sf(trend_BM$mean, trend_BM$sd)
sf_BMP   <- make_sf(trend_BMP$mean, trend_BMP$sd)

world <- ne_countries(scale="medium", returnclass="sf")
world_north <- world[st_coordinates(st_centroid(world))[,2] > 0,]
world_aeqd <- st_transform(world_north, crs=aeqd_proj)

rg_mean <- range(trend_INDEP$mean,
                 trend_BM$mean,
                 trend_BMP$mean)

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
    labs(title = title, color = "") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
}

p1 <- plot_mean(sf_INDEP, "Aggregated Trend - IND")
p2 <- plot_mean(sf_BM, "Aggregated Trend - Weekly BYM")
p3 <- plot_mean(sf_BMP, "Aggregated Trend - Weekly BYM+")

cowplot::plot_grid(p1, p2, p3, nrow=1)


# ------------------------------------------------------------
#  SD MAP (log scale)
# ------------------------------------------------------------

rg_sd <- range(
  log(trend_INDEP$sd),
  log(trend_BM$sd),
  log(trend_BMP$sd)
)

plot_sd <- function(sf_obj, title) {
  ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
    geom_sf(data = sf_obj,
            aes(color = log(sd)),
            size = 2, shape = 18) +
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




make_table_row <- function(trend) {
  
  inc_total  <- sum(trend$mean > 0)
  dec_total  <- sum(trend$mean < 0)
  
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

# LaTeX table (with vertical lines)
kable(
  tb,
  format = "latex",
  booktabs = FALSE,
  align = "c",
  col.names = colnames(tb)
)

colnames(tb) <- c("Decrease","Increase","Mean","Median","SD")
rownames(tb) <- c("IND","BYM","BYM+")

knitr::kable(tb, format="latex", digits=4)

# ------------------------------------------------------------
# 8️⃣ Regression
# ------------------------------------------------------------
M1 <- lm(trend_BM$mean ~ lats + elev)
knitr::kable(summary(M1)$coefficients,
             format="latex", digits=5)
