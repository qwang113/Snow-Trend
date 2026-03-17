rm(list = ls())

library(sf)
library(igraph)
library(ggplot2)
library(rnaturalearth)
library(magick)
library(reticulate)

# ============================================================
# 0️⃣  Use correct Python (CPD env with numpy)
# ============================================================
use_condaenv("CPD", required = TRUE)
np <- import("numpy")

# ============================================================
# CONFIG
# ============================================================
DIST_TH <- 0.22
aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"

no_nbs <- c(57, 170, 236, 269, 343, 685, 946, 947, 989,
            1037, 1084, 1090, 1109, 1118, 1127, 1176, 1203)

# ============================================================
# 1️⃣ Load snow + reproduce Python filtering
#     (two largest connected components)
# ============================================================
snow_full <- readRDS("snow_cleaned_full.Rda")
snow_full <- snow_full[-no_nbs,]

coords_full <- as.matrix(snow_full[,1:2])
y_full <- as.matrix(snow_full[,-c(1,2)])

coords_sf <- st_as_sf(
  data.frame(LON=coords_full[,1], LAT=coords_full[,2]),
  coords=c("LON","LAT"),
  crs=4326
)

coords_aeqd <- st_transform(coords_sf, aeqd_proj)
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

# ============================================================
# 2️⃣ Load BYM weekly + cov (strict version, NPZ)
# ============================================================
bymp_data <- np$load("trend_weekly_bym+cov.npz", allow_pickle=TRUE)

weekly_ini_BMP   <- bymp_data$f[["weekly_ini"]]
weekly_final_BMP <- bymp_data$f[["weekly_final"]]

# ============================================================
# 3️⃣ Weekly difference (final - initial)
# ============================================================
d <- weekly_final_BMP[,,,2] - weekly_ini_BMP[,,,2]

diff_mean <- apply(d, c(2,3), mean)
diff_sd   <- apply(d, c(2,3), sd)

# ============================================================
# 4️⃣ Create SF object (aligned)
# ============================================================
trend_sf <- st_transform(
  st_as_sf(
    data.frame(LON=coords[,1],
               LAT=coords[,2],
               diff_mean,
               diff_sd),
    coords=c("LON","LAT"),
    crs=4326
  ),
  crs=aeqd_proj
)

# ============================================================
# 5️⃣ Background map + equator
# ============================================================
world <- ne_countries(scale="medium", returnclass="sf")
world_north <- world[st_coordinates(st_centroid(world))[,2] > 0,]
world_aeqd <- st_transform(world_north, crs=aeqd_proj)

equator_points <- data.frame(
  lon = seq(-180,180,length.out=200),
  lat = rep(0,200)
)
equator_sf <- st_as_sf(equator_points,
                       coords=c("lon","lat"),
                       crs=4326)
equator_aeqd <- st_transform(equator_sf, crs=aeqd_proj)

# ============================================================
# 6️⃣ Week labels (month annotation)
# ============================================================
wk_names <- substr(colnames(y)[1:52],6,11)

month_raw <- substr(wk_names, 1, 2)
month_num <- as.integer(month_raw)
month_abbr <- month.abb[month_num]

# global color range
rg_mean <- range(diff_mean)

# ============================================================
# 7️⃣ Animated GIF (BYM weekly + cov only)
# ============================================================
for(i in 1:52){
  
  trend_sf$val <- trend_sf[[i]]
  
  p <- ggplot() +
    geom_sf(data=world_aeqd, fill="lightgray", color=NA) +
    geom_sf(data=trend_sf,
            aes(color=val),
            size=3, shape=18) +
    geom_sf(data=world_aeqd, fill=NA, color="black") +
    geom_sf(data=equator_aeqd,
            color="red", linetype="dashed", size=0.1) +
    scale_color_gradient2(
      low="red",
      mid="white",
      high="blue",
      midpoint=0,
      limits=rg_mean,
      guide=guide_colorbar(barwidth=25, barheight=0.5)
    ) +
    theme_minimal() +
    labs(
      title=paste("Trend Mean (BYM Weekly + Covariates) - Week", wk_names[i]),
      color=""
    ) +
    annotate("text",
             x = -7000000,
             y = -1000000,
             label = month_abbr[i],
             size = 14,
             fontface = "bold") +
    theme(
      legend.position="bottom",
      plot.title=element_text(hjust=0.5)
    )
  
  ggsave(paste0("plot_",i,".png"),
         plot=p,
         width=10,height=10,dpi=120)
}

# ============================================================
# 8️⃣ Create GIF
# ============================================================
png_files <- list.files(pattern="plot_\\d+\\.png$")
png_files <- png_files[order(as.numeric(gsub("\\D","",png_files)))]

gif <- image_read(png_files)
gif <- image_animate(gif, fps=5)
image_write(gif, "trend_animation_bym_weekly_cov.gif")




i = 17
trend_sf$val <- trend_sf[[i]]
p1 <- ggplot() +
  geom_sf(data=world_aeqd, fill="lightgray", color=NA) +
  geom_sf(data=trend_sf,
          aes(color=val),
          size=3, shape=18) +
  geom_sf(data=world_aeqd, fill=NA, color="black") +
  geom_sf(data=equator_aeqd,
          color="red", linetype="dashed", size=0.1) +
  scale_color_gradient2(
    low="red",
    mid="white",
    high="blue",
    midpoint=0,
    limits=rg_mean,
    guide=guide_colorbar(barwidth=25, barheight=0.5)
  ) +
  theme_minimal() +
  labs(
    title=paste("Trend Mean (BYM Weekly + Covariates) - Week", wk_names[i]),
    color=""
  ) +
  annotate("text",
           x = -7000000,
           y = -1000000,
           label = month_abbr[i],
           size = 14,
           fontface = "bold") +
  theme(
    legend.position="bottom",
    plot.title=element_text(hjust=0.5)
  )



i = 21
trend_sf$val <- trend_sf[[i]]
p2 <- ggplot() +
  geom_sf(data=world_aeqd, fill="lightgray", color=NA) +
  geom_sf(data=trend_sf,
          aes(color=val),
          size=3, shape=18) +
  geom_sf(data=world_aeqd, fill=NA, color="black") +
  geom_sf(data=equator_aeqd,
          color="red", linetype="dashed", size=0.1) +
  scale_color_gradient2(
    low="red",
    mid="white",
    high="blue",
    midpoint=0,
    limits=rg_mean,
    guide=guide_colorbar(barwidth=25, barheight=0.5)
  ) +
  theme_minimal() +
  labs(
    title=paste("Trend Mean (BYM Weekly + Covariates) - Week", wk_names[i]),
    color=""
  ) +
  annotate("text",
           x = -7000000,
           y = -1000000,
           label = month_abbr[i],
           size = 14,
           fontface = "bold") +
  theme(
    legend.position="bottom",
    plot.title=element_text(hjust=0.5)
  )

cowplot::plot_grid(p1,p2)







i = 38
trend_sf$val <- trend_sf[[i]]
p1 <- ggplot() +
  geom_sf(data=world_aeqd, fill="lightgray", color=NA) +
  geom_sf(data=trend_sf,
          aes(color=val),
          size=3, shape=18) +
  geom_sf(data=world_aeqd, fill=NA, color="black") +
  geom_sf(data=equator_aeqd,
          color="red", linetype="dashed", size=0.1) +
  scale_color_gradient2(
    low="red",
    mid="white",
    high="blue",
    midpoint=0,
    limits=rg_mean,
    guide=guide_colorbar(barwidth=25, barheight=0.5)
  ) +
  theme_minimal() +
  labs(
    title=paste("Trend Mean (BYM Weekly + Covariates) - Week", wk_names[i]),
    color=""
  ) +
  annotate("text",
           x = -7000000,
           y = -1000000,
           label = month_abbr[i],
           size = 14,
           fontface = "bold") +
  theme(
    legend.position="bottom",
    plot.title=element_text(hjust=0.5)
  )



i = 44
trend_sf$val <- trend_sf[[i]]
p2 <- ggplot() +
  geom_sf(data=world_aeqd, fill="lightgray", color=NA) +
  geom_sf(data=trend_sf,
          aes(color=val),
          size=3, shape=18) +
  geom_sf(data=world_aeqd, fill=NA, color="black") +
  geom_sf(data=equator_aeqd,
          color="red", linetype="dashed", size=0.1) +
  scale_color_gradient2(
    low="red",
    mid="white",
    high="blue",
    midpoint=0,
    limits=rg_mean,
    guide=guide_colorbar(barwidth=25, barheight=0.5)
  ) +
  theme_minimal() +
  labs(
    title=paste("Trend Mean (BYM Weekly + Covariates) - Week", wk_names[i]),
    color=""
  ) +
  annotate("text",
           x = -7000000,
           y = -1000000,
           label = month_abbr[i],
           size = 14,
           fontface = "bold") +
  theme(
    legend.position="bottom",
    plot.title=element_text(hjust=0.5)
  )

cowplot::plot_grid(p1,p2)





