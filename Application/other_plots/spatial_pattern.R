rm(list = ls())

library(sf)
library(igraph)
library(ggplot2)
library(rnaturalearth)

# ============================================================
# CONFIG
# ============================================================

DIST_TH <- 0.22
aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"

no_nbs <- c(
  57,170,236,269,343,685,946,947,989,
  1037,1084,1090,1109,1118,1127,1176,1203
)

# ============================================================
# LOAD DATA
# ============================================================

snow <- readRDS("snow_cleaned_full.Rda")

coords <- as.matrix(snow[,1:2])
y_full <- as.matrix(snow[,-c(1,2)])

coords <- coords[-no_nbs,]
y_full <- y_full[-no_nbs,]

# ============================================================
# TWO LARGEST COMPONENTS
# ============================================================

coords_sf <- st_as_sf(
  data.frame(LON=coords[,1], LAT=coords[,2]),
  coords=c("LON","LAT"),
  crs=4326
)

coords_aeqd <- st_transform(coords_sf, aeqd_proj)

xy <- st_coordinates(coords_aeqd) / 1e6

D <- as.matrix(dist(xy))

W_full <- (D <= DIST_TH) * 1
diag(W_full) <- 0

g <- graph_from_adjacency_matrix(W_full, mode="undirected")

comp <- components(g)

sizes <- comp$csize
order_idx <- order(sizes, decreasing=TRUE)

keep <- sort(which(comp$membership %in% order_idx[1:2]))

coords <- coords[keep,]
y <- y_full[keep,]

S <- nrow(coords)
T <- ncol(y)

cat("Using S =", S, "\n")

# ============================================================
# REBUILD ADJACENCY
# ============================================================

coords_sf <- st_as_sf(
  data.frame(LON=coords[,1], LAT=coords[,2]),
  coords=c("LON","LAT"),
  crs=4326
)

coords_aeqd <- st_transform(coords_sf, aeqd_proj)

xy <- st_coordinates(coords_aeqd) / 1e6

D <- as.matrix(dist(xy))

W <- (D <= DIST_TH) * 1
diag(W) <- 0

w_sum <- sum(W)

# ============================================================
# WEEK INDEX (same as Python)
# ============================================================

week <- (0:(T-1)) %% 52

# ============================================================
# OBSERVED WEEKLY MORAN'S I
# ============================================================

I_obs_week <- numeric(52)

for(w in 0:51){
  
  idx <- which(week == w)
  
  y_bar <- rowMeans(y[,idx])
  
  yc <- y_bar - mean(y_bar)
  
  I_obs_week[w+1] <- (S / w_sum) *
    ( t(yc) %*% W %*% yc ) /
    ( t(yc) %*% yc )
  
}

I_obs_week <- as.numeric(I_obs_week)

# ============================================================
# WEEKLY SNOW PROPORTION
# ============================================================

snow_prop <- numeric(52)

for(w in 0:51){
  
  idx <- which(week == w)
  
  snow_prop[w+1] <- mean(y[,idx])
  
}

# ============================================================
# DATA FRAME
# ============================================================

df <- data.frame(
  week = 1:52,
  snow = snow_prop,
  moran = I_obs_week
)

# ============================================================
# SEASON BLOCKS
# Week1 = Aug
# ============================================================
season_df <- data.frame(
  xmin=c(1,5,18,31,44),
  xmax=c(4,17,30,43,52),
  season=c("Summer","Fall","Winter","Spring","Summer")
)

season_colors <- c(
  Winter="#CFE8FF",
  Spring="#DFF5D8",
  Summer="#FFE7A3",
  Fall="#FFD1B3"
)

# ============================================================
# Snow proportion
# ============================================================

p1 <- ggplot(df,aes(week,snow)) +
  
  geom_rect(
    data=season_df,
    aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=season),
    alpha=0.35,
    inherit.aes=FALSE
  ) +
  
  geom_line(size=1.2,color="#2C7FB8") +
  geom_point(size=2,color="#2C7FB8") +
  
  scale_fill_manual(values=season_colors) +
  
  scale_x_continuous(breaks=seq(1,52,4)) +
  
  theme_bw() +
  
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    legend.position="top",
    plot.title=element_text(hjust=0.5)
  ) +
  
  labs(
    x="Week (Week 1 = Aug 1)",
    y="Snow proportion",
    fill="Season",
    title="Weekly Proportion of Snow Presence Across Locations"
  )

# ============================================================
# Moran's I
# ============================================================

p2 <- ggplot(df,aes(week,moran)) +
  
  geom_rect(
    data=season_df,
    aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=season),
    alpha=0.35,
    inherit.aes=FALSE
  ) +
  
  geom_line(size=1.2,color="#D7301F") +
  geom_point(size=2,color="#D7301F") +
  
  scale_fill_manual(values=season_colors) +
  
  scale_x_continuous(breaks=seq(1,52,4)) +
  
  theme_bw() +
  
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    legend.position="top",
    plot.title=element_text(hjust=0.5)
  ) +
  
  labs(
    x="Week (Week 1 = Aug 1)",
    y="Moran's I",
    fill="Season",
    title="Weekly Moran's I of Snow Presence Across Locations"
  )

p1
p2
