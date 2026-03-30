rm(list = ls())

library(sf)
library(igraph)
library(ggplot2)
library(rnaturalearth)
library(cowplot)

# =====================================================
# LOAD p-value CSV (from python)
# =====================================================
df <- read.csv("ppp_cell_week20.csv")   # :contentReference[oaicite:0]{index=0}

coords <- as.matrix(df[,c("lon","lat")])

# =====================================================
# AEQD projection
# =====================================================
aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"

make_sf <- function(pvals){
  st_transform(
    st_as_sf(
      data.frame(
        LON = coords[,1],
        LAT = coords[,2],
        p   = pvals
      ),
      coords = c("LON","LAT"),
      crs = 4326
    ),
    crs = aeqd_proj
  )
}

sf_IID  <- make_sf(df$p_iid)
sf_BYM  <- make_sf(df$p_bym)
sf_BYMp <- make_sf(df$p_bym_cov)

# =====================================================
# BASE MAP
# =====================================================
world <- ne_countries(scale="medium", returnclass="sf")
world_north <- world[st_coordinates(st_centroid(world))[,2] > 0,]
world_aeqd <- st_transform(world_north, crs=aeqd_proj)

# =====================================================
# PLOT
# =====================================================
plot_p <- function(sf_obj, title){
  ggplot() +
    geom_sf(data = world_aeqd, fill = "lightgray", color = NA) +
    geom_sf(data = sf_obj, aes(color = p), size = 2, shape = 18) +
    geom_sf(data = world_aeqd, fill = NA, color = "black") +
    scale_color_viridis_c(
      option = "C",
      limits = c(0,1),
      guide = guide_colorbar(barwidth = 20, barheight = 0.5)
    ) +
    theme_minimal() +
    labs(title = title, color = "p-value") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
}

p1 <- plot_p(sf_IID,  "P-value - IND")
p2 <- plot_p(sf_BYM,  "P-value - Weekly BYM")
p3 <- plot_p(sf_BYMp, "P-value - Weekly BYM+")

cowplot::plot_grid(p1, p2, p3, nrow = 1)