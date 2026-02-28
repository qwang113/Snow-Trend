rm(list = ls())

library(sf)
library(ggplot2)
library(rnaturalearth)
library(cowplot)


df <- read.csv("D:/77/Research/temp/snow/ppc_time_avg_pvalues.csv")


aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"


world <- ne_countries(scale = "medium", returnclass = "sf")

world_north <- world[
  st_coordinates(st_centroid(world))[,2] > 0,
]

world_aeqd <- st_transform(world_north, crs = aeqd_proj)


pts_sf <- st_as_sf(df, coords = c("x","y"), crs = 4326)
pts_aeqd <- st_transform(pts_sf, crs = aeqd_proj)


equator_points <- data.frame(
  lon = seq(-180, 180, length.out = 300),
  lat = rep(0, 300)
)

equator_sf <- st_as_sf(
  equator_points,
  coords = c("lon","lat"),
  crs = 4326
)

equator_aeqd <- st_transform(equator_sf, crs = aeqd_proj)


plot_ppc_map <- function(column_name, plot_title){
  
  df_plot <- pts_aeqd
  df_plot$color_val <- df_plot[[column_name]]
  
  ggplot() +
    
    geom_sf(data = world_aeqd,
            fill = "lightgray",
            color = NA) +
    
    geom_sf(data = df_plot,
            aes(color = color_val),
            size = 2,
            shape = 18) +
    
    geom_sf(data = world_aeqd,
            fill = NA,
            color = "black") +
    
    geom_sf(data = equator_aeqd,
            color = "red",
            linetype = "dashed",
            size = 0.1) +
    
    scale_color_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0.5,
      limits = c(0,1)
    ) +
    
    theme_minimal() +
    
    labs(
      title = plot_title,
      color = ""
    ) +
    
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 20)
    ) +
    
    guides(
      color = guide_colorbar(
        barwidth = 10,
        barheight = 0.5
      )
    )
}



p1 <- plot_ppc_map("p_iid",
                   "IID")

p2 <- plot_ppc_map("p_week",
                   "BYM weekly")

p3 <- plot_ppc_map("p_week_cov",
                   "BYM weekly")


combined_plot <- plot_grid(p1, p2, p3, nrow = 1)

print(combined_plot)


ggsave(
  "ppc_time_avg_maps.png",
  combined_plot,
  width = 18,
  height = 6,
  dpi = 300
)