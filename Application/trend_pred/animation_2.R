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
np <- import("numpy", convert = TRUE)

# ============================================================
# CONFIG
# ============================================================
DIST_TH <- 0.22
aeqd_proj <- "+proj=aeqd +lat_0=90 +lon_0=-100"

# ============================================================
# 1️⃣ Load snow (FULL DATA)
# ============================================================
snow_full <- readRDS("snow_cleaned_full.Rda")

coords <- as.matrix(snow_full[,1:2])
y <- as.matrix(snow_full[,-c(1,2)])

S <- nrow(coords)
cat("FULL S =", S, "\n")

# ============================================================
# 2️⃣ Load BYM weekly + cov (strict version, NPZ)
# ============================================================
bymp_data <- np$load("D:/77/research/temp/snow/trend_weekly_bym+cov.npz", allow_pickle=TRUE)

weekly_ini_BMP   <- bymp_data$f[["weekly_ini"]]
weekly_final_BMP <- bymp_data$f[["weekly_final"]]

# ============================================================
# 3️⃣ Weekly difference (final - initial)
# ============================================================
d_py <- weekly_final_BMP[,,,1] - weekly_ini_BMP[,,,1]

diff_mean <- py_to_r(np$mean(d_py, axis = 0L)$tolist())
diff_sd   <- py_to_r(np$std(d_py, axis = 0L)$tolist())

diff_mean <- do.call(rbind, diff_mean)
diff_sd   <- do.call(rbind, diff_sd)

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

all_vals <- as.numeric(diff_mean)
# q_global <- quantile(all_vals, c(0.025, 0.975), na.rm = TRUE)
q_global <- c(-0.25,0.25)

# ============================================================
# 7️⃣ Plot loop
# ============================================================
for(i in 2:52){
  
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
      limits=q_global,
      oob=scales::squish,
      guide=guide_colorbar(barwidth=25, barheight=0.5)
    ) +
    
    theme_minimal() +
    
    labs(
      title=paste("Trend Mean (bym+) - Week", wk_names[i]),
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
         width=10,
         height=10,
         dpi=150)
}

# ============================================================
# 8️⃣ Create GIF
# ============================================================

png_files <- list.files(pattern="plot_\\d+\\.png$")
png_files <- png_files[order(as.numeric(gsub("\\D","",png_files)))]

png_files <- png_files[as.numeric(gsub("\\D","",png_files)) >= 2]

library(magick)

gif <- image_read(png_files)
gif <- image_animate(gif, fps=5)
image_write(gif, "animation_bym+.gif")

pos_count <- colSums(diff_mean > 0)[-1]
df_plot <- data.frame(
  week = 2:52,
  n_pos = pos_count
)

ggplot(df_plot, aes(x = week, y = n_pos)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  # bounds
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = S, linetype = "dashed", color = "blue") +
  
  scale_x_continuous(breaks = seq(1, 52, by = 4)) +
  labs(
    x = "Week",
    y = "Number of cells with trend > 0",
    title = "Weekly count of cells with positive trend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )




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
    limits=q_global,          # 🔥 用quantile
    oob=scales::squish,
    guide=guide_colorbar(barwidth=25, barheight=0.5)
  ) +
  theme_minimal() +
  labs(
    title=paste("Trend Mean - Week", wk_names[i]),
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
    limits=q_global,          # 🔥 用quantile
    oob=scales::squish,
    guide=guide_colorbar(barwidth=25, barheight=0.5)
  ) +
  theme_minimal() +
  labs(
    title=paste("Trend Mean - Week", wk_names[i]),
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
    limits=q_global,          # 🔥 用quantile
    oob=scales::squish,
    guide=guide_colorbar(barwidth=25, barheight=0.5)
  ) +
  theme_minimal() +
  labs(
    title=paste("Trend Mean - Week", wk_names[i]),
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
    limits=q_global,
    oob=scales::squish,
    guide=guide_colorbar(barwidth=25, barheight=0.5)
  ) +
  theme_minimal() +
  labs(
    title=paste("Trend Mean - Week", wk_names[i]),
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


library(reticulate)
use_condaenv("CPD", required = TRUE)

py_run_string("
import numpy as np

data = np.load('D:/77/research/temp/snow/trend_weekly_bym+cov.npz', allow_pickle=True)

ini   = data['weekly_ini']
final = data['weekly_final']

d = final[:,:,:,0] - ini[:,:,:,0]

weeks = [17,21,38,44]

rows = []

for w in weeks:
    samples = d[:,:,w]
    
    mean = samples.mean(axis=0)
    q025 = np.quantile(samples, 0.025, axis=0)
    q975 = np.quantile(samples, 0.975, axis=0)
    
    inc = int(np.sum(mean > 0))
    dec = int(np.sum(mean < 0))
    
    sig_inc = int(np.sum(q025 > 0))
    sig_dec = int(np.sum(q975 < 0))
    
    rows.append((w, inc, dec, sig_inc, sig_dec, mean.mean(), np.median(mean)))

# ---- build LaTeX ----
latex = \"\"\"\\begin{table}[htbp]
\\centering
\\small
\\caption{Summary of spatial snow trend changes (posterior quantile-based).}
\\begin{tabular}{ccccccc}
\\hline
Week & Increase & Decrease & Sig. Increase & Sig. Decrease & Mean & Median \\\\
\\hline
\"\"\"

for r in rows:
    latex += \"%d & %d & %d & %d & %d & %.3f & %.3f \\\\\\n\" % r

latex += \"\"\"\\hline
\\end{tabular}
\\end{table}\"\"\"

print(latex)
")



