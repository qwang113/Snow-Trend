args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_holdout_scores.R <draw_scores.csv> <output.png>")
}

suppressPackageStartupMessages(library(ggplot2))
library(grid)

scores <- read.csv(args[[1]], check.names = FALSE)
scores$Model <- factor(scores$Model, levels = c("IND", "BYM", "BYM+"))

model_colors <- c(
  "IND" = "#0072B2",
  "BYM" = "#E69F00",
  "BYM+" = "#009E73"
)

make_panel <- function(metric, panel_title) {
  ggplot(scores, aes(x = Model, y = .data[[metric]], fill = Model)) +
    geom_boxplot(
      width = 0.58,
      linewidth = 0.55,
      alpha = 0.9,
      outlier.shape = 16,
      outlier.size = 0.45,
      outlier.alpha = 0.18
    ) +
    scale_fill_manual(values = model_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
    labs(title = panel_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(size = 13, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 12, 8, 12)
    )
}

plots <- list(
  make_panel("Weekly probability MSE", "MSE"),
  make_panel("Log loss", "Log Loss"),
  make_panel("ROC-AUC", "ROC-AUC"),
  make_panel("Accuracy @ 0.5", "Accuracy @ 0.5")
)

png(args[[2]], width = 2400, height = 1900, res = 220, bg = "white")
grid.newpage()
layout <- grid.layout(nrow = 2, ncol = 2)
pushViewport(viewport(layout = layout))
for (i in seq_along(plots)) {
  row <- ((i - 1L) %/% 2L) + 1L
  col <- ((i - 1L) %% 2L) + 1L
  print(plots[[i]], vp = viewport(layout.pos.row = row, layout.pos.col = col))
}
dev.off()
