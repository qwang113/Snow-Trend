rm(list = ls())
library(dplyr)
library(ggplot2)
df = read.csv("D:/Snow-Trend/py_models/PPC/ppp_weekly.csv")
df <- df %>%
  mutate(
    extreme = (p_value < 0.05 | p_value > 0.95)
  )

colors <- c(
  "IID" = "#4C72B0",
  "BYM weekly" = "#DD8452",
  "BYM weekly + cov" = "#55A868"
)

ggplot(df, aes(x = week, y = p_value, color = model)) +
  
  geom_line(linewidth = 0.9, alpha = 0.9) +
  
  geom_point(
    data = df %>% filter(!extreme),
    size = 1.2,
    alpha = 0.35
  ) +
  
  geom_point(
    data = df %>% filter(extreme),
    size = 2.6,
    alpha = 1
  ) +
  
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             color = "gray50") +
  
  geom_hline(yintercept = c(0.05, 0.95),
             linetype = "dashed",
             color = "#C44E52",
             alpha = 0.8) +
  
  scale_color_manual(values = colors) +
  
  labs(
    x = "Week",
    y = "Posterior predictive p-value",
    color = "Model"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 15)
  )
