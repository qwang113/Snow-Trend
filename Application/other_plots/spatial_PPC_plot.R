rm(list = ls())
library(dplyr)
library(ggplot2)

df = read.csv("D:/Snow-Trend/py_models/PPC/ppp_weekly.csv")

df <- df %>%
  mutate(
    model = recode(model,
                   "BYM weekly" = "BYM",
                   "BYM weekly + cov + lon" = "BYM+"),
    extreme = (p_value < 0.05 | p_value > 0.95)
  )

colors <- c(
  "IID" = "#0072B2",
  "BYM" = "#E69F00", 
  "BYM+" = "#009E73"
)

linetypes <- c(
  "IID" = "dotted", 
  "BYM" = "dashed", 
  "BYM+" = "dotdash"
)

ggplot(df, aes(x = week, y = p_value, color = model, linetype = model)) +
  
  geom_line(linewidth = 1.6, alpha = 0.95) +
  
  geom_point(
    data = df %>% filter(!extreme),
    size = 1.2,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  
  geom_point(
    data = df %>% filter(extreme),
    shape = 4,
    stroke = 1.2,
    size = 3,
    show.legend = FALSE
  ) +
  
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             color = "gray50") +
  
  geom_hline(yintercept = c(0.05, 0.95),
             linetype = "dashed",
             color = "#C44E52",
             alpha = 0.8) +
  
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  
  labs(
    x = "Week (WCY)",
    y = "Posterior predictive p-value",
    color = "Model",
    linetype = "Model"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 15),
    legend.key.width = unit(2.5, "cm")
  ) 
