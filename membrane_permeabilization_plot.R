# ============================================================
# Supplemental Figure S19
# Flow cytometry: DAPI and PI staining
#
# Bars: Mean
# Error bars: SEM
# Points: Biological replicates (n = 3)
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1. Experimental data
# ------------------------------------------------------------

df <- data.frame(
  Condition = rep(c(
    "DAPI only",
    "DAPI + PI",
    "i-PrOH permeabilized",
    "NEO + glucose",
    "NEO + arabinose"
  ), each = 3),
  
  Replicate = rep(1:3, 5),
  
  DAPI = c(
    99.71, 99.92, 99.89,
    99.87, 99.98, 99.94,
    99.90, 99.80, 99.78,
    99.90, 99.90, 100.00,
    99.34, 97.80, 99.09
  ),
  
  DAPI_PI = c(
    2.51, 2.42, 2.69,
    5.17, 7.48, 5.74,
    82.20, 89.70, 92.30,
    13.10, 18.00, 20.00,
    7.54, 16.90, 7.19
  )
)

# ------------------------------------------------------------
# 2. Define groups and condition order
# ------------------------------------------------------------

df <- df %>%
  mutate(
    Group = ifelse(
      Condition %in% c(
        "DAPI only",
        "DAPI + PI",
        "i-PrOH permeabilized"
      ),
      "Controls",
      "NEO expression"
    ),
    
    Group = factor(
      Group,
      levels = c("Controls", "NEO expression")
    ),
    
    Condition = factor(
      Condition,
      levels = c(
        "DAPI only",
        "DAPI + PI",
        "i-PrOH permeabilized",
        "NEO + glucose",
        "NEO + arabinose"
      )
    )
  )

# ------------------------------------------------------------
# 3. Convert to long format
# ------------------------------------------------------------

df_long <- df %>%
  pivot_longer(
    cols = c(DAPI, DAPI_PI),
    names_to = "Population",
    values_to = "Percentage"
  ) %>%
  mutate(
    Population = factor(
      Population,
      levels = c("DAPI", "DAPI_PI"),
      labels = c("DAPI+", "DAPI+ / PI+")
    )
  )

# ------------------------------------------------------------
# 4. Calculate mean and SEM
# ------------------------------------------------------------

summary_df <- df_long %>%
  group_by(Group, Condition, Population) %>%
  summarise(
    Mean = mean(Percentage),
    SEM = sd(Percentage) / sqrt(n()),
    .groups = "drop"
  )

print(summary_df)

# ------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------

p <- ggplot(
  summary_df,
  aes(x = Condition, y = Mean, fill = Population)
) +
  
  # Mean bars
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  # SEM error bars
  geom_errorbar(
    aes(
      ymin = Mean - SEM,
      ymax = Mean + SEM
    ),
    position = position_dodge(width = 0.8),
    width = 0.15,
    linewidth = 0.5
  ) +
  
  # Individual biological replicates
  geom_point(
    data = df_long,
    aes(
      x = Condition,
      y = Percentage,
      group = Population
    ),
    inherit.aes = FALSE,
    position = position_jitterdodge(
      jitter.width = 0.05,
      dodge.width = 0.8
    ),
    shape = 21,
    size = 2,
    fill = "white",
    color = "black"
  ) +
  
  # Facets
  facet_grid(
    ~ Group,
    scales = "free_x",
    space = "free_x"
  ) +
  
  # Muted colors
  scale_fill_manual(
    values = c(
      "DAPI+" = "#91B9CC",
      "DAPI+ / PI+" = "#A4C9AC"
    )
  ) +
  
  # Axes
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  coord_cartesian(ylim = c(0, 110)) +
  
  labs(
    x = NULL,
    y = "Cells (%)",
    fill = NULL
  ) +
  
  # Standard ggplot theme
  theme_bw(base_size = 12) +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    ),
    
    legend.position = "top",
    legend.justification = "right"
  )

# ------------------------------------------------------------
# 6. Display in RStudio
# ------------------------------------------------------------

print(p)