library(tidyverse)
library(ggplot2)
library(readr)
# ============================================================
# Read data
# ============================================================

df_old <- read_delim("/Users/muratbuyukyoruk/Documents/DRT2/DRT2_biochem/FACS_plot.txt", delim = "\t", show_col_types = FALSE)
df <- read_delim("/Users/muratbuyukyoruk/Desktop/DRT2_Cell_RTR/Code/input_data_for_codes/20261009_Compiled_percent_depolarized_data.txt", delim = "\t", show_col_types = FALSE)


# ============================================================
# Format labels and set ordering
# ============================================================

df <- df %>%
  mutate(
    Sample = factor(
      Sample,
      levels = c(
        "MG1655",
        "EGFP (Y66A)",
        "KpnNEO",
        "NEOc3.1",
        "NEOc4.1",
        "NEOc7.6"
      ),
      labels = c(
        "MG1655",
        "EGFP (Y66A)",
        "KpnNEO",
        "NEOc3.1",
        "NEOc4.1",
        "NEOc7.6"
      )
    ),
    
    condition = factor(
      condition,
      levels = c(
        "N/A",
        "CCCP",
        "Glucose",
        "Arabinose",
        "Arabinose + CCCP"
      ),
      labels = c(
        "No CCCP",
        "CCCP",
        "Glucose",
        "Arabinose",
        "Arabinose + CCCP"
      )
    )
  )


# ============================================================
# Plot
# ============================================================

p <- ggplot(
  df,
  aes(
    x = condition,
    y = average
  )
) +
  
  # Bars
  geom_col(
    width = 0.72,
    fill = "#8EC5E8",
    color = "black",
    linewidth = 0.35
  ) +
  
  # SEM
  geom_errorbar(
    aes(
      ymin = pmax(average - SEM, 0),
      ymax = average + SEM
    ),
    width = 0.20,
    linewidth = 0.45
  ) +
  
  # Separate panels for MG1655, EGFP, and KpnNeo
  facet_grid(
    . ~ Sample,
    scales = "free_x",
    space = "free_x"
  ) +
  
  # Y axis
  scale_y_continuous(
    limits = c(0, 105),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = NULL,
    y = "% Cells Depolarized"
  ) +
  
  theme_classic(base_size = 10) +
  
  theme(
    # --------------------------------------------------------
    # Facet headers
    # --------------------------------------------------------
    strip.background = element_rect(
      fill = "grey90",
      color = "grey45",
      linewidth = 0.5
    ),
    
    strip.text = element_text(
      size = 9,
      face = "plain",
      color = "black",
      margin = margin(
        t = 3,
        r = 4,
        b = 3,
        l = 4
      )
    ),
    
    # --------------------------------------------------------
    # X axis
    # --------------------------------------------------------
    axis.text.x = element_text(
      size = 8,
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    
    axis.ticks.x = element_blank(),
    
    axis.title.x = element_blank(),
    
    # --------------------------------------------------------
    # Y axis
    # --------------------------------------------------------
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 10,
      face = "bold",
      color = "black",
      margin = margin(r = 7)
    ),
    
    axis.ticks.y = element_line(
      color = "black",
      linewidth = 0.4
    ),
    
    # --------------------------------------------------------
    # Panel borders
    # --------------------------------------------------------
    panel.border = element_rect(
      color = "grey45",
      fill = NA,
      linewidth = 0.5
    ),
    
    panel.spacing.x = unit(
      0.12,
      "cm"
    ),
    
    # --------------------------------------------------------
    # Grid lines
    # --------------------------------------------------------
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    
    panel.grid.minor.y = element_blank(),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.minor.x = element_blank(),
    
    # --------------------------------------------------------
    # Legend / margins
    # --------------------------------------------------------
    legend.position = "none",
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )


# ============================================================
# Display plot in RStudio
# ============================================================

p

