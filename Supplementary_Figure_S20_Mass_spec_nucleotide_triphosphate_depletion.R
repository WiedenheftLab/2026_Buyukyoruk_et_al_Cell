library(tidyverse)
library(stringr)

input_file <- "~/Desktop/Mass_Spec_ATPs_06042026_plot.txt"
input_file <- "~/Desktop/Mass_Spec_GTPs_06042026_plot.txt"
input_file <- "~/Desktop/Mass_Spec_CTPs_06042026_plot.txt"
input_file <- "~/Desktop/Mass_Spec_UTPs_06042026_plot.txt"

df <- read_tsv(input_file, show_col_types = FALSE)

sample_order <- c(
  "EGFP, Gluc.",
  "EGFP, Arab.",
  "Neo, Gluc.",
  "Neo, Arab.",
  "Neo-EGFP, Gluc.",
  "Neo-EGFP, Arab.",
  "Neo (TM)-EGFP, Gluc.",
  "Neo (TM)-EGFP, Arab."
)

ntp_order <- c("ATP", "UTP", "GTP", "CTP")

ntp_colors <- c(
  "ATP" = "#d7d1e7",
  "UTP" = "#ecd2d1",
  "GTP" = "#dae6d9",
  "CTP" = "#d1deec"
)

df <- df %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    NTP = factor(NTP, levels = ntp_order)
  )

p <- ggplot(df, aes(x = sample, y = Average, fill = NTP)) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(
      ymin = Average - SEM,
      ymax = Average + SEM
    ),
    width = 0.18,
    linewidth = 0.45
  ) +
  facet_wrap(~NTP, nrow = 1) +
  scale_fill_manual(values = ntp_colors) +
  scale_y_continuous(
    limits = c(0, 125),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Normalized Peak Area Ratios (NTP/C13ATP Internal Standard) Percentage (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      size = 15,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    panel.spacing = unit(1.1, "lines"),
    plot.margin = margin(10, 12, 10, 10),
    
    panel.grid.major.y = element_line(
      color = "gray80",
      linewidth = 0.5,
      linetype = "solid"
    ),
    panel.grid.minor.y = element_line(
      color = "gray85",
      linewidth = 0.5,
      linetype = "22"
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    x = NULL,
    y = str_wrap(
      "Normalized Peak Area Ratios (NTP/C13ATP Internal Standard) Percentage (%)",
      width = 38
    )
  )

print(p)
