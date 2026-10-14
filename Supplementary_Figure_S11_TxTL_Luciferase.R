#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(readr)
  library(RColorBrewer)
})

# ==== Load data ====
data_file <- "/Users/muratbuyukyoruk/Desktop/RT_projects/DRT2_data/TxTL/20250815_TXTL.txt"  # Replace with actual path to your file

raw <- read_delim(file = data_file,
                  delim = "\t",
                  locale = locale(encoding = "UTF-16LE"),
                  col_types = cols(),
                  trim_ws = TRUE)

names(raw) <- make.names(names(raw), unique = TRUE)

# ==== Split Average and SEM ====
avg_df <- raw %>%
  filter(Type == "Average") %>%
  select(-Name, -Type) %>%
  pivot_longer(cols = -Time, names_to = "Condition", values_to = "Average")

sem_df <- raw %>%
  filter(Type == "SEM") %>%
  select(-Name, -Type) %>%
  pivot_longer(cols = -Time, names_to = "Condition", values_to = "SEM")

# ==== Merge and process ====
df <- left_join(avg_df, sem_df, by = c("Time", "Condition")) %>%
  mutate(
    Time = as.numeric(Time),
    Average = as.numeric(Average),
    SEM = as.numeric(SEM)
  )

# ==== Fix encoding-renamed columns ====
df$Condition <- recode(df$Condition, "X70_EtOH" = "70% EtOH")

# ==== Assign Groups and Labels ====
df <- df %>%
  mutate(
    Group = case_when(
      Condition %in% c("Luciferase", "70% EtOH", "Ampicillin", "Erythromycin") ~ "Control",
      Condition %in% c("EGFP", "Neo") ~ "Sample",
      TRUE ~ "Other"
    ),
    Label = paste(Group, "-", Condition)
  )

# ==== Set Factor Order for Consistency ====
df$Label <- factor(df$Label, levels = c(
  "Control - Luciferase",
  "Control - 70% EtOH",
  "Control - Ampicillin",
  "Control - Erythromycin",
  "Sample - EGFP",
  "Sample - Neo"
))

# ==== Define Manual Color Palette with Swap ====
labels <- levels(df$Label)
palette <- setNames(brewer.pal(6, "Dark2"), labels)

# Swap Ampicillin and Neo colors
tmp <- palette["Control - Ampicillin"]
palette["Control - Ampicillin"] <- palette["Sample - Neo"]
palette["Sample - Neo"] <- tmp

# ==== Plot ====
p <- ggplot(df, aes(x = Time, group = Label)) +
  # SEM ribbon
  geom_ribbon(aes(ymin = Average - SEM, ymax = Average + SEM, fill = Label), alpha = 0.5, color = NA) +
  # SEM border lines
  geom_line(aes(y = Average + SEM, color = Label), size = 0.25, show.legend = FALSE) +
  geom_line(aes(y = Average - SEM, color = Label), size = 0.25, show.legend = FALSE) +
  # Main curve
  geom_line(aes(y = Average, color = Label), size = 0.75, show.legend = TRUE) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(
    title = "Luciferase Signal Over Time",
    x = "Time (minutes)",
    y = expression(bold("Luminescence (RLU)")),
    color = "Condition",
    fill  = "Condition"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x       = element_text(face = "bold"),
    axis.title.y       = element_text(face = "bold"),
    axis.text.x        = element_text(face = "bold"),
    axis.text.y        = element_text(face = "bold"),
    panel.border       = element_rect(colour = "black", fill = NA, size = 1),
    axis.line          = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    legend.position    = "right",
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.key         = element_blank()
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 1)),
    fill  = guide_legend(order = 1, override.aes = list(alpha = 0.5))
  )

# ==== Display ====
print(p)