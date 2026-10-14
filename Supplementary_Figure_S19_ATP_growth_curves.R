library(ggplot2)
library(tidyverse)

# Load data
data <- read.delim("~/Desktop/ATP_assay_growrth_curves_20250626.txt", header = TRUE)

# Simplify condition names
condition_labels <- c(
  Glucose_0.2 = "Glucose",
  Arabinose_0.2 = "Arabinose"
)

# Prepare average and SEM data
df_avg <- data %>%
  filter(Type == "Average") %>%
  pivot_longer(cols = c(Glucose_0.2, Arabinose_0.2),
               names_to = "Condition", values_to = "Value")

df_sem <- data %>%
  filter(Type == "SEM") %>%
  pivot_longer(cols = c(Glucose_0.2, Arabinose_0.2),
               names_to = "Condition", values_to = "SEM")

# Merge and annotate
df_plot <- left_join(df_avg, df_sem, by = c("Sample", "Time", "Condition")) %>%
  mutate(
    ConditionSimple = recode(Condition, !!!condition_labels),
    Label = paste(Sample, ConditionSimple, sep = " - "),
    ConditionSimple = factor(ConditionSimple, levels = c("Glucose", "Arabinose"))
  )

# Set custom label order
df_plot$Label <- factor(df_plot$Label, levels = c(
  "EGFP - Glucose", "EGFP - Arabinose",
  "MS2L - Glucose", "MS2L - Arabinose",
  "NEO - Glucose",  "NEO - Arabinose"
))

# Plot
p <- ggplot(df_plot, aes(x = Time, group = Label)) +
  geom_ribbon(aes(ymin = Value - SEM,
                  ymax = Value + SEM,
                  fill = Sample),
              alpha = 0.5,
              color = NA) +
  geom_line(aes(y = Value + SEM, color = Sample),
            size = 0.25, show.legend = FALSE) +
  geom_line(aes(y = Value - SEM, color = Sample),
            size = 0.25, show.legend = FALSE) +
  geom_line(aes(y = Value,
                color = Sample,
                linetype = ConditionSimple),
            size = 0.75, show.legend = TRUE) +
  
  # NEW POINTS
  geom_point(aes(y = Value,
                 color = Sample,
                 shape = ConditionSimple),
             size = 2.2,
             stroke = 0.3,
             fill = "black",
             show.legend = TRUE) +
  
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("Glucose" = "dashed", "Arabinose" = "solid")) +
  scale_shape_manual(values = c("Glucose" = 21, "Arabinose" = 19)) +   # NEW
  
  labs(
    title = "Growth Curves (ATP Assay)",
    x = "Time (minutes)",
    y = expression(bold(Growth~(OD[600]))),
    color = "Sample",
    fill  = "Sample",
    linetype = "Condition",
    shape   = "Condition"
  ) +
  scale_y_continuous(limits = c(0, 0.3045), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-60, NA), expand = c(0, 0)) +
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
    fill  = guide_legend(order = 1, override.aes = list(alpha = 0.5)),
    linetype = guide_legend(order = 2, override.aes = list(color = "black")),
    shape = guide_legend(order = 2)  # NEW
  )

# Show plot
print(p)