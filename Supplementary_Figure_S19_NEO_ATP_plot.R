library(tidyverse)
library(RColorBrewer)

sci_labeller <- function(x) {
  labs <- sapply(x, function(v) {
    if (v == 0) return("0")
    e <- floor(log10(v))
    c <- v / 10^e
    sprintf("bold(%g %s 10^%d)", c, "%*%", e)
  })
  parse(text = labs)
}

# 1) pivot to long so we can map aesthetic by Protein & Condition
data <- read.delim("~/Desktop/NEO_ATP_plot_data.txt", header = TRUE)
data <- read.delim("~/Desktop/NEO_ATP_plot_data_20250626.txt", header = TRUE)

df_avg <- data %>% 
  filter(Type == "Average") %>% 
  pivot_longer(
    cols      = c(EGFP_Glucose, EGFP_Arabinose, NEO_Glucose, NEO_Arabinose, MS2L_Glucose, MS2L_Arabinose),
    names_to  = c("Expression", "Condition"),
    names_sep = "_",
    values_to = "Value"
  )

# 2) pivot the SEM rows the same way
df_sem <- data %>% 
  filter(Type == "SEM") %>% 
  pivot_longer(
    cols      = c(EGFP_Glucose, EGFP_Arabinose, NEO_Glucose, NEO_Arabinose, MS2L_Glucose, MS2L_Arabinose),
    names_to  = c("Expression", "Condition"),
    names_sep = "_",
    values_to = "SEM"
  )

# 3) join so each row has both Value and SEM
df_plot <- left_join(df_avg, df_sem,
                     by = c("Sample", "Time", "Expression", "Condition"))
# 2) define colors (EGFP=green, NEO=first Dark2) and line types (Arabinose=solid, Glucose=dashed)
dark2cols <- brewer.pal(3, "Dark2")
col_map  <- c(EGFP = dark2cols[1], NEO = dark2cols[3], MS2L = dark2cols[2] )
lt_map   <- c(Arabinose = "solid", Glucose = "dashed")

# 3) plot
ggplot(df_plot, aes(x = Time, y = Value, color = Expression, fill = Expression,
                    linetype = Condition, group = interaction(Expression,Condition))) +
  geom_ribbon(
    aes(ymin = Value - SEM, ymax = Value + SEM, colour = Expression),
    alpha    = 0.2,
    size     = 0.3,
    linetype = "solid"   # force a solid border on every ribbon
  ) +
  geom_line(size = 1) +
  scale_color_manual(values = col_map) +
  scale_fill_manual(values = col_map) +
  scale_linetype_manual(values = lt_map) +
  labs(
    x     = "Time after induction (min)",
    y     = "NEO expression impact on ATP - Luminescence (RLU)",
  ) +
  scale_y_continuous(
    limits = c(0, 1000000),
    expand = c(0, 0),
    labels = sci_labeller    # <— re‑add your 1×10^5 style tick labels
  ) +
  scale_x_continuous( expand = c(0, 0.02)) +
  theme_classic(base_size = 14) +
  theme(
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(face = "bold"),
    panel.border  = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid    = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

# ── Barplot version ───────────────────────────────
# make Time discrete, keep numeric order
df_bar <- df_plot %>%
  mutate(Time = factor(Time, levels = sort(unique(Time))))

# one shared dodger for bars + error bars
pd <- position_dodge(width = 0.70)

ggplot(df_bar, aes(x = Time, y = Value, fill = Expression)) +
  geom_col(position = pd, width = 0.62, colour = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = Value - SEM, ymax = Value + SEM, group = Expression),  # <- key
    position = pd, width = 0.18, linewidth = 0.4, colour = "black"
  ) +
  facet_wrap(~ Condition, nrow = 1) +
  scale_fill_manual(values = col_map, breaks = names(col_map)) +
  labs(x = "Time after induction (min)",
       y = "NEO expression impact on ATP - Luminescence (RLU)") +
  scale_y_continuous(limits = c(0, 1e6), expand = c(0, 0), labels = sci_labeller) +
  theme_classic(base_size = 14) +
  theme(
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing  = unit(10, "pt"),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.6),
    strip.text       = element_text(face = "bold", size = 13),
    axis.title   = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 12),
    axis.ticks.length = unit(2.2, "mm"),
    legend.position = "top",
    legend.title    = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


########

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(colorspace)  # robust lightening while preserving hue
  library(cowplot)     # stacking + shared y-axis layout & legend extraction
  library(grid)        # textGrob for the y-axis title
})

# ==== INPUT ====
infile <- "~/Desktop/NEO_ATP_plot_data_20250626.txt"

# ==== scientific tick labels (bold c × 10^e) ====
sci_labeller <- function(x) {
  labs <- sapply(x, function(v) {
    if (v == 0) return("0")
    e <- floor(log10(v)); c <- v / 10^e
    sprintf("bold(%g %s 10^%d)", c, "%*%", e)
  })
  parse(text = labs)
}

# ==== load & reshape (replicate columns expected) ====
raw <- read.delim(infile, header = TRUE, check.names = FALSE)

keep_cols <- c("Sample", "Time",
               grep("^(EGFP|NEO|MS2L)_(Glucose|Arabinose)", names(raw), value = TRUE))
raw <- raw[, keep_cols, drop = FALSE]

df_raw <- raw %>%
  pivot_longer(
    cols      = -c(Sample, Time),
    names_to  = c("Expression", "Condition", "Replicate"),
    names_sep = "_",
    values_to = "RLU"
  ) %>%
  mutate(
    Replicate  = ifelse(is.na(Replicate) | Replicate == "", "R1", Replicate),
    Expression = factor(Expression, levels = c("EGFP","NEO","MS2L")),
    Condition  = factor(Condition,  levels = c("Glucose","Arabinose")),
    Time       = factor(Time, levels = sort(unique(Time)))
  )

# ==== Dark2 per gene; Glucose = lightened shade, Arabinose = base ====
dark2 <- brewer.pal(3, "Dark2")
expr_base <- c(EGFP = dark2[1],  # green
               NEO  = dark2[3],  # purple (as required)
               MS2L = dark2[2])  # orange

# strong lightening for clear contrast; increase amount for even lighter fills
GLU_LIGHTEN <- 0.65
df_raw <- df_raw %>%
  mutate(
    base_col = expr_base[as.character(Expression)],
    fill_col = ifelse(Condition == "Glucose",
                      colorspace::lighten(base_col, amount = GLU_LIGHTEN),
                      base_col)
  )

# ==== one panel builder (no titles; x shown only on bottom) ====
# one shared dodger (works even if some Time×Condition combos are missing)
pd <- position_dodge2(width = 0.70, preserve = "single")

make_panel <- function(expr, show_x = FALSE) {
  d <- df_raw %>% filter(Expression == expr)
  
  # per-panel colors: Arabinose = base (darker), Glucose = lightened base
  base_col <- expr_base[[as.character(expr)]]
  cond_cols <- c(
    Glucose   = colorspace::lighten(base_col, amount = GLU_LIGHTEN),
    Arabinose = base_col
  )
  
  ggplot(d, aes(x = Time, y = RLU, fill = Condition)) +
    # all replicate dots (including outliers) — correctly aligned with boxes
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.70),
      shape = 21, size = 1.1, stroke = 0.3, colour = "black", alpha = 1
    ) +
    # box summarizing the same data — dodged by Condition at each Time
    geom_boxplot(
      position = pd, width = 0.60,
      outlier.shape = NA, alpha = 1,
      colour = "black", linewidth = 0.5
    ) +
    scale_fill_manual(values = cond_cols, guide = "none") +  # no condition legend
    scale_y_continuous(limits = c(0, 1e6), expand = c(0, 0), labels = sci_labeller) +
    labs(x = if (show_x) "Time after induction (min)" else NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      axis.title   = element_text(face = "bold"),
      axis.text    = element_text(face = "bold", colour = "black"),
      axis.text.x  = if (show_x) element_text(face = "bold") else element_blank(),
      axis.ticks.x = if (show_x) element_line()             else element_blank(),
      plot.margin  = margin(t = 4, r = 10, b = 4, l = 10)
    )
}
genes  <- levels(df_raw$Expression)
panels <- lapply(seq_along(genes), function(i) make_panel(genes[i], show_x = i == length(genes)))

# ==== gene legend (Dark2 bases) ====
legend_df <- tibble(Expression = factor(genes, levels = genes), y = 1)
legend_plot <- ggplot(legend_df, aes(Expression, y, fill = Expression)) +
  geom_col(width = 0.65, colour = "black") +
  scale_fill_manual(values = expr_base, name = "Gene") +
  theme_void(base_size = 14) +
  theme(
    legend.position   = "top",
    legend.title      = element_text(face = "bold"),
    legend.text       = element_text(face = "bold"),
    legend.key.height = unit(10, "pt"),
    legend.key.width  = unit(22, "pt")
  )
gene_legend <- cowplot::get_legend(legend_plot)

# ==== stack panels and add a single centered Y-axis title (reliably) ====
# 1) stack panels (no legends in panels)
stacked <- plot_grid(plotlist = panels, ncol = 1, align = "v", axis = "lr",
                     rel_heights = rep(1, length(panels)))

# 2) put the legend above the stack
with_legend <- plot_grid(gene_legend, stacked, ncol = 1, rel_heights = c(0.14, 1))

# 3) add a real, centered y-axis label using a left column grob
y_grob <- textGrob("NEO expression impact on ATP – Luminescence (RLU)",
                   gp = gpar(fontsize = 14, fontface = "bold"),
                   rot = 90)
final_plot <- plot_grid(y_grob, with_legend, ncol = 2, rel_widths = c(0.06, 1))

# show
print(final_plot)

# ==== export ====
# ggsave("neo_atp_boxpoints_stacked_pub_final.pdf", final_plot,
#        width = 7.6, height = 8.6, units = "in", device = cairo_pdf)
# ggsave("neo_atp_boxpoints_stacked_pub_final.png", final_plot,
#        width = 7.6, height = 8.6, units = "in", dpi = 300)
