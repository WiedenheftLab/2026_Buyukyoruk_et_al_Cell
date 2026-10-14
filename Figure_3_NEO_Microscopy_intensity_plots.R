suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)    # for trans_breaks, trans_format, math_format
  library(tools)     # for file_path_sans_ext()
})

# ─── A) LIST INPUT FILES ───────────────────────────────────────────────────────
setwd("/Users/muratbuyukyoruk/Desktop/RT_projects/DRT2_data/Microscopy_intensity_plots")
files <- list.files(pattern = "\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("No .txt files found in working dir")

# ─── B) PROCESS ONE FILE ──────────────────────────────────────────────────────
process_one <- function(infile) {
  sample <- file_path_sans_ext(basename(infile))
  raw    <- read.delim(infile, check.names = FALSE)
  nc     <- ncol(raw); stopifnot(nc %% 2 == 0)
  ncells <- nc/2
  colnames(raw) <- unlist(lapply(seq_len(ncells), function(i){
    c(paste0("Distance_cell", i), paste0("Gray_cell", i))
  }))
  
  long <- raw %>%
    mutate(rowid = row_number()) %>%
    pivot_longer(
      cols          = -rowid,
      names_to      = c(".value", "cell"),
      names_pattern = "(.*)_cell(\\d+)"
    ) %>%
    select(-rowid) %>%
    group_by(cell) %>%
    mutate(
      mid          = (min(Distance, na.rm = TRUE) + max(Distance, na.rm = TRUE)) / 2,
      CenteredDist = Distance - mid
    ) %>%
    ungroup() %>%
    mutate(
      Bin    = round(CenteredDist / 0.025) * 0.025,
      sample = sample
    )
  
  summary_raw <- long %>%
    group_by(sample, Bin) %>%
    summarize(
      n_cells   = sum(!is.na(Gray)),
      mean_gray = mean(Gray, na.rm = TRUE),
      sd_gray   = sd(Gray,   na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    arrange(sample, Bin)
  
  fit_mean <- loess(mean_gray ~ Bin, data = summary_raw, span = 0.15)
  fit_sd   <- loess(sd_gray   ~ Bin, data = summary_raw, span = 0.15)
  sm       <- summary_raw %>%
    mutate(
      mean_s = predict(fit_mean, newdata = summary_raw),
      sd_s   = predict(fit_sd,   newdata = summary_raw)
    )
  
  list(long       = long,
       summary_raw = summary_raw,
       smooth      = sm)
}

# ─── C) PROCESS ALL FILES ─────────────────────────────────────────────────────
all_data        <- lapply(files, process_one)
long_all        <- bind_rows(lapply(all_data, `[[`, "long"))
summary_raw_all <- bind_rows(lapply(all_data, `[[`, "summary_raw"))
smoothed_all    <- bind_rows(lapply(all_data, `[[`, "smooth"))

# ─── D) COMMON AXIS BREAKS ────────────────────────────────────────────────────
global_half <- max(abs(long_all$Bin), na.rm = TRUE)
x_breaks   <- seq(-global_half, global_half, by = 0.2)

# ─── E) OVERLAID RAW‐BINNED PLOT (log₁₀ y) ────────────────────────────────────
p_ori <- ggplot(summary_raw_all,
                aes(x = Bin, y = mean_gray,
                    color = sample, fill = sample, group = sample)) +
  geom_ribbon(aes(ymin = pmax(mean_gray - sd_gray, 1e-3),  # floor to avoid ≤0
                  ymax = mean_gray + sd_gray),
              alpha = 0.2, size = 0.25) +
  geom_line(size = 0.5) +
  scale_x_continuous(breaks = x_breaks,
                     limits = c(-global_half, global_half),
                     expand = c(0, 0)) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = expression(Distance~"["*mu*m*"]"),
       y = "Mean pixel intensity (log10)") +
  theme_classic(base_size = 12) +
  theme(
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks  = element_line(color = "black"),
    axis.text   = element_text(color = "black"),
    axis.title  = element_text(color = "black")
  )

# ─── F) OVERLAID SMOOTHED PLOT (log₁₀ y) ─────────────────────────────────────
p_smooth <- ggplot(smoothed_all,
                   aes(x = Bin, y = mean_s,
                       color = sample, fill = sample, group = sample)) +
  geom_ribbon(aes(ymin = pmax(mean_s - sd_s, 1e-3),
                  ymax = mean_s + sd_s),
              alpha = 0.2, size = 0.25) +
  geom_line(size = 0.5) +
  scale_x_continuous(breaks = x_breaks,
                     limits = c(-global_half, global_half),
                     expand = c(0, 0)) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = expression(Distance~"["*mu*m*"]"),
       y = "Mean pixel intensity (log10)") +
  theme_classic(base_size = 12) +
  theme(
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks  = element_line(color = "black"),
    axis.text   = element_text(color = "black"),
    axis.title  = element_text(color = "black")
  )

# ─── G) INDIVIDUAL RAW‐VALUES DEBUG PLOTS ────────────────────────────────────
for (smpl in unique(long_all$sample)) {
  df_s <- filter(long_all, sample == smpl)
  p <- ggplot(df_s,
              aes(x = Bin,
                  y = Gray,
                  group = cell,
                  color = cell)) +
    geom_line(size = 0.8, alpha = 0.8) +
    scale_x_continuous(breaks = x_breaks,
                       limits = c(-global_half, global_half),
                       expand = c(0, 0)) +
    labs(
      x     = expression(Distance~"["*mu*m*"]"),
      y     = "Gray value",
      title = paste("Raw values:", smpl)
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black")
    )
  print(p)
}

# ─── H) DRAW THE COMBINED PROFILES ────────────────────────────────────────────
print(p_ori)
print(p_smooth)