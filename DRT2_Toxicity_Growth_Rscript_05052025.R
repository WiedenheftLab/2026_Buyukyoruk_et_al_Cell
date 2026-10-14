library(ggplot2)
library(tidyverse)
library(ggpubr)

# data <- read.delim("~/Desktop/Toxicity_plate_reader_results_TM_muts.txt",header=TRUE)
# data <- read.delim("~/Desktop/Toxicity_plate_reader_results_TM_muts_combined.txt",header=TRUE)
# data <- read.delim("~/Desktop/Toxicity_plate_reader_results.txt",header=TRUE)
data <- read.delim("~/Desktop/Toxicity_plate_reader_results_combined_08192025.txt",header=TRUE) # all NEOs
# data <- read.delim("~/Desktop/Toxicity_plate_reader_results_ORF_w_Stop_combined_08212025.txt",header=TRUE) # ORFs w/ stop
# data <- read.delim("~/Desktop/DRT2_Cell_RTR/EGFP-NEO_tox_assays_for_revision.txt",header=TRUE)
# data <- read.delim("~/Desktop/DRT2_Cell_RTR/NEOc8.2_corrected_growth_curve_data.txt",header=TRUE)
samples <- unique(data$Sample)

for(i in 1:length(samples)){
  this_sample <- samples[i]
  df <- subset(data,Sample==this_sample)
  
  # 1) pivot averages long
  df_avg <- df %>% 
    filter(Type == "Average") %>% 
    pivot_longer(
      cols      = c(Glucose_0.2, Arabinose_0.2, Arabinose_0.002),
      names_to  = "Condition",
      values_to = "Value"
    )
  
  # 2) pivot SEM long
  df_sem <- df %>% 
    filter(Type == "SEM") %>% 
    pivot_longer(
      cols      = c(Glucose_0.2, Arabinose_0.2, Arabinose_0.002),
      names_to  = "Condition",
      values_to = "SEM"
    )
  
  # 3) join so each row has Value + SEM
  df_plot <- left_join(df_avg, df_sem,
                       by = c("Sample", "Time", "Condition"))
  
  # 4) detect the modal ΔTime per condition
  allowed_gaps <- df_plot %>%
    group_by(Condition) %>%
    arrange(Time) %>%
    summarize(
      allowed_gap = as.numeric(
        names(which.max(table(round(diff(Time), 6))))
      ),
      .groups = "drop"
    )
  
  # 5) compute true ΔTime and segment IDs (first delta = NA)
  df_seg <- df_plot %>%
    left_join(allowed_gaps, by = "Condition") %>%
    arrange(Condition, Time) %>%
    group_by(Condition) %>%
    mutate(
      delta   = Time - lag(Time),  
      segment = cumsum(!is.na(delta) & abs(delta - allowed_gap) > 1e-6)
    ) %>%
    ungroup()
  
  # 6) find the long‐gap boundaries
  gap_info <- df_seg %>%
    filter(!is.na(delta) & abs(delta - allowed_gap) > 1e-6) %>%
    mutate(
      start = Time - delta,
      end   = Time
    ) %>%
    slice(1)            # assume one induction event
  
  ind_start <- gap_info$start
  ind_end   <- gap_info$end
  break_times <- c(ind_start, ind_end)
  
  cond_labels <- c(
    Glucose_0.2     = "0.2% Glucose",
    Arabinose_0.002 = "0.002% Arabinose",
    Arabinose_0.2   = "0.2% Arabinose"
  )
  
  
  # 7) plot with shading and v‑lines at both boundaries
  p <- ggplot(df_seg,
         aes(x = Time,
             y = Value,
             color = Condition,
             fill  = Condition,
             group = interaction(Condition, segment))) +
    # shaded “Induced” region between start and end
    geom_rect(
      inherit.aes = FALSE,
      data        = tibble(xmin = ind_start, xmax = ind_end),
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill  = "gray70", alpha = 0.3
    ) +
    # dashed lines at both induction boundaries
    geom_vline(xintercept = break_times,
               linetype   = "dashed",
               color      = "gray40") +
    # ribbon for ±SEM
    geom_ribbon(aes(ymin = Value - SEM, ymax = Value + SEM),
                alpha  = 0.2, colour = NA) +
    # thin ±SEM boundary lines
    geom_line(aes(y = Value + SEM), size = 0.25) +
    geom_line(aes(y = Value - SEM), size = 0.25) +
    # main mean curve
    geom_line(size = 0.75) +
    # label “Induced” placed vertically in the shaded region
    annotate(
      "text",
      x     = (ind_start + ind_end) / 2,
      y     = 0.35, #max(df_seg$Value + df_seg$SEM)/2,
      label = "Induced",
      angle = 90,
      vjust = 0.5,
      fontface = "bold",
      color = "gray40"
    ) +
    # distinct palette, axis labels, theme
    scale_color_brewer(
      palette = "Dark2",
      name    = "Condition",    # legend title
      labels  = cond_labels      # legend entries
    ) +
    scale_fill_brewer(
      palette = "Dark2",
      name    = "Condition",
      labels  = cond_labels
    ) +
    labs(
      title = unique(df_plot$Sample),
      x     = "Time (hours)",
      y     = expression(bold(Growth~(OD[600]))),
      color = "Condition",
      fill  = "Condition"
    ) +   scale_y_continuous(
      limits = c(0, 0.7),  # lower = 0, upper = auto
      expand = c(0, 0)    # no padding below or above
    ) +
    scale_x_continuous(
      limits = c(0, 4),  # lower = 0, upper = auto
      expand = c(0, 0.05)    # no padding below or above
    ) + 
    theme_classic(base_size = 14) +
    theme(
      axis.title.x       = element_text(face = "bold"),
      axis.title.y       = element_text(face = "bold"),
      axis.text.x       = element_text(face = "bold"),
      axis.text.y       = element_text(face = "bold"),
      panel.border       = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90"),
      legend.position    = "right",
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 12)
    )  + guides(
      color = guide_legend(reverse = TRUE),
      fill  = guide_legend(reverse = TRUE)
    ) 
  assign(paste0("p_", i), p, envir = .GlobalEnv)
}

plot_list <- mget(paste0("p_", seq_along(samples)))

# paginate into chunks of up to 12 plots each
per_page   <- 12
plot_chunks <- split(plot_list, ceiling(seq_along(plot_list) / per_page))

for (page in seq_along(plot_chunks)) {
  chunk <- plot_chunks[[page]]
  
  combined <- ggarrange(
    plotlist      = chunk,
    ncol          = 3,
    nrow          = 4,          # FORCE 4 rows on every page
    common.legend = TRUE,
    legend        = "top"
  )
  
  print(combined)
  
  # optionally save:
  # ggsave(
  #   sprintf("page_%02d.pdf", page),
  #   combined,
  #   width  = 3 * 4,  # e.g. each panel 4" wide
  #   height = 3 * 4   # each panel 3" tall
  # )
}

