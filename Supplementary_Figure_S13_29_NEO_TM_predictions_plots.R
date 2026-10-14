library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(gggenes)
library(patchwork)

# read full dataset once
df_all <- read.delim(
  "~/Desktop/NEO_TM_predictions.txt", header = TRUE,
  col.names = c("pos","aa","inner","outer","membrane","signal")
)

all_segments <- tibble(
  molecule    = character(),
  gene        = character(),
  start       = integer(),
  end         = integer(),
  orientation = integer()
)

# loop over each unique sequence identifier
seqs <- unique(df_all$aa)
for (i in seq_along(seqs)) {
  this_seq <- seqs[i]
  df       <- subset(df_all, aa == this_seq)
  
  # 1) melt out the three line–features (inner, outer, signal)
  df_lines <- df %>% 
    select(pos, inner, outer, signal) %>% 
    pivot_longer(-pos, names_to = "feature", values_to = "prob") 
  
  # 2) keep the membrane probabilities for bars
  df_bars  <- df %>% select(pos, membrane)
  
  # 3) compute the bottom track
  df_track <- df %>% 
    mutate(
      maxfeat = pmap_chr(
        list(inner, outer, membrane, signal),
        ~c(inner = ..1, outer = ..2, membrane = ..3, signal = ..4) %>%
          which.max() %>% names
      )
    )
  
  probability_cutoff <- 0.5
  
  tm_segments <- df %>% 
    mutate(is_tm = membrane > probability_cutoff) %>% 
    filter(is_tm) %>% 
    arrange(pos) %>% 
    mutate(gap = c(TRUE, diff(pos) != 1),
           seg = cumsum(gap)) %>% 
    group_by(seg) %>% 
    summarize(start = min(pos), end = max(pos)) %>% 
    ungroup()
  
  tm_segments$gene <- unique(df$aa)
  
  # 4) define your colors
  dark2     <- brewer.pal(4, "Dark2")
  line_cols <- c(inner  = dark2[3],
                 outer  = dark2[1],
                 signal = "black")
  bar_fill  <- "purple3"
  track_cols <- c(membrane = "purple3")
  
  # 5a) now plot main
  main <- ggplot() +
    geom_hline(yintercept = probability_cutoff, linetype = "dashed", color = "black") +
    geom_col(data = df_bars,
             aes(x = pos, y = membrane, fill = "membrane"),
             width = 0.75,
             color = "black", size = 0) +
    geom_line(data = subset(df_lines, feature == "signal"),
              aes(x = pos, y = prob, color = feature),
              size = 0.75) +
    geom_tile(data = df_track,
              aes(x = pos, y = -0.05, fill = maxfeat),
              height = 0.05) +
    scale_color_manual(
      name   = "Feature",
      values = line_cols,
      labels = c(inner = "Inner",
                 outer = "Outer",
                 signal = "Signal peptide")
    ) +
    scale_fill_manual(
      name   = "Feature",
      values = track_cols,
      labels = c(membrane = "Transmembrane")
    ) +
    labs(y = "Probability", title = unique(df$aa)) +
    scale_y_continuous(limits = c(-0.05, 1), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
    theme_classic(base_size = 14) +
    theme(
      text         = element_text(family = "Arial"),
      axis.line      = element_blank(),
      axis.title     = element_text(face = "bold", size = 10),
      axis.title.x   = element_blank(),
      axis.text.x    = element_blank(),
      axis.ticks.x   = element_blank(),
      axis.text      = element_text(face = "bold"),
      panel.border   = element_rect(colour = "black", fill = NA, size = 1),
      panel.grid     = element_blank(),
      plot.title     = element_text(face = "bold", hjust = 0.5, size=14),
      plot.margin    = margin(b = 0, t = 3, r = 5)
    )
  
  # 5b) gene panel
  n <- max(df$pos)
  if(i == 1){
    seg_df <- tibble(
      gene   = unique(df$aa),
      start  = c(0),
      end    = c(n),
      name   = c("MS2 Lysis"),
      strand = "+"
    )
  } else {
    seg_df <- tibble(
      gene   = unique(df$aa),
      start  = c(floor(2*n/3)-(n/20), floor(n/3)-(n/20), 0),
      end    = c(n, floor(2*n/3)+(n/40), floor(n/3)+(n/40)),
      name   = c("Repeat 3","Repeat 2","Repeat 1"),
      strand = "+"
    )
  }
  
  gene <- ggplot() +
    geom_gene_arrow(data = seg_df,
                    aes(xmin = start, xmax = end, y = gene),
                    arrowhead_height  = unit(4, "mm"),
                    arrowhead_width   = unit(2, "mm"),
                    arrow_body_height = unit(4, "mm"),
                    fill  = "gray90", color = "black",
                    size  = 0.5) +
    geom_gene_label(data = seg_df,
                    aes(xmin = start, xmax = end, y = gene, label = name, fontface = "bold"),
                    height = grid::unit(4, "mm")) +
    geom_gene_arrow(data = tm_segments,
                    aes(xmin = start, xmax = end, y = gene, forward = TRUE),
                    arrowhead_height  = unit(0, "mm"),
                    arrowhead_width   = unit(0, "mm"),
                    arrow_body_height = unit(4, "mm"),
                    fill               = bar_fill,
                    color              = "black",
                    alpha              = 1) +
    scale_x_continuous(limits = c(1, n), expand = c(0, 0)) +
    theme(
      text         = element_text(family = "Arial"),
      axis.title    = element_text(face = "bold", size = 10),
      axis.text     = element_text(face = "bold", size = 12),
      panel.border  = element_blank(), # element_rect(colour = "black", fill = NA, size = 1),
      panel.grid    = element_blank(),
      panel.background = element_blank(),
      axis.text.y   = element_blank(),
      axis.title.y  = element_blank(),
      axis.ticks.y  = element_blank(),
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.margin   = margin(t = 0, r = 5)
    ) +
    scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
    labs(x = "Position (aa)")
  
  # 5c) combine
  combined <- main / gene + plot_layout(ncol = 1, heights = c(4, 0.75))
  
  # 6) save
  assign(paste0("main_",     i), main)
  assign(paste0("gene_",     i), gene)
  assign(paste0("combined_", i), combined)
  
  ann_seg <- seg_df %>%
    transmute(
      molecule    = this_seq,
      gene        = name,      # your seg_df$name
      start,
      end,
      orientation = 1
    )
  ann_tm  <- tm_segments %>%
    transmute(
      molecule    = this_seq,
      gene        = "TM",
      start,
      end,
      orientation = 1
    )
  all_segments <- bind_rows(all_segments, ann_seg, ann_tm)
  
}
# collect all combined_i into a list
combined_list <- mget(paste0("combined_", seq_along(seqs)))

# split into pages of up to 3
per_page   <- 15
plot_pages <- split(combined_list, ceiling(seq_along(combined_list)/per_page))

for(page in seq_along(plot_pages)) {
  chunk <- plot_pages[[page]]
  
  # if fewer than 3 plots on this page, pad with blank spaces
  if(length(chunk) < per_page) {
    chunk <- c(chunk, rep(list(plot_spacer()), per_page - length(chunk)))
  }
  
  # stack them in one column (3 rows guaranteed)
  p <- wrap_plots(chunk, ncol = 3) +
    plot_layout(guides = "collect") &    # if you want to collect legends
    theme(legend.position = "bottom")
  
  print(p)
}

all_segments
