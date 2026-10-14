############################################################
## Volcano plot for EGFP vs Neo – fully fixed version
############################################################

library(ggplot2)
library(ggrepel)
library(scales)

## 1) Read and subset data ---------------------------------
dat <- read.delim("~/Desktop/metabolite_data.tsv",
                  header = TRUE,
                  sep = "\t",
                  check.names = FALSE,
                  stringsAsFactors = FALSE)

dat <- dat[dat$Label %in% c("EGFP", "Neo"), ]


## 2) Compute log2FC and p values --------------------------
metabolite_names <- setdiff(colnames(dat), "Label")

calc_stats <- function(colname) {
  egfp <- dat[dat$Label == "EGFP", colname]
  neo  <- dat[dat$Label == "Neo",  colname]
  
  tt <- t.test(neo, egfp)
  
  log2FC <- log2(mean(neo,  na.rm = TRUE) /
                   mean(egfp, na.rm = TRUE))
  pval   <- tt$p.value
  
  c(log2FC = log2FC, pval = pval)
}

stats_mat <- sapply(metabolite_names, calc_stats)
res <- as.data.frame(t(stats_mat))
res$Metabolite <- rownames(res)
res$negLog10P  <- -log10(res$pval)


## 2b) Dynamically determine top p-value bin ---------------
pmax <- max(res$pval, na.rm = TRUE)

# Select a reasonable top cutoff:
if (pmax > 0.5) {
  top_bin <- 1.0
} else if (pmax > 0.1) {
  top_bin <- 0.5
} else if (pmax > 0.05) {
  top_bin <- 0.1
} else {
  top_bin <- 0.05
}

# Dynamic breaks:
p_breaks <- c(0, 0.001, 0.01, 0.05, 0.1, top_bin)

# Corresponding labels:
p_labels <- c("0.001", "0.01", "0.05", "0.1", as.character(top_bin))

# Assign bins
res$p_bin <- cut(
  res$pval,
  breaks = p_breaks,
  labels = p_labels,
  include.lowest = TRUE
)

# Ensure correct legend order
res$p_bin <- factor(res$p_bin,
                    levels = rev(p_labels),
                    ordered = TRUE)


## 3) Dot sizes (scaled larger for visibility) --------------
size_values <- seq(2, 6, length.out = length(p_labels))
size_values <- rev(size_values)        # <-- IMPORTANT LINE
names(size_values) <- p_labels


## 4) Plot thresholds ---------------------------------------
p_cutoff   <- 0.05
fc_cutoff1 <- 1
fc_cutoff2 <- -1


## 5) Volcano plot ------------------------------------------
p <- ggplot(res, aes(x = log2FC, y = negLog10P)) +
  
  # Labels BELOW points
  geom_text_repel(
    aes(label = Metabolite),
    colour        = "black",
    size          = 4,
    max.overlaps  = Inf,
    box.padding   = 0.5,
    point.padding = 0.4,
    segment.size  = 0.2,
    show.legend   = FALSE
  ) +
  
  # Points ABOVE labels
  geom_point(
    aes(fill = log2FC, size = p_bin),
    shape  = 21,
    colour = "black",
    stroke = 0.4
  ) +
  
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_vline(xintercept = fc_cutoff1,       linetype = "dashed") +
  geom_vline(xintercept = fc_cutoff2,       linetype = "dashed") +
  
  scale_fill_gradientn(
    colours = c("#2166AC", "gray", "#B2182B"),
    values  = rescale(c(min(res$log2FC, -6), 0, max(res$log2FC, 2))),
    limits  = c(-6, 2),
    breaks  = c(-6, -4, -2, 0, 2),
    name    = "Log2(FC)",
    guide   = guide_colourbar(order = 2)
  ) +
  
  scale_size_manual(
    name   = "P-value",
    values = size_values,
    breaks = rev(p_labels),
    labels = rev(p_labels),
    guide  = guide_legend(
      order = 1,
      override.aes = list(
        shape  = 21,
        fill   = NA,
        colour = "black",
        stroke = 0.4
      )
    )
  ) +
  
  labs(
    x = "log2(FC) (Neo / EGFP)",
    y = "-log10(p-value)"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(colour = "black"),
    legend.position = "right"
  ) +
  guides(
    fill = guide_colourbar(order = 1),
    size = guide_legend(order = 2,
                        override.aes = list(fill = NA))
  )

print(p)