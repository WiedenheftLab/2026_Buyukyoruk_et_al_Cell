# Bar plot with SEM error bars, colored by condition (Dark2)
suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)

# Read data: either from CLI arg or interactively choose a file
# infile <- if (length(args) >= 1) args[1] else {
#   message("No input file given; opening file chooser...")
#   file.choose()
# }

infile <- "/Users/muratbuyukyoruk/Documents/DRT2/DRT2_biochem/FACS_plot.txt"

# Optional output path (png/pdf); if missing, plot to device
outfile <- if (length(args) >= 2) args[2] else NULL

# Read tab-delimited file (expects columns: Sample, condition, average, SEM)
df <- utils::read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)

# Basic sanity checks and coercions
required_cols <- c("Sample", "condition", "average", "SEM")
missing_cols  <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Input file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

df$Sample    <- factor(df$Sample,    levels = unique(df$Sample))
df$condition <- factor(df$condition, levels = unique(df$condition))
df$average   <- as.numeric(df$average)
df$SEM       <- as.numeric(df$SEM)

# Position settings
pd <- position_dodge(width = 0.55)
bar_width <- 0.5
err_width <- 0.2

p <- ggplot(df, aes(x = Sample, y = average, fill = condition)) +
  geom_col(
    aes(color = condition), # outline mapped to condition
    position = pd, 
    width = bar_width,
    alpha = 0.4,            # transparency of the fill
    linewidth = 0.6         # thickness of outline
  ) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none") + 
  geom_errorbar(
    aes(ymin = pmax(0, average - SEM), ymax = average + SEM),
    position = pd, width = err_width, linewidth = 0.6
  ) + 
  labs(
    x = "Samples",
    y = "% Cells with red/green ratio < 0.3",
    fill = "Condition",
    title = "Membrane Depolarization"
  ) +
  # expand_limits(y = 0) +
  theme_classic(base_size = 16) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  theme(
    axis.ticks.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0)
  )

# Write to file if provided; otherwise print
if (!is.null(outfile)) {
  ext <- tolower(tools::file_ext(outfile))
  if (ext %in% c("png", "jpg", "jpeg", "tiff")) {
    # Raster output
    if (ext == "png") {
      png(outfile, width = 1800, height = 1200, res = 200)
    } else if (ext %in% c("jpg", "jpeg")) {
      jpeg(outfile, width = 1800, height = 1200, res = 200, quality = 95)
    } else if (ext == "tiff") {
      tiff(outfile, width = 1800, height = 1200, res = 200, compression = "lzw")
    }
    print(p)
    dev.off()
  } else if (ext == "pdf") {
    pdf(outfile, width = 8, height = 5)
    print(p)
    dev.off()
  } else {
    warning("Unknown output extension '", ext, "'. Printing to screen instead.")
    print(p)
  }
} else {
  print(p)
}