# ============================================================
# PLOT 14
#
# SAME ALLUVIAL GEOMETRY AS ORIGINAL P14
#
# - Keep ALL clades
# - Keep ALL MMseqs2 clusters
# - Keep stratum widths/heights proportional to sequence counts
# - Keep small_cluster_cutoff
# - Keep EXACT same alluvial arrangement
#
# ONLY CHANGE:
# One facet per clade.
#
# In each facet:
#   selected clade ribbons = visible
#   all other ribbons       = transparent
#
# Because ALL observations remain in EVERY facet,
# ggalluvial preserves the same geometry in every panel.
# ============================================================

library(tidyverse)
library(ggalluvial)

small_cluster_cutoff <- 5

input_file <- "~/Desktop/DRT2_Cell_RTR/DRT2_NEO_MMSeq2_Clusters.txt"

data <- read.delim(
  input_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(data) <- c("sequence", "cluster")

data <- data %>%
  mutate(
    clade = str_extract(sequence, "Clade[0-9]+"),
    clade_num = as.integer(str_remove(clade, "Clade")),
    cluster_num = as.integer(str_remove(cluster, "Cluster"))
  )

# ------------------------------------------------------------
# 1. MMseqs2 cluster sizes
# ------------------------------------------------------------

cluster_sizes <- data %>%
  count(
    cluster,
    cluster_num,
    name = "cluster_size"
  )


# ------------------------------------------------------------
# 2. Collapse small MMseqs2 clusters
#
# SAME behavior as original P14
# ------------------------------------------------------------

small_cluster_label <- paste0(
  "Clusters < ",
  small_cluster_cutoff,
  " NEOs"
)


alluvial_base <- data %>%
  
  left_join(
    cluster_sizes,
    by = c(
      "cluster",
      "cluster_num"
    )
  ) %>%
  
  mutate(
    
    cluster_display = if_else(
      cluster_size < small_cluster_cutoff,
      small_cluster_label,
      cluster
    )
    
  ) %>%
  
  count(
    clade,
    clade_num,
    cluster_display,
    name = "n"
  )


# ------------------------------------------------------------
# 3. Numerical clade order
# ------------------------------------------------------------

clade_order <- alluvial_base %>%
  distinct(
    clade,
    clade_num
  ) %>%
  arrange(clade_num) %>%
  pull(clade)


alluvial_base <- alluvial_base %>%
  mutate(
    clade = factor(
      clade,
      levels = clade_order
    )
  )


# ------------------------------------------------------------
# 4. Cluster order
#
# Keep the same basic ordering as your original P14:
# retained clusters numerically, followed by <5 category.
# ------------------------------------------------------------

retained_clusters <- cluster_sizes %>%
  
  filter(
    cluster_size >= small_cluster_cutoff
  ) %>%
  
  arrange(
    cluster_num
  ) %>%
  
  pull(
    cluster
  )


cluster_order <- c(
  retained_clusters,
  small_cluster_label
)


alluvial_base <- alluvial_base %>%
  mutate(
    cluster_display = factor(
      cluster_display,
      levels = cluster_order
    )
  )


# ------------------------------------------------------------
# 5. DUPLICATE THE COMPLETE DATASET FOR EVERY CLADE
#
# THIS is the critical part.
#
# Every facet receives the COMPLETE alluvial dataset.
#
# Therefore:
#   - Clade widths remain identical
#   - Cluster widths remain identical
#   - All strata remain identical
#   - All ribbon positions remain identical
#
# Only ribbon visibility changes.
# ------------------------------------------------------------

facet_data <- crossing(
  panel_clade = clade_order,
  alluvial_base
) %>%
  
  mutate(
    
    panel_clade = factor(
      panel_clade,
      levels = clade_order
    ),
    
    # Selected clade = fully visible
    #
    # Everything else = essentially invisible
    ribbon_alpha = if_else(
      as.character(clade) ==
        as.character(panel_clade),
      0.85,
      0.00
    )
    
  )


# ------------------------------------------------------------
# 6. PLOT
# ------------------------------------------------------------

p14 <- ggplot(
  facet_data,
  aes(
    axis1 = clade,
    axis2 = cluster_display,
    y = n
  )
) +
  
  # ----------------------------------------------------------
# ALL ribbons are still physically present.
#
# Only selected clade has alpha > 0.
#
# This preserves ggalluvial geometry.
# ----------------------------------------------------------

geom_alluvium(
  aes(
    fill = clade,
    alpha = ribbon_alpha
  ),
  width = 1 / 12
) +
  
  
  # ----------------------------------------------------------
# Original strata
#
# Heights remain proportional to number of sequences.
# ----------------------------------------------------------

geom_stratum(
  width = 1 / 12,
  fill = "grey92",
  color = "grey40",
  linewidth = 0.25
) +
  
  
  # ----------------------------------------------------------
# Original labels
# ----------------------------------------------------------

geom_text(
  stat = "stratum",
  aes(
    label = after_stat(stratum)
  ),
  size = 2.5
) +
  
  
  # ----------------------------------------------------------
# Alpha
#
# identity means use 0 / 0.85 literally.
# ----------------------------------------------------------

scale_alpha_identity() +
  
  
  # ----------------------------------------------------------
# X-axis
# ----------------------------------------------------------

scale_x_discrete(
  limits = c(
    "Previous NEO clade",
    "MMseqs2 cluster"
  ),
  expand = c(
    0.12,
    0.12
  )
) +
  
  
  # ----------------------------------------------------------
# One panel per clade
# ----------------------------------------------------------

facet_wrap(
  ~ panel_clade,
  ncol = 4
) +
  
  
  # ----------------------------------------------------------
# No fill legend needed
# ----------------------------------------------------------

guides(
  fill = "none"
) +
  
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_size = 11
) +
  
  theme(
    
    axis.text.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    axis.line.y = element_blank(),
    
    axis.text.x = element_text(
      size = 8,
      face = "bold"
    ),
    
    axis.ticks.x = element_blank(),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    panel.spacing = unit(
      1.2,
      "lines"
    )
    
  )


print(p14)