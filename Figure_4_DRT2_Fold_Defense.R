# ============================================================
# DRT2 plaque assay replicate-matched fold-defense bar plot
#
# Normalization:
#   Each sample replicate is normalized to the corresponding
#   MG1655 replicate from the same phage.
#
#   Replicate 1 sample / MG1655 replicate 1
#   Replicate 2 sample / MG1655 replicate 2
#   Replicate 3 sample / MG1655 replicate 3
#
# Fold defense:
#   MG1655 PFU/mL / sample PFU/mL
#
# MG1655 therefore equals 1-fold defense.
# On a log10 axis, 1-fold is located at 10^0,
# which corresponds mathematically to log10(1) = 0.
#
# Rows with missing plaque_count or dilution_factor are ignored.
#
# Facets:
#   rows    = phage
#   columns = condition
#
# Plot:
#   y-axis = unchanged log10 fold-defense axis
#   bars   = mean fold defense
#   points = individual replicate-matched values
#   error  = mean ± SEM
# ============================================================


# -----------------------------
# Input file
# -----------------------------

input_file <- "~/Desktop/DRT2_plaque_assay_counts_paper.csv"


# -----------------------------
# Load packages
# -----------------------------

required_packages <- c(
  "tidyverse",
  "scales"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(scales)


# -----------------------------
# Confirm that the file exists
# -----------------------------

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file was not found:\n",
      input_file
    )
  )
}


# -----------------------------
# Read input file
# Automatically detects tab- or comma-separated data
# -----------------------------

first_line <- readLines(
  input_file,
  n = 1,
  warn = FALSE
)

if (grepl("\t", first_line)) {
  
  raw_data <- read.delim(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "N/A", "na", "n/a")
  )
  
} else {
  
  raw_data <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "N/A", "na", "n/a")
  )
}


# -----------------------------
# Standardize column names
#
# spot_volume_uL becomes spot_volume_ul
# -----------------------------

colnames(raw_data) <- colnames(raw_data) %>%
  trimws() %>%
  tolower() %>%
  gsub("[[:space:]]+", "_", .) %>%
  gsub("[^a-z0-9_]", "", .)


# -----------------------------
# Required columns
# -----------------------------

required_columns <- c(
  "phage",
  "sample",
  "replicate",
  "plaque_count",
  "dilution_factor",
  "spot_volume_ul",
  "condition"
)

missing_columns <- setdiff(
  required_columns,
  colnames(raw_data)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing required column(s):\n",
      paste0(
        "  - ",
        missing_columns,
        collapse = "\n"
      ),
      "\n\nColumns found in the input file:\n",
      paste0(
        "  - ",
        colnames(raw_data),
        collapse = "\n"
      )
    )
  )
}


# -----------------------------
# Text-cleaning function
# -----------------------------

clean_text <- function(x) {
  
  x <- as.character(x)
  
  # Non-breaking spaces
  x <- gsub(
    "\u00A0",
    " ",
    x,
    fixed = TRUE
  )
  
  # En dash, em dash, and Unicode minus
  x <- gsub(
    "[\u2013\u2014\u2212]",
    "-",
    x
  )
  
  # Repeated spaces
  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )
  
  x <- trimws(x)
  
  x[x == ""] <- NA_character_
  
  x
}


# -----------------------------
# Remove completely empty spreadsheet rows
# -----------------------------

raw_data <- raw_data %>%
  mutate(
    across(
      everything(),
      ~ {
        value <- trimws(as.character(.x))
        value[value == ""] <- NA_character_
        value
      }
    )
  ) %>%
  filter(
    !if_all(
      everything(),
      is.na
    )
  )


# -----------------------------
# Clean the imported data
# -----------------------------

raw_data <- raw_data %>%
  mutate(
    source_row = row_number() + 1,
    
    phage = clean_text(phage),
    sample = clean_text(sample),
    replicate = clean_text(replicate),
    condition = clean_text(condition),
    
    plaque_count = suppressWarnings(
      as.numeric(plaque_count)
    ),
    
    dilution_factor = suppressWarnings(
      as.numeric(dilution_factor)
    ),
    
    spot_volume_ul = suppressWarnings(
      as.numeric(spot_volume_ul)
    ),
    
    condition = case_when(
      tolower(condition) == "control" ~ "Control",
      tolower(condition) == "sample" ~ "Sample",
      TRUE ~ condition
    )
  )


# -----------------------------
# Validate descriptive columns
#
# Missing plaque count and dilution factor are allowed.
# Those rows represent unfinished replicates and are ignored.
# -----------------------------

invalid_description_rows <- raw_data %>%
  filter(
    is.na(phage) |
      is.na(sample) |
      is.na(replicate) |
      is.na(condition) |
      is.na(spot_volume_ul)
  )

if (nrow(invalid_description_rows) > 0) {
  
  print(
    invalid_description_rows %>%
      select(
        source_row,
        phage,
        sample,
        replicate,
        plaque_count,
        dilution_factor,
        spot_volume_ul,
        condition
      )
  )
  
  stop(
    paste0(
      nrow(invalid_description_rows),
      " row(s) are missing phage, sample, replicate, ",
      "condition, or spot volume."
    )
  )
}


# -----------------------------
# Sample order and displayed labels
# -----------------------------

sample_order <- c(
  "MG1655",
  "DRT2 - WT",
  "DRT2 - YCAA",
  "DRT2 - I40*",
  "DRT2 - TM mut",
  "DRT2 - Wilkinson Control"
)

sample_labels <- c(
  "MG1655" = "MG1655",
  "DRT2 - WT" = "WT",
  "DRT2 - YCAA" = "YCAA",
  "DRT2 - I40*" = "I40Stop",
  "DRT2 - TM mut" = "TM mut",
  "DRT2 - Wilkinson Control" = "Wilkinson\ncontrol"
)


# -----------------------------
# Validate sample names
# -----------------------------

unknown_sample_rows <- raw_data %>%
  filter(
    !sample %in% sample_order
  ) %>%
  distinct(
    source_row,
    sample
  )

if (nrow(unknown_sample_rows) > 0) {
  
  print(unknown_sample_rows)
  
  stop(
    paste0(
      "Unrecognized sample name(s) were found.\n\n",
      "Allowed sample names:\n",
      paste0(
        "  - ",
        sample_order,
        collapse = "\n"
      )
    )
  )
}


# -----------------------------
# Validate condition names
# -----------------------------

unknown_condition_rows <- raw_data %>%
  filter(
    !condition %in% c(
      "Control",
      "Sample"
    )
  ) %>%
  distinct(
    source_row,
    condition
  )

if (nrow(unknown_condition_rows) > 0) {
  
  print(unknown_condition_rows)
  
  stop(
    "Condition values must be Control or Sample."
  )
}


# -----------------------------
# Convert variables to ordered factors
# -----------------------------

raw_data <- raw_data %>%
  mutate(
    sample = factor(
      sample,
      levels = sample_order
    ),
    
    condition = factor(
      condition,
      levels = c(
        "Control",
        "Sample"
      )
    ),
    
    phage = factor(
      phage,
      levels = unique(phage)
    )
  )


# -----------------------------
# Exclude unfinished replicate rows
#
# A row is included only when both plaque_count
# and dilution_factor are available.
# -----------------------------

excluded_rows <- raw_data %>%
  filter(
    is.na(plaque_count) |
      is.na(dilution_factor)
  )

if (nrow(excluded_rows) > 0) {
  message(
    nrow(excluded_rows),
    " unfinished replicate row(s) were ignored."
  )
}

analysis_data <- raw_data %>%
  filter(
    !is.na(plaque_count),
    !is.na(dilution_factor)
  )


# -----------------------------
# Validate completed numerical values
# -----------------------------

invalid_numeric_rows <- analysis_data %>%
  filter(
    plaque_count <= 0 |
      dilution_factor <= 0 |
      spot_volume_ul <= 0
  )

if (nrow(invalid_numeric_rows) > 0) {
  
  print(
    invalid_numeric_rows %>%
      select(
        source_row,
        phage,
        sample,
        replicate,
        plaque_count,
        dilution_factor,
        spot_volume_ul,
        condition
      )
  )
  
  stop(
    paste0(
      "Completed rows must have plaque_count, dilution_factor, ",
      "and spot_volume_ul greater than zero."
    )
  )
}


# -----------------------------
# Check for duplicate rows
#
# There must be only one completed measurement for each:
# phage + sample + replicate
# -----------------------------

duplicate_measurements <- analysis_data %>%
  count(
    phage,
    sample,
    replicate,
    name = "number_of_rows"
  ) %>%
  filter(
    number_of_rows > 1
  )

if (nrow(duplicate_measurements) > 0) {
  
  print(duplicate_measurements)
  
  stop(
    paste0(
      "Duplicate completed measurements were found for the same ",
      "phage, sample, and replicate."
    )
  )
}


# -----------------------------
# Calculate PFU/mL
#
# PFU/mL =
# plaque count /
# dilution factor /
# spot volume in mL
# -----------------------------

pfu_data <- analysis_data %>%
  mutate(
    spot_volume_ml = spot_volume_ul / 1000,
    
    pfu_per_ml =
      plaque_count /
      dilution_factor /
      spot_volume_ml
  )


# -----------------------------
# Build replicate-specific MG1655 reference
#
# The reference retains one PFU/mL value for each:
#   phage + replicate
#
# Condition is deliberately not included because MG1655
# is the common reference for samples in both facet columns.
# -----------------------------

mg1655_reference <- pfu_data %>%
  filter(
    sample == "MG1655"
  ) %>%
  select(
    phage,
    replicate,
    mg1655_pfu_per_ml = pfu_per_ml
  )


# -----------------------------
# Ensure MG1655 has only one measurement
# for each phage and replicate
# -----------------------------

duplicate_mg1655_references <- mg1655_reference %>%
  count(
    phage,
    replicate,
    name = "number_of_mg1655_rows"
  ) %>%
  filter(
    number_of_mg1655_rows > 1
  )

if (nrow(duplicate_mg1655_references) > 0) {
  
  print(duplicate_mg1655_references)
  
  stop(
    paste0(
      "More than one MG1655 measurement exists for the same ",
      "phage and replicate."
    )
  )
}


# -----------------------------
# Match every sample replicate to the corresponding
# MG1655 replicate
# -----------------------------

fold_data <- pfu_data %>%
  left_join(
    mg1655_reference,
    by = c(
      "phage",
      "replicate"
    )
  )


# -----------------------------
# Report completed sample rows that have no matching
# completed MG1655 replicate
# -----------------------------

missing_mg1655_matches <- fold_data %>%
  filter(
    is.na(mg1655_pfu_per_ml)
  ) %>%
  select(
    source_row,
    phage,
    sample,
    replicate,
    condition,
    pfu_per_ml
  )

if (nrow(missing_mg1655_matches) > 0) {
  
  print(missing_mg1655_matches)
  
  stop(
    paste0(
      "One or more completed sample replicates do not have a ",
      "completed MG1655 measurement with the same phage and replicate."
    )
  )
}


# -----------------------------
# Calculate replicate-matched fold defense
#
# MG1655 replicate 1 / sample replicate 1
# MG1655 replicate 2 / sample replicate 2
# etc.
#
# MG1655 divided by itself equals 1.
# -----------------------------

fold_data <- fold_data %>%
  mutate(
    fold_defense =
      mg1655_pfu_per_ml /
      pfu_per_ml
  )


# -----------------------------
# Validate calculated values
# -----------------------------

invalid_fold_rows <- fold_data %>%
  filter(
    !is.finite(fold_defense) |
      fold_defense <= 0
  )

if (nrow(invalid_fold_rows) > 0) {
  
  print(
    invalid_fold_rows %>%
      select(
        source_row,
        phage,
        condition,
        sample,
        replicate,
        pfu_per_ml,
        mg1655_pfu_per_ml,
        fold_defense
      )
  )
  
  stop(
    "Invalid fold-defense values were calculated."
  )
}


# -----------------------------
# Verify that completed MG1655 rows equal 1
# -----------------------------

mg1655_check <- fold_data %>%
  filter(
    sample == "MG1655"
  )

if (
  any(
    abs(mg1655_check$fold_defense - 1) >
    sqrt(.Machine$double.eps)
  )
) {
  stop(
    paste0(
      "MG1655 normalization check failed. ",
      "MG1655 fold-defense values should equal 1."
    )
  )
}


# -----------------------------
# Print replicate-matched calculations
# -----------------------------

replicate_output <- fold_data %>%
  arrange(
    phage,
    condition,
    sample,
    replicate
  ) %>%
  select(
    phage,
    condition,
    sample,
    replicate,
    plaque_count,
    dilution_factor,
    pfu_per_ml,
    mg1655_pfu_per_ml,
    fold_defense
  )

print(replicate_output)


# -----------------------------
# Summarize mean ± SEM
#
# The mean and SEM are calculated after
# replicate-matched normalization.
# -----------------------------

summary_data <- fold_data %>%
  group_by(
    phage,
    condition,
    sample
  ) %>%
  summarise(
    mean_fold_defense = mean(
      fold_defense,
      na.rm = TRUE
    ),
    
    sd_fold_defense = sd(
      fold_defense,
      na.rm = TRUE
    ),
    
    n = sum(
      !is.na(fold_defense)
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    sem_fold_defense = case_when(
      n > 1 ~ sd_fold_defense / sqrt(n),
      TRUE ~ 0
    ),
    
    errorbar_ymin =
      mean_fold_defense -
      sem_fold_defense,
    
    errorbar_ymax =
      mean_fold_defense +
      sem_fold_defense
  )


# -----------------------------
# Build condition-specific x positions
# -----------------------------

x_positions <- raw_data %>%
  distinct(
    condition,
    sample
  ) %>%
  mutate(
    sample_numeric_order = match(
      as.character(sample),
      sample_order
    ),
    
    condition_offset = case_when(
      condition == "Control" ~ 0,
      condition == "Sample" ~ 100,
      TRUE ~ 200
    )
  ) %>%
  arrange(
    condition,
    sample_numeric_order
  ) %>%
  group_by(condition) %>%
  mutate(
    x_within_condition = row_number(),
    
    x_pos =
      condition_offset +
      x_within_condition
  ) %>%
  ungroup() %>%
  mutate(
    x_label = unname(
      sample_labels[
        as.character(sample)
      ]
    )
  )


# -----------------------------
# Add x positions
# -----------------------------

fold_data <- fold_data %>%
  left_join(
    x_positions %>%
      select(
        condition,
        sample,
        x_pos
      ),
    by = c(
      "condition",
      "sample"
    )
  )

summary_data <- summary_data %>%
  left_join(
    x_positions %>%
      select(
        condition,
        sample,
        x_pos
      ),
    by = c(
      "condition",
      "sample"
    )
  )


# -----------------------------
# Keep the original logarithmic y-axis
#
# MG1655 remains numerically equal to 1.
# Its log10 position is therefore 0 = 10^0.
# -----------------------------

plot_lower_limit <- 0.5
bar_baseline <- 0.6

summary_data <- summary_data %>%
  mutate(
    errorbar_ymin = pmax(
      errorbar_ymin,
      plot_lower_limit + 0.000001
    ),
    
    bar_ymin = bar_baseline,
    bar_ymax = mean_fold_defense
  )


# -----------------------------
# Print summary values
# -----------------------------

summary_output <- summary_data %>%
  arrange(
    phage,
    condition,
    sample
  ) %>%
  select(
    phage,
    condition,
    sample,
    n,
    mean_fold_defense,
    sem_fold_defense,
    errorbar_ymin,
    errorbar_ymax
  )

print(summary_output)


# -----------------------------
# Plot
# -----------------------------

set.seed(1)

p <- ggplot() +
  
  # Mean bars
  geom_rect(
    data = summary_data,
    aes(
      xmin = x_pos - 0.34,
      xmax = x_pos + 0.34,
      ymin = bar_ymin,
      ymax = bar_ymax,
      fill = condition
    ),
    color = "black",
    linewidth = 0.55
  ) +
  
  # Individual replicate-matched values
  geom_point(
    data = fold_data,
    aes(
      x = x_pos,
      y = fold_defense
    ),
    size = 1.8,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.45,
    position = position_jitter(
      width = 0.055,
      height = 0,
      seed = 1
    )
  ) +
  
  # Mean ± SEM
  geom_errorbar(
    data = summary_data,
    aes(
      x = x_pos,
      ymin = errorbar_ymin,
      ymax = errorbar_ymax
    ),
    width = 0.22,
    linewidth = 0.65,
    color = "black"
  ) +
  
  # Facets
  facet_grid(
    rows = vars(phage),
    cols = vars(condition),
    scales = "free_x",
    space = "free_x",
    drop = TRUE
  ) +
  
  # X-axis
  scale_x_continuous(
    breaks = x_positions$x_pos,
    labels = x_positions$x_label,
    expand = expansion(
      add = c(
        0.55,
        0.55
      )
    )
  ) +
  
  # Bar colors
  scale_fill_manual(
    values = c(
      "Control" = "#D55E00",
      "Sample" = "#0072B2"
    ),
    drop = FALSE
  ) +
  
  # Original unchanged y-axis
  scale_y_log10(
    limits = c(
      0.5,
      NA
    ),
    
    breaks = c(
      1,
      10,
      100,
      1000,
      10000,
      100000,
      1000000
    ),
    
    labels = c(
      expression(10^0),
      expression(10^1),
      expression(10^2),
      expression(10^3),
      expression(10^4),
      expression(10^5),
      expression(10^6)
    ),
    
    minor_breaks = c(
      0.6, 0.7, 0.8, 0.9,
      
      2, 3, 4, 5, 6, 7, 8, 9,
      
      20, 30, 40, 50, 60, 70, 80, 90,
      
      200, 300, 400, 500, 600, 700, 800, 900,
      
      2000, 3000, 4000, 5000,
      6000, 7000, 8000, 9000,
      
      20000, 30000, 40000, 50000,
      60000, 70000, 80000, 90000,
      
      200000, 300000, 400000, 500000,
      600000, 700000, 800000, 900000
    ),
    
    expand = expansion(
      mult = c(
        0.01,
        0.12
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = "Fold defense relative to MG1655"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    legend.position = "none",
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.65
    ),
    
    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 0.65
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 14,
      color = "black",
      margin = margin(r = 8)
    ),
    
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 1,
      color = "black",
      size = 8.5,
      lineheight = 0.85,
      margin = margin(t = 5)
    ),
    
    axis.text.y = element_text(
      color = "black",
      size = 11
    ),
    
    axis.line = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.45
    ),
    
    axis.ticks.length = unit(
      0.14,
      "cm"
    ),
    
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.28
    ),
    
    panel.grid.minor.y = element_line(
      color = "grey94",
      linewidth = 0.18
    ),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    panel.spacing.x = unit(
      0.45,
      "lines"
    ),
    
    panel.spacing.y = unit(
      0.55,
      "lines"
    ),
    
    plot.margin = margin(
      8,
      12,
      18,
      8
    )
  )


# -----------------------------
# Display plot
# -----------------------------

print(p)