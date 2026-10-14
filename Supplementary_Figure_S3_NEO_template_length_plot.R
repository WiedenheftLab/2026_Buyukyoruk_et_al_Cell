# install.packages(c("tidyverse","tidytext"))
library(tidyverse)
library(tidytext)
library(stringr)  # for str_extract()

# ——— 1) read & prep data ——————————————————————————————————
df <- read.delim(
  "~/Desktop/NEO_all_cargo_sequences_05192025_plot_data.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(df) <- c("ID","Length","Clade")

df <- df %>%
  # make Length numeric and compute mod3 factor
  mutate(
    Length = as.integer(Length),
    mod3    = factor(Length %% 3, levels = c(0,1,2))
  )

# pull out the numeric part of each Clade, dedupe & sort, re‐build the levels
nums         <- unique(as.integer(str_extract(df$Clade, "\\d+")))
clade_levels <- paste0("Clade", sort(nums))

df <- df %>%
  mutate(
    Clade = factor(Clade, levels = clade_levels)
  )

# ——— 2) plot ——————————————————————————————————————————
ggplot(df, aes(
  x    = reorder_within(ID, Length, Clade, .desc = FALSE),
  y    = Length,
  fill = mod3
)) +
  geom_col() +
  
  # ▼ put strips on the bottom using switch="x"
  facet_grid(
    . ~ Clade,
    scales = "free_x",
    space  = "free_x",
    switch = "x"
  ) +
  
  scale_x_reordered(expand = expansion(add = c(10, 10))) +
  scale_fill_manual(
    name   = "nt increments",
    values = c("0" = "lightblue3", "1" = "gray50", "2" = "gray80")
  ) +
  scale_y_continuous(
    breaks       = seq(90, 240, by = 6),
    minor_breaks = seq(93, 240, by = 6),
    expand       = expansion(add = 0)
  ) +
  coord_cartesian(ylim = c(90, 240)) +
  labs(x = NULL, y = "Length (nt)") +
  theme_bw() +
  theme(
    # clean up panel & plot bg
    panel.background   = element_rect(fill = NA),
    plot.background    = element_rect(fill = NA),
    
    # dashed horizontal grid only
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.grid.minor.y = element_line(linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    # drop x‑axis text/ticks
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    
    # move strip labels outside & style them
    strip.placement    = "outside",
    strip.background   = element_blank(),
    strip.text.x       = element_text(
      angle = 60,
      face  = "bold",
      hjust = 0.5
    ),
    
    panel.spacing.x    = unit(0.5, "lines")
  )