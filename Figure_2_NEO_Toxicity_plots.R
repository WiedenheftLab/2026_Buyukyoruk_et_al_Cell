library(ggplot2)
# df <- read.delim("/Users/muratbuyukyoruk/Documents/DRT2//DRT2_biochem/Neo_Toxicity_growth_curve_repeat_10012026.txt",header = TRUE)
df <- read.delim("/Users/muratbuyukyoruk/Documents/DRT2/DRT2_biochem/ORF_w_stop_Toxicity_growth_curve_repeat_10012026.txt",header = TRUE)

df[, 1] <- factor(df[, 1] , levels = unique(df[, 1]))

df$label <- dplyr::case_when(
  is.na(df$Average) ~ "", 
  df$Average >= 0 ~ ""
)

df$Name <- factor(df$Name, levels = rev(levels(df$Name)))

p2 <- ggplot(subset(df,Condition=="0.2% arabinose"),aes(x = Name, y = Average, fill = Repeat)) + 
  geom_bar(stat="identity", size = 0.25, color = "black") +
  geom_errorbar(aes(ymax = Average+SEM, ymin = Average-SEM), position="dodge", width = 0.5, size = 0.25) + 
  facet_grid(Group ~ ., scales="free", space="free_y") +  # <-- Changed here too
  geom_text(aes(y = 0.1, label = label)) +
  scale_fill_manual(values=c('blue1','orange','purple',"darkgreen")) +
  theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "gray80", size = 0.5),   # Only major x gridlines
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),        # Border around plot
    strip.background = element_rect(color = "black", linewidth = 0.7),
    strip.text.y = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10), face = "bold", size = 10),
    axis.text.x = element_text(face = "bold", size =  10),
    legend.text = element_text(face = "bold")
  ) +
  labs(x="NEO", y = "Growth Rate (ΔOD600/hr)", na.value = "No data") +
  coord_flip()
p2
