# ============================================================
# Pairwise protein percent identity matrix + heatmap
# Sequences are hard-coded
# Order in matrix is preserved from input order
# Pairwise global alignments are calculated for unaligned sequences
# ============================================================

# -----------------------------
# Install Biostrings if needed
# -----------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# if (!requireNamespace("Biostrings", quietly = TRUE)) {
#   BiocManager::install("Biostrings")
# }

library(Biostrings)
library(ggplot2)

# -----------------------------
# Hard-coded protein sequences
# -----------------------------

seqs <- c(
  "KpnNEO"    = "MNIIQVLLTVVTVIGFTKTCGRLDNLLDMLYNVISSVESI",
  "NEO_c1.1" = "MEISSLLTQSCICSKWGLISTFLLVTIYLTQISQVCTMKILDSFKPICYNVIIEILLVLRKEII",
  "NEO_c1.2" = "MELIITLILILSNWLIQELLIVLEVNLEQFLTEHVDKVILDSLKLLCYNEVNENLKGTIW",
  "NEO_c3.1" = "MITDLLVDLMATLAANLIINIAELDNAKHIIVYWIVEGLDI",
  "NEO_c3.2" = "MDNVIDIINCIILILVKNVLDTEYDVMYNNIISGILLRIGGNT",
  "NEO_c4.0" = "MFTEQLSTETIFNLRISISRLTTMINYVIIVKEESPIITLVFWLVTI",
  "NEO_c4.1" = "MITIAIQASLNPQAGECADNVWNVAIASFLILTAVLAYARIELEDT",
  "NEO_c4.2" = "MWSACWLMCEFVSTSELCQCVLCVLLLARPLTSMNNYAMIVEEACS",
  "NEO_c5.1" = "MTTLTRIMTAIFMRLVSSLGAVIEINLNTCIFGIIIVYLTYRNKSLVN",
  "NEO_c5.2" = "MKHPIVKVILAVMILLGQLMASASSEINSGIACISEIVIVYLHINSNKSI",
  "NEO_c6.1" = "MTELKILPITLLLMQLIVTGAQICEIGIKFVTLLTPINIIVST",
  "NEO_c6.2" = "MLNKSSQRNNKFLVTIVITLLIGIFFIILAKDIENGIEIAILVTKGLNT",
  "NEO_c7.1" = "MLLLIFKCIFICIDVNEICTTASKLFGGSWINHYNNINLLKSFSI",
  "NEO_c7.2" = "MNVKFNTNIINSLNLIIFVQGLDILPSLQIG",
  "NEO_c7.3" = "MKRNINQSLLKSTMKDIAISILTFILSSYLYISEELVSLIIEVLVIVINVINTY*",
  "NEO_c7.4" = "MNLISLAFSKGKHRISRCSWYDTIFQFLRLSVIKLISLLIIQKLFVYLSIESELVIYIIHKVI",
  "NEO_c7.5" = "MNQYEYIKQKQQMKMKVSIINFILRIANSLLILSKNLVFLEINVIFTTDDFLYVILLTLNSIIT",
  "NEO_c7.6" = "MQLIGNFDLILMKANRTKSHWLLHIIINWVEALTSFLYITENISDILMLCQESYVSLPLVVLL",
  "NEO_c7.7" = "MRNYSLRIIIFVLKMILQHLDGGMLTSFITEIFITGVLEWFKTLFNKHNKANQT*",
  "NEO_c8.1" = "MNQISSDTQSVILAWLLSYLLRTIIFSNLNNYYVIVASKYILCNVITEI",
  "NEO_c8.2" = "MAWAAALKLSLLAIREKINLNRNNVMKSIMNLTLMAILTILQLYSYSQEY*",
  "NEO_c9.1" = "MLSVSFDWVAVVAMLISIVLPAATFEVDSLNARRYNVTVEVRA",
  "NEO_c9.2" = "MLEVSSIRIAAIVLTLVGIVDAGQYEVDSFNARQYSVVVEVRA",
  "NEO_c10.1" = "MRNPVIALLSIVLAITLPVGLDNLLGMTYNVHMLLEVS",
  "NEO_c10.2" = "MKITIRGLLAALLVILGCYTWDAALGLDSLLGSTYNVSVTLGVS",
  "NEO_c11.1" = "MSFIGIVLKSKSSIQIAGIIHMILLETVFPIDSVNCIVYTVIKL",
  "NEO_c11.2" = "MNNCVMILLEQLTTLMVQFFITIFLNKLLDSSNLICYTIGASM",
  "NEO_c12.1" = "MAFATVLAAATRIALAVLNLVCYSSTMCRGDE",
  "NEO_c12.2" = "MIEIKPTRNPLTQIALAVLNLAWYVVTTMCLEGCD"
)

# -----------------------------
# Clean sequences
# -----------------------------
# Removes stop codons (*) because they are not amino acids for protein identity calculation.

seqs <- gsub("\\*", "", seqs)

seq_names <- names(seqs)

# -----------------------------
# Pairwise percent identity function
# -----------------------------

pairwise_percent_identity <- function(seq1, seq2) {
  
  aln <- pairwiseAlignment(
    pattern = AAString(seq1),
    subject = AAString(seq2),
    type = "global",
    substitutionMatrix = "BLOSUM62",
    gapOpening = 10,
    gapExtension = 0.5
  )
  
  aligned_1 <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  aligned_2 <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  
  # Ignore positions where either sequence has a gap
  keep <- aligned_1 != "-" & aligned_2 != "-"
  
  aligned_1_nogap <- aligned_1[keep]
  aligned_2_nogap <- aligned_2[keep]
  
  identical_positions <- sum(aligned_1_nogap == aligned_2_nogap)
  compared_positions <- length(aligned_1_nogap)
  
  percent_identity <- (identical_positions / compared_positions) * 100
  
  return(percent_identity)
}

# -----------------------------
# Create identity matrix
# -----------------------------

identity_matrix <- matrix(
  NA,
  nrow = length(seqs),
  ncol = length(seqs),
  dimnames = list(seq_names, seq_names)
)

for (i in seq_along(seqs)) {
  for (j in seq_along(seqs)) {
    
    if (i == j) {
      identity_matrix[i, j] <- 100
    } else {
      identity_matrix[i, j] <- pairwise_percent_identity(seqs[i], seqs[j])
    }
  }
}

identity_matrix <- round(identity_matrix, 1)

print(identity_matrix)

# -----------------------------
# Save matrix as CSV
# -----------------------------

write.csv(
  identity_matrix,
  file = "NEO_pairwise_percent_identity_matrix.csv",
  quote = FALSE
)

# -----------------------------
# Convert matrix to long format for ggplot
# -----------------------------

df <- as.data.frame(as.table(identity_matrix))
colnames(df) <- c("Seq1", "Seq2", "Percent_identity")

# Preserve input order in both axes
df$Seq1 <- factor(df$Seq1, levels = rev(seq_names))
df$Seq2 <- factor(df$Seq2, levels = seq_names)

# -----------------------------
# Heatmap plot
# -----------------------------

p <- ggplot(df, aes(x = Seq2, y = Seq1, fill = Percent_identity)) +
  
  geom_tile(color = "white", linewidth = 0.25) +
  
  geom_text(
    
    aes(label = paste0(round(Percent_identity, 1), "%")),
    
    size = 1.8
    
  ) +
  
  scale_fill_gradientn(
    
    colors = c(
      
      "white",
      
      "#d9f0ff",
      
      "#91c9f7",
      
      "#4292c6",
      
      "#08519c"
      
    ),
    
    limits = c(0, 100),
    
    breaks = seq(0, 100, 20),
    
    name = "% identity"
    
  ) +
  
  coord_fixed() +
  
  theme_minimal(base_size = 11) +
  
  theme(
    
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      
      angle = 45,
      
      hjust = 1,
      
      vjust = 1,
      
      size = 8
      
    ),
    
    axis.text.y = element_text(size = 8),
    
    panel.grid = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(size = 10),
    
    legend.text = element_text(size = 9),
    
    plot.title = element_text(hjust = 0.5, face = "bold")
    
  ) +
  
  ggtitle("Pairwise NEO protein percent identity")

print(p)
