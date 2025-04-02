##############################################################################
# Figure description: Dotplot displaying combinatorial binding 

# Stage specific ATAC-seq peaks were used for performing tfcomb analysis
# (https://tf-comb.readthedocs.io/en/latest/examples/TFBS_from_motifs.html)

library(ggplot2)
library(data.table)

# Define TFs of interest
foi <- c("NFIB_1","TCF4_1", "MEIS2_1")

# Load predicted binding information from tf comb analysis for e12 stage
cobind_e12 <- read.csv("Processed_Objects/Ce12_TFbinding.csv")
cobind_e12_sub <- cobind_e12[cobind_e12$TF1 %in% foi & cobind_e12$TF2 %in% foi, ]
cobind_e12_sub2 <- cobind_e12_sub[cobind_e12_sub$X %in% c("NFIB_1-TCF4_1", "NFIB_1-MEIS2_1"), ]
cobind_e12_sub3 <- cobind_e12_sub2[c("X", "TF1_TF2_count", "cosine", "zscore")]
cobind_e12_sub3$stage <- "e12"
cobind_e12_sub4 <- melt(cobind_e12_sub3)

# Load predicted binding information from tf comb analysis for e16 stage
cobind_e16 <- read.csv("Processed_Objects/Ce16_TFbinding.csv")
cobind_e16_sub <- cobind_e16[cobind_e16$TF1 %in% foi & cobind_e16$TF2 %in% foi, ]
cobind_e16_sub2 <- cobind_e16_sub[cobind_e16_sub$X %in% c("NFIB_1-TCF4_1", "NFIB_1-MEIS2_1"), ]
cobind_e16_sub3 <- cobind_e16_sub2[c("X", "TF1_TF2_count", "cosine", "zscore")]
cobind_e16_sub3$stage <- "e16"
cobind_e16_sub4 <- melt(cobind_e16_sub3)

# Merge data frames for e12 and e16 stages
merged_df <- rbind(cobind_e12_sub3, cobind_e16_sub3)

write.table(merged_df, "Results/source_data/edf5d.csv", sep = ",", quote = F, row.names = F)

# Create a dot plot using ggplot2 with the merged_df
ggplot(merged_df, aes(x = X, y = cosine, color = stage, size = TF1_TF2_count)) +
  geom_point() +
  labs(title = "Combinatorial Binding",
       x = "TF1_TF2",
       y = "cosine",
       color = "Stage",
       size = "TF1_TF2 Count") +
  scale_color_manual(values = c("e12" = "#E16B8F", "e16" = "#5776A7")) +
  theme_minimal()
