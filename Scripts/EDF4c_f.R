## EDF 4 c-g

library(ArchR)
library(parallel)
library(ggplot2)

## load archr project:
archr_proj <- loadArchRProject("/data/mayerlab/neuhaus/dorsal_ventral_comp/cfse_network/results/archr_proj_fAnn_wLabelTransfer/", force = T)

## EDF 4c:
marker_genes <- c("Fabp7","Ccnd2","Dlx5","Gad2","Maf","Ebf1")
p <- plotEmbedding(
  ArchRProj = archr_proj, 
  colorBy = "GeneScoreMatrix", 
  name = marker_genes, 
  embedding = "UMAPHarmony_S",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

## EDF 4d:

load("Processed_Objects/Inhibitory_datasets.Rdata")
cfse_seurat_sub <- subset(Inhibitory_datasets, subset = Experiment %in% c("E12+6hCFSE","E16+6hCFSE"))

## do label transfer:
archr_proj2 <- addGeneIntegrationMatrix(
  ArchRProj = archr_proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = cfse_seurat_sub,
  addToArrow = FALSE,
  groupRNA = "Broad_annotation",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
archr_proj2$stage_broad_annotation <- paste0(archr_proj2$predictedGroup_Un, "_" ,archr_proj2$dataset)
dps <- plotEmbedding(archr_proj2, embedding = "UMAPHarmony", colorBy = "cellColData", name = c("dataset","predictedGroup_Un","stage_broad_annotation"))
wrap_plots(dps$dataset, dps$predictedGroup_Un, dps$stage_broad_annotation)

## EDF 4e:

grange_e12 <- readRDS("Processed_Objects/e12-reproduciblePeaks.gr.rds")
length(grange_e12)
grange_e16 <- readRDS("Processed_Objects/e16-reproduciblePeaks.gr.rds")
length(grange_e16)

overlap_e12_e16 <- findOverlaps(grange_e12, grange_e16, minoverlap = 100)
overlap_e12_e16
overlapping_peaks <- grange_e12[queryHits(overlap_e12_e16)]
length(overlapping_peaks)

non_overlapping_e12 <- subsetByOverlaps(grange_e12, grange_e16, invert = TRUE)
length(non_overlapping_e12)
non_overlapping_e16 <- subsetByOverlaps(grange_e16, grange_e12, invert = TRUE)
length(non_overlapping_e16)

non_overlapping_e12_df <- as.data.frame(non_overlapping_e12); non_overlapping_e12_df$class <- "e12_specific"
non_overlapping_e16_df <- as.data.frame(non_overlapping_e16); non_overlapping_e16_df$class <- "e16_specific"
overlapping_df <- as.data.frame(overlapping_peaks); overlapping_df$class <- "overlapping"

merge_df <- rbind(non_overlapping_e12_df, non_overlapping_e16_df, overlapping_df)
write.table(merge_df, "Results/source_data/edf4e.csv", sep = ",", quote = F, row.names = F)
ggplot(merge_df, aes(x = class, fill = peakType)) +
  geom_bar(stat = "count", position = "dodge") +
  theme_bw()


## EDF 4f:
markerTest_stage <- getMarkerFeatures(
  ArchRProj = archr_proj, 
  useMatrix = "PeakMatrix",
  groupBy = "dataset",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "e12",
  bgdGroups = "e16"
)
pv <- plotMarkers(seMarker = markerTest_stage, name = "e12", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")

log2fc <- assay(markerTest_stage, "Log2FC")
Pval <- assay(markerTest_stage, "Pval") # If FDR is not available, MeanBGD could represent similar info
row_data <- as.data.frame(rowData(markerTest_stage))

# Combine the data into a single data frame for plotting
plot_df <- cbind(row_data, Log2FC = log2fc, Pval = Pval) # Adjust 'FDR' if it's in a different assay
colnames(plot_df) <- c("seqnames", "idx" , "start", "end", "Log2FC", "Pval")
head(plot_df)


# Create a new column to define up/down-regulated and non-significant points
plot_df$Peaks <- ifelse(plot_df$Pval <= 0.1 & plot_df$Log2FC >= 1, "e12.5 enriched",
                        ifelse(plot_df$Pval <= 0.1 & plot_df$Log2FC <= -1, "e16.5 enriched", 
                               "Non-significant"))

# Count the number of enriched peaks
num_e12_enriched <- sum(plot_df$Peaks == "e12.5 enriched")
num_e16_enriched <- sum(plot_df$Peaks == "e16.5 enriched")
num_non_significant <- sum(plot_df$Peaks == "Non-significant")

# Create custom labels for the legend with the counts included
custom_labels <- c(paste0("e12.5 enriched (", num_e12_enriched, ")"),
                   paste0("e16.5 enriched (", num_e16_enriched, ")"),
                   paste0("Non-significant (", num_non_significant, ")"))

write.table(plot_df, "Results/source_data/edf4f.csv", sep = ",", quote = F, row.names = F)

# Create the volcano plot with specific colors and enriched peak numbers in the legend
ggplot(plot_df, aes(x = Log2FC, y = -log10(Pval))) +
  geom_point(aes(color = Peaks), size = 2) +  # Use Regulation for color mapping
  scale_color_manual(values = c("e12.5 enriched" = "#E16B8F", "e16.5 enriched" = "#5776A7", "Non-significant" = "gray"),
                     labels = custom_labels) +  # Use custom labels with counts
  theme_minimal() +
  labs(title = "Volcano Plot for Differential Peaks",
       x = "Log2 Fold Change",
       y = "-log10 Pval") +
  theme(legend.position = "right")


