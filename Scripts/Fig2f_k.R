#Scripts to generate Figure 2f to 2k
library(ArchR)
library(parallel)
library(GenomicRanges)
library(ggplot2)

##############################################################################
# Figure number: 2f
# Figure description: scATAC-seq UMAP based on stage.

ArchRProject <- loadArchRProject("/path/to/archr_proj/", force = T)


# Plotting figure
plotEmbedding(
  ArchRProj = ArchRProject, 
  colorBy = "cellColData", 
  name = "dataset", 
  embedding = "UMAPHarmony_S", 
  size = 0.01, 
  pal = c("e12" = "#E16B8F", "e16" = "#5776A7"), 
  labelMeans = FALSE
) +
  theme(
    legend.text = element_text(size = 10), 
    legend.position = "right", 
    legend.direction = "vertical"
  ) + 
  guides(color = guide_legend(override.aes = list(size = 5)))

##############################################################################



##############################################################################
# Figure number: 2g
# Figure description: scATAC-seq temporal dynamics

# Metacell scripts: https://jdblischak.github.io/singleCellSeq/analysis/coverage-endogenous.html

# Loading peak files

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

## edf 4e:
non_overlapping_e12_df <- as.data.frame(non_overlapping_e12); non_overlapping_e12_df$class <- "e12_specific"
non_overlapping_e16_df <- as.data.frame(non_overlapping_e16); non_overlapping_e16_df$class <- "e16_specific"
overlapping_df <- as.data.frame(overlapping_peaks); overlapping_df$class <- "overlapping"

merge_df <- rbind(non_overlapping_e12_df, non_overlapping_e16_df, overlapping_df)
write.table(merge_df, "Results/source_data/edf4e.csv", sep = ",", quote = F, row.names = F)
ggplot(merge_df, aes(x = class, fill = peakType)) +
  geom_bar(stat = "count", position = "dodge") +
  theme_bw()


# Load required libraries
library("genomation")
library("Rsamtools")
library("plyr")
library("tidyverse")
library("ggplot2")
theme_set(theme_bw(base_size = 14))
theme_update(panel.grid.minor.x = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.major.y = element_blank())

# Loading BAM files
bam <- c("Processed_Objects/filtered_E12.bam",
         "Processed_Objects/filtered_E16.bam")


# Read BED file
bed_file <- readBed("Processed_Objects/non_overlapping_e12.bed")
bed_file


# ScoreMatrixList for E16 peaks
peak_score <- ScoreMatrixList(target = bam, windows = non_overlapping_e16, type = "bam",
                              rpm = TRUE, strand.aware = TRUE)
names(peak_score) <- c("E12_ATAC", "E16_ATAC")
peak_score_sm_df <- ldply(peak_score, colMeans, .id = "sample_id")
colnames(peak_score_sm_df)[-1] <- paste0("p", 1:(ncol(peak_score_sm_df) - 1))
peak_score_sm_df$feature <- "E16_Peak"
peak_score_sm_df_long <- gather(peak_score_sm_df, key = "pos", value = "rpm", p1:p500)

# ScoreMatrixList for overlapping peaks
peak_score_overlap <- ScoreMatrixList(target = bam, windows = overlapping_peaks, type = "bam",
                                      rpm = TRUE, strand.aware = TRUE)
names(peak_score_overlap) <- c("E12_ATAC", "E16_ATAC")
peak_score_overlap_df <- ldply(peak_score_overlap, colMeans, .id = "sample_id")
colnames(peak_score_overlap_df)[-1] <- paste0("p", 1:(ncol(peak_score_overlap_df) - 1))
peak_score_overlap_df$feature <- "Overlapping_Peak"
peak_score_overlap_df_long <- gather(peak_score_overlap_df, key = "pos", value = "rpm", p1:p500)

# ScoreMatrixList for E12 peaks
peak_score_e12 <- ScoreMatrixList(target = bam, windows = non_overlapping_e12, type = "bam",
                                  rpm = TRUE, strand.aware = TRUE)
names(peak_score_e12)
names(peak_score_e12) <- c("E12_ATAC", "E16_ATAC")
head(peak_score_e12)
peakscoree12_sm_df <- ldply(peak_score_e12, colMeans, .id = "sample_id")
head(peakscoree12_sm_df)
colnames(peakscoree12_sm_df)[-1] <- paste0("p", 1:(ncol(peakscoree12_sm_df) - 1))
peakscoree12_sm_df$feature <- "E12_Peak"
peakscoree12_sm_df_long <- gather(peakscoree12_sm_df, key = "pos", value = "rpm", p1:p500)

# Plotting figures
features <- do.call("rbind", list(peakscoree12_sm_df_long, peak_score_sm_df_long, peak_score_overlap_df_long))
# Convert base position back to integer value
features$pos <- sub("p", "", features$pos)
features$pos <- as.numeric(features$pos)
# Subtract 1001 to recalibrate as +/- 1 kb
features$pos <- features$pos - 250
features$feature <- factor(features$feature, levels = c("E12_Peak", "E16_Peak", "Overlapping_Peak"),
                           labels = c("E12 Peak","E16 Peak", "Overlapping Peak"))

write.table(features, "Results/source_data/fig2g_upper.csv", sep = ",", quote = F, row.names = F)
ggplot(features, aes(x = pos, y = rpm, color = sample_id)) +
  geom_line() +
  facet_wrap(~ feature, scales = "free_y") +
  scale_color_manual(name = "Sample", values = c("#E16B8F", "#5776A7")) +
  labs(x = "Relative position (bp)",
       y = "Counts per million (mean)",
         title = "Peak Signal Coverage") +
  theme_classic()



# Define paths to the H3K4me1 BigWig files for different developmental stages
bw_me1 <- c( "Processed_Objects/H3K4me1_e12.5.bw",
             "Processed_Objects/H3K4me1_e16.5.bw") 

# Resize peaks by extending 2000 base pairs on each side while maintaining the center position
non_overlapping_e12b <- resize(non_overlapping_e12, width = width(non_overlapping_e12) + 2000, fix = "center")
non_overlapping_e16b <- resize(non_overlapping_e16, width = width(non_overlapping_e16) + 2000, fix = "center")
overlapping_peaks_b <- resize(overlapping_peaks, width = width(overlapping_peaks) + 2000, fix = "center")

# Calculate peak scores for E12 peaks using H3K4me1 BigWig files
me1_score = ScoreMatrixList(target = bw_me1, windows = non_overlapping_e12b, type = "bigWig", strand.aware = T)
names(me1_score) <- c("E12.5_H3K4me1", "E16.5_H3K4me1")
me1_score_sm_df <- ldply(me1_score, colMeans, .id = "sample_id")
colnames(me1_score_sm_df)[-1] <- paste0("p", 1:(ncol(me1_score_sm_df) - 1))
me1_score_sm_df$feature = "E12_Peak"
me1_score_sm_df_long <- gather(me1_score_sm_df, key = "pos", value = "rpm", p1:p2501)


# Calculate peak scores for E16 peaks using H3K4me1 BigWig files
me1_score_e16 = ScoreMatrixList(target = bw_me1, windows = non_overlapping_e16b, type = "bigWig", strand.aware = T)
names(me1_score_e16) <- c("E12.5_H3K4me1", "E16.5_H3K4me1")
me1_score_e16_df <- ldply(me1_score_e16, colMeans, .id = "sample_id")
colnames(me1_score_e16_df)[-1] <- paste0("p", 1:(ncol(me1_score_e16_df) - 1))
me1_score_e16_df$feature = "E16_Peak"
me1_score_e16_long <- gather(me1_score_e16_df, key = "pos", value = "rpm", p1:p2501)

# Calculate peak scores for overlapping peaks using H3K4me1 BigWig files
me1_score_overlap = ScoreMatrixList(target = bw_me1, windows = overlapping_peaks_b, type = "bigWig", strand.aware = T)
names(me1_score_overlap) <- c("E12.5_H3K4me1","E16.5_H3K4me1")
me1_score_overlap_df <- ldply(me1_score_overlap, colMeans, .id = "sample_id")
colnames(me1_score_overlap_df)[-1] <- paste0("p", 1:(ncol(me1_score_overlap_df) - 1))
me1_score_overlap_df$feature = "Overlapping_Peak"
me1_score_overlap_df_long <- gather(me1_score_overlap_df, key = "pos", value = "rpm", p1:p2501)

# Merge peak score data frames for E12, E16, and overlapping peaks
features <- do.call("rbind", list(me1_score_sm_df_long, me1_score_e16_long, me1_score_overlap_df_long))
# Convert base position back to integer value
features$pos <- sub("p", "", features$pos)
features$pos <- as.numeric(features$pos)
# Subtract 1001 to recalibrate as +/- 1 kb
features$pos <- features$pos - 1250

# Order factor so that TSS is displayed left of TES
features$feature <- factor(features$feature, levels = c("E12_Peak", "E16_Peak", "Overlapping_Peak"),
                           labels = c("E12_Peak",
                                      "E16_Peak", "Overlapping_Peak"))


write.table(features, "Results/source_data/fig2g_bottom.csv", sep = ",", quote = F, row.names = F)
# Plot peak signal coverage
ggplot(features, aes(x = pos, y = rpm, color = sample_id)) +
  geom_line() +
  facet_wrap(~ feature, scales = "free_y") +
  scale_color_manual(name = "Sample", values = c("#638B66", "#F47942")) +
  labs(x = "Relative position (bp)",
       y = "Counts per million (mean)",
       title = "Peak Signal Coverage") +
  theme_classic()
##############################################################################



##############################################################################
#Figure number: 2h
#Figure description: heatmap of peaks along trajectory

# get trajectories from project:
e12_traj_matrix  <- getTrajectory(ArchRProj = ArchRProject, name = "e12_Trajectory", useMatrix = "PeakMatrix")
e16_traj_matrix  <- getTrajectory(ArchRProj = ArchRProject, name = "e16_Trajectory", useMatrix = "PeakMatrix")

# Extracting trajectory information for e12
e12_trajectory_heatmap <- plotTrajectoryHeatmap(e12_traj_matrix, pal = paletteContinuous(set = "solarExtra"), returnMatrix = TRUE)
rownames(e12_trajectory_heatmap) <- gsub("_", "-", rownames(e12_trajectory_heatmap))
head(e12_trajectory_heatmap)

# Extracting trajectory information for e16
e16_trajectory_heatmap <- plotTrajectoryHeatmap(e16_traj_matrix, pal = paletteContinuous(set = "solarExtra"), returnMatrix = TRUE)
rownames(e16_trajectory_heatmap) <- gsub("_", "-", rownames(e16_trajectory_heatmap))
head(e16_trajectory_heatmap)

# Find the common names
common_peaks <- intersect(rownames(e12_trajectory_heatmap), rownames(e16_trajectory_heatmap))
length(common_peaks)

# Find the unique names in e12
unique_names_e12 <- setdiff(rownames(e12_trajectory_heatmap), rownames(e16_trajectory_heatmap))
length(unique_names_e12)

# Find the unique names in e16
unique_names_e16 <- setdiff(rownames(e16_trajectory_heatmap), rownames(e12_trajectory_heatmap))
length(unique_names_e16)

# Combine all names into a single data frame
all_names <- unique(c(common_peaks, unique_names_e12, unique_names_e16))

# Making annotation of E12 heatmap 
table_common <- cbind(unique(common_peaks), annotation = "Common")
head(table_common)
table_e12 <- cbind(unique(unique_names_e12), annotation = "E12")
head(table_e12)
table_combined <- as.data.frame(rbind(table_common, table_e12))
head(table_combined)
table_final <- table_combined[-1]
rownames(table_final) <- table_combined$V1
annotation_colors_e12 <- list(annotation = c(Common = "grey", E12 = "#E16B8F"))
heatmap_e12 <- pheatmap(as.matrix(e12_trajectory_heatmap), cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, color = h6, 
                        use_raster = FALSE, annotation_row = table_final, annotation_colors = annotation_colors_e12)
heatmap_e12


# Making annotation of E16 heatmap 
table_common <- cbind(unique(common_peaks), annotation = "Common")
table_e16 <- cbind(unique(unique_names_e16), annotation = "E16")
table_combined <- as.data.frame(rbind(table_common, table_e16))
table_final <- table_combined[-1]
rownames(table_final) <- table_combined$V1
annotation_colors_e16 <- list(annotation = c(Common = "grey", E16 = "#5776A7"))

heatmap_e16 <- pheatmap(as.matrix(e16_trajectory_heatmap), cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, color = h6, 
                        use_raster = FALSE, annotation_row = table_final, annotation_colors = annotation_colors_e16)
heatmap_e16



##############################################################################
# Figure Number: 2i
# Volcano Plot

library(ggrepel)

changelimit = 0.19
pvallimit = 130

topgenes = read.table("Processed_Objects/bindetect_results.txt", stringsAsFactor = FALSE, h=TRUE)
topgenes$diffexpressed <- rep("NO",length(topgenes$E12_E16_change))

min(topgenes$E12_E16_change)
max(topgenes$E12_E16_change)
View(topgenes)

topgenes$diffexpressed[topgenes$E12_E16_change > 0.1 & -log10(topgenes$E12_E16_pvalue) > 100] <- "UP"
topgenes$diffexpressed[topgenes$E12_E16_change < -0.1 & -log10(topgenes$E12_E16_pvalue) > 100] <- "DOWN"
topgenes$labels <- rep(NA,length(topgenes$E12_E16_change))
topgenes$labels[topgenes$diffexpressed != "NO"] <- topgenes$name[topgenes$diffexpressed != "NO"]
table(topgenes$diffexpressed)

write.table(topgenes, "Results/source_data/fig2i.csv", sep = ",", quote = F, row.names = F)
topgenes$E12_E16_change <- -topgenes$E12_E16_change
ggplot(data=topgenes, aes(x=E12_E16_change, y=-log10(E12_E16_pvalue), col = diffexpressed), label = labels) +
  geom_point(size = 1) +
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  geom_text_repel(label = topgenes$labels, size = 3, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#5776A7","#bdbdbd", "#E16B8F"),labels=c("Higher scores in E16", "No changes", "Higher scores in E12")) +
  geom_vline(xintercept=c(-0.09, 0.09),linetype = "dashed",  alpha=0.4) +
  geom_hline(yintercept=0.09, linetype = "dashed", alpha=0.4)  +
  labs(x ="Change", y = expression("Significance (-log1])")) + 
  expand_limits(x=c(-0.2,0.2)) + 
  guides(colour=guide_legend(title=NULL))



##############################################################################
# Figure number: 2j
# Figure description: Generating aggregate footprint plots for TFs of interest using TOBIAS pipeline

# TOBIAS pipeline command for generating aggregate footprint plots

# Replace placeholders with actual file paths or command parameters as needed.

# /home/brightann/anaconda3/bin/TOBIAS PlotAggregate \
# --TFBS path/to/TFBS_file.bed \  # Replace path/to/TFBS_file.bed with the actual path to the TFBS BED file.
# --signals path/to/filtered_E12_corrected.bw path/to/filtered_E16_corrected.bw \  # Replace path/to/filtered_E12_corrected.bw and path/to/filtered_E16_corrected.bw with the actual paths to the signal files.
# --output NFIB_all.pdf \
# --share_y both \
# --plot_boundaries \
# --signal-on-x

##############################################################################



##############################################################################
# Figure number: 2k
# Figure description: NFIB Coverage Plot 
library(ChIPpeakAnno)

tobias_nfib <- read.table("Processed_Objects/NFIB_1_MA1643.1_all.bed",
                          header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
colnames(tobias_nfib) <- c("seqnames","start","end","name","score","strand","peak_seqnames","peak_start","peak_end","peak_name","V11","V12","ens_id","symbol","V15","V16")
tobias_nfib_gr <- toGRanges(tobias_nfib)

tobias_nfib_grb <- resize(tobias_nfib_gr, width = width(tobias_nfib_gr) + 478, fix = "center")
width(tobias_nfib_grb)[1:5]
nfib_score = ScoreMatrixList(target = bam, windows = tobias_nfib_grb, type = "bam",
                             rpm = TRUE, strand.aware = TRUE)
names(nfib_score) <- c("E12_ATAC", "E16_ATAC")
nfib_score_df <- ldply(nfib_score, colMeans, .id = "sample_id")
colnames(nfib_score_df)[-1] <- paste0("p", 1:(ncol(nfib_score_df) - 1))
nfib_score_df$feature = "NFIB_all"
nfib_score_df_long <- gather(nfib_score_df, key = "pos", value = "rpm", p1:p500)

features <- nfib_score_df_long
# Convert base position back to integer value
features$pos <- sub("p", "", features$pos)
features$pos <- as.numeric(features$pos)
# Subtract 1001 to recalibrate as +/- 1 kb
features$pos <- features$pos - 250

features$feature <- factor(features$feature, levels = c("E12_Peak", "E16_Peak", "Overlapping_Peak"),
                           labels = c("E12 Peak","E16 Peak", "Overlapping Peak"))
write.table(features, "Results/source_data/fig2k.csv",sep = ",", row.names = F, quote = F)
ggplot(features, aes(x = pos, y = rpm, color = sample_id)) +
  geom_line() +
  facet_wrap(~ feature, scales = "free_y") +
  scale_color_manual(name = "Sample", values = c("#E16B8F", "#5776A7")) +
  labs(x = "Relative position (bp)",
       y = "Counts per million (mean)",
       title = "Peak Signal Coverage") +
  theme_classic()



