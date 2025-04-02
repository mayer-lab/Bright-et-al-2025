library(vioplot)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(Seurat)
library(ggrepel)
library(VennDiagram)


# Analysis of FlashTag transcriptome datasets

#Fig2
# Subset FlashTag data from the common pool ----
load("Processed_Objects/Inhibitory_datasets.Rdata")
Inhibitory_datasets <- SetIdent(Inhibitory_datasets, value = "Experiment" )
CFSE <- subset(Inhibitory_datasets, ident = c("E16+6hCFSE", "E12+6hCFSE","E12+4daysCFSE"))

#Make a dataframe ----
#CFSE_DF <-data.frame(c (t(as.matrix(CFSE@assays$RNA@data))),CFSE@meta.data )
CFSE_DF <- CFSE@meta.data
CFSE_DF <- CFSE_DF %>%
  mutate(Pseudotime = as.numeric(Pseudotime),
         Dataset = as.factor(Dataset))

CFSE_DF$Stage <- as.character(CFSE_DF$Dataset)
CFSE_DF$Stage <- factor(CFSE_DF$Stage, levels = unique(CFSE_DF$Stage))
levels(CFSE_DF$Stage)
levels(CFSE_DF$Stage)[levels(CFSE_DF$Stage) == "E028"] = "E12+4days"
levels(CFSE_DF$Stage)[levels(CFSE_DF$Stage) == "E007"] = "E12"
levels(CFSE_DF$Stage)[levels(CFSE_DF$Stage) == "E008"] = "E16"
levels(CFSE_DF$Stage)[levels(CFSE_DF$Stage) == "KB5"] = "E12"
levels(CFSE_DF$Stage)[levels(CFSE_DF$Stage) == "KB6"] = "E16"

write.table(CFSE_DF, "Results/source_data/fig2c.csv", sep = ",", quote = F)
#Make a violin plot of the pseudotime change () ----
vioplot(CFSE_DF$Pseudotime ~ CFSE_DF$Stage,
        horizontal = F,col = c("#E16B8F", "#5776A7", "orange")) 


#Calculate significance  ----

downsampled_df <- CFSE_DF %>%
  group_by(Stage) %>%
  sample_n(1000, replace = TRUE) %>%
  ungroup()

wilcox_test_result <- downsampled_df %>% 
  wilcox_test(Pseudotime ~ Stage, conf.level = 0.95) %>%
  add_significance()
wilcox_test_result
downsampled_df %>% wilcox_effsize(Pseudotime ~ Stage)


#Subset the postmitotic cells ----
CFSE <- SetIdent(CFSE, value = "Fine_annotation")
CFSE_postmitotic <- subset(CFSE, ident = c("Fabp7", "Fabp7_Ccnd2", "Top2a", "Ube2c"), invert = T)
#Find marker genes in postmitotic cells e12.5+6h and e16.5+6h ----
CFSE_postmitotic <- SetIdent(CFSE_postmitotic, value = "Experiment")
markers.df <- FindMarkers(CFSE_postmitotic, ident.1 = "E16+6hCFSE", ident.2 = "E12+6hCFSE", ) # make E16 as first group
#Differential gene expression volcano plot (Fig2d)
markers.df <- as.data.frame(markers.df)
markers.df$diffexpressed <- "NO"
markers.df$diffexpressed[markers.df$avg_log2FC > 1 & markers.df$p_val_adj < 0.05] <- "UP"
markers.df$diffexpressed[markers.df$avg_log2FC< -1 & markers.df$p_val_adj < 0.05] <- "DOWN"
markers.df$delabel <- NA
markers.df$delabel[markers.df$diffexpressed != "NO"] <- rownames(markers.df)[markers.df$diffexpressed != "NO"]

write.table(markers.df, "Results/source_data/fig2d.csv", quote = F, sep = ",")

ggplot(data=markers.df, aes(x=avg_log2FC, y=-log10(p_val_adj), col = diffexpressed), label = common) +
  geom_point(size = 1) +
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  geom_text_repel(label = markers.df$delabel, size = 2.5, color = "black") +
  scale_color_manual(values = c("#E16B8F","#bdbdbd", "#5776A7"))

#Differential genes between FlashTag 6h datasets visualised in all the FlashTag datasets (Fig.2e + SFig.8b)
upregulated <- markers.df$delabel[markers.df$avg_log2FC > 1 & markers.df$p_val_adj < 0.05] 
downregulated <- markers.df$delabel[markers.df$avg_log2FC < -1 & markers.df$p_val_adj < 0.05] 
e12_e16_de_markers <- union(upregulated,downregulated)

CFSE_postmitotic.avg <- AverageExpression(CFSE_postmitotic, return.seurat = T)

write.table(CFSE_postmitotic.avg@assays$RNA@scale.data[e12_e16_de_markers, ], "Results/source_data/fig2e.csv", sep = ",", quote = F)

DoHeatmap(CFSE_postmitotic.avg, features = c(e12_e16_de_markers), raster = F) +
  scale_fill_gradientn(colors = c('#9970ab','#f7f7f7', '#001F1F')) 

## EDF3 a ##
CFSE$cell_state <- "post-mitotic"; CFSE$cell_state[CFSE$Mitotic_and_tips == "Mitotic"] <- "mitotic"
table(CFSE$cell_state)
plot_df <- CFSE@meta.data[, c("Experiment","cell_state")]
plot_df %>%
  group_by(Experiment, cell_state) %>%
  dplyr::count() -> plot_df
plot_df %>%
  group_by(Experiment) %>%
  dplyr::mutate(cells_per_exp = sum(n)) -> plot_df
plot_df$frac_cells <- plot_df$n / plot_df$cells_per_exp * 100

write.table(plot_df, "Results/source_data/edf3a.csv", sep = ",", quote = F, row.names = F)

ggplot(plot_df, aes(x = Experiment, y = frac_cells, fill = cell_state)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw()

## EDF 3b ##
Inhibitory_datasets$broad_state <- Inhibitory_datasets$Mitotic_and_tips
Inhibitory_datasets$broad_state[Inhibitory_datasets$broad_state %in% c(
  "Ebf1_Isl1","Gucy1a3","Nr2f2","Maf","Snhg11"
)] <- "tip"
Inhibitory_datasets$broad_state[Inhibitory_datasets$broad_state == "No"] <- "trunk"

# e12.5+6:
tip_cells <- colnames(Inhibitory_datasets)[Inhibitory_datasets$broad_state == "tip" & Inhibitory_datasets$Experiment == "E12+6hCFSE"]
trunk_cells <- colnames(Inhibitory_datasets)[Inhibitory_datasets$broad_state == "trunk" & Inhibitory_datasets$Experiment == "E12+6hCFSE"]
DimPlot(Inhibitory_datasets, reduction = "umap2", group.by = "broad_state", cells.highlight = list("tip" = tip_cells, "trunk" = trunk_cells), 
        cols.highlight = c("blue","darkblue"), cols = "grey")
# e12.5+96:
tip_cells <- colnames(Inhibitory_datasets)[Inhibitory_datasets$broad_state == "tip" & Inhibitory_datasets$Experiment == "E12+4daysCFSE"]
trunk_cells <- colnames(Inhibitory_datasets)[Inhibitory_datasets$broad_state == "trunk" & Inhibitory_datasets$Experiment == "E12+4daysCFSE"]
DimPlot(Inhibitory_datasets, reduction = "umap2", group.by = "broad_state", cells.highlight = list("tip" = tip_cells, "trunk" = trunk_cells), 
        cols.highlight = c("blue","darkblue"), cols = "grey")
# e16.5+6:
tip_cells <- colnames(Inhibitory_datasets)[Inhibitory_datasets$broad_state == "tip" & Inhibitory_datasets$Experiment == "E16+6hCFSE"]
trunk_cells <- colnames(Inhibitory_datasets)[Inhibitory_datasets$broad_state == "trunk" & Inhibitory_datasets$Experiment == "E16+6hCFSE"]
DimPlot(Inhibitory_datasets, reduction = "umap2", group.by = "broad_state", cells.highlight = list("tip" = tip_cells, "trunk" = trunk_cells), 
        cols.highlight = c("blue","darkblue"), cols = "grey")


## EDF 3c:
CFSE <- SetIdent(Inhibitory_datasets, value = "Experiment" )
E12_6h <- subset(CFSE, ident = c("E12+6hCFSE"))
E16_6h <- subset(CFSE, ident = c("E16+6hCFSE"))
E12_96h <- subset(CFSE, ident = c("E12+4daysCFSE"))

E12_6h <- SetIdent(E12_6h, value = "Fine_annotation" )
E16_6h <- SetIdent(E16_6h, value = "Fine_annotation" )
E12_96h <- SetIdent(E12_96h, value = "Fine_annotation" )

E12_6h <- AverageExpression(E12_6h, return.seurat = T)
E16_6h <- AverageExpression(E16_6h, return.seurat = T)
E12_96h <- AverageExpression(E12_96h, return.seurat = T)

load("Processed_Objects/top5_Inhibitory_datasets.Rdata")

## edf3c was assembled in illustrator from these 3 heatmaps:
DoHeatmap(E12_6h, features = c(top5_Inhibitory_datasets$gene), raster = F,draw.lines = F) + NoLegend()
DoHeatmap(E16_6h, features = c(top5_Inhibitory_datasets$gene), raster = F,draw.lines = F) + NoLegend()
DoHeatmap(E12_96h, features = c(top5_Inhibitory_datasets$gene), raster = F,draw.lines = F) + NoLegend()


## CFSE DE gene overlap (EDF3 d,e)

CFSE <- SetIdent(CFSE, value = "Experiment")

E12_6_markers<- FindMarkers(CFSE, ident.1 = "E12+6hCFSE", only.pos = T)
E12_6 <- rownames(E12_6_markers)[E12_6_markers$avg_log2FC > 0.25 & E12_6_markers$p_val_adj < 0.05 ] 

E12_96_markers<- FindMarkers(CFSE, ident.1 = "E12+4daysCFSE", only.pos = T)
E12_96 <- rownames(E12_96_markers)[E12_96_markers$avg_log2FC > 0.25 & E12_96_markers$p_val_adj < 0.05 ] 

E16_6_markers<- FindMarkers(CFSE, ident.1 = "E16+6hCFSE", only.pos = T)
E16_6 <- rownames(E16_6_markers)[E16_6_markers$avg_log2FC > 0.25 & E16_6_markers$p_val_adj < 0.05 ] 

venn.plot <- venn.diagram(
  x = list(E16_6,E12_96, E12_6),
  category.names = c("E16_6", "E12_96", "E12_6"), # Labels for the lists
  filename = NULL, # Don't save to file
  fill = c("#5776A7","#5AAC56","#E16B8F"), # Colors for the circles
  alpha = 0.5, # Transparency level
  cex = 2, # Font size for the intersection numbers
  cat.cex = 1, # Font size for the category labels
  cat.pos = c(-45, 45, 135), # Positions for category labels
  cat.dist = c(0.05, 0.05, 0.05) # Distance of labels from circles
)

######
dev.off()
grid.draw(venn.plot)

## create data.frame for source data
list_df <- data.frame(
  "experiment" = c("E12+6hCFSE", "E12+4daysCFSE", "E16+6hCFSE"),
  "marker_genes" = c(paste(E12_6, collapse = "/"), paste(E12_96, collapse = "/"), paste(E16_6, collapse = "/"))
)
write.table(list_df, "Results/source_data/edf3e.csv", sep = ",", quote = F, row.names = F)

