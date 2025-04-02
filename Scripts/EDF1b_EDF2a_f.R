library(Seurat)
library(dplyr)

load("Processed_Objects/Inhibitory_datasets.Rdata")

################################################################################
## EDF1b

VlnPlot(Inhibitory_datasets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Dataset", pt.size = 0)

df <- data.frame(
  "dataset" = Inhibitory_datasets$Dataset,
  "nFeature_RNA" = Inhibitory_datasets$nFeature_RNA, 
  "nCount_RNA" = Inhibitory_datasets$nCount_RNA,
  "percent.mt" = Inhibitory_datasets$percent.mt
)
write.table(df, "Results/source_data/edf1b.csv", sep = ",", quote = F, row.names = T)

################################################################################
## EDF2 ##

##
Inhibitory_datasets$Fine_annotation[Inhibitory_datasets$Fine_annotation == "Snhg11_Lhx8"] <- "Snhg11"

## umaps (x axis is flipped in final panels)
DimPlot(Inhibitory_datasets, group.by = "Fine_annotation", reduction = "umap2", label = T)

FeaturePlot(Inhibitory_datasets, "Pseudotime", reduction = "umap2")

FeaturePlot(Inhibitory_datasets, reduction = "umap2", features = c(
  "Nes", "Ascl1", "Meis2", "Tcf4", "Fabp7", "Ccnd2", "Gad2", "Dlx6os1"
), ncol = 4)

bc_vec <- c("Fabp7"="AP","Fabp7_Ccnd2"="BP","Top2a"="BP","Ube2c"="BP","Abracl"="PN_precursor","Tshz1"="PN_precursor",
            "Six3_Gucy1a3"="PN_precursor","Gucy1a3"="PN_precursor","Ebf1_Isl1"="PN_precursor","Nkx2_1"="IN_precursor",
            "Npy"="IN_precursor","Tcf4_Nr2f2"="IN_precursor","Maf_Sst"="IN_precursor","Snhg11_Lhx8"="IN_precursor",
            "Snhg11"="IN_precursor")
Inhibitory_datasets$broad_cell_state <- bc_vec[Inhibitory_datasets$Fine_annotation]
DimPlot(Inhibitory_datasets, reduction = "umap2", group.by = "broad_cell_state")

DimPlot(Inhibitory_datasets, reduction = "umap2", group.by = "Phase")

## heatmap of marker genes:
Inhibitory_datasets <- SetIdent(Inhibitory_datasets, value = "Fine_annotation")
cluster_markers <- FindAllMarkers(Inhibitory_datasets, only.pos = T, logfc.threshold = 1)
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
#Inhibitory_datasets <- ScaleData(Inhibitory_datasets, features = top5$gene, assay = "RNA")
inh_avg <- AverageExpression(Inhibitory_datasets, return.seurat = T, assay = "RNA", slot = "scale.data", group.by = "Fine_annotation", verbose = T)
inh_avg$cluster <- rownames(inh_avg@meta.data)
inh_avg$cluster <- factor(inh_avg$cluster, levels = c("Fabp7","Fabp7_Ccnd2","Top2a","Ube2c","Abracl","Tshz1","Tcf4_Nr2f2","Six3_Gucy1a3","Ebf1_Isl1","Gucy1a3","Nkx2_1","Npy","Maf_Sst","Snhg11"))

plot_df <- inh_avg@assays$RNA@scale.data[top5$gene, ]
write.table(plot_df, "Results/source_data/edf2a.csv", sep = ",", quote = F, row.names = F)
DoHeatmap(inh_avg, features = top5$gene, draw.lines = F, group.by = "cluster")
