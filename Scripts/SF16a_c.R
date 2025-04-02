## MILO Excitatory gNfib/x ##

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(miloR)
library(patchwork)
library(SingleCellExperiment)
library(pals)

## load pre-processed seurat:
gNFI_seurat <- readRDS("Processed_Objects/gNFI_merged_seurat_wLabel2.rds")


## -----------------------------------------------------------------------------
## subset for excitatory clusters:
gNFI_exc_sub <- subset(gNFI_seurat, subset = cluster_gene_annot %in% c("Gm29260_Hist1h1b","Unc5d_Nrg1","Satb2_9130024F11Rik","Tafa1_Adgrl3","Kcnip4_Nrg3"))

DimPlot(gNFI_exc_sub, group.by = c("cluster_gene_annot","predicted_EI_annotation"), label = T)

## re-process:
gNFI_exc_sub <- NormalizeData(gNFI_exc_sub, normalization.method = "LogNormalize", scale.factor = 10000)
gNFI_exc_sub <- FindVariableFeatures(gNFI_exc_sub, selection.method = "vst", nfeatures = 2000)
gNFI_exc_sub <- ScaleData(gNFI_exc_sub, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
gNFI_exc_sub <- RunPCA(gNFI_exc_sub, features = VariableFeatures(object = gNFI_exc_sub))
ElbowPlot(gNFI_exc_sub)

gNFI_exc_sub <- FindNeighbors(gNFI_exc_sub, dims = 1:15)
gNFI_exc_sub <- FindClusters(gNFI_exc_sub, resolution = 0.5)
gNFI_exc_sub$seurat_clusters_r_0_5 <- gNFI_exc_sub$seurat_clusters

gNFI_exc_sub <- RunUMAP(gNFI_exc_sub, dims = 1:15)

u1 <- DimPlot(gNFI_exc_sub, reduction = "umap", label = T, group.by = "seurat_clusters_r_0_5") + NoLegend()
u2 <- DimPlot(gNFI_exc_sub, reduction = "umap", label = T, group.by = "replicate") + NoLegend()
u3 <- DimPlot(gNFI_exc_sub, reduction = "umap", label = T, group.by = "predicted_broad_annotation") + NoLegend()
u1+u2+u3

FeaturePlot(gNFI_exc_sub, features = c("Nes","Fabp7","Meis2","Gad2", "Dlx5","Lhx6","Sst","Eomes", "Neurog2","Neurod2", "Neurod6"))

## run DE analysis:
Idents(gNFI_exc_sub) <- gNFI_exc_sub$seurat_clusters_r_0_5
gNFI_exc_cluster_markers <- FindAllMarkers(gNFI_exc_sub, assay = "RNA", slot = "data", logfc.threshold = 1, only.pos = TRUE)

## annotate clusters:
cluster_annot <- c()
for(ct in unique(gNFI_exc_sub$seurat_clusters_r_0_5)) {
  df_sub <- gNFI_exc_cluster_markers[gNFI_exc_cluster_markers$cluster == ct, ]
  df_sub <- df_sub[order(df_sub$avg_log2FC, decreasing = T), ]
  cluster_annot[as.character(ct)] <- paste0(df_sub$gene[1], "_", df_sub$gene[2])
}

gNFI_exc_sub$cluster_gene_annot <- cluster_annot[as.character(gNFI_exc_sub$seurat_clusters_r_0_5)]
DimPlot(gNFI_exc_sub, group.by = "cluster_gene_annot", label = T) + NoLegend() + ggtitle("Cluster annotated by top markers")


## plot Heatmap for excitatory clusters:
gNFI_exc_cluster_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top_cluster_markers
cluster_ID_annot_map <- sapply(unique(gNFI_exc_sub$seurat_clusters_r_0_5), function(cid) {
  unique(gNFI_exc_sub$cluster_gene_annot[gNFI_exc_sub$seurat_clusters_r_0_5 == cid])
})
names(cluster_ID_annot_map) <- unique(gNFI_exc_sub$seurat_clusters_r_0_5)
top_cluster_markers$cluster <- cluster_ID_annot_map[top_cluster_markers$cluster]

Idents(gNFI_exc_sub) <- "cluster_gene_annot"
DoHeatmap(gNFI_exc_sub, features = top_cluster_markers$gene) + NoLegend()

## -----------------------------------------------------------------------------

## subset for cells that contain a guide:
gNFI_exc_sub <- subset(gNFI_exc_sub, subset = guide_annot != "no guide")

## re-run pre-processing:
gNFI_exc_sub <- NormalizeData(gNFI_exc_sub, normalization.method = "LogNormalize", scale.factor = 10000)
gNFI_exc_sub <- FindVariableFeatures(gNFI_exc_sub, selection.method = "vst", nfeatures = 2000)
gNFI_exc_sub <- ScaleData(gNFI_exc_sub, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
gNFI_exc_sub <- RunPCA(gNFI_exc_sub, features = VariableFeatures(object = gNFI_exc_sub))
ElbowPlot(gNFI_exc_sub, ndims = 40)

gNFI_exc_sub <- FindNeighbors(gNFI_exc_sub, dims = 1:15)
gNFI_exc_sub <- FindClusters(gNFI_exc_sub, resolution = 0.5)
gNFI_exc_sub <- RunUMAP(gNFI_exc_sub, dims = 1:15)


## convert object:
gNFI_exc_sce <- as.SingleCellExperiment(gNFI_exc_sub)
gNFI_exc_sce
gNFI_exc_milo <- Milo(gNFI_exc_sce)


## build KNN graph:
## findNeighbours uses k=20
gNFI_exc_milo <- buildGraph(gNFI_exc_milo, k = 40, d = 20)
gNFI_exc_milo

gNFI_exc_milo <- makeNhoods(gNFI_exc_milo, prop = 0.1, k = 40, d=20, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(gNFI_exc_milo)

nhood_mtx <- nhoods(gNFI_exc_milo)
nhood_mtx_gNFI <- nhood_mtx[gNFI_exc_milo@colData$guide_annot == "gNFI", ]
nhood_mtx_gLacZ <- nhood_mtx[gNFI_exc_milo@colData$guide_annot == "gLacZ", ]

hist(colSums(nhood_mtx_gNFI), breaks = 100)
hist(colSums(nhood_mtx_gLacZ), breaks = 100)

## do counting:
colData(gNFI_exc_milo)$sample_annot <- paste(colData(gNFI_exc_milo)$guide_annot, colData(gNFI_exc_milo)$replicate, sep = "__")

gNFI_exc_milo <- countCells(gNFI_exc_milo, meta.data = as.data.frame(colData(gNFI_exc_milo)), sample="sample_annot")
head(nhoodCounts(gNFI_exc_milo))

gNFI_design <- data.frame(colData(gNFI_exc_milo))[,c("sample_annot", "guide_annot", "replicate")]

## Convert batch info from integer to factor
gNFI_design$replicate <- as.factor(gNFI_design$replicate) 
gNFI_design <- distinct(gNFI_design)
rownames(gNFI_design) <- gNFI_design$sample_annot

gNFI_design

## compute neighborhood connectivity:
gNFI_exc_milo <- calcNhoodDistance(gNFI_exc_milo, d=20, reduced.dim = "PCA")

## testing:
da_results <- testNhoods(gNFI_exc_milo, design = ~ replicate + guide_annot, design.df = gNFI_design)
head(da_results)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

## look at DA results:
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)


## build graph:
gNFI_exc_milo <- buildNhoodGraph(gNFI_exc_milo)
## Plot neighbourhood graph
plotNhoodGraphDA(gNFI_exc_milo, da_results, layout="UMAP",alpha=0.1)  +
  theme(legend.position = "none")
