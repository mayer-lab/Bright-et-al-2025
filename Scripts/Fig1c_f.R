## Fig1 UMAP plots ##

source("Scripts/lib.R")

## load pre-processed objects:
load("Processed_Objects/Inhibitory_datasets.Rdata")
EI_seurat <- readRDS("Processed_Objects/EXCIT_INHIBIT_cleaned_sub.rds")

## for the final plots in the paper the first UMAP axis was flipped

## F1c:
DimPlot(Inhibitory_datasets, reduction = "umap2", group.by = "Experiment2")

## F1d:
DimPlot(Inhibitory_datasets, reduction = "umap2", group.by = "Fine_annotation")

## F1e:
Inhibitory_datasets_sub <- subset(Inhibitory_datasets, subset = Experiment2 == "WT")
DimPlot(Inhibitory_datasets_sub, reduction = "umap2", group.by = "Collection_stage")

## F1f:
DimPlot(EI_seurat, reduction = "umap2", group.by = "Stage_DV2")
