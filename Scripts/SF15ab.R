## MILO Inhibitory gNfib/x ##


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(miloR)
library(patchwork)
library(SingleCellExperiment)
library(pals)

gNFI_inh_sub <- readRDS("Processed_Objects/gNFI_inhibitory_sub_wLabels_wPseudotime.rds")

## subset for cells that contain a guide:
gNFI_inh_sub <- subset(gNFI_inh_sub, subset = guide_annot != "no guide")

## re-run pre-processing:
gNFI_inh_sub <- NormalizeData(gNFI_inh_sub, normalization.method = "LogNormalize", scale.factor = 10000)
gNFI_inh_sub <- FindVariableFeatures(gNFI_inh_sub, selection.method = "vst", nfeatures = 2000)
gNFI_inh_sub <- ScaleData(gNFI_inh_sub, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
gNFI_inh_sub <- RunPCA(gNFI_inh_sub, features = VariableFeatures(object = gNFI_inh_sub))
ElbowPlot(gNFI_inh_sub, ndims = 40)

gNFI_inh_sub <- FindNeighbors(gNFI_inh_sub, dims = 1:20)
gNFI_inh_sub <- FindClusters(gNFI_inh_sub, resolution = 0.5)
gNFI_inh_sub <- RunUMAP(gNFI_inh_sub, dims = 1:20)

## plotting:
cluster_state_map <- c("Fabp7_Dbi"="mt","Hist1h1b_Top2a"="mt","Ube2c_Cenpf"="mt","Ccnd2_Abracl"="mt","Adarb2_Npas3"="in","Nxph1_Sst"="in",
                       "Tshz1_Erbb4"="in","Six3_Six3os1"="i-pn","Ptprm_Gucy1a1"="i-pn","Cntnap5b_Dcc"="i-pn","Ebf1_Pou3f1"="i-pn","Cntn5_Cdh8"="i-pn","Gm26917_Grm5"="na")
cell_state_vec <- cluster_state_map[as.character(gNFI_inh_sub$cluster_gene_annot)]
gNFI_inh_sub$cluster_gene_annot2 <- str_replace_all(gNFI_inh_sub$cluster_gene_annot,"_","/")
gNFI_inh_sub$state_cluster_gene_annot <- paste(cell_state_vec, gNFI_inh_sub$cluster_gene_annot2, sep = ":")
table(gNFI_inh_sub$state_cluster_gene_annot)

gNFI_inh_sub$state_cluster_gene_annot <- factor(gNFI_inh_sub$state_cluster_gene_annot, levels = c(
  "mt:Fabp7/Dbi","mt:Hist1h1b/Top2a","mt:Ube2c/Cenpf","mt:Ccnd2/Abracl","in:Adarb2/Npas3","in:Nxph1/Sst",
  "in:Tshz1/Erbb4","i-pn:Six3/Six3os1","i-pn:Ptprm/Gucy1a1","i-pn:Cntnap5b/Dcc","i-pn:Ebf1/Pou3f1","i-pn:Cntn5/Cdh8","na:Gm26917/Grm5"
))
col_vec <- stepped3(length(levels(gNFI_inh_sub$state_cluster_gene_annot)))
names(col_vec) <- levels(gNFI_inh_sub$state_cluster_gene_annot)

col_vec["i-pn:Six3/Six3os1"] <- "#DADAEB"
col_vec["i-pn:Ptprm/Gucy1a1"] <- "#C7E9C0"
col_vec["i-pn:Cntn5/Cdh8"] <- "#31A354"

DimPlot(gNFI_inh_sub, reduction = "umap", group.by = "state_cluster_gene_annot", label = T, cols = col_vec)


## convert object:
gNFI_inh_sce <- as.SingleCellExperiment(gNFI_inh_sub)
gNFI_inh_sce
gNFI_inh_milo <- Milo(gNFI_inh_sce)


gNFI_inh_milo <- buildGraph(gNFI_inh_milo, k = 40, d = 20)

gNFI_inh_milo <- makeNhoods(gNFI_inh_milo, prop = 0.1, k = 40, d=20, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(gNFI_inh_milo)

nhood_mtx <- nhoods(gNFI_inh_milo)
nhood_mtx_gNFI <- nhood_mtx[gNFI_inh_milo@colData$guide_annot == "gNFI", ]
nhood_mtx_gLacZ <- nhood_mtx[gNFI_inh_milo@colData$guide_annot == "gLacZ", ]

hist(colSums(nhood_mtx_gNFI), breaks = 100)
hist(colSums(nhood_mtx_gLacZ), breaks = 100)

## do counting:
colData(gNFI_inh_milo)$sample_annot <- paste(colData(gNFI_inh_milo)$guide_annot, colData(gNFI_inh_milo)$replicate, sep = "__")

gNFI_inh_milo <- countCells(gNFI_inh_milo, meta.data = as.data.frame(colData(gNFI_inh_milo)), sample="sample_annot")
head(nhoodCounts(gNFI_inh_milo))

gNFI_design <- data.frame(colData(gNFI_inh_milo))[,c("sample_annot", "guide_annot", "replicate")]

## Convert batch info from integer to factor
gNFI_design$replicate <- as.factor(gNFI_design$replicate) 
gNFI_design <- distinct(gNFI_design)
rownames(gNFI_design) <- gNFI_design$sample_annot

gNFI_design

## compute neighbourhood connectivity:
gNFI_inh_milo <- calcNhoodDistance(gNFI_inh_milo, d=20, reduced.dim = "PCA")

## testing:
da_results <- testNhoods(gNFI_inh_milo, design = ~ replicate + guide_annot, design.df = gNFI_design)
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
gNFI_inh_milo <- buildNhoodGraph(gNFI_inh_milo)
## Plot neighbourhood graph
plotNhoodGraphDA(gNFI_inh_milo, da_results, layout="UMAP",alpha=0.1) + 
  theme(legend.position="none")


da_results <- annotateNhoods(gNFI_inh_milo, da_results, coldata_col = "cluster_gene_annot")
head(da_results)

da_results$cluster_gene_annot <- ifelse(da_results$cluster_gene_annot_fraction < 0.7, "Mixed", da_results$cluster_gene_annot)
plotDAbeeswarm(da_results, group.by = "cluster_gene_annot")
