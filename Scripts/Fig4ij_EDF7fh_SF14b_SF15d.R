## Nfib OE ##

library(Seurat)
library(ggplot2)
library(pals)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggrepel)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(ggpubr)

NFI_OE_seurat <- readRDS("Processed_Objects/nfib_oe_s3_harmony2b.rds")

cell_state_vec <- as.character(NFI_OE_seurat$cluster_annotation_broad)
cell_state_vec[cell_state_vec %in% c("AP", "exc_BP","inh_BP")] <- "mt"
cell_state_vec[cell_state_vec == "Interneurons"] <- "in"; cell_state_vec[cell_state_vec == "Projection neurons"] <- "i-pn"
cell_state_vec[cell_state_vec == "exc_precursor"] <- "e-pn"; cell_state_vec[cell_state_vec == "unknown"] <- "na"

NFI_OE_seurat$cluster_annotation <- str_replace_all(NFI_OE_seurat$cluster_annotation, "_","/")
NFI_OE_seurat$state_cluster_annotation <- paste(cell_state_vec, NFI_OE_seurat$cluster_annotation, sep = ":")

## fig 4f,g: UMAPs

DimPlot(NFI_OE_seurat, reduction = "umap", group.by = "predicted_EI_annotation", label = T)
FeaturePlot(NFI_OE_seurat, reduction = "umap", features = "pseudotime")


## fig 4i,j: cluster abundance, pseudotime

CellComp_Poisson<-function(seurat_obj,celltype="seurat_clusters_r_0_5",perturbations="plasmid_group2",batch="Replicate",cutoff=10) {
  res <- list()
  meta=seurat_obj@meta.data[,c(celltype, perturbations, batch)]
  colnames(meta)=c("CellType","Pert","Batch")
  
  ## remove cells with no guide:
  meta <- meta[meta$Pert != "others", ]
  
  meta <- meta %>% group_by(CellType,Pert,Batch) %>% summarise(Num=length(CellType)) %>% as.data.frame()
  meta2<-meta %>% group_by(Pert,Batch) %>% summarise(Tot=sum(Num)) %>% as.data.frame()
  meta=left_join(meta,meta2)
  
  meta=meta[meta[,"Tot"]>cutoff,]
  
  meta["Pert"]=relevel(factor(meta[,"Pert"]),ref="egfp")
  
  lst=list()
  for(i in unique(meta[,"CellType"])){lst[[i]]=meta[meta[,"CellType"]==i,]}
  
  out <- lapply(lst,function(cur){
    celltype <- cur[1,"CellType"]
    
    ## check if celltype is contained across all replicates:
    if(nrow(cur) < 2*length(unique(meta$Batch))) {
      print(paste0("Skipping celltype: ", celltype))
      return(NULL)
    } else {
      print(as.character(celltype))
      cur["logTot"]=log(cur[,"Tot"])
      fit<-glm(Num~offset(logTot)+Batch+Pert,data=cur,family="poisson")
      tab=summary(fit)
      tab=tab$coefficients
      tab=data.frame(tab)
      tab=tab[grep("Pert",rownames(tab)),]
      tab["Gene"]=sub("Pert","",rownames(tab))
      tab["CellType"]=celltype
      tab=tab[,c(5,6,4,1,2,3)]
      colnames(tab)[3]="pval"
      return(tab)
    }
  })
  tab=do.call(rbind,out)
  rownames(tab)=NULL
  tab=tab[order(tab[,"pval"]),]
  tab["padj"]=p.adjust(tab[,"pval"],"fdr")
  tab$CellType <- as.character(tab$CellType)
  res$tab <- tab
  
  meta <- meta[meta$CellType %in% tab$CellType, ]
  res$meta <- meta
  
  ## also create proportion dataframe:
  ## xNfib and egfp are hard-coded!
  prop_df <- left_join(
    x = meta[meta$Pert == "xNfib", ],
    y = meta[meta$Pert == "egfp", ],
    by = c("CellType", "Batch")
  )
  prop_df$prop_xNfib <- prop_df$Num.x / prop_df$Tot.x
  prop_df$prop_egfp <- prop_df$Num.y / prop_df$Tot.y
  
  prop_df$prop_change <- log10(prop_df$prop_xNfib / prop_df$prop_egfp)
  res$prop <- prop_df
  
  return(res)
}

nfib_oe_sub <- subset(NFI_OE_seurat, subset = predicted_EI_annotation != "unknown")

## PROPORTION CHANGE ##
nfib_oe_EI_label_composition <- CellComp_Poisson(
  seurat_obj = nfib_oe_sub,
  celltype = "predicted_EI_annotation",
  perturbations = "plasmid_group2",
  batch = "Replicate",
  cutoff = 10
)

prop_df <- nfib_oe_EI_label_composition$prop

cell_state_map <- c(
  "Hes1_Fabp7"="mt","Fabp7_Mt3"="mt","Gas1_Ldha"="mt","Hist1h1b_Top2a"="mt","Ccnd2_Nudt4"="i-pn","Isl1_Zfp503"="i-pn","Ebf1_Foxp1"="i-pn",
  "Foxp1_Gucy1a3"="i-pn","Nr2f2_Nr2f1"="in","Nkx2-1_Lhx8"="in","Npy_Nxph1"="in","Sst_Maf"="in","Neurog2_Rrm2"="mt","Neurog2_Eomes"="mt",
  "Neurod2_Neurod6"="e-pn","Neurod6_Mef2c"="e-pn","na"
)

prop_df$CellType <- paste(cell_state_map[as.character(prop_df$CellType)], prop_df$CellType, sep = ":")
prop_df$CellType <- stringr::str_replace_all(prop_df$CellType, "_", "/")
prop_df$CellType <- factor(prop_df$CellType, levels = c(
  "mt:Hes1/Fabp7","mt:Fabp7/Mt3","mt:Gas1/Ldha","mt:Hist1h1b/Top2a","in:Nr2f2/Nr2f1","in:Nkx2-1/Lhx8","in:Npy/Nxph1","in:Sst/Maf",
  "i-pn:Ccnd2/Nudt4","i-pn:Isl1/Zfp503","i-pn:Ebf1/Foxp1","i-pn:Foxp1/Gucy1a3",
  "mt:Neurog2/Rrm2","mt:Neurog2/Eomes","e-pn:Neurod2/Neurod6","e-pn:Neurod6/Mef2c"
))

prop_df_sub <- prop_df[!prop_df$CellType %in% c("mt:Neurog2/Rrm2","mt:Neurog2/Eomes","e-pn:Neurod2/Neurod6","e-pn:Neurod6/Mef2c"), ]

write.table(prop_df_sub, "Results/source_data/fig4i.csv", sep = ",", quote = F, row.names = F)

ggplot(prop_df_sub, aes(x = CellType, y = prop_change)) +
  geom_point(size = 2) +
  theme_bw() +
  ylab("log10(xNfib-fraction / eGFP-fraction)") +
  xlab("Predicted label") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_hline(yintercept = 0, linetype = 'dotted')

## PSEUDOTIME ##
pt_plot_df <- data.frame(
  "plasmid" = NFI_OE_seurat$plasmid_group2,
  "pseudotime" = NFI_OE_seurat$pseudotime,
  "broad_label" = NFI_OE_seurat$predicted_broad_annotation2
)

pt_plot_df <- pt_plot_df[pt_plot_df$plasmid != "others",]
pt_plot_df <- pt_plot_df[pt_plot_df$broad_label %in% c("AP","inh_BP","Interneurons","Projection neurons"), ]
pt_plot_df$plasmid <- factor(pt_plot_df$plasmid, levels = c("egfp","xNfib"))

bc_col_vec <- c("AP" = "#6BAED6", "inh_BP" = "#C6DBEF", "Interneurons" = "#FD8D3C", "Projection neurons" = "#74C476")

write.table(pt_plot_df, "Results/source_data/fig4j_lower.csv", sep = ",", quote = F, row.names = T)

ggplot(pt_plot_df, aes(x = plasmid, y = pseudotime)) +
  geom_violin(aes(fill = broad_label)) +
  theme_bw() +
  facet_wrap(~broad_label) +
  geom_signif(comparisons = list(c("egfp","xNfib")), test = "wilcox.test", map_signif_level = TRUE) +
  stat_summary(fun.y=median, geom="point", size=2, color="black") +
  scale_fill_manual(values = bc_col_vec)

## EDF 7f,h: broad abundance, num DE genes
nfib_oe_sub$predicted_broad_annotation3 <- as.character(nfib_oe_sub$predicted_broad_annotation2)
nfib_oe_sub$predicted_broad_annotation3[nfib_oe_sub$predicted_broad_annotation2 %in% c("AP","exc_BP","inh_BP")] <- "mitotic"

nfib_oe_broad_label_composition <- CellComp_Poisson(
  seurat_obj = nfib_oe_sub,
  celltype = "predicted_broad_annotation3",
  perturbations = "plasmid_group2",
  batch = "Replicate",
  cutoff = 10
)

prop_df <- nfib_oe_broad_label_composition$prop
prop_df$CellType <- factor(prop_df$CellType, levels = c(
  "mitotic","Interneurons","Projection neurons","exc_precursor"
))

write.table(prop_df, "Results/source_data/edf7f.csv", sep = ",", quote = F, row.names = F)
ggplot(prop_df, aes(x = CellType, y = prop_change)) +
  geom_point(size = 2) +
  theme_bw() +
  ylab("log10(xNfib-fraction / eGFP-fraction)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 'dotted')

## DE genes:
DE_gene_list <- list()
for(cluster_name in unique(NFI_OE_seurat$cluster_annotation)) {
  print(cluster_name)
  cluster_sub <- subset(NFI_OE_seurat, subset = cluster_annotation == cluster_name)
  Idents(cluster_sub) <- "plasmid_group2"
  cluster_markers <- FindMarkers(cluster_sub, ident.1 = "xNfib", ident.2 = "egfp", test.use = "wilcox")
  ## more annot:
  cluster_markers$cluster <- cluster_name
  cluster_markers$gene <- rownames(cluster_markers)
  cluster_markers$upregulated <- cluster_markers$avg_log2FC > 0
  ## filter:
  cluster_markers <- cluster_markers[cluster_markers$p_val_adj < 0.01, ]
  DE_gene_list[[cluster_name]] <- cluster_markers
}

all_marker_df <- do.call(rbind.data.frame, DE_gene_list)

all_marker_df %>%
  group_by(cluster, upregulated) %>%
  summarize(number_of_genes = n()) -> plot_df

write.table(plot_df, "Results/source_data/edf7h.csv", quote = F, sep = ",", row.names = F)

ggplot(plot_df, aes(cluster, number_of_genes, fill = upregulated)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c("FALSE" = "#fc8d62", "TRUE" = "#8da0cb")) +
  ylab("number of DE genes") + xlab("cluster")


## SF 14b: clusters per replicate
plot_df <- NFI_OE_seurat@meta.data[, c("Replicate", "state_cluster_annotation")]
plot_df$state_cluster_annotation <- factor(plot_df$state_cluster_annotation, levels = c(
  "mt:Fabp7/Dbi","mt:Creb5/Gm29260","mt:Top2a/Sox1ot","i-pn:Gm38505/Sox2ot","in:Erbb4/Nxph1","i-pn:Ebf1/Mgat4c","na:Tuba1a/Tubb2b","na:Gm42418/Camk1d",
  "na:Galntl6/Adarb2","na:Reln/Clstn2","na:Apoe/C1qb","mt:Gm29260/Nrg1", "mt:Plcb1/Unc5d","e-pn:Adgrl3/Pcdh9","e-pn:Tiam2/Satb2","na:Opcml/Lrfn5"
))

write.table(plot_df, "Results/source_data/sf14b.csv", sep = ",", quote = F, row.names = F)
ggplot(plot_df, aes(Replicate, fill = state_cluster_annotation)) +
  geom_bar(stat = "count", position = "fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = stepped3(n = 16)) +
  ylab("percent cells per replicate")


## SF 15d: DE genes in C&R
## TODO !!
promoter_sub <- readRDS("Processed_Objects/cutnrun_promoter_ranges.rds")
Nfib_OE_DE_genes <- readRDS("Processed_Objects/NFI_OE_plasmid_markers_per_cluster_wilcox.rds")


Nfib_OE_DE_gene_vec <- unique(Nfib_OE_DE_genes$gene)
Nfib_OE_bound <- Nfib_OE_DE_gene_vec %in% promoter_sub$gene_symbol

NFI_OE_up <- unique(Nfib_OE_DE_genes$gene[Nfib_OE_DE_genes$upregulated == T])
NFI_OE_down <- unique(Nfib_OE_DE_genes$gene[Nfib_OE_DE_genes$upregulated == F])

ggvenn(list("Nfib_OE_UP" = NFI_OE_up, "Nfib_OE_DOWN" = NFI_OE_down, "Nfib_CnR" = unique(promoter_sub$gene_symbol)), fill_color = c("#8da0cb","#fc8d62","#bd0026"))
