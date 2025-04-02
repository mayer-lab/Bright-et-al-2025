## gNFI KO ##

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
library(ggsignif)
library(ggpubr)

## load objects:
gNFI_seurat <- readRDS("Processed_Objects/gNFI_merged_seurat_wLabel2.rds")
gNFI_inh_sub <- readRDS("Processed_Objects/gNFI_inhibitory_sub_wLabels_wPseudotime.rds")
NFI_OE_seurat <- readRDS("Processed_Objects/nfib_oe_s3_harmony2b.rds")
EI_seurat <- readRDS("Processed_Objects/EXCIT_INHIBIT_cleaned_sub.rds")

## fig 4d,e: UMAPS
DimPlot(gNFI_seurat, reduction = "umap", group.by = "cluster_gene_annot", label = T)
FeaturePlot(gNFI_inh_sub, reduction = "umap", features = "pseudotime")

## fig 4h,j,k: cluster abundance, pseudotime & heatmap

## ABUNDANCE CHANGE ##
x <- c("replicate 1" = "experiment 1", "replicate 2" = "experiment 1", "replicate 3" = "experiment 2", "replicate 4" = "experiment 2")
gNFI_seurat$biological_replicate <- x[gNFI_seurat$replicate]

CellComp_Poisson<-function(seurat_obj,celltype="cluster_gene_annot",perturbations="guide_annot",batch="biological_replicate",cutoff=10) {
  res <- list()
  meta=seurat_obj@meta.data[,c(celltype, perturbations, batch)]
  colnames(meta)=c("CellType","Pert","Batch")
  
  ## remove cells with no guide:
  meta <- meta[meta$Pert != "no guide", ]
  
  meta <- meta %>% group_by(CellType,Pert,Batch) %>% summarise(Num=length(CellType)) %>% as.data.frame()
  meta2<-meta %>% group_by(Pert,Batch) %>% summarise(Tot=sum(Num)) %>% as.data.frame()
  meta=left_join(meta,meta2)

  meta=meta[meta[,"Tot"]>cutoff,]
  
  meta["Pert"]=relevel(factor(meta[,"Pert"]),ref="gLacZ")
  
  lst=list()
  for(i in unique(meta[,"CellType"])){lst[[i]]=meta[meta[,"CellType"]==i,]}
  
  out <- lapply(lst,function(cur){
    celltype=cur[1,"CellType"]
    
    print(celltype)
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
  ## gNFI and gLAcZ are hard-coded!
  prop_df <- left_join(
    x = meta[meta$Pert == "gNFI", ],
    y = meta[meta$Pert == "gLacZ", ],
    by = c("CellType", "Batch")
  )
  prop_df$prop_gNFI <- prop_df$Num.x / prop_df$Tot.x
  prop_df$prop_gLacZ <- prop_df$Num.y / prop_df$Tot.y
  
  prop_df$prop_change <- log10(prop_df$prop_gNFI / prop_df$prop_gLacZ)
  res$prop <- prop_df
  
  return(res)
}

gNFI_cluster_composition <- CellComp_Poisson(
  seurat_obj = gNFI_seurat,
  celltype = "cluster_gene_annot",
  perturbations = "guide_annot",
  batch = "biological_replicate",
  cutoff = 10
)

prop_df <- gNFI_cluster_composition$prop
cell_state_map <- c(
  "Fabp7_Dbi"="mt","Top2a_Hmgb2"="mt","Ccnd2_Abracl"="mt","Adarb2_Erbb4"="in","Nxph1_Sst"="in","Tshz1_Erbb4"="in","Six3_Six3os1"="i-pn",
  "Ptprm_Gucy1a1"="i-pn","Cntnap5b_Isl1"="i-pn","Gm26917_Lingo2"="i-pn","Ebf1_Cntn5"="i-pn",
  "Gm29260_Hist1h1b"="mt","Unc5d_Nrg1"="mt","Satb2_9130024F11Rik"="e-pn","Tafa1_Adgrl3"="e-pn","Kcnip4_Nrg3"="e-pn"
)
prop_df$CellType <- paste(cell_state_map[prop_df$CellType], prop_df$CellType, sep = ":")
prop_df$CellType <- stringr::str_replace_all(prop_df$CellType, "_", "/")
prop_df$CellType <- factor(prop_df$CellType, levels = c(
  "mt:Fabp7/Dbi","mt:Top2a/Hmgb2","mt:Ccnd2/Abracl","in:Adarb2/Erbb4","in:Nxph1/Sst","in:Tshz1/Erbb4","i-pn:Six3/Six3os1","i-pn:Ptprm/Gucy1a1",
  "i-pn:Cntnap5b/Isl1","i-pn:Gm26917/Lingo2","i-pn:Ebf1/Cntn5", "mt:Gm29260/Hist1h1b","mt:Unc5d/Nrg1","e-pn:Satb2/9130024F11Rik","e-pn:Tafa1/Adgrl3","e-pn:Kcnip4/Nrg3"
))

prop_df_sub <- prop_df[! prop_df$CellType %in% c("mt:Gm29260/Hist1h1b","mt:Unc5d/Nrg1","e-pn:Satb2/9130024F11Rik","e-pn:Tafa1/Adgrl3","e-pn:Kcnip4/Nrg3"), ]

write.table(prop_df_sub, "Results/source_data/fig4h.csv", sep = ",", quote = F, row.names = F)

ggplot(prop_df_sub, aes(x = CellType, y = prop_change)) +
  geom_point(size = 2) +
  theme_bw() +
  ylab("log10(gNFI-fraction / gLacZ-fraction)") +
  xlab("Cluster") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 'dotted')


## PSEUDOTIME ##
borad_cluster_vec <- c(
  "Fabp7_Dbi" = "AP","Hist1h1b_Top2a" = "BP","Ube2c_Cenpf" = "BP","Ccnd2_Abracl" = "BP","Adarb2_Npas3" = "IN","Nxph1_Sst" = "IN","Tshz1_Erbb4" = "IN",
  "Six3_Six3os1" = "PN","Ptprm_Gucy1a1" = "PN","Cntnap5b_Dcc" = "PN","Ebf1_Pou3f1" = "PN","Gm26917_Grm5" = "unknown","Cntn5_Cdh8" = "PN"
)
gNFI_inh_sub$cluster_broad_annot <- borad_cluster_vec[gNFI_inh_sub$cluster_gene_annot]

pt_plot_df <- data.frame(
  "guide" = gNFI_inh_sub$guide_annot,
  "pseudotime" = gNFI_inh_sub$pseudotime,
  "broad_cluster" = gNFI_inh_sub$cluster_broad_annot
)

pt_plot_df <- pt_plot_df[pt_plot_df$guide != "no guide",]
pt_plot_df$guide <- factor(pt_plot_df$guide, levels = c("gLacZ","gNFI"))

write.table(pt_plot_df, "Results/source_data/fig4j_upper.csv", quote = F, sep = ",", row.names = T)

ggplot(pt_plot_df, aes(guide, pseudotime)) +
  geom_violin() +
  theme_bw() +
  facet_wrap(~broad_cluster) +
  geom_signif(comparisons = list(c("gLacZ","gNFI")), test = "wilcox.test", map_signif_level = TRUE) +
  stat_summary(fun.y=median, geom="point", size=2, color="red")

## HEATMAP ##
plot_difference_heatmap <- function(gNFI_sub, NFI_OE_sub, gene_list, merged_DE_df, cluster_cols = FALSE) {
  gene_list <- lapply(gene_list, function(l) {
    l[l %in% rownames(gNFI_sub) & l %in% rownames(NFI_OE_sub)]
  })
  
  all_genes <- unlist(gene_list)
  ## get average expression
  gNFI_avg_expr <- AverageExpression(gNFI_sub, slot = "data", assays = "RNA", group.by = "condition_cluster")
  NFI_OE_avg_expr <- AverageExpression(NFI_OE_sub, slot = "data", assays = "RNA", group.by = "condition_cluster")
  
  ## select for genes:
  gNFI_gene_expr <- as.data.frame(gNFI_avg_expr$RNA[all_genes, ])
  NFI_OE_gene_expr <- as.data.frame(NFI_OE_avg_expr$RNA[all_genes, ])
  
  ## transform:
  gNFI_mtx <- as.data.frame(t(as.matrix(gNFI_gene_expr)))
  NFI_OE_mtx <- as.data.frame(t(as.matrix(NFI_OE_gene_expr)))
  
  ## add some more metadata:
  gNFI_mtx$condition <- sapply(rownames(gNFI_mtx), function(el) {strsplit(el,"__")[[1]][1]})
  gNFI_mtx$cluster <- sapply(rownames(gNFI_mtx), function(el) {strsplit(el,"__")[[1]][2]})
  
  NFI_OE_mtx$condition <- sapply(rownames(NFI_OE_mtx), function(el) {strsplit(el,"__")[[1]][1]})
  NFI_OE_mtx$cluster <- sapply(rownames(NFI_OE_mtx), function(el) {strsplit(el,"__")[[1]][2]})
  
  ## create merged dataframe for plotting:
  gene_dfs <- list()
  
  for(gene_name in all_genes) {
    #print(gene_name)
    ## gNFI top
    upper_df <- data.frame(
      "expression_diff" = gNFI_mtx[gNFI_mtx$condition == "perturb", gene_name] - gNFI_mtx[gNFI_mtx$condition == "control", gene_name],
      row.names = gNFI_mtx$cluster[gNFI_mtx$condition == "control"]
    )
    colnames(upper_df) <- gene_name
    ## OE bottom
    lower_df <- data.frame(
      "expression_diff" = NFI_OE_mtx[NFI_OE_mtx$condition == "perturb", gene_name] - NFI_OE_mtx[NFI_OE_mtx$condition == "control", gene_name],
      row.names = NFI_OE_mtx$cluster[NFI_OE_mtx$condition == "control"]
    )
    colnames(lower_df) <- gene_name
    ## merge:
    df <- rbind(upper_df, lower_df)
    gene_dfs[[gene_name]] <- df
  }
  
  merged_plot_df <- data.frame(do.call('cbind', gene_dfs))
  
  ## order rows:
  gNFI_levels <- c("Fabp7_Dbi_1","Top2a_Hmgb2","Ccnd2_Abracl","Adarb2_Erbb4","Nxph1_Sst","Tshz1_Erbb4","Six3_Six3os1","Ptprm_Gucy1a1","Cntnap5b_Isl1","Ebf1_Cntn5")
  NFI_OE_levels <- c("Fabp7_Dbi_2","Creb5_Gm29260","Top2a_Sox1ot","Gm38505_Sox2ot","Erbb4_Nxph1","Ebf1_Mgat4c")
  
  merged_plot_df <- merged_plot_df[c(gNFI_levels, NFI_OE_levels), ]
  
  ## annotation
  row_annot <- rowAnnotation(
    broad_cluster = c("AP","BP","BP","IN","IN","IN","PN","PN","PN","PN",
                      "AP","BP","BP","BP","IN","PN"),
    experiment = c(rep("gNFI", 10), rep("Nfib_OE",6)),
    col = list(
      broad_cluster = c("AP" = "#6BAED6", "BP" = "#C6DBEF", "IN" = "#FD8D3C", "PN" = "#74C476"),
      experiment = c("gNFI" = "#bebada", "Nfib_OE" = "#fb8072")
    )
  )
  
  col_annot <- columnAnnotation(
    gene_list = unlist(sapply(1:length(gene_list), function(i) {rep(names(gene_list)[i], length(gene_list[[i]]))})),
    col = list(
      gene_list = c("STEM-STUFF" = "#80b1d3", "DLX"= "#ffff99", "MIGRATION" = "#e41a1c", "INTERNEURON" = "#FD8D3C", "PN" = "#74C476")
    )
  )
  
  ## column breaks:
  col_breaks <- unlist(sapply(1:length(gene_list), function(i) {rep(names(gene_list)[i], length(gene_list[[i]]))}))
  col_breaks <- factor(col_breaks, levels = names(gene_list))
  ## row breaks:
  row_breaks <- c(rep("gNFI", 10), rep("Nfib_OE",6))
  
  ## color function:
  max_plot_value <- 5
  col_fun = colorRamp2(c(-max_plot_value,0,max_plot_value), c("#2166ac", "#f7f7f7", "#b2182b"))
  
  merged_plot_mtx <- as.matrix(merged_plot_df)
  #hist(as.numeric(merged_plot_mtx), breaks = 100)
  
  write.table(merged_plot_mtx, "Results/source_data/fig4k.csv", sep = ",", quote = F)
  
  if(!cluster_cols) {
    Heatmap(
      merged_plot_mtx,
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      left_annotation = row_annot,
      row_split = row_breaks,
      top_annotation = col_annot,
      column_split = col_breaks,
      show_column_names = TRUE,
      cell_fun = function(j, i, x, y, width, height, fill) {
        cluster_name <- rownames(merged_plot_mtx)[i]
        gene_name <- colnames(merged_plot_mtx)[j]
        cluster_de_genes <- merged_DE_df$gene[merged_DE_df$cluster == cluster_name]
        if(gene_name %in% cluster_de_genes) {
          grid.text("*", x, y, gp = gpar(fontsize = 10))
        } else {
          grid.text("", x, y, gp = gpar(fontsize = 10))
        }
        
      }
    )
  } else {
    ## or with clustered columns:
    Heatmap(
      merged_plot_mtx,
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      left_annotation = row_annot,
      row_split = row_breaks,
      top_annotation = col_annot,
      show_column_names = TRUE,
      column_names_gp = grid::gpar(fontsize = 4),
      #row_names_gp = grid::gpar(fontsize = 8),
      cell_fun = function(j, i, x, y, width, height, fill) {
        cluster_name <- rownames(merged_plot_mtx)[i]
        gene_name <- colnames(merged_plot_mtx)[j]
        cluster_de_genes <- merged_DE_df$gene[merged_DE_df$cluster == cluster_name]
        if(gene_name %in% cluster_de_genes) {
          grid.text("*", x, y, gp = gpar(fontsize = 10))
        } else {
          grid.text("", x, y, gp = gpar(fontsize = 10))
        }
        
      }
    )
  }
}

gNFI_DE_gene_df <- readRDS("Processed_Objects/gNFI_guide_markers_per_cluster_wilcox.rds")
NFI_OE_DE_gene_df <- readRDS("Processed_Objects/NFI_OE_plasmid_markers_per_cluster_wilcox.rds")

gNFI_DE_gene_df$experiment <- "gNFI"
NFI_OE_DE_gene_df$experiment <- "NFI_OE"
gNFI_DE_gene_df$cluster[gNFI_DE_gene_df$cluster == "Fabp7_Dbi"] <- "Fabp7_Dbi_1"
NFI_OE_DE_gene_df$cluster[NFI_OE_DE_gene_df$cluster == "Fabp7_Dbi"] <- "Fabp7_Dbi_2"

merged_DE_df <- rbind(gNFI_DE_gene_df, NFI_OE_DE_gene_df)

gNFI_sub <- subset(gNFI_seurat, subset = guide_annot != "no guide")
NFI_OE_sub <- subset(NFI_OE_seurat, subset = plasmid_group2 != "others")
# inhibitory clusters:
gNFI_sub <- subset(gNFI_sub, subset = cluster_gene_annot %in% c(
  "Fabp7_Dbi","Top2a_Hmgb2","Ccnd2_Abracl","Adarb2_Erbb4","Nxph1_Sst","Tshz1_Erbb4",
  "Six3_Six3os1","Ptprm_Gucy1a1","Cntnap5b_Isl1","Ebf1_Cntn5"
))
NFI_OE_sub <- subset(NFI_OE_sub, subset = cluster_annotation %in% c(
  "Fabp7_Dbi","Creb5_Gm29260","Top2a_Sox1ot","Gm38505_Sox2ot","Erbb4_Nxph1","Ebf1_Mgat4c"
))
## update metadata:
gNFI_cond_vec <- c("gNFI" = "perturb", "gLacZ" = "control")
gNFI_sub$condition <- gNFI_cond_vec[gNFI_sub$guide_annot]

NFI_OE_cond_vec <- c("xNfib" = "perturb", "egfp" = "control")
NFI_OE_sub$condition <- NFI_OE_cond_vec[as.vector(NFI_OE_sub$plasmid_group2)]

gNFI_sub$cluster_gene_annot[gNFI_sub$cluster_gene_annot == "Fabp7_Dbi"] <- "Fabp7_Dbi_1"
NFI_OE_sub$cluster_annotation[NFI_OE_sub$cluster_annotation == "Fabp7_Dbi"] <- "Fabp7_Dbi_2"

gNFI_sub$condition_cluster <- paste(gNFI_sub$condition, gNFI_sub$cluster_gene_annot, sep = "__")
NFI_OE_sub$condition_cluster <- paste(NFI_OE_sub$condition, NFI_OE_sub$cluster_annotation, sep = "__")

## plot gene list:
slim_gene_list <- list(
  "STEM-STUFF" = c("Nes","Vim","Tuba1b","Tubb2a","Sox2","Ezh2","Ccnd2","Ascl1"),
  "DLX" = c("Dlx1","Dlx2","Dlx5"),
  "MIGRATION" = c("Epha5","Robo2","Slit1","Pak3","Ctnna2","Shtn1"),
  "INTERNEURON" = c("Tcf4","Erbb4","Maf","Prox1"),
  "PN" = c("Meis2","Gucy1a3","Sp9","Ebf1","Isl1")
)

plot_difference_heatmap(gNFI_sub, NFI_OE_sub, slim_gene_list, merged_DE_df, cluster_cols = F)


## EDF 7e,g: broad abundance, num DE genes

gNFI_EI_broad_cluster_composition <- CellComp_Poisson(
  seurat_obj = gNFI_seurat,
  celltype = "cluster_broad_annot2",
  perturbations = "guide_annot",
  batch = "biological_replicate",
  cutoff = 10
)

prop_df <- gNFI_EI_broad_cluster_composition$prop
prop_df$CellType <- factor(prop_df$CellType, levels = c(
  "mitotic","interneuron","projection_neuron","exc_PN"
))

write.table(prop_df, "Results/source_data/edf7e.csv", sep = ",", quote = F, row.names = F)

ggplot(prop_df, aes(x = CellType, y = prop_change)) +
  geom_point(size = 2) +
  theme_bw() +
  ylab("log10(gNFI-fraction / gLacZ-fraction)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  xlab("Broad Cluster") +
  geom_hline(yintercept = 0, linetype = 'dotted')


DE_gene_list <- list()
for(cluster_name in unique(gNFI_seurat$cluster_gene_annot)) {
  print(cluster_name)
  cluster_sub <- subset(gNFI_seurat, subset = cluster_gene_annot == cluster_name)
  Idents(cluster_sub) <- "guide_annot"
  cluster_markers <- FindMarkers(cluster_sub, ident.1 = "gNFI", ident.2 = "gLacZ", test.use = "wilcox")
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

write.table(plot_df, "Results/source_data/edf7g.csv", quote = F, sep = ",", row.names = F)

ggplot(plot_df, aes(cluster, number_of_genes, fill = upregulated)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c("FALSE" = "#fc8d62", "TRUE" = "#8da0cb")) +
  ylab("number of DE genes") + xlab("cluster")



## SF 14a: clusters per replicate

plot_df <- gNFI_seurat@meta.data[, c("replicate", "cluster_gene_annot")]
plot_df$cluster_gene_annot <- factor(plot_df$cluster_gene_annot, levels = c(
  "Fabp7_Dbi","Top2a_Hmgb2","Ccnd2_Abracl","Adarb2_Erbb4","Nxph1_Sst","Tshz1_Erbb4","Gm26917_Lingo2","Six3_Six3os1","Ptprm_Gucy1a1","Cntnap5b_Isl1","Ebf1_Cntn5",
  "Gm29260_Hist1h1b","Unc5d_Nrg1","Satb2_9130024F11Rik","Kcnip4_Nrg3","Tafa1_Adgrl3"
))
fill_vec <- stepped3(n = 16)
names(fill_vec) <- levels(plot_df$cluster_gene_annot)

write.table(plot_df, "Results/source_data/sf14a.csv", quote = F, sep = ",", row.names = F)
ggplot(plot_df, aes(replicate, fill = cluster_gene_annot)) +
  geom_bar(stat = "count", position = "fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = fill_vec) +
  ylab("percent cells per replicate")

## SF 15c: DE genes in C&R
## TODO!!
promoter_sub <- readRDS("Processed_Objects/cutnrun_promoter_ranges.rds")
gNFI_DE_genes <- readRDS("Processed_Objects/gNFI_guide_markers_per_cluster_wilcox.rds")

gNFI_DE_gene_vec <- unique(gNFI_DE_genes$gene)

gNFI_up <- unique(gNFI_DE_genes$gene[gNFI_DE_genes$upregulated == T])
gNFI_down <- unique(gNFI_DE_genes$gene[gNFI_DE_genes$upregulated == F])

ggvenn(list("gNFI_UP" = gNFI_up, "gNFI_DOWN" = gNFI_down, "Nfib_CnR" = unique(promoter_sub$gene_symbol)), fill_color = c("#8da0cb","#fc8d62","#bd0026"))

