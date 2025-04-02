## cell type abundance of post-mitotic cells across stages ##

source(file = "Scripts/lib.R")
library(data.table)
library(ggalluvial)


## load data:
EXCIT_INHIBIT_cleaned_sub <- readRDS(file = "Processed_Objects/EXCIT_INHIBIT_cleaned_sub.rds")

plot_df <- EXCIT_INHIBIT_cleaned_sub@meta.data[,c("Celltype2","Gene_Annotation","General2","Stage_DV1")]

## 1: ventral telencephalon:
ventral_list_of_tables <- list(
  "12" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E12_I"]),
  "13" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E13_I"]),
  "14" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E14_I"]),
  "15" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E15_I"]),
  "16" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E16_I"])
)
v_blacklist <- c("Apical progenitors", "Basal progenitors (inhibitory)", "Six3_Mpped2", "Nkx2-1_Lhx8", "Ccnd2_Nudt4")
v_celltypes <- unique(unlist(lapply(ventral_list_of_tables, names)))
v_celltypes <- v_celltypes[!v_celltypes %in% v_blacklist]

ventral_list_of_tables <- lapply(ventral_list_of_tables, function(vec) {
  vec <- vec[!names(vec) %in% v_blacklist]
  for(ct in v_celltypes) {
    if(!ct %in% names(vec)) {vec[ct] <- 0}
  }
  return(vec)
})
ventral_list_of_tables <- lapply(ventral_list_of_tables, function(vec) {
  vec <- vec / sum(vec) * 100
  return(vec)
})

tmp_l <- unlist(ventral_list_of_tables)
pd1 <- data.frame("celltype" = sapply(names(tmp_l), function(el) {strsplit(el, "[.]")[[1]][2]}),
                  "relative_number_of_cells" = tmp_l,
                  "stage" = sapply(names(tmp_l), function(el) {strsplit(el, "[.]")[[1]][1]}))

write.table(pd1, "source_data/fig1h_right.csv", sep = ",", quote = FALSE, row.names = FALSE)

g1 <- ggplot(pd1, aes(stage, relative_number_of_cells, col = celltype, group = celltype)) +
  geom_point(size = 3) +
  geom_line() +ylim(0,100) +
  ylab("Relative number of cells") + xlab("Stage") + ggtitle("Inhibitory Trajectory") +
  scale_color_manual(values=c("#54990F", "#78B33E", "#0F8299","#3E9FB3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



## 2: dorsal telencephalon:
dorsal_list_of_tables <- list(
  "12" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E12_E"]),
  "13" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E13_E"]),
  "14" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E14_E"]),
  "15" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E15_E"]),
  "16" = table(EXCIT_INHIBIT_cleaned_sub$Celltype2[EXCIT_INHIBIT_cleaned_sub$Stage_DV2 == "E16_E"])
)
d_blacklist <- c("Apical progenitors", "Basal Progenitors (excitatory)","Immature neurons","Excitatory post-mitotic neurons",
                 "Intermediate progenitors","unknown")
d_celltypes <- unique(unlist(lapply(dorsal_list_of_tables, names)))
d_celltypes <- d_celltypes[!d_celltypes %in% d_blacklist]

dorsal_list_of_tables <- lapply(dorsal_list_of_tables, function(vec) {
  vec <- vec[!names(vec) %in% d_blacklist]
  for(ct in d_celltypes) {
    if(!ct %in% names(vec)) {vec[ct] <- 0}
  }
  return(vec)
})
dorsal_list_of_tables <- lapply(dorsal_list_of_tables, function(vec) {
  vec <- vec / sum(vec) * 100
  return(vec)
})

tmp_l2 <- unlist(dorsal_list_of_tables)
pd2 <- data.frame("celltype" = sapply(names(tmp_l2), function(el) {strsplit(el, "[.]")[[1]][2]}),
                  "relative_number_of_cells" = tmp_l2,
                  "stage" = sapply(names(tmp_l2), function(el) {strsplit(el, "[.]")[[1]][1]}))
write.table(pd2, "source_data/fig1h_left.csv", sep = ",", quote = FALSE, row.names = FALSE)

g2 <- ggplot(pd2, aes(stage, relative_number_of_cells, col = celltype, group = celltype)) +
  geom_point(size = 3) +
  geom_line() +ylim(0,100) +
  ylab("Relative number of cells") + xlab("Stage") + ggtitle("Excitatory Trajectory") +
  scale_color_manual(values = c("#C7B8E6","#3D0F99","#CC7A88","#653EB3","#967ACC")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

wrap_plots(g2,g1)


###########################################################################################
## label transfer ##

load("Processed_Objects//Inhibitory_datasets.Rdata")

## function for label-transfer from Mayer et al.: ##
## https://github.com/ChristophH/in-lineage/blob/master/R/lib.R ##

load("Processed_Objects//transcriptome.integrated.TIS2.Rdata")

P10 <- transcriptome.integrated
remove(transcriptome.integrated)

selection.P10 <- c("2-Inhib Neuron OB Meis2", "6-Inhib Neuron OB Synpr","7a-D2 SPNs",
                   "7b-D1 SPNs","8-Inhib ITC Amygdala",
                   "34-Inhib PN Ventral Striatum/Eac","13a-MGE IN Snhg11",
                   "19b-CGE Neurogliaform IN","19a-CGE VIP IN")


Idents(object = P10) <- "refined_COUP_clust"
P10 <- subset(x = P10, idents = selection.P10, invert = F)
Idents(P10) <- droplevels(Idents(P10))
DefaultAssay(object = P10) <- "RNA"
P10 <- NormalizeData(P10, normalization.method = "LogNormalize", scale.factor = 10000)
P10 <-ScaleData(P10)

## also re-scale inhibitory WT tips data
inh_tip_seurat <- subset(Inhibitory_datasets, subset = Tips != "No")
inh_tip_seurat <- NormalizeData(inh_tip_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
inh_tip_seurat <- ScaleData(inh_tip_seurat)

## function from CH:
transfer.label.cor <- function(expr.a, labels.a, expr.b) {
  # average expr by label
  expr.a.agg <- t(apply(expr.a, 1, function(x) aggregate(x, by=list(labels.a), FUN=mean)$x))
  cmat <- cor(expr.b, expr.a.agg, method='pearson', use='pairwise.complete.obs')
  colnames(cmat) <- paste0('cor.', sort(unique(labels.a)))
  res.df <- data.frame(label=sort(unique(labels.a))[apply(cmat, 1, which.max)], cor=apply(cmat, 1, max))
  # get background correlation based on randomized expr.b
  cor.bg <- sapply(1:100, function(x) apply(cor(matrix(sample(expr.b), nrow(expr.b)), expr.a.agg, method='pearson', use='pairwise.complete.obs'), 1, max))
  res.df$cor.p <- sapply(res.df$cor, function(x) mean(cor.bg >= x))
  res.df$cor.p.adjust <- p.adjust(res.df$cor.p, method = 'fdr')
  res.df <- cbind(res.df, cmat)
  return(res.df)
}

## get average expression per cluster in sticr dataset:
sticr_average_expression <- sapply(unique(P10$refined_COUP_clust), function(cluster_label) {
  apply(P10@assays$RNA@scale.data[, P10$refined_COUP_clust == cluster_label], 1, mean)
})
head(sticr_average_expression)
dim(sticr_average_expression)

## get overlapping HVGs ##
hvg_genes <- rownames(sticr_average_expression)[rownames(sticr_average_expression) %in% VariableFeatures(inh_tip_seurat)]
length(hvg_genes)

## run label transfer ##
inh_map <- transfer.label.cor(
  expr.a = sticr_average_expression[hvg_genes, ],
  labels.a = colnames(sticr_average_expression),
  expr.b = inh_tip_seurat@assays$RNA@scale.data[hvg_genes, ]
)

head(inh_map)
table(inh_map$label, inh_map$cor.p.adjust < 0.1)

## also map back to original ##
inh_map_filtered <- inh_map[inh_map$cor.p.adjust < 0.1, ]
predicted_id_CH <- sapply(colnames(Inhibitory_datasets), function(cell_id) {
  if(cell_id %in% rownames(inh_map_filtered)) {
    inh_map_filtered[cell_id, "label"]
  } else {
    NA
  }
})
## add to object:
Inhibitory_datasets$predicted_id_adult_CH <- predicted_id_CH

## SF2a:
plot_df <- data.frame(
  "cluster" = Inhibitory_datasets$Tips,
  "label" = Inhibitory_datasets$predicted_id_adult_CH
)
plot_df <- plot_df[plot_df$cluster != "No", ]
plot_df$label[is.na(plot_df$label)] <- "uncertain"

plot_df$label <- factor(plot_df$label, levels = c("2-Inhib Neuron OB Meis2", "6-Inhib Neuron OB Synpr","7a-D2 SPNs","7b-D1 SPNs","8-Inhib ITC Amygdala","13a-MGE IN Snhg11","19a-CGE VIP IN","19b-CGE Neurogliaform IN", "34-Inhib PN Ventral Striatum/Eac", "uncertain"))
col_vec <- c("#c7e9c0","#a1d99b","#74c476","#31a354","#e5f5e0","#9ecae1","#6baed6","#4292c6","#006d2c","#d9d9d9")

write.table(plot_df, "../Results/source_data/sf2a.csv", sep = ",", quote = F, row.names = T)
ggplot(plot_df, aes(x = cluster, fill = label)) +
  geom_bar(stat = "count", position = "fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = col_vec)


###########################################################################################
## compare Monocle3 to Seurat clusters ##

load("Processed_Objects/Inhibitory_datasets.Rdata")

ElbowPlot(Inhibitory_datasets)
Inhibitory_datasets <- FindNeighbors(Inhibitory_datasets, dims = 1:15)
Inhibitory_datasets <- FindClusters(Inhibitory_datasets, resolution = 0.25)
DimPlot(Inhibitory_datasets, reduction = "umap2")
Inhibitory_datasets$seurat_clusters_r0_25 <- Idents(Inhibitory_datasets)

cluster_df <- data.frame(
  "monocle_clusters" = Inhibitory_datasets$Fine_annotation,
  "seurat_clusters" = Inhibitory_datasets$seurat_clusters_r0_25
)

cluster_df <- cluster_df[cluster_df$monocle_clusters %in% c("Ebf1_Isl1"," Gucy1a3","Tcf4_Nr2f2","Maf_Sst","Snhg11","Snhg11_Lhx8"), ]
## re-format
cluster_df %>%
  dplyr::count(monocle_clusters, seurat_clusters) -> cluster_df_allu

is_alluvia_form(cluster_df_allu)
write.table(cluster_df_allu, "../Results/source_data/sf2b.csv", sep = ",", quote = F, row.names = F)
ggplot(cluster_df_allu, aes(y = n, axis1 = seurat_clusters, axis2 = monocle_clusters)) +
  geom_alluvium(aes(fill = monocle_clusters), width = 1/12) +
  geom_stratum(width = 1/12, fill = "beige", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Seurat cluster", "Monocle3 cluster"), expand = c(.1, .1)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Relationship between Monocle3 and Seurat clusters") +
  theme_bw()
