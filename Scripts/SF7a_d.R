## lineage z-score analysis of GE e12/e16 data ##

library(Seurat)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)

## load data:
load("E12_lineage.Rdata")
load("E16_lineage.Rdata")

## clonal distribution per cluster:
E12_lineage$has_cloneID <- !is.na(E12_lineage$cloneID.umi6)
E16_lineage$has_cloneID <- !is.na(E16_lineage$cloneID.LINEAGE)

plot_df <- data.frame(
  "has_clone_id" = c(E12_lineage$has_cloneID, E16_lineage$has_cloneID),
  "stage" = c(rep("e12", ncol(E12_lineage)), rep("e16", ncol(E16_lineage))),
  "fine_annotation" = c(E12_lineage$Fine_annotation, E16_lineage$Fine_annotation)
)
plot_df <- plot_df[plot_df$has_clone_id, ]
plot_df$fine_annotation <- factor(plot_df$fine_annotation, levels = row_order)
ggplot(plot_df, aes(fine_annotation)) +
  geom_bar(stat = "count", position = "dodge2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~stage)


## write data-frames for z-score analysis:
e12_lin_subset_df <- data.frame(
  "cellID" = colnames(E12_lineage),
  "cloneID" = E12_lineage$cloneID.umi6,
  "ident" = E12_lineage$Fine_annotation
)
rownames(e12_lin_subset_df) <- 1:nrow(e12_lin_subset_df)
e12_lin_subset_df <- e12_lin_subset_df[!is.na(e12_lin_subset_df$cloneID), ]
clone_dist_e12 <- table(e12_lin_subset_df$cloneID)
e12_lin_subset_df <- e12_lin_subset_df[clone_dist_e12[e12_lin_subset_df$cloneID] > 1, ]
write.table(e12_lin_subset_df, "e12_lineage_table_fine_annot.csv", sep = ",", row.names = F)


e16_lin_subset_df <- data.frame(
  "cellID" = colnames(E16_lineage),
  "cloneID" = E16_lineage$cloneID.LINEAGE,
  "ident" = E16_lineage$Fine_annotation
)
rownames(e16_lin_subset_df) <- 1:nrow(e16_lin_subset_df)
e16_lin_subset_df <- e16_lin_subset_df[!is.na(e16_lin_subset_df$cloneID), ]
clone_dist_e16 <- table(e16_lin_subset_df$cloneID)
e16_lin_subset_df <- e16_lin_subset_df[clone_dist_e16[e16_lin_subset_df$cloneID] > 1, ]
write.table(e16_lin_subset_df, "e16_lineage_table_fine_annot.csv", sep = ",", row.names = F)


## run lineage coupling analysis using standard method and updated scripts (incorporating p-value calculation) ##
## code repo for this steps can be found under: https://github.com/mayer-lab/Bandler-et-al_lineage
#system("Rscript lineage_coupling/1_generate_input_for_lineage_coupling_analysis.R")
#system("bash lineage_coupling/2.2_lineage_coupling_analysis.sh")

## visualize results ##
row_order <- c("Fabp7","Fabp7_Ccnd2","Top2a","Ube2c","Abracl","Nkx2_1","Npy","Maf_Sst","Snhg11_Lhx8","Snhg11","Tcf4_Nr2f2","Tshz1","Six3_Gucy1a3","Gucy1a3","Ebf1_Isl1")

e12_mtx <- read.table("lineage_coupling_out/GE_e12/lineage_coupling_scores_matrix.csv", h=T, sep=",",row.names = 1)
e12_mtx <- e12_mtx[row_order, row_order]

e16_mtx <- read.table("lineage_coupling_out/GE_e16/lineage_coupling_scores_matrix.csv", h=T, sep=",",row.names = 1)
e16_mtx <- e16_mtx[row_order, row_order]

pal_length <- 50
color_vec <- colorRampPalette(c("blue", "grey", "red"))(pal_length)
min_value <- min(c(min(e12_mtx), min(e16_mtx))); max_value <- max(c(max(e12_mtx), max(e16_mtx)))
breaks_vec <- c(seq(min_value, 0, length.out=ceiling(pal_length/2) + 1), 
              seq(max_value/pal_length,max_value, length.out=floor(pal_length/2)))


p11 <- pheatmap(e12_mtx, cluster_rows = F, cluster_cols = F, display_numbers = T, na_col="white", color = color_vec, breaks = breaks_vec)[[4]]
p12 <- pheatmap(e16_mtx, cluster_rows = F, cluster_cols = F, display_numbers = T, na_col="white", color = color_vec, breaks = breaks_vec)[[4]]

grid.arrange(p11,p12, ncol = 2, left = "e12", right = "e16")

e12_pvalue_mtx <- read.table("lineage_coupling_out/GE_e12/lineage_coupling_p_values_matrix.csv", h=T, sep=",",row.names = 1)
e12_pvalue_mtx <- e12_pvalue_mtx[row_order, row_order]

e16_pvalue_mtx <- read.table("lineage_coupling_out/GE_e16/lineage_coupling_p_values_matrix.csv", h=T, sep=",",row.names = 1)
e16_pvalue_mtx <- e16_pvalue_mtx[row_order, row_order]

## substitute p-values with symbols:
e12_pvalue_mtx[e12_pvalue_mtx < 0.001] <- "***"
e12_pvalue_mtx[e12_pvalue_mtx < 0.01 & e12_pvalue_mtx >= 0.001] <- "**"
e12_pvalue_mtx[e12_pvalue_mtx < 0.05 & e12_pvalue_mtx >= 0.01] <- "*"
e12_pvalue_mtx[e12_pvalue_mtx >= 0.05] <- ""

e16_pvalue_mtx[e16_pvalue_mtx < 0.001] <- "***"
e16_pvalue_mtx[e16_pvalue_mtx < 0.01 & e16_pvalue_mtx >= 0.001] <- "**"
e16_pvalue_mtx[e16_pvalue_mtx < 0.05 & e16_pvalue_mtx >= 0.01] <- "*"
e16_pvalue_mtx[e16_pvalue_mtx >= 0.05] <- ""

## only show upper triangle:
e12_mtx_triangle <- e12_mtx;  e12_mtx_triangle[lower.tri(e12_mtx_triangle)] <- NA
e16_mtx_triangle <- e16_mtx;  e16_mtx_triangle[lower.tri(e16_mtx_triangle)] <- NA
e12_pvalue_mtx[lower.tri(e12_pvalue_mtx)] <- ""
e16_pvalue_mtx[lower.tri(e16_pvalue_mtx)] <- ""

p21 <- pheatmap(e12_mtx_triangle, cluster_rows = F, cluster_cols = F, 
                display_numbers = e12_pvalue_mtx, na_col = "white", border_color ="white", color = color_vec, breaks = breaks_vec)[[4]]
p22 <- pheatmap(e16_mtx_triangle, cluster_rows = F, cluster_cols = F, 
                display_numbers = e16_pvalue_mtx, na_col = "white", border_color ="white", color = color_vec, breaks = breaks_vec)[[4]]

grid.arrange(p21,p22, ncol = 2, left = "e12", right = "e16")


## cluster rows and cols:
p31 <- pheatmap(e12_mtx, cluster_rows = T, cluster_cols = T, display_numbers = e12_pvalue_mtx, color = color_vec, breaks = breaks_vec)[[4]]
p32 <- pheatmap(e16_mtx, cluster_rows = T, cluster_cols = T, display_numbers = e16_pvalue_mtx, color = color_vec, breaks = breaks_vec)[[4]]

grid.arrange(p31,p32, ncol = 2, left = "e12", right = "e16")

