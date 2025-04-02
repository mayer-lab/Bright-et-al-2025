## get subnetworks from scenicplus network ##

library(ggrepel)
library(igraph)

source(file = "Scripts/lib.R")

## load data:
eRegulon_md <- read.table("Processed_Objects/eRegulon_metadata_filtered.tsv", h = T, sep = "\t")
eRegulon_AUC <- read.table("Processed_Objects/cfse_network_wArchR_peaks_eRegulon_AUC_filtered_gene_based.csv", h=T, sep = ",", row.names = 1)
eRegulon_AUC_thresholds <- read.table("Processed_Objects/cfse_network_wArchR_peaks_eRegulon_AUC_filtered_gene_based_thresholds.csv", h=T, sep = ",", row.names = 1)
cell_metadata <- read.table("Processed_Objects/cfse_network_wArchR_peaks_cell_metadata.csv", h=T, sep = ",", row.names = 1)
cfse_seurat <- readRDS("Processed_Objects/CFSE_sub.rds")


## get subnetwork per stage and cell state
stage_broad_vec <- c("AP_e12", "AP_e16","BP_e12","BP_e16","precursor_e12","precursor_e16")
eRegulon_list <- lapply(stage_broad_vec, get_sub_network)
names(eRegulon_list) <- stage_broad_vec


## get complete networks:
AP_e12_dg <- get_complete_network(eRegulon_list$AP_e12, sba = "AP_e12")
AP_e16_dg <- get_complete_network(eRegulon_list$AP_e16, sba = "AP_e16")
BP_e12_dg <- get_complete_network(eRegulon_list$BP_e12, sba = "BP_e12")
BP_e16_dg <- get_complete_network(eRegulon_list$BP_e16, sba = "BP_e16")
precursor_e12_dg <- get_complete_network(eRegulon_list$precursor_e12, sba = "precursor_e12")
precursor_e16_dg <- get_complete_network(eRegulon_list$precursor_e16, sba = "precursor_e16")

cell_state_dgs <- list("AP_e12" = AP_e12_dg, "AP_e16" = AP_e16_dg, "BP_e12" = BP_e12_dg, "BP_e16" = BP_e16_dg,
                       "precursor_e12" = precursor_e12_dg, "precursor_e16" = precursor_e16_dg)

## network with only TFs as nodes ##
mm10_tfs <- read.table("Processed_Objects/mm_mgi_tfs.txt"); mm10_tfs <- mm10_tfs$V1

AP_e12_tf_dg <- get_complete_network(eRegulon_list$AP_e12, sba = "AP_e12", only_TFs = TRUE)
AP_e16_tf_dg <- get_complete_network(eRegulon_list$AP_e16, sba = "AP_e16", only_TFs = TRUE)
BP_e12_tf_dg <- get_complete_network(eRegulon_list$BP_e12, sba = "BP_e12", only_TFs = TRUE)
BP_e16_tf_dg <- get_complete_network(eRegulon_list$BP_e16, sba = "BP_e16", only_TFs = TRUE)
precursor_e12_tf_dg <- get_complete_network(eRegulon_list$precursor_e12, sba = "precursor_e12", only_TFs = TRUE)
precursor_e16_tf_dg <- get_complete_network(eRegulon_list$precursor_e16, sba = "precursor_e16", only_TFs = TRUE)

cell_state_tf_tf_dgs <- list("AP_e12" = AP_e12_tf_dg, "AP_e16" = AP_e16_tf_dg, "BP_e12" = BP_e12_tf_dg, "BP_e16" = BP_e16_tf_dg,
                       "precursor_e12" = precursor_e12_tf_dg, "precursor_e16" = precursor_e16_tf_dg)


## networks across stages per celltype:
## AP:
AP_merged_tf_dg <- combine_networks(e12_dg = AP_e12_tf_dg, e16_dg = AP_e16_tf_dg,
                                    celltype = "AP", differential_edge_colors = FALSE)
##additional plot aesthetics:
E(AP_merged_tf_dg)$width <- 0.5
alpha_vec <- E(AP_merged_tf_dg)$TF2G_rho*255
alpha_vec[alpha_vec < 100] <- 100
E(AP_merged_tf_dg)$color <- rgb(105,105,105, alpha = alpha_vec, maxColorValue = 255)
V(AP_merged_tf_dg)$size <- (igraph::degree(AP_merged_tf_dg, mode = "out") + 5)*0.3
E(AP_merged_tf_dg)$arrow.size <- 0.4

AP_merged_tf_df <- as_data_frame(AP_merged_tf_dg, what = "both")
write.table(AP_merged_tf_df$vertices, "Results/source_data/fig3a_nodes.csv", sep = ",", quote = T, row.names = F)
write.table(AP_merged_tf_df$edges, "Results/source_data/fig3a_edges.csv", sep = ",", quote = T, row.names = F)

set.seed(1)
plot(AP_merged_tf_dg, edge.lengths = 20, edge.color="grey",
     vertex.label.color = "black", vertex.frame.color="grey", edge.arrow.size = 0.3)


## merged BP & precursor network:
BP_merged_tf_dg <- combine_networks(e12_dg = BP_e12_tf_dg, e16_dg = BP_e16_tf_dg,
                                    celltype = "BP", differential_edge_colors = FALSE)
##additional plot aesthetics:
V(BP_merged_tf_dg)$size <- (igraph::degree(BP_merged_tf_dg, mode = "out") + 5)*0.6
E(BP_merged_tf_dg)$width <- 0.5
E(BP_merged_tf_dg)$arrow.size <- 0.5
alpha_vec <- E(BP_merged_tf_dg)$TF2G_rho*255
alpha_vec[alpha_vec < 100] <- 100
E(BP_merged_tf_dg)$color <- rgb(105,105,105, alpha = alpha_vec, maxColorValue = 255)

BP_merged_tf_df <- as_data_frame(BP_merged_tf_dg, what = "both")
write.table(BP_merged_tf_df$vertices, "Results/source_data/edf5b_nodes.csv", sep = ",", quote = T, row.names = F)
write.table(BP_merged_tf_df$edges, "Results/source_data/edf5b_edges.csv", sep = ",", quote = T, row.names = F)

set.seed(1)
plot(BP_merged_tf_dg, edge.lengths = 20,
     vertex.label.color = "black", vertex.frame.color="grey", edge.arrow.size = 0.3, edge.color="grey")

##
precursor_merged_tf_dg <- combine_networks(e12_dg = precursor_e12_tf_dg, e16_dg = precursor_e16_tf_dg,
                                    celltype = "precursor", differential_edge_colors = FALSE)
##additional plot aesthetics:
V(precursor_merged_tf_dg)$size <- (igraph::degree(precursor_merged_tf_dg, mode = "out") + 5)*0.6
E(precursor_merged_tf_dg)$width <- 0.5
E(precursor_merged_tf_dg)$arrow.size <- 0.5
alpha_vec <- E(precursor_merged_tf_dg)$TF2G_rho*255
alpha_vec[alpha_vec < 100] <- 100
E(precursor_merged_tf_dg)$color <- rgb(105,105,105, alpha = alpha_vec, maxColorValue = 255)

precursor_merged_tf_df <- as_data_frame(precursor_merged_tf_dg, what = "both")
write.table(precursor_merged_tf_df$vertices, "Results/source_data/edf5c_nodes.csv", sep = ",", quote = T, row.names = F)
write.table(precursor_merged_tf_df$edges, "Results/source_data/edf5c_edges.csv", sep = ",", quote = T, row.names = F)

set.seed(66)
plot(precursor_merged_tf_dg, edge.lengths = 20,
     vertex.label.color = "black", vertex.frame.color="grey", edge.arrow.size = 0.3, edge.color="grey")


## STAGE SPECIFIC SUBNETWORKS

maturation_TF_list <- readRDS("Processed_Objects/maturation_n2_TF_list.rds")
maturation_gene_list <- readRDS("Processed_Objects/maturation_n2_gene_list.rds")

## combine cell state networks per stage:

## e12 TF-TF
e12_regulon_df <- rbind(eRegulon_list$AP_e12, eRegulon_list$BP_e12, eRegulon_list$precursor_e12)
e12_complete_tf_tf_gene_dgs <- get_complete_network_per_stage(e12_regulon_df, stage = "E12+6hCFSE", only_TFs = TRUE)
e12_ventral_maturation_l1_tf_tf_dgs <- subset_network(e12_complete_tf_tf_gene_dgs, foi_vec = maturation_TF_list$inhibitory, lev = 1, just_upstream = TRUE)
E(e12_ventral_maturation_l1_tf_tf_dgs)$width <- E(e12_ventral_maturation_l1_tf_tf_dgs)$width / 2

e12_ventral_maturation_l1_tf_tf_df <- as_data_frame(e12_ventral_maturation_l1_tf_tf_dgs, what = "both")
write.table(e12_ventral_maturation_l1_tf_tf_df$vertices, "Results/source_data/edf5f_nodes.csv", sep = ",", quote = T, row.names = F)
write.table(e12_ventral_maturation_l1_tf_tf_df$edges, "Results/source_data/edf5f_edges.csv", sep = ",", quote = T, row.names = F)

set.seed(1)
plot(e12_ventral_maturation_l1_tf_tf_dgs, edge.arrow.size = .4, edge.lengths = 20, vertex.size = 7,
     vertex.frame.color = "grey", edge.color = "grey", vertex.label.color = "black")

# same for e16
e16_regulon_df <- rbind(eRegulon_list$AP_e16, eRegulon_list$BP_e16, eRegulon_list$precursor_e16)
e16_complete_tf_tf_gene_dgs <- get_complete_network_per_stage(e16_regulon_df, stage = "E16+6hCFSE", only_TFs = TRUE)
e16_ventral_maturation_l1_tf_tf_dgs <- subset_network(e16_complete_tf_tf_gene_dgs, foi_vec = maturation_TF_list$inhibitory, lev = 1, just_upstream = TRUE)
E(e16_ventral_maturation_l1_tf_tf_dgs)$width <- E(e16_ventral_maturation_l1_tf_tf_dgs)$width / 2

e16_ventral_maturation_l1_tf_tf_df <- as_data_frame(e16_ventral_maturation_l1_tf_tf_dgs, what = "both")
write.table(e16_ventral_maturation_l1_tf_tf_df$vertices, "Results/source_data/edf5g_nodes.csv", sep = ",", quote = T, row.names = F)
write.table(e16_ventral_maturation_l1_tf_tf_df$edges, "Results/source_data/edf5g_edges.csv", sep = ",", quote = T, row.names = F)

set.seed(1)
plot(e16_ventral_maturation_l1_tf_tf_dgs, edge.arrow.size = .4, edge.lengths = 20, vertex.size = 7,
     vertex.frame.color = "grey", edge.color = "grey", vertex.label.color = "black")

## quantify binding events:
e12_complete_tf_gene_dgs <- get_complete_network_per_stage(e12_regulon_df, stage = "E12+6hCFSE", only_TFs = FALSE)
e12_ventral_maturation_l1_tf_gene_dgs <- subset_network(e12_complete_tf_gene_dgs, foi_vec = maturation_gene_list$inhibitory, lev = 1, just_upstream = TRUE)

e16_complete_tf_gene_dgs <- get_complete_network_per_stage(e16_regulon_df, stage = "E16+6hCFSE", only_TFs = FALSE)
e16_ventral_maturation_l1_tf_gene_dgs <- subset_network(e16_complete_tf_gene_dgs, foi_vec = maturation_gene_list$inhibitory, lev = 1, just_upstream = TRUE)

e12_ventral_diff_degrees <- degree(e12_ventral_maturation_l1_tf_gene_dgs, mode = "out")
e16_ventral_diff_degrees <- degree(e16_ventral_maturation_l1_tf_gene_dgs, mode = "out")

e12_ventral_diff_degrees_df <- data.frame("gene" = names(e12_ventral_diff_degrees), "stage" = "e12", "degree" = as.numeric(e12_ventral_diff_degrees))
e12_ventral_diff_degrees_df <- e12_ventral_diff_degrees_df[order(e12_ventral_diff_degrees_df$degree, decreasing = F), ]
e12_ventral_diff_degrees_df$idx <- 1:nrow(e12_ventral_diff_degrees_df)
e16_ventral_diff_degrees_df <- data.frame("gene" = names(e16_ventral_diff_degrees), "stage" = "e16", "degree" = as.numeric(e16_ventral_diff_degrees))
e16_ventral_diff_degrees_df <- e16_ventral_diff_degrees_df[order(e16_ventral_diff_degrees_df$degree, decreasing = F), ]
e16_ventral_diff_degrees_df$idx <- 1:nrow(e16_ventral_diff_degrees_df)

diff_TFs_degree_df <- rbind(e12_ventral_diff_degrees_df, e16_ventral_diff_degrees_df)
diff_TFs_degree_df$gene[diff_TFs_degree_df$degree < 8] <- NA

write.table(diff_TFs_degree_df, "Results/source_data/edf5h.csv", sep = ",", quote = F, row.names = F)

ggplot(diff_TFs_degree_df, aes(idx, degree, color = stage, label = gene)) +
  geom_point() +
  facet_wrap(~stage, ncol = 1) +
  geom_text_repel(max.overlaps = 25, size = 8) +
  scale_color_manual(values = c("#E16B8F", "#5776A7")) +
  ylab("Number of maturation genes regulated by TF") + xlab("Index") +
  theme(legend.position = 'none', axis.text=element_text(size=14),
        axis.title=element_text(size=14), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=14))



## plot Nfib-Tcf4-Meis2 neighbor tf-gene network in e16:
e16_regulon_sub_df <- e16_regulon_df[e16_regulon_df$TF %in% c("Nfib","Tcf4","Meis2"), ]
FOI_e16_tf_gene_dg <- get_complete_network_per_stage(e16_regulon_sub_df, stage = "E16+6hCFSE", only_TFs = F)

V(FOI_e16_tf_gene_dg)$label.cex <- 1
V(FOI_e16_tf_gene_dg)$label.cex[V(FOI_e16_tf_gene_dg)$name %in% c(c("Nfib","Tcf4","Meis2"))] <- 2

expr_vec <- V(FOI_e16_tf_gene_dg)$avg_expression
expr_vec <- round(expr_vec, digits = 1)
expr_vec_levels <- unique(expr_vec)[order(unique(expr_vec))]
colfunc <- colorRampPalette(c("white", "#5776A7"))
col_vec <- colfunc(length(expr_vec_levels))
names(col_vec) <- expr_vec_levels
V(FOI_e16_tf_gene_dg)$color <- col_vec[as.character(expr_vec)]
edge_df <- as_long_data_frame(FOI_e16_tf_gene_dg)
E(FOI_e16_tf_gene_dg)$color <- "grey"
E(FOI_e16_tf_gene_dg)$color[edge_df$from_name %in% c("Nfib","Tcf4","Meis2") &
                              edge_df$to_name %in% c("Nfib","Tcf4","Meis2")] <- "blue"
E(FOI_e16_tf_gene_dg)$width <- E(FOI_e16_tf_gene_dg)$width / 2

FOI_e16_tf_gene_df <- as_data_frame(FOI_e16_tf_gene_dg, what = "both")
write.table(FOI_e16_tf_gene_df$vertices, "Results/source_data/fig3b_nodes.csv", sep = ",", quote = T, row.names = F)
write.table(FOI_e16_tf_gene_df$edges, "Results/source_data/fig3b_edges.csv", sep = ",", quote = T, row.names = F)

set.seed(6)
plot(FOI_e16_tf_gene_dg, vertex.size = 5, edge.lengths = 8, edge.arrow.size = 0.5,
     vertex.label.color = "black", vertex.frame.color = "grey",
     label.font = "Helvetica", vertex.frame.width = 0.7)

