## Co-Binding according to scenicplus network ##

library(ggplot2)
library(Seurat)
library(UpSetR)
library(igraph)
library(rcartocolor)
library(ggvenn)

## load data:
eRegulon_md <- read.table("Processed_Objects/eRegulon_metadata_filtered.tsv", h = T, sep = "\t")
cfse_seurat <- readRDS("Processed_Objects/CFSE_sub.rds")
cell_state_dgs <- readRDS(file = "Processed_Objects/cell_states_tf_gene_subnetworks.rds")
eRegulon_list <- readRDS("Processed_Objects/cfse_network_eRegulon_stage_broad_list.rds")

## region aware co-binding over all states ##

## get all TFs binding to regions:
regulon_enh_list <- list()
for(i in 1:nrow(eRegulon_md)) {
  if(eRegulon_md$Region[i] %in% names(regulon_enh_list)) {
    regulon_enh_list[[eRegulon_md$Region[i]]] <- c(regulon_enh_list[[eRegulon_md$Region[i]]], eRegulon_md$TF[i])
  } else {
    regulon_enh_list[[eRegulon_md$Region[i]]] <- eRegulon_md$TF[i]
  }
}

## check out target gene overlap per cell state ##
nfib_target_list <- lapply(eRegulon_list, function(df) {unique(df$Gene[df$TF == "Nfib"])})
nfix_target_list <- lapply(eRegulon_list, function(df) {unique(df$Gene[df$TF == "Nfix"])})
tcf4_target_list <- lapply(eRegulon_list, function(df) {unique(df$Gene[df$TF == "Tcf4"])})
meis2_target_list <- lapply(eRegulon_list, function(df) {unique(df$Gene[df$TF == "Meis2"])})

plot_list <- list()
for(cell_state in names(eRegulon_list)) {
  plot_list[[cell_state]] <- ggvenn(list("Tcf4"=tcf4_target_list[[cell_state]], "Nfib"=nfib_target_list[[cell_state]], "Nfix"=nfix_target_list[[cell_state]], "Meis2"=meis2_target_list[[cell_state]]), 
                                    auto_scale = F, fill_color = carto_pal(10, "Prism")[c(1,3,6,8)])
}
patchwork::wrap_plots(plot_list, ncol = 2, byrow = T)


## check out common target genes of Nfib, Tcf4 and Meis2 at e16:
nfib_e16_targets <- unique(c(nfib_target_list$AP_e16, nfib_target_list$BP_e16, nfib_target_list$precursor_e16))
tcf4_e16_targets <- unique(c(tcf4_target_list$AP_e16, tcf4_target_list$BP_e16, tcf4_target_list$precursor_e16))
meis2_e16_targets <- unique(c(meis2_target_list$AP_e16, meis2_target_list$BP_e16, meis2_target_list$precursor_e16))

## shared genes as in bound by Nfib & Tcf4 OR Nfib & Meis2:
shared_genes3 <- unique(nfib_e16_targets[nfib_e16_targets %in% tcf4_e16_targets | nfib_e16_targets %in% meis2_e16_targets])
paste(shared_genes3, collapse = ", ")

go_res3 <- read.table("nfib_tcf4_meis2_e16_shared_targets_BP_GO_enrich.txt.txt", sep = "\t", h=T)
go_res3$Term <- factor(go_res3$Term, levels = go_res3$Term[order(go_res3$Bonferroni, decreasing = T)])
ggplot(go_res3[go_res3$Benjamini < 0.05, ], aes(x = Term, y = -log10(Benjamini))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(legend.position = 'none', axis.text=element_text(size=6),
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=6))
