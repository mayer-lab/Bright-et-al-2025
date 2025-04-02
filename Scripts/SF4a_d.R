## genes linked to maturation (maturation genes) ##

source(file = "Scripts/lib.R")

## load data:
EXCIT_INHIBIT_cleaned_sub <- readRDS(file ="Processed_Objects/EXCIT_INHIBIT_cleaned_sub.rds")
traj_list_merged <- readRDS(file = "Processed_Objects/traj_list_cellIDs_merged2.rds")

## create subsets:
inhibitory_clusters <- c("Gas1_Ldha","Hes1_Fabp7","Fabp7_Mt3","Hist1h1b_Top2a","Ccnd2_Nudt4","Nkx2-1_Lhx8","Maf_Npy", "Sst_Maf","Nr2f2_Nr2f1","Foxp1_Gucy1a3", "Ebf1_Foxp1","Isl1_Zfp503")
excitatory_clusters <- c("Gas1_Ldha","Hes1_Fabp7","Fabp7_Mt3","Neurog2_Rrm2","Neurog2_Eomes","Neurod2_Neurod6","Neurod6_Mef2c")

## subset seurat for type and stage:
ventral_sub <- subset(EXCIT_INHIBIT_cleaned_sub, subset = General2 == "ventral")
ventral_sub <- add_smooth_expression(ventral_sub)
dorsal_sub <- subset(EXCIT_INHIBIT_cleaned_sub, subset = General2 == "dorsal")
dorsal_sub <- add_smooth_expression(dorsal_sub)


e12_ventral_sub <- subset(ventral_sub, subset = Dataset_DV %in% c("INHIBITORY_E12_1","INHIBITORY_E12_2"))
e12_ventral_sub <- add_smooth_expression(e12_ventral_sub)

e13_ventral_sub <- subset(ventral_sub, subset = Dataset_DV %in% c("INHIBITORY_E13"))
e13_ventral_sub <- add_smooth_expression(e13_ventral_sub)

e14_ventral_sub <- subset(ventral_sub, subset = Dataset_DV %in% c("INHIBITORY_E14"))
e14_ventral_sub <- add_smooth_expression(e14_ventral_sub)

e15_ventral_sub <- subset(ventral_sub, subset = Dataset_DV %in% c("INHIBITORY_E15"))
e15_ventral_sub <- add_smooth_expression(e15_ventral_sub)

e16_ventral_sub <- subset(ventral_sub, subset = Dataset_DV %in% c("INHIBITORY_E16_1", "INHIBITORY_E16_2", "INHIBITORY_E16_3"))
e16_ventral_sub <- add_smooth_expression(e16_ventral_sub)


e12_dorsal_sub <- subset(dorsal_sub, subset = Dataset_DV %in% c("EXCITATORY_E12"))
e12_dorsal_sub <- add_smooth_expression(e12_dorsal_sub)

e13_dorsal_sub <- subset(dorsal_sub, subset = Dataset_DV %in% c("EXCITATORY_E13"))
e13_dorsal_sub <- add_smooth_expression(e13_dorsal_sub)

e14_dorsal_sub <- subset(dorsal_sub, subset = Dataset_DV %in% c("EXCITATORY_E14"))
e14_dorsal_sub <- add_smooth_expression(e14_dorsal_sub)

e15_dorsal_sub <- subset(dorsal_sub, subset = Dataset_DV %in% c("EXCITATORY_E15"))
e15_dorsal_sub <- add_smooth_expression(e15_dorsal_sub)

e16_dorsal_sub <- subset(dorsal_sub, subset = Dataset_DV %in% c("EXCITATORY_E16"))
e16_dorsal_sub <- add_smooth_expression(e16_dorsal_sub)


## get maturation genes:
e12_dorsal_maturation_genes <- get_maturation_genes(e12_dorsal_sub, nbins = 10)
e13_dorsal_maturation_genes <- get_maturation_genes(e13_dorsal_sub, nbins = 10)
e14_dorsal_maturation_genes <- get_maturation_genes(e14_dorsal_sub, nbins = 10)
e15_dorsal_maturation_genes <- get_maturation_genes(e15_dorsal_sub, nbins = 10)
e16_dorsal_maturation_genes <- get_maturation_genes(e16_dorsal_sub, nbins = 10)

excitatory_maturation_genes <- e12_dorsal_maturation_genes[(e12_dorsal_maturation_genes %in% e13_dorsal_maturation_genes) &
                                                              (e12_dorsal_maturation_genes %in% e14_dorsal_maturation_genes) &
                                                              (e12_dorsal_maturation_genes %in% e15_dorsal_maturation_genes) &
                                                              (e12_dorsal_maturation_genes %in% e16_dorsal_maturation_genes)]

## interneuron and PN maturation genes (num_trajectories = 2) ##
e12_n2_ventral_maturation_genes <- get_maturation_genes(e12_ventral_sub, traj_list = traj_list_merged, nbins = 10, num_trajectories = 2)
e13_n2_ventral_maturation_genes <- get_maturation_genes(e13_ventral_sub, traj_list = traj_list_merged, nbins = 10, num_trajectories = 2)
e14_n2_ventral_maturation_genes <- get_maturation_genes(e14_ventral_sub, traj_list = traj_list_merged, nbins = 10, num_trajectories = 2)
e15_n2_ventral_maturation_genes <- get_maturation_genes(e15_ventral_sub, traj_list = traj_list_merged, nbins = 10, num_trajectories = 2)
e16_n2_ventral_maturation_genes <- get_maturation_genes(e16_ventral_sub, traj_list = traj_list_merged, nbins = 10, num_trajectories = 2)

inhibitory_n2_maturation_genes <- e12_n2_ventral_maturation_genes[(e12_n2_ventral_maturation_genes %in% e13_n2_ventral_maturation_genes) &
                                                              (e12_n2_ventral_maturation_genes %in% e14_n2_ventral_maturation_genes) &
                                                              (e12_n2_ventral_maturation_genes %in% e15_n2_ventral_maturation_genes) &
                                                              (e12_n2_ventral_maturation_genes %in% e16_n2_ventral_maturation_genes)]


maturation_n2_gene_list <- list(
  "ventral_e12_n2" = e12_n2_ventral_maturation_genes,
  "ventral_e13_n2" = e13_n2_ventral_maturation_genes,
  "ventral_e14_n2" = e14_n2_ventral_maturation_genes,
  "ventral_e15_n2" = e15_n2_ventral_maturation_genes,
  "ventral_e16_n2" = e16_n2_ventral_maturation_genes,
  "dorsal_e12" = e12_dorsal_maturation_genes,
  "dorsal_e13" = e13_dorsal_maturation_genes,
  "dorsal_e14" = e14_dorsal_maturation_genes,
  "dorsal_e15" = e15_dorsal_maturation_genes,
  "dorsal_e16" = e16_dorsal_maturation_genes,
  "inhibitory" = inhibitory_n2_maturation_genes,
  "excitatory" = excitatory_maturation_genes
)

## subset TFs:
mm10_tfs <- read.table("Processed_Objects/mm_mgi_tfs.txt"); mm10_tfs <- mm10_tfs$V1
mm10_tfs <- mm10_tfs[mm10_tfs %in% VariableFeatures(EXCIT_INHIBIT_cleaned_sub)]

maturation_tf_list <- lapply(maturation_gene_list, function(maturation_genes) {
  return(maturation_genes[maturation_genes %in% mm10_tfs])
})

## plotting ##
## plot inhibitory maturation genes:
plot_mtx_with_ggplot_v2(maturation_gene_list$inhibitory, s = ventral_sub, mpm_thres = -1, legend = F, ysize = 5, sort = T, annotate_only_TFs = F)

plot_mtx_with_ggplot_v2(maturation_gene_list$excitatory, s = dorsal_sub, mpm_thres = -1, legend = F, ysize = 7, sort = T)

## overlap between exc/ inh genes:
ggvenn::ggvenn(list("inhibitory" = maturation_gene_list$inhibitory, "excitatory" = maturation_gene_list$excitatory), auto_scale = T, 
               fill_color = c("#24796C","#DAA51B"))

## plot expression of overlapping genes:
ei_overlap_genes <- maturation_gene_list$inhibitory[maturation_gene_list$inhibitory %in% maturation_gene_list$excitatory]
plot_mtx_with_ggplot_v2(ei_overlap_genes, s = ventral_sub, mpm_thres = -1, legend = F, ysize = 13, sort = F, annotate_only_TFs = F)
plot_mtx_with_ggplot_v2(ei_overlap_genes, s = dorsal_sub, mpm_thres = -1, legend = F, ysize = 13, sort = F, annotate_only_TFs = F)
