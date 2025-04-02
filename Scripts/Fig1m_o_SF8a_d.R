#Upset plots ----
source("Scripts/lib.R")
library(UpSetR)
library(tidyr)

#Prep

load("Processed_Objects/Seurat_objects_All.Rdata")

##
exp2_vec <- c("FlashTag" = "CFSE", "Reference" = "WT", "TrackerSeq" = "LINEAGE")
tip_vec <- c("Fabp7"="No","Fabp7_Ccnd2"="No","Snhg11"="Snhg11","Npy"="No","Tshz1"="No", "Ube2c"="No", "Ebf1_Isl1"="Ebf1_Isl1","Nkx2_1"="No",
             "Top2a"="No","Tcf4_Nr2f2"="Nr2f2", "Abracl"="No", "Maf_Sst"="Maf",  "Six3_Gucy1a3"="No", "Gucy1a3"="Gucy1a3")
mt_vec <- c("Fabp7"="Mitotic","Fabp7_Ccnd2"="Mitotic","Snhg11"="Snhg11","Npy"="No","Tshz1"="No", "Ube2c"="Mitotic", "Ebf1_Isl1"="Ebf1_Isl1","Nkx2_1"="No",
            "Top2a"="Mitotic","Tcf4_Nr2f2"="Nr2f2", "Abracl"="No", "Maf_Sst"="Maf",  "Six3_Gucy1a3"="No", "Gucy1a3"="Gucy1a3")


Seurat_objects_TempGO$Fine_annotation <- Seurat_objects_TempGO$Cluster_annotation
Seurat_objects_TempGO$Experiment2 <- exp2_vec[as.character(Seurat_objects_TempGO$Dataset_type)]
Seurat_objects_TempGO$Tips <- tip_vec[as.character(Seurat_objects_TempGO$Cluster_annotation)]
Seurat_objects_TempGO$Mitotic_and_tips <- mt_vec[as.character(Seurat_objects_TempGO$Cluster_annotation)]

Seurat_objects_TempGO$Experiment <- Seurat_objects_TempGO$Experiment2
Seurat_objects_TempGO$Experiment[Seurat_objects_TempGO$Experiment == "LINEAGE" & Seurat_objects_TempGO$Injection_stage == "E12"] <- "E12+4daysLINEAGE"
Seurat_objects_TempGO$Experiment[Seurat_objects_TempGO$Experiment == "LINEAGE" & Seurat_objects_TempGO$Injection_stage == "E16"] <- "E16+4daysLINEAGE"
Seurat_objects_TempGO$Experiment[Seurat_objects_TempGO$Experiment == "CFSE" & Seurat_objects_TempGO$Injection_stage == "E16"] <- "E16+6hCFSE"
Seurat_objects_TempGO$Experiment[Seurat_objects_TempGO$Experiment == "CFSE" & Seurat_objects_TempGO$Injection_stage == "E12" & Seurat_objects_TempGO$Collection_stage == "E12"] <- "E12+6hCFSE"
Seurat_objects_TempGO$Experiment[Seurat_objects_TempGO$Experiment == "CFSE" & Seurat_objects_TempGO$Injection_stage == "E12" & Seurat_objects_TempGO$Collection_stage == "E16"] <- "E12+4daysCFSE"
table(Seurat_objects_TempGO$Experiment)

Inhibitory_datasets <- Seurat_objects_TempGO


# subset:
Inhibitory_datasets <- SetIdent(Inhibitory_datasets, value = "Experiment2")
Lineage_datasets <- subset(Inhibitory_datasets, ident = c("LINEAGE"))
Lineage_datasets <- SetIdent(Lineage_datasets, value = "Tips")
Tips <- subset(Lineage_datasets, ident = c("No"), invert = T)
Tips <- SetIdent(Tips, value = "Experiment")
E12_lineage <- subset(Tips, ident = "E12+4daysLINEAGE")
E16_lineage <- subset(Tips, ident = "E16+4daysLINEAGE")
E12_lineage <- SetIdent(E12_lineage, value = "Tips")
E16_lineage <- SetIdent(E16_lineage, value = "Tips")

# UpsetPlot E12 
pool_upset.E12 <- FetchData(E12_lineage,c("cloneID.umi6","Dataset","Tips"))
pool_upset.E12 <- pool_upset.E12 %>% drop_na(cloneID.umi6)
pool_upset.E12 <-  pool_upset.E12 %>% ungroup() %>% group_by(cloneID.umi6) %>% filter(n() != 1)# Take everything that is not single clone

write.table(fromList(split(pool_upset.E12$cloneID.umi6, pool_upset.E12$Tips)), "Results/source_data/fig1m.csv", sep = ",", quote = F, row.names = F)

upset(fromList(split(pool_upset.E12$cloneID.umi6, pool_upset.E12$Tips)),
      nintersects =40,
      nsets = 10,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.3, 0.5)
)

# UpsetPlot E16
E16_lineage$cloneID.LINEAGE <- E16_lineage$Clones_e16
pool_upset.E16<- FetchData(E16_lineage, c("cloneID.LINEAGE","Dataset","Tips"))
pool_upset.E16<- pool_upset.E16%>% drop_na(cloneID.LINEAGE)# Drop rows containing missing values
pool_upset.E16<-  pool_upset.E16%>% ungroup() %>% group_by(cloneID.LINEAGE) %>% filter(n() != 1)# Take everything that is not single clone

write.table(fromList(split(pool_upset.E16$cloneID.LINEAGE, pool_upset.E16$Tips)), "Results/source_data/fig1n.csv", sep = ",", quote = F, row.names = F)
upset(fromList(split(pool_upset.E16$cloneID.LINEAGE, pool_upset.E16$Tips)),
      nintersects = 40,
      nsets = 20,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.3, 0.7)
)


# Correlation analysis (example: Ebf1_Isl1 cell state) ----

#Prep 
Inhibitory_datasets <- SetIdent(Inhibitory_datasets, value = "Experiment")
Inhibitory_datasets_WT <- subset(Inhibitory_datasets, ident = "WT")
Inhibitory_datasets_WT <- SetIdent(Inhibitory_datasets_WT, value = "Mitotic_and_tips")

Ebf1_Isl1_WT <- subset(Inhibitory_datasets_WT, ident = "Ebf1_Isl1")
mit_WT <- subset(Inhibitory_datasets_WT, ident = c("Mitotic"))
Ebf1_Isl1_WT <- subset(Ebf1_Isl1_WT, downsample = 1000)
mit_WT <- subset(mit_WT, downsample = 1000)

#Ebf1 ----
Inhibitory_datasets <- SetIdent(Inhibitory_datasets, value = "clones_both_stages") #all clones from e12.5 and e16.5
Ebf1_Isl1_clone <- subset(Inhibitory_datasets, ident = c("Mitotic/Ebf1_Isl1", "Ebf1_Isl1", "Ebf1_Isl1/No", "Mitotic/Ebf1_Isl1/No")) # all the clones with cells in Ebf1 state
Ebf1_Isl1_clone <- SetIdent(Ebf1_Isl1_clone, value = "Mitotic_and_tips")
Ebf1_Isl1_clone_mit <- subset(Ebf1_Isl1_clone, ident = c("Mitotic"))
Ebf1_Isl1_clone_post <- subset(Ebf1_Isl1_clone, ident = c("Ebf1_Isl1"))
#dataframe Ebf1 reference cluster
scaled_matrix_Ebf1_Isl1_WT <- as(as.matrix(GetAssayData(Ebf1_Isl1_WT, assay = "RNA", slot = "scale.data")), 'sparseMatrix')
Ebf1_Isl1_WT.df <- as.data.frame(scaled_matrix_Ebf1_Isl1_WT)
#dataframe random mitotic cells 
scaled_matrix_mit_WT <- as(as.matrix(GetAssayData(mit_WT, assay = "RNA", slot = "scale.data")), 'sparseMatrix')
mit_WT.df <- as.data.frame(scaled_matrix_mit_WT)
#dataframe mitotic cells in Ebf1 clones
scaled_matrix_Ebf1_Isl1_clone_mit <- as(as.matrix(GetAssayData(Ebf1_Isl1_clone_mit, assay = "RNA", slot = "scale.data")), 'sparseMatrix')
Ebf1_Isl1_clone_mit.df <- as.data.frame(scaled_matrix_Ebf1_Isl1_clone_mit)
#dataframe postmitotic cells in Ebf1 clones
scaled_matrix_Ebf1_Isl1_clone_post <- as(as.matrix(GetAssayData(Ebf1_Isl1_clone_post, assay = "RNA", slot = "scale.data")), 'sparseMatrix')
Ebf1_Isl1_clone_post.df <- as.data.frame(scaled_matrix_Ebf1_Isl1_clone_post)

#Correlation: all the possible cell combinations
cor_vals_post_clone_Ebf1 <- combn(1:ncol(Ebf1_Isl1_clone_post.df), 2, function(x) cor(Ebf1_Isl1_WT.df[, x[1]], Ebf1_Isl1_clone_post.df[, x[2]]))
cor_vals_mit_clone_Ebf1 <- combn(1:ncol(Ebf1_Isl1_clone_mit.df), 2, function(x) cor(Ebf1_Isl1_WT.df[, x[1]], Ebf1_Isl1_clone_mit.df[, x[2]]))
cor_vals_mit_Wt <- combn(1:ncol(mit_WT.df), 2, function(x) cor(Ebf1_Isl1_WT.df[, x[1]], mit_WT.df[, x[2]]))
cor_vals_post_Ebf1 <- combn(1:ncol(Ebf1_Isl1_WT.df), 2, function(x) cor(Ebf1_Isl1_WT.df[, x[1]], Ebf1_Isl1_WT.df[, x[2]]))
hist(cor_vals_post_clone_Ebf1, xlim=range(-1,1),ylim=range(0,10), freq = F)
hist(cor_vals_mit_clone_Ebf1, xlim=range(-1,1),ylim=range(0,10), freq = F)
hist(cor_vals_mit_Wt, xlim=range(-1,1),ylim=range(0,10), freq = F)
hist(cor_vals_post_Ebf1, xlim=range(-1,1),ylim=range(0,10), freq = F)

