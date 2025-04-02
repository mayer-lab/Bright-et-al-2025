library(ggplot2)
library(ggrepel)
library(biomaRt)
library(dplyr)
library(countToFPKM)
library(expss)
library(tibble)
library(Seurat)
library(vioplot)
library(rstatix)


load("Processed_Objects/E12_transplant_matrix.Rdata")
load("Processed_Objects/E16_transplanted.Rdata")


load("Processed_Objects/Inhibitory_datasets.Rdata") #reference scRNAseq datasets
cell_states <- c("Mitotic","Abracl","Tshz1","Tcf4_Nr2f2","Six3_Gucy1a3","Ebf1_Isl1","Gucy1a3","Nkx2_1","Npy","Maf_Sst","Snhg11")
cell_states_col <- c('#555D55','#e6f598','#fee08b','#fdae61','#f46d43','#9e0142','#d53e4f','#abdda4','#66c2a5','#3288bd','#5e4fa2')
levels_mtx <- c("Counts_WT1","Counts_WT2","Counts_WT3","Counts_WT4","Counts_WT5","Counts_TR1","Counts_TR2","Counts_TR3","Counts_TR4", "Counts_TR5")

# Transplantation data analysis: e12.5 cells  transplanted to e12.5 or e16.5

# fpkm ----
e12_transplant_matrix <- na.omit(e12_transplant_matrix)
e12_transplant_matrix[, -1] <- lapply(e12_transplant_matrix[, -1], as.numeric)
# Calculate mean length and convert count data to FPKM
mean_length <- mean(e12_transplant_matrix$LENGTH)
fpkm_e12 <- countToFPKM::fpkm(e12_transplant_matrix[, -1], e12_transplant_matrix$LENGTH, rep(mean_length, ncol(e12_transplant_matrix) - 1))
# Calculate log10 of FPKM+1 and create a boxplot (SFig.13)
log10_fpkmData12 <- log10(fpkm_e12 + 1)
boxplot(data.frame(log10_fpkmData12), ylab = "log10(fpkmData12)", main = "E12 transplant")
# Bisque mapping ----      
Inhibitory_datasets <- SetIdent(Inhibitory_datasets, value = "PT_parts")
sc.eset <- BisqueRNA::SeuratToExpressionSet(Inhibitory_datasets, delimiter="_", position=2, version="v3")
bulk.eset <- Biobase::ExpressionSet(assayData = fpkm_e12)
e12_Bisque <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, use.overlap=F)

estimates_E12fpkm <- as.data.frame(e12_Bisque$bulk.props)
estimates_E12fpkm.df.1 <- estimates_E12fpkm %>% 
  rownames_to_column(var = "Cell_types") %>% 
  tidyr::pivot_longer(cols = -Cell_types, names_to = "variable", values_to = "value")
#plot (EDF 6e)

write.table(estimates_E12fpkm.df.1, "Results/source_data/edf6e_left.csv", sep = ",", quote = F, row.names = F)

ggplot(estimates_E12fpkm.df.1, aes(
  x = factor(variable, levels = c(levels_mtx)),
  y = factor(Cell_types, levels = c(cell_states)),
  size = value, color = Cell_types
)) +
  scale_color_manual(values = cell_states_col, breaks = cell_states)+
  geom_point() +
  theme_classic() +
  scale_size_area(max_size = 10) + NoLegend()

# Calculate the average PT from the reference dataset---- 
Inhibitory_datasets@meta.data$Pseudotime <- as.numeric(Inhibitory_datasets@meta.data$Pseudotime)
list_fine_split <- SplitObject(Inhibitory_datasets, split.by = "PT_parts")
# Initialize an empty vector to store the median values
Avg_score_list <- numeric(length(cell_states))
# Loop through each cell state and calculate the median
for (i in 1:length(cell_states)) {
  Avg_score_list[i] <- median(list_fine_split[[cell_states[i]]]@meta.data$Pseudotime)
}

# Weight PT based on Bisque mapping ----
estimates_E12fpkm_ordered <- estimates_E12fpkm[match(cell_states, rownames(estimates_E12fpkm)), ] #reorder the fpkm bisque estimates
avgPTdf <- data.frame(Dataset = c(rownames(estimates_E12fpkm_ordered)), 
                      avgPT = c(Avg_score_list))#create the average score dataframe 
E12fpkm_assigned2 <- cbind(estimates_E12fpkm_ordered, avgPTdf) #unite the sverage scores with the E12 fpkm normalised bisque estimates
E12fpkm_assigned2 <- E12fpkm_assigned2[ -c(11) ] #remove reduntant names

num_columns <- 10  # Number of new columns to create
# Create new columns
for (i in 1:num_columns) {
  col_name <- ifelse(i <= 5, paste0("W", i), paste0("T", i - 5))
  avgPTdf[[col_name]] <- avgPTdf$avgPT
}

avgPTdf2 <- avgPTdf[ -c(1,2) ]
rownames(avgPTdf2) <- rownames(E12fpkm_assigned2)
avgPTdf3  <- reshape2::melt(avgPTdf2)
avgPTdf3 <- as.data.frame(avgPTdf3)
weights <-  E12fpkm_assigned2[ -c(11) ]
weights1 <- as.data.frame(reshape2::melt(weights[-11] * 100))
tr12_weighted.df <-  as.data.frame(expss::weight_by(avgPTdf3, weight = weights1$value))
tr12_weighted.df$W_T <- ifelse(
  tr12_weighted.df$variable %in% c("W1", "W2", "W3", "W4", "W5"),
  "WT",
  ifelse(
    tr12_weighted.df$variable %in% c("T1", "T2", "T3", "T4", "T5"),
    "TR",
    tr12_weighted.df$variable  
  )
)

#plot (Fig 4b)
vioplot(tr12_weighted.df$value ~ tr12_weighted.df$W_T,
        horizontal = F,col = c("#5776A7","#E16B8F"))
#Wilcoxon test ----
stat.test_tr12 <- tr12_weighted.df %>% 
  wilcox_test(value ~ W_T) %>%
  add_significance()


#Prep metadata for DESeq----
e12_transplant_matrix2<-e12_transplant_matrix[ ,-c(1)]
df_metadata_12 <- data.frame(Dataset = c(colnames(e12_transplant_matrix2)) , 
                             Condition= c("control", "control", "control", "control", "control", "treated","treated", "treated","treated", "treated"),
                             Batch=c("1", "2", "3","4","5", "1","2", "3","4","5")) #control = e12->e12 transplant, treated = e12->e16 transplant
rownames(df_metadata_12) <- df_metadata_12$Dataset
df_metadata_12<-df_metadata_12[ ,-c(1)]
mtx_counts12 <-(data.matrix(frame = e12_transplant_matrix2, rownames.force = T))
# DESeq ----
neurogenes2000 <- VariableFeatures(Inhibitory_datasets)
E12_matrix_neuro <- e12_transplant_matrix2[ c(neurogenes2000),]
E12_matrix_neuro <- na.omit(E12_matrix_neuro)
E12_matrix_neuro2 <-(data.matrix(frame = E12_matrix_neuro, rownames.force = T))


#The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input.
library("DESeq2")
dds_12_neurogenes <- DESeqDataSetFromMatrix(countData = E12_matrix_neuro2, colData = df_metadata_12, design = ~ Condition)
dds_12_neurogenes_deseq <- DESeq(dds_12_neurogenes)
res_12_neurogenes <- results(dds_12_neurogenes_deseq)
res_12_neurogenes <- res_12_neurogenes[!is.na(res_12_neurogenes$padj),]

upregulated <- rownames(res_12_neurogenes[res_12_neurogenes$log2FoldChange > 1 & res_12_neurogenes$pvalue < 0.05,])
downregulated <- rownames(res_12_neurogenes[res_12_neurogenes$log2FoldChange < -1 & res_12_neurogenes$pvalue < 0.05,])

#loading genes targeted by Nfib, Tcf4, Meis2
Target_genes_1 <- readRDS("/datastore_share/Users/neuhaus/for_Yana/e16_FOI_target_genes.rds")
Target_genes <- c(Target_genes_1$Nfib,Target_genes_1$Tcf4, Target_genes_1$Meis2,
                  Target_genes_1$Nfib_and_Tcf4_and_Meis2, Target_genes_1$`(Nfib_and_Tcf4)_or_(Nfib_and_Meis2)` )

#Plotting up- and down-regulated targets of Nfib, Tcf4 and Meis2 
res12_df <- as.data.frame(res_12_neurogenes)
res12_df$diffexpressed <- "NO"
res12_df$diffexpressed[res12_df$log2FoldChange > 1 & res12_df$pvalue < 0.05] <- "UP"
res12_df$diffexpressed[res12_df$log2FoldChange < -1 & res12_df$pvalue < 0.05] <- "DOWN"
res12_df$downstream[rownames(res12_df) %in% c(Target_genes)] <- rownames(res12_df)[rownames(res12_df) %in% c(Target_genes)]
res12_df$intersection <- ifelse(res12_df$diffexpressed %in% c("UP", "DOWN") & !is.na(res12_df$downstream), "UpDown", "None")

ggplot(data = res12_df, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_point() +
  geom_text_repel(aes(label = ifelse(intersection == "UpDown", as.character(downstream), "")), size = 3, color = "black", max.overlaps = 100) +
  #scale_color_manual(values = c("#efedf5", "#A4AEC8")) +
  theme_minimal()


################################################################################
## code for EDF 6 f,g (e16 transplant)
head(E16_transplanted)

# fpkm ----
load(file = "Processed_Objects/fpkm_16.Rdata")
head(fpkm_16)

# Calculate log10 of FPKM+1 and create a boxplot (SFig.13)
log10_fpkmData16 <- log10(fpkm_16 + 1)
boxplot(data.frame(log10_fpkmData16), ylab = "log10(fpkmData16)", main = "E16 transplant")
# Bisque mapping ----      
Inhibitory_datasets <- SetIdent(Inhibitory_datasets, value = "PT_parts")
sc.eset <- BisqueRNA::SeuratToExpressionSet(Inhibitory_datasets, delimiter="_", position=2, version="v3")
bulk.eset <- Biobase::ExpressionSet(assayData = fpkm_16)
e16_Bisque <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, use.overlap=F)

estimates_E16fpkm <- as.data.frame(e16_Bisque$bulk.props)
estimates_E16fpkm.df.1 <- estimates_E16fpkm %>% 
  rownames_to_column(var = "Cell_types") %>% 
  tidyr::pivot_longer(cols = -Cell_types, names_to = "variable", values_to = "value")

## update levels vector:
levels_mtx <- c("Counts_WT1_16","Counts_WT2_16", "Counts_WT3_16", "Counts_WT4_16", "Counts_TR1_16", "Counts_TR2_16", "Counts_TR3_16", "Counts_TR4_16")

write.table(estimates_E16fpkm.df.1, "Results/source_data/edf6e_right.csv", sep = ",", quote = F, row.names = F)

ggplot(estimates_E16fpkm.df.1, aes(
  x = factor(variable, levels = c(levels_mtx)),
  y = factor(Cell_types, levels = c(cell_states)),
  size = value, color = Cell_types
)) +
  scale_color_manual(values = cell_states_col, breaks = cell_states)+
  geom_point() +
  theme_classic() +
  scale_size_area(max_size = 10) + NoLegend()

# Calculate the average PT from the reference dataset---- 
Inhibitory_datasets@meta.data$Pseudotime <- as.numeric(Inhibitory_datasets@meta.data$Pseudotime)
list_fine_split <- SplitObject(Inhibitory_datasets, split.by = "PT_parts")
# Initialize an empty vector to store the median values
Avg_score_list <- numeric(length(cell_states))
# Loop through each cell state and calculate the median
for (i in 1:length(cell_states)) {
  Avg_score_list[i] <- median(list_fine_split[[cell_states[i]]]@meta.data$Pseudotime)
}

# Weight PT based on Bisque mapping ----
estimates_E16fpkm_ordered <- estimates_E16fpkm[match(cell_states, rownames(estimates_E16fpkm)), ] #reorder the fpkm bisque estimates
avgPTdf <- data.frame(Dataset = c(rownames(estimates_E16fpkm_ordered)), 
                      avgPT = c(Avg_score_list))#create the average score dataframe 
E16fpkm_assigned2 <- cbind(estimates_E16fpkm_ordered, avgPTdf) #unite the sverage scores with the E12 fpkm normalised bisque estimates
E16fpkm_assigned2 <- E16fpkm_assigned2[ -c(9) ] #remove reduntant names

num_columns <- 8  # Number of new columns to create
# Create new columns
for (i in 1:num_columns) {
  col_name <- ifelse(i <= 4, paste0("W", i), paste0("T", i - 4))
  avgPTdf[[col_name]] <- avgPTdf$avgPT
}

avgPTdf2 <- avgPTdf[ -c(1,2) ]
rownames(avgPTdf2) <- rownames(E16fpkm_assigned2)
avgPTdf3  <- reshape2::melt(avgPTdf2)
avgPTdf3 <- as.data.frame(avgPTdf3)
weights <-  E16fpkm_assigned2[ -c(9) ]
weights1 <- as.data.frame(reshape2::melt(weights[-9] * 100))
tr16_weighted.df <-  as.data.frame(expss::weight_by(avgPTdf3, weight = weights1$value))
tr16_weighted.df$W_T <- ifelse(
  tr16_weighted.df$variable %in% c("W1", "W2", "W3", "W4"),
  "WT",
  ifelse(
    tr16_weighted.df$variable %in% c("T1", "T2", "T3", "T4"),
    "TR",
    tr16_weighted.df$variable  
  )
)

#plot
vioplot(tr16_weighted.df$value ~ tr16_weighted.df$W_T,
        horizontal = F,col = c("#5776A7","#E16B8F"))
#Wilcoxon test ----
stat.test_tr16 <- tr16_weighted.df %>% 
  wilcox_test(value ~ W_T) %>%
  add_significance()


#Prep metadata for DESeq----
load("Processed_Objects/Transpl_E16.Rdata")
E16_transplanted <- Transpl_E16
load("Processed_Objects/neurogenes2000.Rdata")
load("Processed_Objects/df_metadata_16.Rdata")
E16_variablegenes <- E16_transplanted[ c(neurogenes2000),]
E16_variablegenes <- na.omit(E16_variablegenes)
mtx_counts16v <-(data.matrix(frame = E16_variablegenes, rownames.force = T))

library("DESeq2")
dds_16_neurogenes <- DESeqDataSetFromMatrix(countData = mtx_counts16v, colData = df_metadata_16, design = ~ Condition)
dds_16_neurogenes_des <- DESeq(dds_16_neurogenes)
res16_neurogenes <- results(dds_16_neurogenes_des)
res16_neurogenes <- res16_neurogenes[!is.na(res16_neurogenes$padj),]
rownames(res16_neurogenes)
min(res16_neurogenes$padj)

up16_neurogenes <- rownames(res16_neurogenes[res16_neurogenes$log2FoldChange > 1 & res16_neurogenes$pvalue < 0.05,])

down16_neurogenes <- rownames(res16_neurogenes[res16_neurogenes$log2FoldChange < - 1 & res16_neurogenes$pvalue < 0.05,])
e16_target_genes <- readRDS("Processed_Objects/e16_FOI_target_genes.rds")

Target_genes <- c(e16_target_genes$Nfib,e16_target_genes$Tcf4, e16_target_genes$Meis2,
                  e16_target_genes$Nfib_and_Tcf4_and_Meis2, e16_target_genes$`(Nfib_and_Tcf4)_or_(Nfib_and_Meis2)` )


list_E16_all1 <- intersect(up16_neurogenes, Target_genes)
list_E12_all1 <- intersect(down16_neurogenes, Target_genes)



res16_df <- as.data.frame(res16_neurogenes)
res16_df$diffexpressed <- "NO"
res16_df$diffexpressed[res16_neurogenes$log2FoldChange > 1 & res16_neurogenes$pvalue < 0.05] <- "UP"
res16_df$diffexpressed[res16_neurogenes$log2FoldChange < -1 & res16_neurogenes$pvalue < 0.05] <- "DOWN"
res16_df$delabel <- NA
res16_df$delabel[res16_df$diffexpressed != "NO"] <- rownames(res16_df)[res16_df$diffexpressed != "NO"]
res16_df$common <- NA
res16_df$common[rownames(res16_df) %in% c(Target_genes)] <- rownames(res16_df)[rownames(res16_df) %in% c(Target_genes)]

ggplot(data=res16_df, aes(x=log2FoldChange, y=-log10(pvalue), col = diffexpressed)) +
  geom_point() +
  theme_minimal() +
  geom_text(label = res16_df$common, color = "black") +
  scale_color_manual(values = c("#807dba","#efedf5", "#e7298a"))

write.table(res16_df, "Results/source_data/edf6g.csv", sep = ",", quote = F, row.names = T)
