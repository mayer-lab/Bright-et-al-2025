## correlation analysis ##
library(Seurat)
library(corrplot)

source("Scripts/lib.R")

## load excitatory inhibitory data:
load("Processed_Objects/EXCIT_INHIBIT_cleaned.Rdata")
EXCIT_INHIBIT_cleaned <- SetIdent(EXCIT_INHIBIT_cleaned, value = "Annotated2")
## subset clusters containing apical progenitors:
AP_DV <- subset(EXCIT_INHIBIT_cleaned, ident = c("Apical_progenitors"))

AP_DV <- SetIdent(AP_DV, value = "Stage_DV2")
AP_DV_same_stages <- subset(AP_DV, ident = c("E12_E", "E14_E", "E16_E", "E12_I", "E14_I", "E16_I"))
AP_DV_same_stages@meta.data$Stage_DV2 <- factor(AP_DV_same_stages@meta.data$Stage_DV2, levels = c("E12_E", "E14_E", "E16_E", "E12_I", "E14_I", "E16_I"))
AP_DV_same_stages <- SetIdent(AP_DV_same_stages, value = "Stage_DV2")

s_genes <- convert_to_mouse(cc.genes$s.genes)
g2m_genes <- convert_to_mouse(cc.genes$g2m.genes)

################################################################################
## fig 1g:

table(AP_DV_same_stages@meta.data$Stage_DV2)

## re-do scaling:
AP_DV_same_stages[["percent.mt"]] <- PercentageFeatureSet(AP_DV_same_stages, pattern = "^mt-")
AP_DV_same_stages <- NormalizeData(AP_DV_same_stages)
AP_DV_same_stages <- FindVariableFeatures(AP_DV_same_stages)
variable_features_before <- VariableFeatures(AP_DV_same_stages)
## remove confounding variable features:
variable_features<- variable_features_before[! variable_features_before %in% c ("Hbb-bs","Hbb-bt", "Hist1h1b", "Hba-a1" ,"Hist1h2ae", "Hba-a2", "mt-Atp6", "mt-Co2",  "mt-Co3",  "Rspo3",   "mt-Nd2" ,"mt-Co1", "Rps12", "Rpl9", "Rpl10", "Rpl6", "Rps26",
                                                                                "Rpl3","Rpl4","Rpl5","Rpl6", "Rpl7","Rpl7a","Rpl8","Rpl9","Rpl10",
                                                                                "Rpl10a", "Rpl11","Rpl12", "Rpl13","Rpl13a","Rpl14","Rpl15","Rpl17", "Rpl18",
                                                                                "Rpl18a","Rpl19","Rpl21","Rpl22","Rpl23","Rpl23a","Rpl24","Rpl26","Rpl27", "Rpl27a",
                                                                                "Rpl28","Rpl29","Rpl30","Rpl31","Rpl32", "Rpl34","Rpl35", "Rpl35a","Rpl36","Rpl44",
                                                                                "Rpl37","Rpl37a","Rpl38","Rpl39","Rpl40","Rpl41","Rplp0", "Rplp1", "Rplp2" )]   
AP_DV_same_stages <- CellCycleScoring(AP_DV_same_stages, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
AP_DV_same_stages$CC.Difference <- AP_DV_same_stages$S.Score - AP_DV_same_stages$G2M.Score
AP_DV_same_stages<- ScaleData(AP_DV_same_stages, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt","CC.Difference"), verbose = T, features = variable_features)

AP_DV_same_stages <- SetIdent(AP_DV_same_stages, value = "Stage_DV2")
av.exp_same <- AverageExpression(AP_DV_same_stages, slot = "scale.data")
ap_cor_same <- cor(av.exp_same$RNA)
testRes = cor.mtest(ap_cor_same, conf.level = 0.95)

## write table:
write.table(ap_cor_same, "Results/source_data/fig1g.csv", sep = ",", quote = F, row.names = T)
corrplot(ap_cor_same, p.mat = testRes$p, method = 'color', type = "upper", diag = T,
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 2, 
         insig = 'label_sig', pch.col = 'grey20', title = "",tl.col = "black")




################################################################################
## sfig 3e:

VlnPlot(object = AP_DV_same_stages, features = c( "nCount_RNA"))

AP_DV_subset<- subset(x = AP_DV_same_stages, subset = nCount_RNA < 10000)
VlnPlot(object = AP_DV_subset, features = c("nCount_RNA"))


#AP_DV_subset <- subset(AP_DV, cells = (WhichCells(AP_DV, downsample = 200)))
table(AP_DV_subset@meta.data$Stage_DV2)
AP_DV_subset <- SetIdent(AP_DV_subset, value = "Stage_DV2")

AP_DV_subset[["percent.mt"]] <- PercentageFeatureSet(AP_DV_subset, pattern = "^mt-")
AP_DV_subset <- NormalizeData(AP_DV_subset)
AP_DV_subset <- FindVariableFeatures(AP_DV_subset)
variable_features_before <- VariableFeatures(AP_DV_same_stages)
## remove confounding variable features:
variable_features<- variable_features_before[! variable_features_before %in% c ("Hbb-bs","Hbb-bt", "Hist1h1b", "Hba-a1" ,"Hist1h2ae", "Hba-a2", "mt-Atp6", "mt-Co2",  "mt-Co3",  "Rspo3",   "mt-Nd2" ,"mt-Co1", "Rps12", "Rpl9", "Rpl10", "Rpl6", "Rps26",
                                                                                "Rpl3","Rpl4","Rpl5","Rpl6", "Rpl7","Rpl7a","Rpl8","Rpl9","Rpl10",
                                                                                "Rpl10a", "Rpl11","Rpl12", "Rpl13","Rpl13a","Rpl14","Rpl15","Rpl17", "Rpl18",
                                                                                "Rpl18a","Rpl19","Rpl21","Rpl22","Rpl23","Rpl23a","Rpl24","Rpl26","Rpl27", "Rpl27a",
                                                                                "Rpl28","Rpl29","Rpl30","Rpl31","Rpl32", "Rpl34","Rpl35", "Rpl35a","Rpl36","Rpl44",
                                                                                "Rpl37","Rpl37a","Rpl38","Rpl39","Rpl40","Rpl41","Rplp0", "Rplp1", "Rplp2" )] 
AP_DV_subset <- CellCycleScoring(AP_DV_subset, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
AP_DV_subset$CC.Difference <- AP_DV_subset$S.Score - AP_DV_subset$G2M.Score
AP_DV_subset<- ScaleData(AP_DV_subset, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt","CC.Difference"), verbose = T, features = variable_features)

AP_DV_subset <- SetIdent(AP_DV_subset, value = "Stage_DV2")
av.exp_same <- AverageExpression(AP_DV_subset, slot = "scale.data")
ap_cor_same <- cor(av.exp_same$RNA)
testRes = cor.mtest(ap_cor_same, conf.level = 0.95)

write.table(ap_cor_same, "Results/source_data/sf3e.csv", sep = ",", quote = F, row.names = T)
corrplot(ap_cor_same, p.mat = testRes$p, method = 'color', type = "upper", diag = T,
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 2, 
         insig = 'label_sig', pch.col = 'grey20', title = "",tl.col = "black")


