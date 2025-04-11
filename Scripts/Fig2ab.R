library(Seurat)
library(dplyr)
library(ggplot2)

## CFSE Transcriptome: UMAP visulaization and cell state abundance ##

## load data and subset:
load("Processed_Objects/Inhibitory_datasets.Rdata")
Inhibitory_datasets <- SetIdent(Inhibitory_datasets, value = "Experiment" )
CFSE <- subset(Inhibitory_datasets, ident = c("E16+6hCFSE", "E12+6hCFSE","E12+4daysCFSE"))

CFSE$Fine_annotation[CFSE$Fine_annotation == "Snhg11_Lhx8"] <- "Snhg11"
CFSE <- SetIdent(CFSE, value = "Fine_annotation")

CFSE_tips <- subset(CFSE, subset = Fine_annotation %in% c("Ebf1_Isl1","Gucy1a3","Maf_Sst","Tcf4_Nr2f2","Snhg11"))

## UMAP:
DimPlot(Inhibitory_datasets, reduction = "umap2", cells.highlight = list(
  "e12+6" = colnames(Inhibitory_datasets)[Inhibitory_datasets$Experiment == "E12+6hCFSE"],
  "e16+6" = colnames(Inhibitory_datasets)[Inhibitory_datasets$Experiment == "E16+6hCFSE"],
  "e12+96" = colnames(Inhibitory_datasets)[Inhibitory_datasets$Experiment == "E12+4daysCFSE"]
), cols.highlight = c("#5776A7","orange","#E16B8F"))


## get abundance:
cell_df <- CFSE_tips@meta.data[, c("Fine_annotation","Experiment")]
cell_df %>%
  group_by(Experiment,Fine_annotation) %>%
  dplyr::count() -> cell_df
cell_df %>%
  group_by(Experiment) %>%
  mutate(cells_per_exp = sum(n)) -> cell_df
cell_df$frac_cells <- cell_df$n / cell_df$cells_per_exp * 100

write.table(cell_df, "Results/source_data/fig2b.csv", sep = ",", quote = F, row.names = F)

ggplot(cell_df, aes(x = Experiment, y = frac_cells, fill = Fine_annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw()
