## GRN description ##

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(EnsDb.Mmusculus.v79)
edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"
library(ggvenn)

## SF10 b-d
eRegulon_md <- read.table("Processed_Objects/eRegulon_metadata_filtered.tsv", h = T, sep = "\t")

num_regions_per_TF <- sapply(unique(eRegulon_md$TF), function(tf) {
  length(unique(eRegulon_md$Region[eRegulon_md$TF == tf]))
})
num_regions_per_TF_df <- data.frame(
  "tf" = names(num_regions_per_TF),
  "number_of_regions" = num_regions_per_TF
)

write.table(num_regions_per_TF_df, "Results/source_data/sf10b.csv", sep = ",", quote = F, row.names = F)
ggplot(num_regions_per_TF_df, aes(number_of_regions)) +
  geom_histogram() +
  geom_vline(xintercept = median(num_regions_per_TF), color = "red") +
  xlab("Number of bound regions per TF") +
  theme(legend.position = 'none', axis.text=element_text(size=14),
        axis.title=element_text(size=14), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=14))

num_genes_per_regions <- sapply(unique(eRegulon_md$Region), function(region) {
  length(unique(eRegulon_md$Gene[eRegulon_md$Region == region]))
})
num_genes_per_regions_df <- data.frame(
  "region" = names(num_genes_per_regions),
  "number_of_genes" = num_genes_per_regions
)
write.table(num_genes_per_regions_df, "Results/source_data/sf10c.csv", sep = ",", quote = F, row.names = F)
ggplot(num_genes_per_regions_df, aes(number_of_genes)) +
  geom_histogram(binwidth = 1) +
  #geom_vline(xintercept = median(num_genes_per_regions), color = "red") +
  xlab("Number of regulated genes per region") +
  theme(legend.position = 'none', axis.text=element_text(size=14),
        axis.title=element_text(size=14), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=14))

num_genes_per_TF <- sapply(unique(eRegulon_md$TF), function(tf) {
  length(unique(eRegulon_md$Gene[eRegulon_md$TF == tf]))
})
num_genes_per_TF_df <- data.frame(
  "TF" = names(num_genes_per_TF),
  "number_of_genes" = num_genes_per_TF
)
write.table(num_genes_per_TF_df, "Results/source_data/sf10d.csv", sep = ",", quote = F, row.names = F)
ggplot(num_genes_per_TF_df, aes(number_of_genes)) +
  geom_histogram() +
  geom_vline(xintercept = median(num_genes_per_TF), color = "red") +
  xlab("Number of regulated genes per TF") +
  theme(legend.position = 'none', axis.text=element_text(size=14),
        axis.title=element_text(size=14), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=14))

