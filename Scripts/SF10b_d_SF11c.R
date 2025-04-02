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


## SF 11c: pct stage-enriched regions contained in network
promoter_sub <- readRDS("Processed_Objects/cutnrun_promoter_ranges.rds")

Nfib_eRegulon_sub <- eRegulon_md[eRegulon_md$TF == "Nfib", ]
nfib_target_peaks <- unique(Nfib_eRegulon_sub$Region)

nfib_target_gr <- GRanges(
  seqnames = sapply(nfib_target_peaks, function(el) {strsplit(el,":")[[1]][1]}),
  ranges = IRanges(
    start = as.numeric(sapply(nfib_target_peaks, function(el) {strsplit(el,"[:-]")[[1]][2]})),
    end = as.numeric(sapply(nfib_target_peaks, function(el) {strsplit(el,"[:-]")[[1]][3]})),
    names = nfib_target_peaks
  ),
  gene = sapply(nfib_target_peaks, function(el) {paste(Nfib_eRegulon_sub$Gene[Nfib_eRegulon_sub$Region == el], collapse = ",")})
)

nfib_target_gr <- annotatePeak(nfib_target_gr, tssRegion = c(-3000,3000),
                               TxDb = edb)
plotAnnoPie(nfib_target_gr)
nfib_promoter_target_gr <- nfib_target_gr@anno[nfib_target_gr@anno$annotation == "Promoter (<=1kb)"]

nfib_target_sub_gr <- subsetByOverlaps(nfib_promoter_target_gr, promoter_sub, maxgap = 100)
print(length(nfib_promoter_target_gr))
print(length(nfib_target_sub_gr))