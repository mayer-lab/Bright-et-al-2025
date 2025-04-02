## NFIB CUT&RUN ##

## heatmaps created using "fluff-heatmap" library:

##################################################
## Figure 3c:

#fluff heatmap -f peakCalling/MACS_set2/NFIB_R1R2_pvalue_peaks.narrowPeak -d $project/NFIB_R1R2.bw $project/H3K4me3_R1R2.bw $project2/filtered_E12_sort.bw $project2/filtered_E16_sort.bw -C k -k 3 -g -M Pearson -r -c "#991616,#666547,#c03e82,#366f9a" -o profiles/heatmap_nfib_pearson_r.png

#################################################
# Figure 3d:

#fluff heatmap -f $peak/non_overlapping_e12.bed -d $project/NFIB_R1R2.bw $project/H3K4me3_R1R2.bw $project2/filtered_E12_sort.bw $project2/filtered_E16_sort.bw -C k -k 3 -g -M Pearson -r -c "#991616,#666547,#c03e82,#366f9a" -o profiles/heatmap_e12peaks_pearson_r.png
#fluff heatmap -f $peak/non_overlapping_e16.bed -d $project/NFIB_R1R2.bw $project/H3K4me3_R1R2.bw $project2/filtered_E12_sort.bw $project2/filtered_E16_sort.bw -C k -k 3 -g -M Pearson -r -c "#991616,#666547,#c03e82,#366f9a" -o profiles/heatmap_e16peaks_pearson_r.png

#################################################
# Figure 3e:

## the following shows code in bash:

# #!/bin/bash
# 
# regions_file="extended_5kb_regions.txt"
# 
# # Create output directory
# mkdir -p profile5kb
# 
# # Loop through regions
# while read -r region; do
# output_file="profile5kb/profile_${region//:/_}.pdf"
# 
# fluff profile \
# -i "$region" \
# -d $project/NFIB_R1.sorted.bam $project2/filtered_E12_sort.bam $project2/filtered_E16_sort.bam peakCalling/MACS_set2/NFIB_R1R2_pvalue_peaks.bed eRegulon_nfib_tm_regions.bed \
# -s 1,2:3,4:5 \
# -c "#991616,#c03e82,#366f9a,#991616,#2F2F5B" \
# -b white \
# -o "$output_file"
# 
# echo "Profile created for $region: $output_file"
# done < "$regions_file"

#################################################
# Figure 3f:

library(JunJunZai)
hr <- loadHomerRes(homerDir = "/datastore_share/Users/ann/Papers/Manuscript_1/Revision_nature/NFIB_CUTNRUN/analysis/homer/summits/", novo = F,motifIndex = 1:249)

View(hr@`known res`)

expressed_tf <- readRDS("Processed_Objects/mm10_tf_list_filtered_by_expression.rds")
motif_names_to_keep <- expressed_tf$TF_th0.1

# Subset the 'known res' data frame
subset_res <- hr@`known res`[tolower(rownames(hr@`known res`)) %in% tolower(motif_names_to_keep), ]

subset_indices <- which(tolower(rownames(hr@`known res`)) %in% tolower(motif_names_to_keep), )

# Subset the PWM and PFM matrices
subset_pwm <- hr@`known motif PWM matrix`[subset_indices]
subset_pfm <- hr@`known motif PFM matrix`[subset_indices]


# Update the homerResult object with the subsets
hr_subset <- hr
hr_subset@`known res` <- subset_res
hr_subset@`known motif PWM matrix` <- subset_pwm
hr_subset@`known motif PFM matrix` <- subset_pfm

# Check if dimensions align
cat("Subset rows in 'res':", nrow(hr_subset@`known res`), "\n")
cat("Subset length of PWM:", length(hr_subset@`known motif PWM matrix`), "\n")

write.table(hr_subset@`known res`, "Results/source_data/fig3f.csv", sep = ",", quote = F, row.names = F)

set.seed(10)
plotMotifHeatmap(object = hr_subset, type = "known", barWidth = 1)




