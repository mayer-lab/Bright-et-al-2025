# Bright-et-al-2025

This repository contains the code associated with the study "Temporal control of progenitor competence shapes maturation in GABAergic neuron development in mice".


## Authors

Ann Rose Bright, Yana Kotlyarenko, Florian Neuhaus, Diana Rodrigues, Chao Feng, Christian Peters, Ilaria Vitali, Elif Dönmez, Michael H. Myoga, Elena Dvoretskova, Christian Mayer

**Correspondence:** [christian.mayer@bi.mpg.de](mailto:christian.mayer@bi.mpg.de)

## Data

#### scRNA-seq
GEO accession number GSE255455: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255455

#### scATAC-seq
GEO accession number GSE255104: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255104

#### TotalRNA-seq
GEO accession number GSE255103: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255103

#### CUT&RUN:
GEO accession number GSE255103: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255103

## Scripts

There are individual scripts for different parts of the analysis.

#### scRNA-seq ($\star$) 
This directory contains code used to process and analyze scRNA-seq data. Below is a list of the files and a brief explanation of what each script does:

- `Fig1c_f.R`: Create initial UMAP visualizations.
- `Fig1g_SF3e.R`: Perform correlation analysis of progenitors across stages and pallium/subpallium.
- `Fig1h_SF2ab_SF5e.R`: This R script calculates and plots the abundance of post-mitotic cell states across stages. For further characterization of post-mitotic precursor states we also perform label transfer, mapping to post-natal cell types and compare different clustering algorithms.
- `EDF1b_EDF2a_f.R`: Characterization of cells in joined inhibitory object.
- `SF4a_d.R`: Inference of genes that are linked to maturation of inhibitory or excitatory neurons. After sorting single cells in each stage by their pseudotime, the script infers genes variable along the trajectory. Following this expression of genes linked to maturation is plotted in heatmaps.

#### Barcode lineage tracing ($\blacktriangle$) 
- `Fig1m_o_SF8a_d.R`: This R script performs the analysis of scRNA-seq TrackerSeq datasets (e12.5+96h and e16.5+96h). It subsets the TrackerSeq datasets from the common scRNA-seq datasets pool and then selects the cells at the tip of the developmental branches. It generates an UpSet plot of cell state intesections within the clone groups. The script also creates cell subsets based on their cell state (mitotic or branch tip) and the experiment type (TrackerSeq or reference scRNAseq without lineage tracing). It calculates the Pearson correlation coefficients between cells in the different groups and the reference group.
- `SF7a_d.R`: Lineage coupling analysis (Bandler et al. 2022). Input data for lineage coupling analysis is generated and results are plotted.

#### Fluorescent birthdating ($\blacksquare$)
- `Fig2ab.R`: Plot the initial UMAP and cell state abundance of cells from FlashTag datasets.
- `Fig2c_e_EDF3a_e_SF5a.R`: This R script performs the analysis of scRNA-seq FlashTag datasets (e12.5+6h, e16.5+6h, e12.5+96h). It subsets the FlashTag datasets from the common scRNA-seq datasets pool and compares the pseudotime scores of different experimental conditions. It proceeds by subsetting the postmitotic cells of the FlashTag datasets, and identifies the genes specific to each 6h cohort (e12.5+6h and e16.5+6h). Finally, it generates a plot of the identified genes in each FlashTag condition.   

#### scATAC-seq
- `Fig2f_k.R`: Analyze scATAC-seq FlashTag datasets from two developmental stages (e12.5+6h, e16.5+6h). Script contains code to identify and plot peak categories such as e12-sites, e16-sites, and overlapping sites from scATAC-seq data. Additionally, it demonstrates how fragment distribution can be calculated and plotted for both scATAC-seq and H3K4me1 Chip-seq data. The script also calculates and plots peak accessibility across pseudotime, categorizing regions into 'initial' and 'intermediate' stages.The script contains codes on how to calculate the aggregate footprint plot for TF of interest based on the analysis performed by TOBIAS on scATAC-seq. Finally, the script generates a coverage plot for NFIB. 
- `EDF4c_f.R`: Create results from extended data figure 4.
- `EDF5d.R`: Plot results from combinatorial binding analysis.

#### Gene Regulatory Network (eGRN) analysis
- `Fig3ab_EDF5a_c_EDF5f_h_SF11a.R`: This R script takes as input the results from eGRN-inference using Scenicplus. To infer gene-regulatory interactions that are specific to cell states and stages, the script uses cluster-specific module enrichment values to subset the large network into subnetworks. These subnetworks are consequently combined to infer modules that are dynamic across either cell state or stage. Subnetworks are represented and plotted using igraph library. Finally the script infers up-stream regulators of genes, that are linked to maturation genes, plots the corresponding network interactions and quantifies top upstream regulators.
- `SF10b_d.R`: Create basic plots to describe the eGRN.
- `SF11ab.R`: This script focuses on the co-binding of Nfib with Meis2 or Tcf4, as predicted by Scenicplus. The number of co-bound genes per cell state and stage is plotted and enriched GO-terms of co-bound genes in e16.5 are plotted.

#### CUT&RUN analysis
- `Fig3cdf_SF12a_d.R`: Analyze results from CUT&RUN: plot heatmaps and generate quantifications of dynamic read coverage.

#### Transplantation datasets analysis
- `Fig4bc_EDF6e_g.R`: This R script performs the analysis of transplantation total RNA-seq datasets obtained by transplanting e12.5 apical progenitors into e12.5 or e16.5 environment (APe12.5->e12.5,APe12.5->e16.5). It normalises the count matrix with the Fragments Per Kilobase of transcript per Million mapped reads (fpkm) method. It proceeds by employing Bisque for cell composition estimation and utilises the common pool of single-cell RNA datasets as reference. This method is adapted from Jew et al., 2020 (Nature Communications), doi.org/10.1038/s41467-020-15816-6. The script also subsets the count matrix by the 2000 highly variable genes obtained from the common pool of single-cell RNA datasets. I then proceeds by identifying the genes specific to each condition (APe12.5->e12.5 or APe12.5->e16.5) and intersecting them with the list of Nfib, Tcf4, Meis2 target genes obtained from the network analysis. The analysis of datasets obtained by tansplanting e16.5 apical progenitors into e16.5 or e12.5 environment (APe16.5->e16.5,APe16.5->e12.5) was performed in the same way and is not included in the script. 

#### Nfib/x KO
- `Fig4hjk_EDF7eg_DF14a_SF15c.R`: Taking the processed seurat object for Nfib/x KO, this script generates plots for differential abundance between guides and control for clusters and croad cell states. Furthermore it plots changes in pseudotime across conditions and differential expression of genes of interest. Additionally this script performs DE analysis and compares DE genes with genes bound by NFIB (according to CUT&RUN).
- `SF15ab.R`: MILO analysis for inhibitory subset of gNfib/x KO
- `SF16a_c.R`: MILO analysis for excitatory subset of gNfib/x KO

#### Nfib OE
- `Fig4ij_EDF7fh_SF14b_SF15d.R`: Taking the processed seurat object for Nfib OE, this script generates plots for differential abundance between guides and control for clusters and croad cell states. Furthermore it plots changes in pseudotime across conditions. Additionally this script performs DE analysis and compares DE genes with genes bound by NFIB (according to CUT&RUN).


Please note that additional descriptions and usage details can be found within each script.
