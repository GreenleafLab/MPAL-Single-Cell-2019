# MPAL-Single-Cell-2019 Analysis Code (Granja JM*, Klemm SK*, McGinnis LM*, et al. 2019)

**Link** : https://www.nature.com/articles/s41587-019-0332-7

## Please cite : Granja JM et al., Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. Nature Biotechnology (2019)

# Brief Descriptions of Analysis Scripts

## scATAC Analyses

**scATAC_01** - Script for reading in 10x scATAC-seq fragments identify cells using number of fragments and TSS enrichment scores and saving fitlered fragments.

**scATAC_02** - Script for pre-clustering using large windows genome-wide and then calling peaks on putative clusters and create a master peak set

**scATAC_03** - LSI-Clustering + UMAP of scATAC-seq data with visualization and demonstration of how to properly save
umap for projection.

**scATAC_04** - Computing Gene Activity Scores using an adapted form of Cicero (Pliner et al 2018).

**scATAC_05** - Identifying potential disease cells by clustering disease w/ healthy reference, and then projecting these
cells onto healthy hematopoiesis.

## scRNA Analyses

**scRNA_01** - LSI-Clustering + UMAP of scRNA-seq data with visualization and demonstration of how to properly save
umap for projection.

**scRNA_02** - Identifying potential disease cells by clustering disease w/ healthy reference, and then projecting these
cells onto healthy hematopoiesis.

## Integration (scATAC + scRNA) Analyses

**scRNA_scATAC_Integration_01** - Alignment of scRNA and scATAC-seq data using Seurat CCA and identifcation of nearest
neighbors across modalities.

**scRNA_scATAC_Integration_02** - Aggregate scRNA + scATAC-seq data for correlation focused analysis.

**scRNA_scATAC_Integration_03** - Identify putative Peak-To-Gene Links with aligned scATAC and scRNA-seq data aggregates.

**scRNA_scATAC_Integration_04** - Link TFs to putative target genes that are differential in both mRNA and nearby accessibility peaks containing motifs of the TFs.

# Additional Data Download Links

### These links may be moved if we can find a better host for better download speed

## Notes

**.rds** file is an R binarized object to read into R use readRDS(filename)

**SummarizedExperiment** is a class in R see : https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html

**deviations** (TF chromVAR) is a class in R see : https://bioconductor.org/packages/release/bioc/html/chromVAR.html

## Healthy Hematopoiesis

**scATAC-seq Hematopoeisis cell x peak Summarized Experiment** :

**scATAC-seq Hematopoeisis cell x TF chromVAR Summarized Experiment** :

**scRNA-seq Hematopoeisis cell x gene Summarized Experiment** :

**scADT-seq Hematopoeisis cell x antibody Summarized Experiment** :

## Healthy + MPAL Data Sets

**scATAC-seq Hematopoeisis + MPAL cell x peak Summarized Experiment** :

**scATAC-seq Hematopoeisis + MPAL cell x TF chromVAR Summarized Experiment** :

**scRNA-seq Hematopoeisis + MPAL cell x gene Summarized Experiment** :

**scADT-seq Hematopoeisis + MPAL cell x antibody Summarized Experiment** :

## LSI-Projection 

**scATAC-seq saved UMAP embedding** :

**scRNA-seq saved UMAP embedding** :

## Integration

**Peak-To-Gene Linkages** :

## Other

**MPAL Clinical FACS Data** : https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data_Revision_MPAL_FACS_FCS.zip



