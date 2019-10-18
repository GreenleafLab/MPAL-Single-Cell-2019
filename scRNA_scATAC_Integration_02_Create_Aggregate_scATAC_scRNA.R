#Clustering and scATAC-seq UMAP for Hematopoiesis data
#06/02/19
#Cite Granja*, Klemm*, Mcginnis* et al. 
#A single cell framework for multi-omic analysis of disease identifies 
#malignant regulatory signatures in mixed phenotype acute leukemia (2019)
#Created by Jeffrey Granja
library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
set.seed(1)

####################################################
#Functions
####################################################

getGeneGTF <- function(file){
  #Import
  message("Reading in GTF...")
  importGTF <- rtracklayer::import(file)
  #Exon Info
  message("Computing Effective Exon Lengths...")
  exonGTF <- importGTF[importGTF$type=="exon",]
  exonList <- reduce(split(exonGTF, mcols(exonGTF)$gene_id))
  exonReduced <- unlist(exonList, use.names=TRUE)
  mcols(exonReduced)$gene_id <- names(exonReduced)
  mcols(exonReduced)$widths <- width(exonReduced)
  exonSplit <- split(exonReduced$widths, mcols(exonReduced)$gene_id)
  exonLengths <- lapply(seq_along(exonSplit), function(x) sum(exonSplit[[x]])) %>% 
    unlist %>% data.frame(row.names=names(exonSplit), effLength=.)
  #Gene Info
  message("Constructing gene GTF...")
  geneGTF1 <- importGTF[importGTF$type=="gene",]
  geneGTF2 <- GRanges(
      seqnames=paste0("chr",seqnames(geneGTF1)),
      ranges=ranges(geneGTF1),
      strand=strand(geneGTF1),
      gene_name=geneGTF1$gene_name,
      gene_id=geneGTF1$gene_id
    ) %>% keepFilteredChromosomes %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
  mcols(geneGTF2)$exonLength <- exonLengths[geneGTF2$gene_id,]
  return(geneGTF2)
}

featureToGR <- function(feature){
  featureSplit <- stringr::str_split(paste0(feature), pattern = "_", n = 3, simplify = TRUE)
  GRanges(
    seqnames = featureSplit[,1],
    ranges = IRanges(
      start = as.integer(featureSplit[,2]),
      end = as.integer(featureSplit[,3]))
  )
}

grToFeature <- function(gr){
  paste(seqnames(gr),start(gr),end(gr),sep="_")
}

####################################################
#Input Data
####################################################

#Input Files
seATAC_file <- "../../../scATAC/data/scATAC_All_Healthy_MPAL_10x-Rename-UMAP-190505.rds"
ciceroKNN_file <- "../../../scATAC/post/heme_mpal/output/Cicero/Cicero_KNN.rds"
ciceroATAC_file <- "../../../scATAC/post/heme_mpal/output/Cicero/Cicero_Obj.rds"
seRNA_file <- "../../../scRNA/data/scRNA_All_Healthy_MPAL_10x-Rename-Filter10k-190505.rds"
CCA_file <- "../Alignment/output/ATAC-RNA-Alignment/CCA-25-Variable-2000-FeatureSelection-All/Aligned-CCA-Summary.rds"
gtf_file <- "../../../../Healthy_10x/Post3/genes.gtf"

#Get Clusters information for each KNN Group top group/exp wins!
scATAC <- readRDS(seATAC_file)
KNN <- data.frame(readRDS(ciceroKNN_file), stringsAsFactors=FALSE)
KNN <- apply(KNN,2,paste0)
KNNClusters <- apply(KNN, 2, function(x) colData(scATAC)[x,"Clusters"])
KNNGroups <- apply(KNN, 2, function(x) colData(scATAC)[x,"Group"])
KNN_Highest_Cluster <- lapply(seq_len(nrow(KNN)), function(x) names(sort(table(KNNClusters[x,]), decreasing=TRUE))[1]) %>% unlist
KNN_Highest_Experiment <- lapply(seq_len(nrow(KNN)), function(x) names(sort(table(KNNGroups[x,]), decreasing=TRUE))[1]) %>% unlist

####################################################
# scATAC-seq Merging from Cicero
####################################################
ciceroObj <- readRDS(ciceroATAC_file)
se <- SummarizedExperiment(
  assays = SimpleList(counts = assayData(ciceroObj)$exprs),
  rowRanges = featureToGR(featureData(ciceroObj)[[1]]),
  colData = DataFrame(
    row.names = colnames(assayData(ciceroObj)$exprs), 
    clustATAC = KNN_Highest_Cluster, 
    groupATAC = KNN_Highest_Experiment
  )
metadata(se)$knn <- KNN
metadata(se)$knnClust <- KNNClusters
metadata(se)$knnGroup <- KNNGroups
rownames(se) <- grToFeature(rowRanges(se))
saveRDS(se, "results/Save-scATAC-Merged-KNN-SVD.rds")

####################################################
# scRNA-seq Merging from Cicero
####################################################
scRNA <- readRDS(seRNA_file)
CCA <- readRDS(CCA_file)
CCA <- CCA[CCA$corCCA >= 0.4,] #Filter low correlation to CCA R > 0.4

#Get scRNA Matrix
scMat <- assay(scRNA)

#Create aggregated Matched RNA Matrix exclude non-mappings
matRNA <- matrix(NA, ncol = nrow(KNN), nrow = nrow(scMat))
for(i in seq_len(nrow(KNN))){
  if(i%%100==0) print(i)
  knnIdx <- paste0(t(KNN[i,]))
  rnaIdx <- CCA$y[match(knnIdx, CCA$x)]
  rnaIdx <- na.omit(rnaIdx)
  matRNA[, i] <- Matrix::rowSums(scMat[,rnaIdx])
}
colnames(matRNA) <- colnames(se)

#Create aggregated summarized experiment
seRNA <- SummarizedExperiment(
    assays = SimpleList(counts = matRNA)
  )
rownames(seRNA) <- rownames(scMat)
gtf <- getGeneGTF(gtf_file)
gtfMatch <- gtf[match(rownames(seRNA),gtf$gene_name)]
names(gtfMatch) <- rownames(seRNA)
rowRanges(seRNA) <- gtfMatch

#Use ATAC cluster info
colData(seRNA) <- DataFrame(row.names = colnames(seRNA), 
  clustATAC = KNN_Highest_Cluster, 
  groupATAC = KNN_Highest_Experiment
  )

saveRDS(seRNA, "results/Save-scRNA-Merged-KNN-SVD.rds")







