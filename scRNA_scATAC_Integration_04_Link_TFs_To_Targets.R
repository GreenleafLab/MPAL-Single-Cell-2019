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

featureName <- function(gr){
  paste(seqnames(gr),start(gr),end(gr),sep="_")
}

####################################################
#Input Data
####################################################

diff_scATAC_file <- "Save-scATAC-PromoterNorm-Differentials.rds"
diff_scRNA_file <- "output/scATAC_scRNA/Differential/Differential_Summary.rds"
p2gLinks_file <- "../../../Integration/heme_mpal/Integration/Save_MPAL_P2G_Links-190511.rds"
motif_matches_file <- "../../../scATAC/post/heme_mpal/output/Motif_ArchRMatches.rds"

#Read Inputs
diffObj <- readRDS(diff_scRNA_file)
p2gLinks <- readRDS(p2gLinks_file)$linksSig
matches <- readRDS(motif_matches_file)
diffATAC <- readRDS(diff_scATAC_file)
rownames(matches) <- paste(seqnames(matches),start(matches),end(matches), sep = "_")
rownames(diffATAC) <- paste(seqnames(diffATAC),start(diffATAC),end(diffATAC), sep = "_")

#Make P2G Mats
names(p2gLinks) <- paste0("l",seq_along(p2gLinks))
dAp2g <- diffATAC[featureName(p2gLinks),colnames(diffObj$diffRNA)]
dRp2g <- diffObj$diffRNA[mcols(p2gLinks)$gene_name,] 
rownames(dAp2g) <- names(p2gLinks)
rownames(dRp2g) <- names(p2gLinks)

#Identify Significant Peaks and Genes per MPAL Subpopulation
sigMatP2G <- assays(dAp2g)[[1]] >= 0.5 & assays(dAp2g)[[2]] <= 0.05 & assays(dRp2g)[[1]] >= 0.5 & assays(dRp2g)[[2]] <= 0.01

#Which have at least 1 MPAL subpopulation with a linked (within p2glinks) diff peak and diff gene
sigP2G <- p2gLinks[which(rowSums(sigMatP2G) > 0)]

#List of TFs to identify Targerts
tfs <-  c("RUNX1_1182")

#Identify Targets
t2gDF <- lapply(seq_along(tfs), function(x){

  message(x)  
  #Subset Linked Peaks that contain TF Motif
  peaksx <- names(which(assay(matches[,tfs[x]])[,1]))

  #Subset sigMat by linked peaks with TF Motif
  sigMatP2Gx <- sigMatP2G[rownames(sigMatP2G) %in% peaksx,]

  #Subset links by linked peaks with TF Motif
  linksx <- sigP2G[featureName(sigP2G) %in% peaksx]

  #Figure out which samples are true for each link
  sigMatP2Gx <- sigMatP2G[names(linksx),]

  #Compute sum of MPAL subpopulations that are diff peak and diff peak for every peak connected to the gene
  sigDFListx <- split(data.frame(sigMatP2Gx), mcols(linksx)$gene_name) %>% lapply(., colSums)

  #Convert to Data Frame
  sigDFx <- data.frame(Reduce("rbind",sigDFListx))
  rownames(sigDFx) <- names(sigDFListx)

  #Set to maximum of 1 for each MPAL subpop for all diff peak gene combo
  sigDFx[sigDFx > 1] <- 1

  #Compute number of positive MPAL subpopulations
  nSubpop <- rowSums(sigDFx)
  
  #Summed R2 linkage for each linked positive diff peak gene combo (Max of Cor Healthy Cor Cancer then Squared for each Gene)
  maxCor <- pmax(mcols(linksx)$CorrelationHealthy,mcols(linksx)$CorrelationCancer,na.rm=TRUE)
  linkageScore <- split(maxCor^2, mcols(linksx)$gene_name) %>% 
    lapply(., sum) %>% 
    unlist(use.names=TRUE)

  #Return Summary Data Frame
  data.frame(TF = tfs[x], Gene = names(nSubpop), N = nSubpop, P = nSubpop / ncol(sigMatP2Gx), linkageScore = linkageScore[names(nSubpop)])

}) %>% Reduce("rbind",.)

#Time to Plot

#Plot TF idx
idx <- "RUNX1_1182"

#Subset TF idx
df <- t2gDF[t2gDF[,1]==idx,]

#Filter Genes with 2 or Less Subpop Diff
df <- df[df$N > 2,]

#Compute Zscore for N Subpop
df$ZN <- scale(df$N)

#Compute Zscore for Linkage Score
df$ZLS <- scale(df$linkageScore)

#Compute Zscore sum for N Subpop and Linkage Score
df$SummedZ <- df$ZN + df$ZLS

#GGPlot2
pdf("results/RUNX1-Targets.pdf")
ggplot(df, aes(N, linkageScore)) + geom_point() + theme_bw() +
  ggrepel::geom_label_repel(
      data = df[head(order(df$SummedZ,decreasing=TRUE),25),], size = 5,
      aes(x=P,y=linkageScore,color=NULL,label=Gene)
    )
dev.off()

saveRDS(df, "results/Save-RUNX1-Targets.rds")


