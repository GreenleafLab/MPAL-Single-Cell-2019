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

Rcpp::sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }

    if(max(idxX)-1 > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }

    if(max(idxY)-1 > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }

    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);
    
    const int nx = X.ncol();
    const int ny = Y.ncol();

    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }

    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }

    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }

    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }

    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);
    } 

    return(cor);

  }'
)

nSample <- function(x, n, type="v"){
  if(type=="v"){
    if(length(x) > n){
      s <- x[sample(seq_along(x),n)]
    }else{
      s <- x
    }
  }else if(type=="c"){
    if(ncol(x) > n){
      s <- x[,sample(seq_len(ncol(x)),n),drop=F]
    }else{
      s <- x
    }
  }else if(type=="r"){
    if(nrow(x) > n){
      s <- x[sample(seq_len(nrow(x)),n),,drop=F]
    }else{
      s <- x
    }
  }else{
    stop(paste0("type ",type," unrecognized..."))
  }
  return(s)
}

getNullCorrelations <- function(seA, seB, o, n){

  set.seed(1)
  o$seq <- seqnames(seA)[o$A]

  nullCor <- lapply(seq_along(unique(o$seq)), function(i){

    #Get chr from olist
    chri <- unique(o$seq)[i]
    cat(paste0(chri), "\n")

    #Randomly get n seA
    transAidx <- nSample( which(as.character(seqnames(seA)) != chri), n, "v")

    #Calculate Correlations
    grid <- expand.grid(transAidx, unique(o[o$seq==chri,]$B))

    idxA <- unique(grid[,1])
    idxB <- unique(grid[,2])

    seSubA <- seA[idxA]
    seSubB <- seB[idxB]

    grid[,3] <- match(grid[,1], idxA)
    grid[,4] <- match(grid[,2], idxB)

    colnames(grid) <- c("A", "B")
    out <- rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
    out <- na.omit(out)

    return(out)

  }) %>% SimpleList

  summaryDF <- lapply(nullCor, function(x){
    data.frame(mean = mean(x), sd = sd(x), median = median(x), n = length(x))
  }) %>% Reduce("rbind",.)

  return(list(summaryDF, unlist(nullCor)))

}

getQuantiles <- function(v, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}


####################################################
#Input Data
####################################################

#Input Files
scATAC_file <- ""
scRNA_file <- ""
gtf_file <- ""

#Params
fixA <- "center"
fixB <- "start"
associationWindow <- 2 * 250*10^3 + 1 #+-250 Kb
corCutOff <- 0.35 #Pearson Correlation Cutoff
fdrCutOff <- 0.1 #FDR Cutoff
distCutOff <- 2500 #Min Dist to TSS

#Input Summarized Experiments Log2 Normalize
seA <- readRDS(scATAC_file) #Aggregated scATAC Summarized Experiment
seB <- readRDS(scRNA_file) #Aggregated scRNA Summarized Experiment
assay(seA) <- log2(edgeR::cpm(assay(seA))/100+1)
assay(seB) <- log2(edgeR::cpm(assay(seB))/100+1)

#Resize B to association Window
seBWindow <- resize(rowRanges(seB), width = 1, fix = fixB) %>%
  {suppressWarnings(resize(., width = associationWindow, fix = "center"))} %>% trim(.)

#Keep only seA within association window
seA <- seA[unique(queryHits(findOverlaps(resize(seA,1,fixA), seBWindow, ignore.strand = TRUE)))]

#Getting distances
message("Getting Distances...")
o <- findOverlaps(seBWindow, resize(rowRanges(seA),1,fixA), ignore.strand = TRUE)

#Get Distance from Fixed point A B correct for minus stranded
mcols(o)$distance <- start(resize(rowRanges(seA),1,fixA))[subjectHits(o)] - start(resize(rowRanges(seB),1,fixB))[queryHits(o)]
mcols(o)$distance[which(as.character(strand(rowRanges(seB)))[queryHits(o)]=="-")] <- -1*mcols(o)$distance[which(as.character(strand(rowRanges(seB)))[queryHits(o)]=="-")]

#Add other info
o <- DataFrame(o)
colnames(o) <- c("B","A","distance")
o <- o[,c("A","B","distance")]

#If you want Peak-To-Gene links from MPAL SEE Below because we split into cells from MPALs and Healthy!
#If you want to use this code for your data just continue
nullCor <- getNullCorrelations(seA, seB, o, 1000)
o$Correlation <- rowCorCpp(as.integer(o[,1]), as.integer(o[,2]), assay(seA), assay(seB))
o$VarAssayA <- getQuantiles(matrixStats::rowVars(assay(seA)))[o$A]
o$VarAssayB <- getQuantiles(matrixStats::rowVars(assay(seB)))[o$B]
o$Pval <- 2*pnorm(-abs(((o$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
o$FDR <- p.adjust(o$Pval, method = "fdr")

#Get GTF
gtf <- getGeneGTF(gtf_file)
tssRNA <- resize(gtf, 1, "start")
strand(tssRNA) <- "*"
peakLinks <- rowRanges(seA)[o[,1]]
geneLinks <- rowRanges(seB) %>% resize(1, "start") %>% {.[o[,2]]}
mcolsLinks <- data.frame(geneLinks)[,c("seqnames","start","strand","gene_name","gene_id","exonLength")]
colnames(mcolsLinks) <- c("gene_chr","gene_start","gene_strand","gene_name","gene_id","exonLength")
mcolsLinks <- cbind(mcolsLinks, data.frame(o))
mcolsLinks$nearestGene <- tssRNA$gene_name[subjectHits(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))]
mcolsLinks$nearestGeneDistance <- mcols(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))$distance
mcols(peakLinks) <- mcolsLinks
peakLinks$peakName <- paste(seqnames(peakLinks), start(peakLinks), end(peakLinks), sep = "_")
peakLinks$sigCorrelation <- peakLinks$Correlation >= corCutOff & peakLinks$FDRCancer <= fdrCutOff & abs(peakLinks$distance) >= distCutOff
linksSig <- peakLinks[which(peakLinks$sigCorrelation)]

outMatch <- list(
  seA = seA[unique(mcols(linksSig)$peakName),], 
  seB = seB[unique(mcols(linksSig)$gene_name),], 
  linksSig = linksSig,
  linksAll = peakLinks
  )

saveRDS(outMatch, "results/Save-P2G-Links.rds")

###################################################
# MPAL Specific P2G Links
###################################################

#MPAL split into healthy and cancer cells
idxCancer <- which(grepl("MPAL|RM",colData(seA)$groupATAC)) #MPAL = MPAL Cells, RM = MPAL Cells from Ravi Majeti
idxHealthy <- which(!grepl("MPAL|RM",colData(seA)$groupATAC)) #MPAL = MPAL Cells, RM = MPAL Cells from Ravi Majeti

nullCorHealthy <- getNullCorrelations(seA[,idxHealthy], seB[,idxHealthy], o, 1000)
gc()
nullCorCancer <- getNullCorrelations(seA[,idxCancer], seB[,idxCancer], o, 1000)
gc()

o$CorrelationCancer <- rowCorCpp(as.integer(o[,1]), as.integer(o[,2]), assay(seA[,idxCancer]), assay(seB[,idxCancer]))
gc()
o$CorrelationHealthy <- rowCorCpp(as.integer(o[,1]), as.integer(o[,2]), assay(seA[,idxHealthy]), assay(seB[,idxHealthy]))
gc()

o$VarCancerA <- getQuantiles(matrixStats::rowVars(assay(seA[,idxCancer])))[o$A]
o$VarCancerB <- getQuantiles(matrixStats::rowVars(assay(seB[,idxCancer])))[o$B]

o$VarHealthyA <- getQuantiles(matrixStats::rowVars(assay(seA[,idxHealthy])))[o$A]
o$VarHealthyB <- getQuantiles(matrixStats::rowVars(assay(seB[,idxHealthy])))[o$B]

o$PvalCancer <- zToPval((o$CorrelationCancer - mean(nullCorCancer[[2]])) / sd(nullCorCancer[[2]]))
o$PvalHealthy <- zToPval((o$CorrelationHealthy - mean(nullCorHealthy[[2]])) / sd(nullCorHealthy[[2]]))

o$FDRCancer <- p.adjust(o$PvalCancer, method = "fdr")
o$FDRHealthy <- p.adjust(o$PvalHealthy, method = "fdr")

#Get GTF
gtf <- getGeneGTF(gtf_file)
tssRNA <- resize(gtf, 1, "start")
strand(tssRNA) <- "*"
peakLinks <- rowRanges(seA)[o[,1]]
geneLinks <- rowRanges(seB) %>% resize(1, "start") %>% {.[o[,2]]}
mcolsLinks <- data.frame(geneLinks)[,c("seqnames","start","strand","gene_name","gene_id","exonLength")]
colnames(mcolsLinks) <- c("gene_chr","gene_start","gene_strand","gene_name","gene_id","exonLength")
mcolsLinks <- cbind(mcolsLinks, data.frame(o))
mcolsLinks$nearestGene <- tssRNA$gene_name[subjectHits(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))]
mcolsLinks$nearestGeneDistance <- mcols(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))$distance
mcols(peakLinks) <- mcolsLinks
peakLinks$peakName <- paste(seqnames(peakLinks), start(peakLinks), end(peakLinks), sep = "_")
peakLinks$sigCancer <- peakLinks$CorrelationCancer >= corCutOff & peakLinks$FDRCancer <= fdrCutOff & abs(peakLinks$distance) >= distCutOff
peakLinks$sigHealthy <- peakLinks$CorrelationHealthy >= corCutOff & peakLinks$FDRHealthy <= fdrCutOff & abs(peakLinks$distance) >= distCutOff

linksSig <- peakLinks[which(peakLinks$sigCancer | peakLinks$sigHealthy)]
linksCancer <- peakLinks[which(peakLinks$sigCancer & !peakLinks$sigHealthy)]
linksHealthy <- peakLinks[which(!peakLinks$sigCancer & peakLinks$sigHealthy)]
linksShared <- peakLinks[which(peakLinks$sigCancer & peakLinks$sigHealthy)]

outMatch <- list(
  seA = seA[unique(mcols(linksSig)$peakName),], 
  seB = seB[unique(mcols(linksSig)$gene_name),], 
  linksCancer = linksCancer,
  linksHealthy = linksHealthy,
  linksShared = linksShared,
  linksSig = linksSig,
  linksAll = peakLinks
  )

saveRDS(outMatch, "Save_MPAL_P2G_Links.rds")

