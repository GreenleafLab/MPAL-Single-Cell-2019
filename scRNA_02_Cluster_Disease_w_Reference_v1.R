#Clustering and scRNA-seq UMAP for Hematopoiesis data
#06/02/19
#Cite Granja*, Klemm*, Mcginnis* et al. 
#A single cell framework for multi-omic analysis of disease identifies 
#malignant regulatory signatures in mixed phenotype acute leukemia (2019)
#Created by Jeffrey Granja
library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(uwot)
library(edgeR)
library(FNN)
library(matrixStats)
library(Rcpp)
set.seed(1)

####################################################
#Functions
####################################################

#Binarize Sparse Matrix
binarizeMat <- function(mat){
    mat@x[mat@x > 0] <- 1
    mat
}

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){

    set.seed(1)

    #TF IDF LSI adapted from flyATAC
    if(binarize){
        message(paste0("Binarizing matrix..."))
        mat@x[mat@x > 0] <- 1 
    }

    if(!is.null(nFeatures)){
        message(paste0("Getting top ", nFeatures, " features..."))
        idx <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
        mat <- mat[idx,] 
    }else{
        idx <- which(Matrix::rowSums(mat) > 0)
        mat <- mat[idx,]
    }

    #Calc RowSums and ColSums
    colSm <- Matrix::colSums(mat)
    rowSm <- Matrix::rowSums(mat)

    #Calc TF IDF
    message("Computing Term Frequency IDF...")
    freqs <- t(t(mat)/colSm)
    idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs

    #Calc SVD then LSI
    message("Computing SVD using irlba...")
    svd <- irlba::irlba(tfidf, nComponents, nComponents)
    svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))

    #Return Object
    out <- list(
        matSVD = matSVD, 
        rowSm = rowSm, 
        colSm = colSm, 
        idx = idx, 
        svd = svd, 
        binarize = binarize, 
        nComponents = nComponents,
        date = Sys.Date(),
        seed = 1)

    out

}

#Sparse Variances Rcpp
sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // [[Rcpp::export]]
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
    const int nv = j.size();
    const int nm = rm.size();
    Rcpp::NumericVector rv(nm);
    Rcpp::NumericVector rit(nm);
    int current;
    // Calculate RowVars Initial
    for (int i = 0; i < nv; ++i) {
      current = j(i) - 1;
      rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
      rit(current) = rit(current) + 1;
    }
    // Calculate Remainder Variance
    for (int i = 0; i < nm; ++i) {
      rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
    }
    rv = rv / (n - 1);
    return(rv);
  }'
)

#Compute Fast Sparse Row Variances
sparseRowVariances <- function (m){
    rM <- Matrix::rowMeans(m)
    rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
    return(rV)
}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

#Seurat SNN
seuratSNN <- function(matSVD, dims.use = 1:50, print.output = TRUE, ...){
  set.seed(1)
  message("Making Seurat Object...")
  mat <- matrix(rnorm(nrow(matSVD) * 3, 1000), ncol = nrow(matSVD), nrow = 3)
  colnames(mat) <- rownames(matSVD)
  obj <- Seurat::CreateSeuratObject(mat, project='scATAC', min.cells=0, min.genes=0)
  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "cell.embeddings", new.data = matSVD)
  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "key", new.data = "PC")
  obj <- Seurat::FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, print.output = print.output, ...)
  clust <- obj@meta.data[,ncol(obj@meta.data)]
  paste0("Cluster",match(clust, unique(clust)))
}

#Optimized LSI for scRNA-seq analysis
optimizeLSI <- function(mat, scaleTo = 10000, priorCount = 3, pcsUse = 1:25, 
    resolution = c(0.2, 0.4, 0.8), varFeatures = c(2500, 2500, 2500), seed = 1){

    set.seed(seed)
    stopifnot(length(resolution) > 1)
    stopifnot(length(resolution) == length(varFeatures))

    #Initialize List
    lsiOut <- list()

    #Initial LSI uses variances that are across all single cells and will have larger batch relationships
    i <- 1
    message("Initial LSI...")
    matNorm <- t(t(mat)/Matrix::colSums(mat)) * scaleTo
    matNorm@x <- log2(matNorm@x + 1)
    idVarFeatures <- head(order(sparseRowVariances(matNorm),decreasing=TRUE), varFeatures[i])
    lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
    clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)

    #Store
    lsiOut[[paste0("iter", i)]] <- list(
        lsiMat = lsiObj$matSVD, 
        varFeatures = idVarFeatures, 
        clusters = clusters
        )

    for(i in seq(2, length(varFeatures))){

       message(sprintf("Additional LSI %s...", i))

        #Run LSI
        clusterMat <- edgeR::cpm(groupSums(mat, clusters, sparse = TRUE), log=TRUE, prior.count = priorCount)
        idVarFeatures <- head(order(rowVars(clusterMat), decreasing=TRUE), varFeatures[i])
        lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
        clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)

        if(i == length(varFeatures)){
            #Save All Information from LSI Attempt
            lsiOut[[paste0("iter", i)]] <- list(
                lsiObj = lsiObj, 
                varFeatures = idVarFeatures, 
                clusters = clusters,
                matNorm = matNorm
                )
        }else{
            lsiOut[[paste0("iter", i)]] <- list(
                lsiMat = lsiObj$matSVD, 
                varFeatures = idVarFeatures, 
                clusters = clusters
                )
        }

    }

    return(lsiOut)

}

####################################################
#Input Data
####################################################
#Read in Summarized Experiment
#Please Note Code here has been modified to work with finalized summarized experiment

#Reference Summarized Experiment
seReference <- readRDS("data/Supplementary_Data_All_Hematopoiesis_MPAL/scRNA-All-Hematopoiesis-MPAL-190429.rds")

#SE Disease Cells
id <- "MPAL1"
seDisease <- seReference[,grep(id, colData(seReference)$Group)]

#SE Healthy Cells
seReference <- seReference[,grep("BMMC|CD34|PBMC", colData(seReference)$Group)]

#Identify Gene Universe
gU <- intersect(rownames(seReference), rownames(seDisease))
gU <- gU[!grepl("^MT", gU)]

#Set Clustering Parameters
resolution <- c(0.2,0.8,0.8) #clustering resolution
varGenesToUse <- c(1000,1000,1000) #number of variable genes

#Optimize LSI Features
matAll <- cbind(assay(seReference[gU,]), assay(seDisease[gU,]))
lsiObj <- optimizeLSI(matAll, resolution = resolution, varFeatures = varGenesToUse)

#UMAP
set.seed(1)
umap <- uwot::umap(
    lsiObj[[length(lsiObj)]]$lsiObj$matSVD[,1:25], 
    n_neighbors = 30, 
    min_dist = 0.5, 
    metric = "euclidean", 
    n_threads = 5, 
    verbose = TRUE, 
    ret_model = FALSE
    )

#Plot Info
cells <- c(rep("reference", ncol(seReference)),rep("disease",ncol(seDisease)))
splitCells <- split(cells,lsiObj[[length(lsiObj)]]$clusters)
df <- data.frame(
    clusters = names(splitCells),
    proportion = unlist(lapply(seq_along(splitCells), function(x) sum(splitCells[[x]]=="disease") / length(splitCells[[x]])))
    )

#Plot UMAP Data Frame
plotDF <- data.frame(umap)
rownames(plotDF) <- c(colnames(seReference), colnames(seDisease))
plotDF$type <- cells
plotDF$clusters <- lsiObj[[length(lsiObj)]]$clusters
plotDF$classification <- 0
#If disease cells are clustered with healthy cluster (proportion > 0.8) we will classify these as healthy-like
plotDF$classification[plotDF$type == "disease" & plotDF$clusters %in% paste0(df$clusters[df[,2] > 0.8])] <- 1
plotDF$classification[plotDF$type == "disease"] <- plotDF$classification[plotDF$type == "disease"] + 1
plotDF <- plotDF[order(plotDF$classification), ]

#Formal Classification
plotDF$classificationSTR <- "reference"
plotDF$classificationSTR[plotDF$classification==1] <- "healthy-like"
plotDF$classificationSTR[plotDF$classification==2] <- "disease-like"

#Plot PDFs
plotDir <- paste0("results/scRNA/classification/")
dir.create(plotDir,recursive=TRUE)
pdf(paste0(plotDir,id,"-Classification-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)

ggplot(plotDF, aes(X1,X2,color=classificationSTR)) + 
    geom_point() +
    theme_bw() +
    xlab("UMAP Dimension 1") + 
    ylab("UMAP Dimension 2") +
    scale_color_manual(values=c("reference"="lightgrey","healthy-like"="dodgerblue3","disease-like"="firebrick3"))

dev.off()

####################################################
#Project Into LSI UMAP
####################################################

#Previous Reference Summarized Experiment
se <- readRDS("data/Supplementary_Data_Hematopoiesis/scRNA-Healthy-Hematopoiesis-190429.rds")

#Load Saved UMAP Manifold
umapManifold <- uwot::load_uwot("data/Supplementary_Data_LSI_Projection/scRNA-Hematopoiesis-UMAP-model.190505.uwot.tar")

#LSI Projection Matrix
lsiGenes <- metadata(se)$variableGenes
matProjectLSI <- assay(seDisease[lsiGenes,])

#LSI Project
lsiReference <- metadata(se)$optimizeLSI[[length(metadata(se)$optimizeLSI)]]$lsiObj
lsiProjection <- projectLSI(matProjectLSI, lsiReference)

#UMAP Projection
#Set Seed Prior to umap_transform (see uwot github)
set.seed(1)
umapProjection <- uwot::umap_transform(as.matrix(lsiProjection$matSVD)[,1:25], umapManifold, verbose = TRUE)

#Plot Projection
refDF <- data.frame(row.names = colnames(se), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "reference")
proDF <- data.frame(row.names = colnames(seDisease), X1 = umapProjection[,1], X2 = umapProjection[,2], Type = plotDF[colnames(seDisease),]$classificationSTR)
projectionDF <- rbind(refDF, proDF)

plotDir <- paste0("results/scRNA/classification/")
dir.create(plotDir,recursive=TRUE)
pdf(paste0(plotDir,id,"-Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)

ggplot(projectionDF, aes(X1,-X2,color=Type)) + 
    geom_point() +
    theme_bw() +
    xlab("UMAP Dimension 1") + 
    ylab("UMAP Dimension 2") +
    scale_color_manual(values=c("reference"="lightgrey","healthy-like"="dodgerblue3","disease-like"="firebrick3"))

dev.off()

####################################################
#Differential Analysis Into LSI UMAP
####################################################

#Input Parameters
input_knn <- 25
scaleTo <- 10000
nMax <- 500

#LSI-SVD
svdReference <- as.data.frame(lsiReference$matSVD)
svdDisease <- as.data.frame(as.matrix(lsiProjection$matSVD))

#Differential Seed
set.seed(1)

#Cells that we are testing of disease
idxDisease <- rownames(plotDF)[plotDF$classificationSTR=="disease-like"]

#If the number of cells exceeds the max downsample to max
if(length(idxDisease) > nMax){
    idxDisease <- sample(idxDisease, nMax)
}

#If the number of cells is greater than 5 continue
stopifnot(length(idxDisease) > 5)

#KNN Nearest Neighbor using FNN
knnDisease <- get.knnx(
    data = svdReference,
    query = svdDisease[idxDisease, ], #Subset by idxDisease 
    k = input_knn)

#Determine the minimum KNN where reference cells are less than 1.25x disease cells
i <- 0
uniqueIdx <- unique(as.vector(knnDisease$nn.index))
while(length(uniqueIdx) > 1.25 * length(idxDisease)){
    i <- i + 1
    uniqueIdx <- unique(as.vector(knnDisease$nn.index[,seq_len(input_knn-i)]))
}

#Reference cells for testing
idxReference <- rownames(svdReference)[uniqueIdx]

#If there are more healthy cells downsample healthy cells
#If there are more disease cells downasmple disease cells
if(length(idxReference) > length(idxDisease)){
    idxReference <- sample(idxReference, length(idxDisease))
}else{
    idxDisease <- sample(idxDisease, length(idxReference))
}
message(sprintf("nDisease = %s\nnHealthy = %s", length(idxDisease), length(idxReference)))

#Disease and Reference Matrix
matHealthy <- assay(seReference[,idxReference])
matDisease <- assay(seDisease[,idxDisease])

#Normalize to scaleTo
matNormDisease <- t(t(matDisease)/Matrix::colSums(matDisease)) * scaleTo
matNormHealthy <- t(t(matHealthy)/Matrix::colSums(matHealthy)) * scaleTo

#T-Test Comparisons
dfTT <- sparseMatTTest(matNormDisease, matNormHealthy)
dfTT$feature <- rownames(matNormDisease)
dfTT$log2Mean <- log2(rowMeans(cbind(dfTT$mean1, dfTT$mean2)) + 10^-4)
dfTT$log2FC <- log2((dfTT$mean1 + 10^-4)/(dfTT$mean2 + 10^-4))

plotDiff <- data.frame(row.names=row.names(dfTT),log2Mean=dfTT$log2Mean,log2FC=dfTT$log2FC,FDR=dfTT$fdr)
plotDiff <- plotDiff[complete.cases(plotDiff),]
plotDiff$type <- "not-differential"
plotDiff$type[plotDiff$log2FC > 0.5 & plotDiff$FDR < 0.01] <- "up-regulated"
plotDiff$type[plotDiff$log2FC < -0.5 & plotDiff$FDR < 0.01] <- "do-regulated"

plotDir <- paste0("results/scRNA/classification/")
dir.create(plotDir,recursive=TRUE)
pdf(paste0(plotDir,id,"-Differential-MA-Plot.pdf"), width = 8, height = 6, useDingbats = FALSE)

ggplot(plotDiff, aes(log2Mean,log2FC,color=type)) + 
    geom_point(size=0.5) +
    theme_bw() +
    xlab("log2 Mean") + 
    ylab("log2 Fold Change") +
    scale_color_manual(values=c("not-differential"="lightgrey", "do-regulated"="dodgerblue3", "up-regulated"="firebrick3"))

dev.off()

#Save Output
readr::write_tsv(dfTT, paste0(plotDir,id,"-Differential-Results.tsv"))











