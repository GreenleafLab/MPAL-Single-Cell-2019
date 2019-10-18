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

#Optimized LSI for scRNA-seq analysis
optimizeLSI <- function(mat, scaleTo = 10000, priorCount = 3, pcsUse = 1:25, 
    resolution = c(0.2, 0.4, 0.8), varFeatures = c(2500, 2500, 2500), seed = 1){

    set.seed(seed)
    stopifnot(length(resolution) > 1)

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
se <- readRDS("data/Supplementary_Data_Hematopoiesis/scRNA-Healthy-Hematopoiesis-190429.rds")

####################################################
#For Clustering Analysis Start Here
####################################################
nPCs <- 1:25 #Number of PCs for clustering
nTop <- 3000 #Choose a higher number of variable peaks
resolution <- c(0.2,0.6,1.0) #Clustering resolutions for Seurat SNN

#Optimize LSI Features
lsiObj <- optimizeLSI(assay(se), 
  resolution = resolution, 
  pcsUse = nPCs,
  varFeatures = nTop)

metadata(se)$optimizeLSI <- lsiObj
metadata(se)$matSVD <- lsiObj[[length(lsiObj)]][[1]][[1]] #Last one
metadata(se)$variableGenes <- rownames(se)[lsi[[length(lsi)]]$varFeatures] #Variable genes

####################################################
#For Creating UMAP Start Here
####################################################
matSVD <- metadata(se)$matSVD
clusters <- colData(se)$Clusters

#Set Seed and perform UMAP on LSI-SVD Matrix
set.seed(1)
uwotUmap <- uwot::umap(
    matSVD, 
    n_neighbors = 35, 
    min_dist = 0.45, 
    metric = "euclidean", 
    n_threads = 1, 
    verbose = TRUE, 
    ret_nn = TRUE,
    ret_model = TRUE
    )

pdf("Plot_UMAP-NN-35-MD-45.pdf", width = 12, height = 12, useDingbats = FALSE)
df <- data.frame(
    x = uwotUmap[[1]][,1],
    y = -uwotUmap[[1]][,2], 
    color = clusters
    )
ggplot(df,aes(x,y,color=color)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(values=metadata(se)$colorMap$Clusters) +
    xlab("UMAP Dimension 1") + 
    ylab("UMAP Dimension 2")
dev.off()

#Add UMAP coordinates to column data in summarized experiment
colData(se)$UMAP1 <- uwotUmap[[1]][,1]
colData(se)$UMAP2 <- uwotUmap[[1]][,2]

#Save Summarized Experiment
#Add UMAP Params
metadata(se)$UMAP_Params <- list(NN = 35, MD = 0.45, PCs = 1:25, VarGenes = 3000, Res = "2.6.10")
saveRDS(se, "results/scRNA-Healthy-Hematopoiesis.rds")

#Save UMAP embedding
save_uwot(uwotUmap, "results/scRNA-Hematopoiesis-UMAP-model.uwot")

#If the above code does not work because tarring doesnt work for some reason on Stanford's compute server
#The following code will do a similar job assumming system commands work
#Adapted from save_uwot
model <- uwotUmap
file <- "results/scRNA-Hematopoiesis-UMAP-model.uwot.tar"
mod_dir <- tempfile(pattern = "dir")
dir.create(mod_dir)
uwot_dir <- file.path(mod_dir, "uwot")
dir.create(uwot_dir)
model_tmpfname <- file.path(uwot_dir, "model")
saveRDS(model, file = model_tmpfname)
metrics <- names(model$metric)
n_metrics <- length(metrics)
for (i in 1:n_metrics) {
    nn_tmpfname <- file.path(uwot_dir, paste0("nn", i))
    if (n_metrics == 1) {
        model$nn_index$save(nn_tmpfname)
        model$nn_index$unload()
        model$nn_index$load(nn_tmpfname)
    }
    else {
        model$nn_index[[i]]$save(nn_tmpfname)
        model$nn_index[[i]]$unload()
        model$nn_index[[i]]$load(nn_tmpfname)
    }
}
setwd(mod_dir)
system(sprintf("tar -cvf %s%s %s", wd, file, "uwot/*"))
setwd(wd)



