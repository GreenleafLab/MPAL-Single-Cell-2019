#Clustering and scATAC-seq UMAP for Hematopoiesis data
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

####################################################
#Functions
####################################################

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){

    set.seed(1)

    #TF IDF LSI adapted from flyATAC
    if(binarize){
        message(paste0("Binarizing matrix..."))
        mat@x[mat@x > 0] <- 1 
    }

    #Filter 0 Sum Peaks
    rowSm <- Matrix::rowSums(mat)
    if(!is.null(nFeatures)){
        message(paste0("Getting top ", nFeatures, " features..."))
        idx1 <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
        idx2 <- which(rowSm>0)
        idx <- intersect(idx1,idx2)
        mat <- mat[idx,,drop=FALSE] 
    }else{
        idx <- which(rowSm>0)
        mat <- mat[idx,,drop=FALSE]
    }

    #Filter 0 Sum Cells
    colSm <- Matrix::colSums(mat)
    if(length(which(colSm==0))>0){
        message("Filtering Cells with 0 ColSums...")
        mat <- mat[,which(colSm>0),drop=FALSE]
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

#Clustering function using seurat SNN (Seurat v2.3.4)
seuratSNN <- function(matSVD, dims.use = 1:50, ...){
  set.seed(1)
  message("Making Seurat Object...")
  mat <- matrix(rnorm(nrow(matSVD) * 3, 1000), ncol = nrow(matSVD), nrow = 3)
  colnames(mat) <- rownames(matSVD)
  obj <- Seurat::CreateSeuratObject(mat, project='scATAC', min.cells=0, min.genes=0)
  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "cell.embeddings", new.data = matSVD)
  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "key", new.data = "PC")
  obj <- Seurat::FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, print.output = TRUE, ...)
  clust <- obj@meta.data[,ncol(obj@meta.data)]
  paste0("Cluster",match(clust, unique(clust)))
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

####################################################
#Input Data
####################################################
#Read in Summarized Experiment
#Please Note Code here has been modified to work with finalized summarized experiment
se <- readRDS("data/Supplementary_Data_Hematopoiesis/scATAC-Healthy-Hematopoiesis-190429.rds")

####################################################
#For Clustering Analysis Start Here
####################################################
nPCs1 <- 1:50 #Number of PCs in first analysis using all peaks
nTop <- 50000 #Choose a higher number of variable peaks across clusters (25-50k) to mitigate batch effects
nPCs2 <- 1:50 #Number of PCs in second analysis using variable peaks across clusters
resolution <- 1.5 #Clustering resolution for Seurat SNN

#RUN LSI 1
message("Running LSI 1...")
mat <- assay(se)
lsi1 <- calcLSI(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL)
clust1 <- seuratSNN(lsi1[[1]], dims.use = nPCs1, resolution = resolution)

#Make Pseudo Bulk Library
message("Making PseudoBulk...")
mat <- mat[,rownames(lsi1[[1]]), drop = FALSE] #sometimes cells are filtered
mat@x[mat@x > 0] <- 1 #binarize
clusterSums <- groupSums(mat = mat, groups = clust1, sparse = TRUE) #Group Sums
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3) #log CPM matrix
varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), nTop) #Top variable peaks

#RUN LSI 2
message("Running LSI 2...")
lsi2 <- calcLSI(mat[varPeaks,,drop=FALSE], nComponents = 50, binarize = TRUE, nFeatures = NULL)
clust2 <- seuratSNN(lsi2[[1]], dims.use = nPCs2, resolution = resolution)

#Append Summarized Experiment
se <- se[,rownames(lsi2[[1]])]
colData(se)$Clusters <- clust2
metadata(se)$LSI <- lsi2
metadata(se)$LSIPeaks <- varPeaks
metadata(se)$matSVD <- lsi2$matSVD

####################################################
#For Creating UMAP Start Here
####################################################
matSVD <- metadata(se)$matSVD
clusters <- colData(se)$Clusters

#Set Seed and perform UMAP on LSI-SVD Matrix
set.seed(1)
uwotUmap <- uwot::umap(
    matSVD[,1:50], 
    n_neighbors = 55, 
    min_dist = 0.45, 
    metric = "euclidean", 
    n_threads = 1, 
    verbose = TRUE, 
    ret_nn = TRUE,
    ret_model = TRUE
    )

pdf("results/Plot_UMAP-NN-55-MD-45.pdf", width = 12, height = 12, useDingbats = FALSE)
df <- data.frame(
    x = uwotUmap[[1]][,1],
    y = uwotUmap[[1]][,2], 
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
metadata(se)$UMAP_Params <- list(NN = 55, MD = 0.45, PCs = 1:50, VarPeaks = 50000, Res = "1.5")
saveRDS(se, "results/scATAC-Healthy-Hematopoiesis.rds")

#Save UMAP embedding
save_uwot(uwotUmap, "results/scATAC-Hematopoiesis-UMAP-model.uwot")

#If the above code does not work because tarring doesnt work for some reason on Stanford's compute server
#The following code will do a similar job assumming system commands work
#Adapted from save_uwot
model <- uwotUmap
file <- "results/scATAC-Hematopoiesis-UMAP-model.uwot.tar"
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



