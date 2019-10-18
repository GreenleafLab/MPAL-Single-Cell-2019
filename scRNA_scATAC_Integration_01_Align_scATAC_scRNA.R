#Clustering and scATAC-seq UMAP for Hematopoiesis data
#06/02/19
#Cite Granja*, Klemm*, Mcginnis* et al. 
#A single cell framework for multi-omic analysis of disease identifies 
#malignant regulatory signatures in mixed phenotype acute leukemia (2019)
#Created by Jeffrey Granja
library(Seurat)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(Rcpp)
set.seed(1)

####################################################
#Functions
####################################################

#Nearest Neighbor differential
findNN <- function(query, reference, method = "euclidean"){
    findClosest <- function(x, m, method = "euclidean"){
        if(method=="euclidean"){
            which.min(sqrt(colSums((t(m) - x) * (t(m) - x))))
        }else if(method=="pearson"){
            which.max(cor(t(m),x,method = method)[,1])
        }else if(method=="spearman"){
            which.max(cor(t(m),x,method = method)[,1])
        }
    }
    pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
    mat <- data.frame(matrix(ncol = 4, nrow = nrow(query)))
    colnames(mat) <- c("x", "i", "y", "j")
    for(i in seq_len(nrow(query))){
        setTxtProgressBar(pb,round(i*100/nrow(query),0))
      j <- findClosest(query[i,], reference, method)
      mat[i,] <- c(x = rownames(query)[i], i = i, y = rownames(reference)[j], j = j)
    }
    return(mat)
}

sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }

    if(max(idxX) > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }

    if(max(idxY) > Y.nrow()){
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

####################################################
#Input Data
####################################################

#Read in Summarized Experiment
#Please Note Code here has been modified to work with finalized summarized experiment

#Prep RNA Matrix from Summarized Experiment
se <- readRDS(opt$input_RNA)
matRNA <- assay(se)

#Prep Gene Score Matrix from Summarized Experiment
seGS <- readRDS(opt$input_GS)
matGS <- assay(seGS)

#Parameters
nCCA <- 20
nVarGenes <- 2500
selectMethod <- "all"

#Gene Universe
geneUniverse <- intersect(rownames(matGS),rownames(matRNA))

#Remove Mito RNA
geneUniverse <- geneUniverse[geneUniverse %ni% grep("^MT", c(rownames(seGS),rownames(se)), value = TRUE)]

#Subset By Gene Universe
matRNA <- matRNA[geneUniverse, ,drop = FALSE]
matGS <- matGS[geneUniverse, ,drop = FALSE]

#Create RNA Seurat
objRNA <- CreateSeuratObject(raw.data = matRNA, project = "RNA")
objRNA <- NormalizeData(object = objRNA)
objRNA <- ScaleData(object = objRNA)
objRNA <- FindVariableGenes(object = objRNA, do.plot = FALSE, selection.method = "dispersion", top.genes = as.integer(nVarGenes))
objRNA@meta.data[, "protocol"] <- "RNA"

#Create GS Seurat
objGS <- CreateSeuratObject(raw.data = matGS, project = "ATAC")
objGS <- NormalizeData(object = objGS)
objGS <- ScaleData(object = objGS)
objGS <- FindVariableGenes(object = objGS, do.plot = FALSE, selection.method = "dispersion", top.genes = as.integer(nVarGenes))
objGS@meta.data[, "protocol"] <- "ATAC"

#Intersect Variable Genes
if(tolower(selectMethod) == "genescores"){
  varGenes <- objGS@var.genes
}else if(tolower(selectMethod) == "rna"){
  varGenes <- objRNA@var.genes
}else if(tolower(selectMethod) == "intersect"){
  varGenes <- intersect(objRNA@var.genes, objGS@var.genes)
}else if(tolower(selectMethod) == "all"){
  varGenes <- unique(c(objRNA@var.genes, objGS@var.genes))
}

#Run CCA Seurat v2.3.4
CCA <- RunCCA(object = objRNA, object2 = objGS, genes.use = varGenes, num.cc = as.integer(nCCA))

#Variance Expectation Ration Seurat v2.3.4
CCA <- CalcVarExpRatio(object = CCA, reduction.type = "pca", grouping.var = "protocol", dims.use = seq_len(as.integer(nCCA)))

#Filter Seurat v2.3.4
CCA <- SubsetData(object = CCA, subset.name = "var.ratio.pca", accept.low = 0.5)

#Align Subspace Seurat v2.3.4
CCA <- AlignSubspace(object = CCA, reduction.type = "cca", grouping.var = "protocol", dims.align = seq_len(as.integer(nCCA)))

saveRDS(CCA, "results/Save-CCA-Alignment-scATAC-scRNA.rds")

#Get CCA Matrix
alignedCCA <- GetCellEmbeddings(CCA, reduction.type = "cca.aligned")

#KNN Search
#Alternatively for speed FNN::getknnx(query, reference, k = 1)
#We just used a simple function
matchedCells <- findNN(
  query = alignedCCA[CCA@meta.data$protocol=="ATAC",],
  reference = alignedCCA[CCA@meta.data$protocol=="RNA",], 
  method = "euclidean")

matchedCells$corCCA <- rowCorCpp(
  match(matchedCells$x, colnames(CCA@data)), 
  match(matchedCells$y, colnames(CCA@data)), 
  alignedCCA, alignedCCA)

matchedCells$corVarRNA <- rowCorCpp(
  match(matchedCells$x, colnames(CCA@data)), 
  match(matchedCells$y, colnames(CCA@data)), 
  t(as.matrix(CCA@data[CCA@var.genes,])), 
  t(as.matrix(CCA@data[CCA@var.genes,])))

matchx <- match(matchedCells$x, colnames(CCA@data))
matchy <- match(matchedCells$y, colnames(CCA@data))
mat <- as.matrix(CCA@data[CCA@var.genes,])

#-------------------------------------------------------
#UMAP
#-------------------------------------------------------
set.seed(1)
umap <- uwot::umap(
    alignedCCA, 
    n_neighbors = 50, 
    min_dist = 0.5, 
    metric = "euclidean", 
    n_threads = 5, 
    verbose = TRUE, 
    ret_model = FALSE)

#Plot DF
plotDF <- data.frame(umap)
rownames(plotDF) <- rownames(alignedCCA)
plotDF[rownames(CCA@meta.data[rownames(plotDF),]),"protocol"] <- CCA@meta.data[rownames(plotDF),]$protocol
plotDF <- plotDF[sample(seq_len(nrow(plotDF)), nrow(plotDF)),, drop = FALSE]

saveRDS(list(plotDF = plotDF, matchedCells = matchedCells), "results/Save-CCA-KNN-UMAP.rds")

