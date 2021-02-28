source("/home/tdeng/SingleCell/FeatureSelection/R/R code/clustering methods/methods_config.R")

run_cluster_methods <- function(method,data){
  switch(method,
         sc3 = cluster.sc3(data),
         liger = cluster.liger(data),
         seurat = cluster.seurat(data),
         cidr = cluster.cidr(data),
         tscan = cluster.tscan(data),
         stop("No such cluster method")
  )
}


### clustering using Seurat
cluster.seurat <- function(data) {
  require(Seurat)
  stopifnot(is(data,"SingleCellExperiment"))
  cnts <- counts(data)
  m_config <- methods.config.seurat
  ## PCA dimention redcution
  seuset <- CreateSeuratObject(cnts, project='simple_accuracy')
  seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize")
  #seuset <- FindVariableFeatures(object = seuset,nfeatures = m_config[['nfeatures']])
  seuset <- ScaleData(object = seuset)
  seuset <- RunPCA(object = seuset, features = rownames(seuset))
  seuset <- FindNeighbors(object = seuset,dims = 1:m_config[['pc_dims']])
  seuset <- FindClusters(object = seuset, resolution = m_config[['resolution']])
  return(as.integer(unname(seuset$seurat_clusters)))
}

### TSCAN
cluster.tscan <- function(data) {
  require(TSCAN)
  stopifnot(is(data,"SingleCellExperiment"))
  num_true_cluster <- length(unique(data@colData@listData[["label"]]))
  cnts <- counts(data)
  m_config <- methods.config.tscan
  procdata <- preprocess(as.matrix(cnts),cvcutoff=m_config[['cvcutoff']])
  if(!purrr::is_null(m_config[['k']]) )
    lpsmclust <- exprmclust(procdata, clusternum=num_true_cluster)
  else
    lpsmclust <- exprmclust(procdata)
  return(as.integer(unname(lpsmclust$clusterid)))
}

### SC3
cluster.sc3 <- function(data) {
  require(SC3)
  stopifnot(is(data,"SingleCellExperiment"))
  num_true_cluster <- length(unique(data@colData@listData[["label"]]))
  counts(data) <- as.matrix(counts(data))
  data <- scater::logNormCounts(data)
  m_config <- methods.config.sc3
  k <- num_true_cluster
  gene_filter <- m_config$gene_filter
  print(stringr::str_glue("gene filter for sc3 is {gene_filter}"))
  rowData(data)$feature_symbol <- rownames(data)
  data <- data[!duplicated(rowData(data)$feature_symbol),]
  data <- sc3_prepare(data,gene_filter = gene_filter)
  if(purrr::is_null(k)){
    data <- sc3_estimate_k(data)## estimate number of clusters
    k <- metadata(sce)$sc3$k_estimation
  }
  data <- sc3_calc_dists(data)
  data <- sc3_calc_transfs(data)
  data <- sc3_kmeans(data, ks = k)
  data <- sc3_calc_consens(data)
  col_data <- colData(data)
  colTb = as.vector(col_data[ , grep("sc3_", colnames(col_data))])
  return(as.integer(colTb))
}

#####liger
cluster.liger <- function(data){
  require(liger)
  m_config <- methods.config.liger
  stopifnot(is(data,"SingleCellExperiment"))
  data_list <- list(counts(data))
  names(data_list) <- metadata(data)$study[[1]]
  liger_data <- createLiger(data_list)
  liger_data <- liger::normalize(liger_data)
  liger_data <- selectGenes(liger_data, var.thresh = c(0.3, 0.875), do.plot = F)
  liger_data <- scaleNotCenter(liger_data)
  if(m_config$suggestK==TRUE){
    print("suggesting K")
    k.suggest <- suggestK(liger_data, num.cores = 5, gen.new = T, plot.log2 = F,nrep = 5)
  }else{
    k.suggest <- if(purrr::is_null(m_config$k.suggest))  25 else m_config$k.suggest
  }
  thresh <- if(purrr::is_null(m_config$thresh)) 5e-5 else m_config$thresh
  lambda <- if(purrr::is_null(m_config$lambda)) 5 else m_config$lambda
  resolution <- if(purrr::is_null(m_config$resolution)) 1.0 else m_config$resolution
  liger_data <- optimizeALS(liger_data, k=k.suggest, thresh = thresh, lambda=lambda,nrep = 3)
  liger_data <- quantileAlignSNF(liger_data,resolution = resolution) #SNF clustering and quantile alignment
  return(as.integer(unname(liger_data@clusters)))
}







### Monocle
cluster.monocle <- function(counts, K) {
  require(monocle)
  dat <- newCellDataSet(counts)
  dat <- estimateSizeFactors(dat)
  dat <- estimateDispersions(dat)
  
  ## pick marker genes for cell clustering
  disp_table <- dispersionTable(dat)
  ix <- order(disp_table[,"mean_expression"], decreasing=TRUE)
  unsup_clustering_genes <- disp_table[ix[50:1000], "gene_id"]
  ## unsup_clustering_genes <- subset(disp_table, mean_expression >= 1)
  dat <- setOrderingFilter(dat, unsup_clustering_genes)
  ## the follwoing step can be slow. Need to keep marker genes number low.
  dat <- reduceDimension(dat, reduction_method="tSNE")
  
  ## clustering
  if( !missing(K) )
    dat <- clusterCells(dat, num_clusters = K)
  else
    dat <- clusterCells(dat)
  
  pData(dat)$Cluster
}

### CIDR
cluster.cidr <- function(counts, K) {
  require(cidr)
  sData <- scDataConstructor(counts)
  sData <- determineDropoutCandidates(sData)
  sData <- wThreshold(sData)
  sData <- scDissim(sData)
  sData <- scPCA(sData, plotPC=FALSE)
  sData <- nPC(sData)
  if(missing(K))
    sData <- scCluster(sData)
  else
    sData <- scCluster(sData, nCluster=K)
  
  return(sData@clusters)
}


