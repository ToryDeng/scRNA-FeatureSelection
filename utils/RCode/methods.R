
run_cluster_methods <- function(method,data){
  switch(method,
         sc3 = cluster.sc3(data),
         seurat = cluster.seurat(data),
         stop("No such cluster method")
  )
}


run_assign_methods <- function(method, train_data, test_data){
  switch(method,
         singleR = assign.singleR(train_data,test_data),
         scPred = assign.scPred(train_data,test_data),
         stop("No such assigning method")
  )
  
}

#### assigning using scPred
assign.scPred <- function(train_data, test_data){
  require(scPred)
  require(Seurat)
  require(magrittr)
  require(doParallel)
  cl <- makePSOCKcluster(detectCores(logical = F) - 1)
  registerDoParallel(cl)
  
  reference <- as.Seurat(train_data, data=NULL)
  query <- as.Seurat(test_data, data=NULL)
  
  reference <- NormalizeData(object = reference, normalization.method = "LogNormalize")
  reference <- ScaleData(object = reference)
  reference <- RunPCA(object = reference, features = rownames(reference))
  reference <- RunUMAP(reference, dims = 1:30)
  
  
  reference <- getFeatureSpace(reference, "label")
  reference <- trainModel(reference, allowParallel = TRUE)
  stopCluster(cl)
  
  query <- NormalizeData(query)
  query <- scPredict(query, reference)
  
  return(query@meta.data[["scpred_prediction"]])
}


#### assigning using singleR
assign.singleR <- function(train_data, test_data, exp_config=NULL){
  require(SingleR)
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  train_data <- scater::logNormCounts(train_data) 
  test_data <- scater::logNormCounts(test_data) 
  preds <- SingleR(test=test_data, ref=train_data, labels=train_data$label, de.n = ncol(train_data))  # disable feature selection
  return(preds$labels)
}


### clustering using Seurat
cluster.seurat <- function(data) {
  require(Seurat)
  require(future)
  require(doParallel)
  stopifnot(is(data,"SingleCellExperiment"))
  cnts <- counts(data)
  
  # Enable parallelization
  plan("multisession", workers = detectCores(logical = F) - 1)
  
  ## PCA dimension reduction
  seuset <- CreateSeuratObject(cnts, project='simple_accuracy')
  seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize")
  seuset <- ScaleData(object = seuset)
  seuset <- RunPCA(object = seuset, features = rownames(seuset))
  seuset <- FindNeighbors(object = seuset)
  seuset <- FindClusters(object = seuset)
  return(as.integer(unname(seuset$seurat_clusters)))
}


### SC3
cluster.sc3 <- function(data) {
  require(SC3)
  stopifnot(is(data,"SingleCellExperiment"))
  num_true_cluster <- dim(unique(colData(data)['label']))[1]
  counts(data) <- as.matrix(counts(data))
  data <- scater::logNormCounts(data)
  k <- num_true_cluster

  rowData(data)$feature_symbol <- rownames(data)
  data <- data[!duplicated(rowData(data)$feature_symbol),]
  data <- sc3_prepare(data,gene_filter = FALSE)
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





  