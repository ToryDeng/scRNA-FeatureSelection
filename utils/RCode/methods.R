source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/methods_config.R")

run_cluster_methods <- function(method,data){
  switch(method,
         sc3 = cluster.sc3(data),
         liger = cluster.liger(data),
         seurat = cluster.seurat(data),
         cidr = cluster.cidr(data),
         tscan = cluster.tscan(data),
         scmap = cluster.scmap(data),
         stop("No such cluster method")
  )
}


run_assign_methods <- function(method,train_data, test_data){
  switch(method,
         scmap_cluster = assign.scmap_cluster(train_data, test_data),
         scmap_cell = assign.scmap_cell(train_data, test_data),
         singlecellnet = assign.singlecellnet(train_data,test_data),
         stop("No such assigning method")
  )
  
}


####assigning using scmap cluster
assign.scmap_cluster <- function(train_data, test_data){
  require(scmap)
  require(scater)
  require(dplyr)
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  counts(train_data) <- as.matrix(counts(train_data))
  counts(test_data) <- as.matrix(counts(test_data))
  train_data <- logNormCounts(train_data)
  rowData(train_data)$feature_symbol <- rownames(train_data)
  train_data <- selectFeatures(train_data,n_features=methods.config.scmap[['nfeatures']]) %>%
    indexCluster(cluster_col = "label")
  test_data <- logNormCounts(test_data)
  rowData(test_data)$feature_symbol <- rownames(test_data)
  results <- scmapCluster(projection=test_data, index_list=list(metadata(train_data)$scmap_cluster_index))
  return(results$combined_labs)
}

####assigning using scmap cell
assign.scmap_cell <- function(train_data, test_data){
  require(scmap)
  require(scater)
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  counts(train_data) <- as.matrix(counts(train_data))
  counts(test_data) <- as.matrix(counts(test_data))
  train_data <- logNormCounts(train_data)
  rowData(train_data)$feature_symbol <- rownames(train_data)
  if(purrr::is_null(methods.config.scmap[['seed']]))
    set.seed(1)
  else
    set.seed(methods.config.scmap[['seed']])
  train_data <- selectFeatures(train_data,n_features=methods.config.scmap[['nfeatures']]) %>%
    indexCell()
  test_data <- logNormCounts(test_data)
  rowData(test_data)$feature_symbol <- rownames(test_data)
  cell_results <- scmapCell(projection=test_data, index_list = list(metadata(train_data)$scmap_cell_index))
  clusters_results <- scmapCell2Cluster(
    cell_results, 
    list(
      as.character(colData(train_data)$label)
    ),threshold = methods.config.scmap[['threshold']]
  )
  return(clusters_results$combined_labs)
}

######assigning using singlecellnet
assign.singlecellnet <- function(train_data, test_data){
  require(singleCellNet)
  stopifnot(is(train_data,"SingleCellExperiment"))
  stopifnot(is(test_data,"SingleCellExperiment"))
  #' extract sampTab and expDat sce object into regular S3 objects
  #' @param sce_object
  #' @param exp_type
  #' @param list
  #' @export
  extractSCE <- function(sce_object, exp_type = "counts"){
    #extract metadata
    sampTab = as.data.frame(colData(sce_object, internal = TRUE))
    sampTab$sample_name = rownames(sampTab)
    
    #extract expression matrix
    if(exp_type == "counts"){
      expDat = counts(sce_object)
    }
    
    if(exp_type == "normcounts"){
      expDat = normcounts(sce_object)
    }
    
    if(exp_type == "logcounts"){
      expDat = logcounts(sce_object)
    }
    
    return(list(sampTab = sampTab, expDat = expDat))
  }
  
  m_config <- methods.config.singlecellnet 
  set.seed(100) #can be any random seed number
  train_scefile <- extractSCE(train_data, exp_type = "counts") 
  train_metadata <- train_scefile$sampTab
  train_expdata <- train_scefile$expDat
  test_scefile <- extractSCE(test_data, exp_type = "counts") 
  test_metadata <- test_scefile$sampTab
  test_expdata <- test_scefile$expDat
  
  if(m_config$cross_species){
    common_gene_file <- str_glue("{data_home}{m_config$common_gene_file}")
    oTab <- utils_loadObject(fname = "../data/human_mouse_genes_Jul_24_2018.rda")
    aa <- csRenameOrth(expQuery = test_expdata, expTrain = train_expdata, orthTable = oTab)
    test_expdata <- aa[['expQuery']]
    train_expdata <- aa[['expTrain']]
  }else{
    commonGenes<-intersect(rownames(train_expdata), rownames(test_expdata))
    train_expdata <- train_expdata[commonGenes, ]
    test_expdata <- test_expdata[commonGenes, ]
  }
  ncells <- if(purrr::is_null(m_config$ncells)) 100 else m_config$ncells
  nTopGenes <- if(purrr::is_null(m_config$nTopGenes)) 10 else m_config$nTopGenes
  nRand <- if(purrr::is_null(m_config$nRand)) 70 else m_config$nRand
  nTrees <- if(purrr::is_null(m_config$nTrees)) 1000 else m_config$nTrees
  nTopGenePairs <- if(purrr::is_null(m_config$nTopGenePairs)) 25 else m_config$nTopGenePairs
  
  class_info<-scn_train(stTrain = train_metadata, expTrain = train_expdata, nTopGenes = nTopGenes, nRand = nRand, nTrees = nTrees, nTopGenePairs = nTopGenePairs, 
                        dLevel = "label", colName_samp = "sample_name")
  #predict
  pred_results <- scn_predict(cnProc=class_info[['cnProc']], expDat=test_expdata, nrand = 0)
  pred_labels <- assign_cate(classRes = pred_results, sampTab = test_metadata, cThresh = 0.5)
  return(pred_labels$category)
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
  require(rliger)
  m_config <- methods.config.liger
  stopifnot(is(data,"SingleCellExperiment"))
  data_list <- list(counts(data))
  names(data_list) <- metadata(data)$study[[1]]
  liger_data <- createLiger(data_list)
  liger_data <- rliger::normalize(liger_data)
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


### CIDR
cluster.cidr <- function(counts, K) {
  require(cidr)
  sData <- scDataConstructor(counts(counts))
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


### scmap
cluster.scmap <- function(data){
  require(scmap)
  normcounts(data) <- counts(data)
  logcounts(data) <- log2(normcounts(data) + 1)
  rowData(data)$feature_symbol <- rownames(data)
  isSpike(data, 'ERCC') <- grepl('^ERCC-', rownames(data))
  # remove features with duplicated names
  data <- data[!duplicated(rownames(data)), ]
  data <- selectFeatures(data)

  data <- indexCluster(data, cluster_col = 'x')
  data <- scmapCluster(data, list(metadata(data)$scmap_cluster_index))
  return(data$scmap_cluster_labs)
}


  