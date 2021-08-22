highly_variable_gene.generate_HVGs <- function(method, data, data_name, n_genes){
  switch(method,
         deviance = highly_variable_gene.deviance(data, data_name),
         m3drop = highly_variable_gene.m3drop(data, data_name),
         scmap = highly_variable_gene.scmap(data, data_name, n_genes),
         stop("No such marker gene generating method")
         )
}

highly_variable_gene.deviance <- function(data, data_name){
  source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/deviance.R")
  count_matrix <- data@assays@data@listData[["counts"]]
  t1 <- proc.time()
  deviance_result <- sort(compute_deviance(count_matrix), decreasing = TRUE)
  t2 <- proc.time()
  # save time and markers
  time_file_name <- paste(data_name, 'time', 'deviance', sep = '_')
  write.csv((t2 - t1)[3], stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{time_file_name}.csv"))
  marker_file_name <- paste(data_name,'markers', 'deviance', sep = '_')
  write.csv(deviance_result, stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{marker_file_name}.csv"))
  
}

highly_variable_gene.m3drop <- function(data, data_name){
  require("M3Drop")
  count_matrix <- data@assays@data@listData[["counts"]]
  t1 <- proc.time()
  norm <- M3DropConvertData(count_matrix, is.counts=TRUE)
  M3Drop_genes <- M3Drop::M3DropFeatureSelection(
    norm,
    mt_method = "fdr",
    mt_threshold = 1, # do not filter genes
    suppress.plot = TRUE
  )
  t2 <- proc.time()
  # save time and markers
  time_file_name <- paste(data_name, 'time', 'm3drop', sep = '_')
  write.csv((t2 - t1)[3], stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{time_file_name}.csv"))
  marker_file_name <- paste(data_name,'markers', 'm3drop', sep = '_')
  write.csv(M3Drop_genes, stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{marker_file_name}.csv"))
}

highly_variable_gene.scmap <- function(data, data_name, n_gene){
  require(scmap)
  require(scater)
  require(dplyr)
  stopifnot(is(data,"SingleCellExperiment"))
  counts(data) <- as.matrix(counts(data))
  data <- logNormCounts(data)
  t1 <- proc.time()
  rowData(data)$feature_symbol <- rownames(data)
  data <- selectFeatures(data, n_gene, suppress_plot = TRUE)
  t2 <- proc.time()
  time_file_name <- paste(data_name, 'time', 'scmap', sep = '_')
  write.csv((t2 - t1)[3], stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{time_file_name}.csv"))
  marker_file_name <- paste(data_name,'markers', 'scmap', sep = '_')
  write.csv(rowData(data), stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{marker_file_name}.csv"), row.names = FALSE)
}


source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/load_sce.R")
setwd('/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData')

args <- commandArgs(T)
data_name <- args[1]
method <- args[2]
n_genes <- as.numeric(args[3])

if("temp_X.csv" %in% list.files() & "temp_y.csv" %in% list.files()){
  sce <- load_sce("all")
} else if("temp_X_train.csv" %in% list.files() & "temp_y_train.csv" %in% list.files()){
  sce <- load_sce("train")
} else {
  stop("ERROR: There are no generating files.")
}


highly_variable_gene.generate_HVGs(method, sce, data_name, n_genes)



