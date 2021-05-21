load_sce <- function(type){
  if (type == "all"){
    count_matrix <- as.matrix(t(read.csv('temp_X.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('temp_y.csv', row.names = 1))
  }
  else if (type == 'train'){
    count_matrix <- as.matrix(t(read.csv('temp_X_train.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('temp_y_train.csv', row.names = 1))
  }
  else if (type == 'test'){
    count_matrix <- as.matrix(t(read.csv('temp_X_test.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('temp_y_test.csv', row.names = 1))
  }
  sce <- SingleCellExperiment(assays = list(counts = count_matrix))
  colData(sce)$label <- make.names(labels)
  metadata(sce) <- list(study=c('temp'))
  return(sce)
}


main <- function(){
  
  source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/methods.R")
  require(SingleCellExperiment)
  require(Hmisc)
  
  setwd('/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData')
  sce_all <- load_sce("all")
  sce_train <- load_sce("train")
  sce_test <- load_sce("test")
  
  cluster_methods <- c('seurat','sc3')#  'tscan', , 'cidr', 'liger', 'scmap'
  assign_methods <- c('scmap_cluster', 'scmap_cell', 'singleR')
  
  for (i in 1:length(cluster_methods)){
    cluster_result <- run_cluster_methods(cluster_methods[i], sce_all)
    write.csv(cluster_result, stringr::str_glue("temp_{cluster_methods[i]}.csv"), row.names = FALSE)
  }
  for (j in 1:length(assign_methods)){
    assign_result <- run_assign_methods(assign_methods[j], sce_train, sce_test)
    write.csv(assign_result, stringr::str_glue("temp_{assign_methods[j]}.csv"), row.names = FALSE)
  }
  
}


main()
