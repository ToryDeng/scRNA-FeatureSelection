

load_sce2 <- function(type){
  require("SingleCellExperiment",character.only = TRUE)
  if (type == "all"){
    count_matrix <- as.matrix(t(read.csv('tempData/temp_X.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('tempData/temp_y.csv', row.names = 1))
  }
  else if (type == 'train'){
    count_matrix <- as.matrix(t(read.csv('tempData/temp_X_train.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('tempData/temp_y_train.csv', row.names = 1))
  }
  else if (type == 'test'){
    count_matrix <- as.matrix(t(read.csv('tempData/temp_X_test.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('tempData/temp_y_test.csv', row.names = 1))
  }
  sce <- SingleCellExperiment(assays = list(counts = count_matrix))
  colData(sce)$label <- make.names(labels)
  metadata(sce) <- list(study=c('temp'))
  return(sce)
}

load_sce <- function(type){
  require("SingleCellExperiment",character.only = TRUE)
  require("arrow")
  if (type == "all"){
    count_matrix <- read_feather('tempData/temp_X.feather')
    labels <- read_feather('tempData/temp_y.feather')
  }
  else if (type == 'train'){
    count_matrix <- read_feather('tempData/temp_X_train.feather')
    labels <- read_feather('tempData/temp_y_train.feather')
  }
  else if (type == 'test'){
    count_matrix <- read_feather('tempData/temp_X_test.feather')
    labels <- read_feather('tempData/temp_y_test.feather')
  }
  
  count_matrix <- as.data.frame(count_matrix)
  rownames(count_matrix) <- count_matrix[,1]
  count_matrix <- count_matrix[,-1]
  
  labels <- as.data.frame(labels)
  rownames(labels) <- labels[,1]
  labels <- labels[,-1]
  
  sce <- SingleCellExperiment(assays = list(counts = t(count_matrix)))
  colData(sce)$label <- make.names(labels)
  metadata(sce) <- list(study=c('temp'))
  return(sce)
}




