

load_sce <- function(type){
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
