# return SingleCellExpieriment object(row:genes, col:cells)

get_data <- function(data_name){
  
  require(SingleCellExperiment)
  require(Hmisc)
  
  if (substr(data_name, 1, 4) == 'PBMC'){
    setwd('/home/tdeng/SingleCell/data/PBMC/integrated data')
    percent <- substr(data_name, 5, 6)
    count_matrix <- as.matrix(t(read.csv(paste0('raw_features_sample', percent, '.csv'), row.names = 1, header = T)))
    labels <- as.matrix(read.csv(paste0('raw_labels_sample', percent, '.csv'), row.names = 1))
  }
  else if(data_name%in%c('muraro','segerstolpe','xin')){
    setwd('/home/tdeng/SingleCell/data/pancreas/separated data/')
    count_matrix <- as.matrix(t(read.csv(paste0(capitalize(data_name), '_pancreas_filtered.csv'), row.names = 1, header = T)))
    labels <- as.matrix(read.csv(paste0(capitalize(data_name), '_trueID_filtered.csv'), row.names = 1))
  }
  else if(data_name == 'all_pancreas'){
    setwd('/home/tdeng/SingleCell/data/pancreas/integrated data/')
    count_matrix <- as.matrix(t(read.csv('features.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('labels.csv', row.names = 1))
  }
  else{
    print("parameter 'data_name' is wrong!")
    return(NULL)
  }
  setwd('/home/tdeng/SingleCell/')
  sce <- SingleCellExperiment(assays = list(counts = count_matrix[which(rowSums(count_matrix) > 0),which(colSums(count_matrix) > 0)]))
  colData(sce)$label <- make.names(labels)
  return(sce)
}