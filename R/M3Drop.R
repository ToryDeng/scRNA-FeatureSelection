library("M3Drop")


get_sce <- function(data_name){
  require(SingleCellExperiment)
  require(Hmisc)
  
  if (data_name == 'PBMC10'){
    setwd('/home/tdeng/SingleCell/data/PBMC/integrated data/')
    count_matrix <- as.matrix(t(read.csv('raw_features_sample10.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('raw_labels_sample10.csv', row.names = 1))
  }
  else if (data_name == 'PBMC20'){
    setwd('/home/tdeng/SingleCell/data/PBMC/integrated data/')
    count_matrix <- as.matrix(t(read.csv('raw_features_sample20.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('raw_labels_sample20.csv', row.names = 1))
  }
  else if (data_name == 'PBMC50'){
    setwd('/home/tdeng/SingleCell/data/PBMC/integrated data/')
    count_matrix <- as.matrix(t(read.csv('raw_features_sample50.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('raw_labels_sample50.csv', row.names = 1))
  }
  else if (data_name == 'PBMC'){
    setwd('/home/tdeng/SingleCell/data/PBMC/integrated data/')
    count_matrix <- as.matrix(t(read.csv('raw_features.csv', row.names = 1, header = T)))
    labels <- as.matrix(read.csv('raw_labels.csv', row.names = 1))
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
  sce <- SingleCellExperiment(assays = list(counts = count_matrix))
  colData(sce)$label <- make.names(labels)
  return(sce)
}


for (data_name in c('muraro','segerstolpe', 'xin','all_pancreas')){
  print(data_name)
  scex <- get_sce(data_name)
  setwd('/home/tdeng/R')
  count_matrix = scex@assays@data@listData[["counts"]]
  norm = M3DropConvertData(count_matrix, is.counts=TRUE)
  M3Drop_genes <- M3Drop::M3DropFeatureSelection(
    norm,
    mt_method = "fdr",
    mt_threshold = 1
  )
  file_name <- paste(data_name,'markers', sep = '_')
  write.csv(M3Drop_genes, stringr::str_glue("{file_name}.csv"))
}

