source("get_data.R")
library("M3Drop")


for (data_name in c('muraro','segerstolpe', 'xin','all_pancreas')){
  print(data_name)
  scex <- get_data(data_name)
  count_matrix = scex@assays@data@listData[["counts"]]
  norm = M3DropConvertData(count_matrix, is.counts=TRUE)
  M3Drop_genes <- M3Drop::M3DropFeatureSelection(
    norm,
    mt_method = "fdr",
    mt_threshold = 1 # do not filter genes
  )
  file_name <- paste(data_name,'markers', sep = '_')
  write.csv(M3Drop_genes, stringr::str_glue("/home/tdeng/R/M3drop/{file_name}.csv"))
}

