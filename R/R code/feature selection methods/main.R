source("/home/tdeng/R/get_data.R")
source("/home/tdeng/R/CellAssign.R")
source("/home/tdeng/R/Deviance.R")
library("M3Drop")

for (data_name in c('PBMC10', 'PBMC20', 'PBMC50', 'muraro','segerstolpe', 'xin', 'all_pancreas')){
  print(data_name)
  sce <- get_data(data_name)
  count_matrix = sce@assays@data@listData[["counts"]]
  file_name <- paste(data_name,'markers', sep = '_')
  # CellAssgin
  cellassign_result <- cellassign(sce, 1000)
  write.csv(cellassign_result, stringr::str_glue("/home/tdeng/R/CellAssign/{file_name}.csv"))
  
  # M3Drop
  norm <- M3DropConvertData(count_matrix, is.counts=TRUE)
  M3Drop_result <- M3Drop::M3DropFeatureSelection(
    norm,
    mt_method = "fdr",
    mt_threshold = 1 # do not filter genes
  )
  write.csv(M3Drop_result, stringr::str_glue("/home/tdeng/R/M3Drop/{file_name}.csv"))
  
  # Deviance
  deviance_result <- sort(compute_deviance(count_matrix), decreasing = TRUE)[1:1000]
  write.csv(deviance_result, stringr::str_glue("/home/tdeng/R/Deviance/{file_name}.csv"))
  
}