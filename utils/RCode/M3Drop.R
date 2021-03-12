source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/get_data.R")
library("M3Drop")


data_name<-commandArgs(T)
scex <- get_data(data_name)
count_matrix = scex@assays@data@listData[["counts"]]
norm = M3DropConvertData(count_matrix, is.counts=TRUE)
M3Drop_genes <- M3Drop::M3DropFeatureSelection(
  norm,
  mt_method = "fdr",
  mt_threshold = 1 # do not filter genes
)
file_name <- paste(data_name,'markers_m3drop', sep = '_')
write.csv(M3Drop_genes, stringr::str_glue("scRNA-FeatureSelection/tempData/{file_name}.csv"))


