main <- function(){
  
  source("/home/tdeng/SingleCell/FeatureSelection/R/R code/clustering methods/methods.R")
  source("/home/tdeng/SingleCell/FeatureSelection/R/R code/get_data.R")
  
  require(SingleCellExperiment)
  require(Hmisc)
  
  setwd('/home/tdeng/SingleCell/FeatureSelection/R/temp filtered data')
  count_matrix <- as.matrix(t(read.csv('temp_X.csv', row.names = 1, header = T)))
  sce <- SingleCellExperiment(assays = list(logcounts = count_matrix))

  cluster_result <- run_cluster_methods('sc3', sce)
  print(cluster_result)
}

main()

"GenomeInfoDb" %in% installed.packages()

