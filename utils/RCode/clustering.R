source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/methods.R")
source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/load_sce.R")
require(SingleCellExperiment)
require(Hmisc)

setwd('/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData')
sce_all <- load_sce("all")
cluster_methods <- c('seurat','sc3', 'cidr')#  'tscan', , 'cidr', 'liger', 'scmap'

for (i in 1:length(cluster_methods)){
  cluster_result <- run_cluster_methods(cluster_methods[i], sce_all)
  write.csv(cluster_result, stringr::str_glue("temp_{cluster_methods[i]}.csv"), row.names = FALSE)
}

#remove.packages(grep("spatstat", installed.packages(), value = T))
#.rs.restartR()
#devtools::install_version("spatstat", version = "1.64-1")