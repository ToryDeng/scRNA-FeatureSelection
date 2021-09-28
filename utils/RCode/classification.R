source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/methods.R")
source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/load_sce.R")
require(SingleCellExperiment)
require(Hmisc)


setwd('/home/tdeng/SingleCell/scRNA-FeatureSelection/')

sce_train <- load_sce("train")
sce_test <- load_sce("test")
assign_methods <- c('singlecellnet', 'singleR')

for (j in 1:length(assign_methods)){
  assign_result <- run_assign_methods(assign_methods[j], sce_train, sce_test)

  write.csv(assign_result, stringr::str_glue("tempData/temp_{assign_methods[j]}.csv"), row.names = FALSE)
}