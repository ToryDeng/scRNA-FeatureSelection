library(kBET)
library(lisi)
source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/load_sce.R")

setwd('/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData')
sce <- load_sce("all")

mtx<-t(counts(sce))
batch_info <- sce@colData@listData[["label"]]


k_neighbors <- round(c(0.05, 0.1, 0.15, 0.2, 0.25) * dim(mtx)[1], 0)

kBET_table <- matrix(NA, length(k_neighbors), 100)
for (i in 1:length(k_neighbors)){
  batch.estimate <- kBET(mtx, batch_info, plot=FALSE, do.pca=FALSE, k0=k_neighbors[i])
  kBET_table[i,] <- batch.estimate$stats$kBET.observed
}
cell_LISI <- lisi::compute_lisi(mtx, data.frame(batch=batch_info), 'batch')
norm_LISI <- (cell_LISI - min(cell_LISI)) / (max(cell_LISI) - min(cell_LISI))

res <- list(kBET=median(na.omit(as.vector(kBET_table))), iLISI=median(unlist(norm_LISI)))
write.csv(res, "temp_correction.csv")

