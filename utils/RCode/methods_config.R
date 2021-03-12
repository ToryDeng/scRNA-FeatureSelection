## cluster method config
methods.config.seurat <- list(nfeatures=1000,pc_dims=10,resolution=0.5)
methods.config.tscan <- list(cvcutoff=0.01,k=8)
methods.config.sc3 <- list(nfeatures=1000,k=8,gene_filter = F)
methods.config.sc3.batch_free <- list(nfeatures=1000,k=8,gene_filter = FALSE)
methods.config.liger <- list(suggestK=F,k.suggest=25,lambda=NULL,resolution=NULL,thresh=NULL)
## assign method config
methods.config.scmap <- list(nfeatures=1000,threshold=0.5,seed=1)
methods.config.singlecellnet <- list(cross_species=FALSE,common_gene_file=NULL,ncells=50,nRand=70,nTrees=1000,nTopGenes=10,nTopGenePairs=25)
