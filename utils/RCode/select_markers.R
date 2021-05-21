markergene.generate_marker_genes <- function(method, data, data_name, n_genes){
  switch(method,
         cellassign = markergene.cellassign(data, data_name, n_genes),
         deviance = markergene.deviance(data, data_name),
         m3drop = markergene.m3drop(data, data_name),
         monocle3 = markergene.monocle3(data, data_name),
         stop("No such marker gene generating method")
         )
}

markergene.cellassign <- function(data, data_name, n_genes){
  require(magrittr)
  require(limma)
  require(org.Hs.eg.db)
  require(edgeR)
  require(matrixStats)
  require(pheatmap)
  require(cellassign)
  t1 = proc.time()
  marker_config <- list(base_cell_type=NULL)
  rowData(data)$geneName <- rownames(data)
  ensembl_map <- dplyr::transmute(tibble::as_tibble(rowData(data)),SYMBOL=geneName)
  #> 'select()' returned 1:1 mapping between keys and columns
  gene_annotations <- ensembl_map %>%
    dplyr::rename(Symbol=SYMBOL)
  cell_lines <- unique(colData(data)$label)
  dge <- DGEList(counts = counts(data), 
                 group = colData(data)$label, 
                 genes = gene_annotations, 
                 remove.zeros = TRUE)
  genes_to_keep <- rowSums(edgeR::cpm(dge$counts) > 0.5) >= 2
  dge_filt <- dge[genes_to_keep,]
  dge_filt <- calcNormFactors(dge_filt, method="TMM")
  
  design <- model.matrix(~ 0+dge_filt$samples$group)
  colnames(design) <- levels(dge_filt$samples$group)
  v <- voom(dge_filt, design)
  fit <- lmFit(v, design)
  base_cell_type <- if(purrr::is_null(marker_config$base_cell_type)) unique(colData(data)$label)[[1]] else marker_config$base_cell_type
  args <- purrr::map(1:ncol(combn(cell_lines,2)), ~{ if(combn(cell_lines,2)[,.][[1]]==base_cell_type)
    return(stringr::str_glue("{combn(cell_lines,2)[,.][[2]]} - {combn(cell_lines,2)[,.][[1]]}"))
    stringr::str_glue("{combn(cell_lines,2)[,.][[1]]} - {combn(cell_lines,2)[,.][[2]]}")
  })
  args$levels <- design
  contrast.matrix <- do.call('makeContrasts',args)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, n=Inf)
  tt_sig <- tt %>%
    dplyr::filter(adj.P.Val < 0.05)
  ##########
  
  
  #####how to determine the base cell type?????
  
  diff_baseline_celltypes <- unlist(purrr::map(cell_lines[which(cell_lines!=base_cell_type)],~{stringr::str_glue("{.}...{base_cell_type}")}))
  lfc_table <- tt_sig[,diff_baseline_celltypes]
  colnames(lfc_table) <- purrr::map(diff_baseline_celltypes,~{stringr::str_split(.,"\\.\\.\\.")[[1]][1]})
  lfc_table[[base_cell_type]] <- 0
  lfc_table <- as.matrix(lfc_table)
  lfc_table <- lfc_table - rowMins(lfc_table)
  lfc_table <- as.data.frame(lfc_table)
  binarize <- function(x, threshold) {
    x[x <= threshold] <- -Inf
    x[x > -Inf] <- 1
    x[x == -Inf] <- 0
    return(x)
  }
  # Find the biggest difference
  maxdiffs <- apply(lfc_table, 1, function(x) max(diff(sort(x))))
  thres_vals <- apply(lfc_table, 1, function(x) sort(x)[which.max(diff(sort(x)))])
  expr_mat_thres <- plyr::rbind.fill(lapply(1:nrow(lfc_table), function(i) {
    binarize(lfc_table[i,], thres_vals[i])
  }))
  rownames(expr_mat_thres) <- rownames(lfc_table)
  
  nth_diff <- sort(maxdiffs, decreasing = T)[n_genes]
  marker_gene_mat <- expr_mat_thres[(maxdiffs >= nth_diff) ,] %>%
    as.matrix
  t2 <- proc.time()
  # save time and markers
  time_file_name <- paste(data_name, 'time', 'cellassign', sep = '_')
  write.csv((t2 - t1)[3], stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{time_file_name}.csv"))
  marker_file_name <- paste(data_name,'markers', 'cellassign', sep = '_')
  write.csv(marker_gene_mat, stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{marker_file_name}.csv"))
  
}

markergene.deviance <- function(data, data_name){
  source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/deviance.R")
  count_matrix <- data@assays@data@listData[["counts"]]
  t1 <- proc.time()
  deviance_result <- sort(compute_deviance(count_matrix), decreasing = TRUE)[1:1000]
  t2 <- proc.time()
  # save time and markers
  time_file_name <- paste(data_name, 'time', 'deviance', sep = '_')
  write.csv((t2 - t1)[3], stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{time_file_name}.csv"))
  marker_file_name <- paste(data_name,'markers', 'deviance', sep = '_')
  write.csv(deviance_result, stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{marker_file_name}.csv"))
  
}

markergene.m3drop <- function(data, data_name){
  require("M3Drop")
  count_matrix <- data@assays@data@listData[["counts"]]
  t1 <- proc.time()
  norm <- M3DropConvertData(count_matrix, is.counts=TRUE)
  M3Drop_genes <- M3Drop::M3DropFeatureSelection(
    norm,
    mt_method = "fdr",
    mt_threshold = 1, # do not filter genes
    suppress.plot = TRUE
  )
  t2 <- proc.time()
  # save time and markers
  time_file_name <- paste(data_name, 'time', 'm3drop', sep = '_')
  write.csv((t2 - t1)[3], stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{time_file_name}.csv"))
  marker_file_name <- paste(data_name,'markers', 'm3drop', sep = '_')
  write.csv(M3Drop_genes, stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{marker_file_name}.csv"))
}

markergene.monocle3 <- function(data, data_name){
  require(monocle3)
  require(magrittr)
  require(dplyr)
  exp_matrix <- data@assays@data@listData[["counts"]]
  t1 <- proc.time()
  cds <- new_cell_data_set(expression_data = exp_matrix)
  ## Step 1: Normalize and pre-process the data
  cds <- preprocess_cds(cds, num_dim = 100)
  ## Step 2 (Remove batch effects with cell alignment) is optional
  ## Step 3: Reduce the dimensions using UMAP
  cds <- reduce_dimension(cds)
  ## Step 4: Cluster the cells
  cds <- cluster_cells(cds)
  marker_test_res <- top_markers(cds,cores=10, genes_to_test_per_group=250, reference_cells = 1000)
  results <- marker_test_res[order(marker_test_res['pseudo_R2'], decreasing = T),]
  drop_dp_result <- results[!duplicated(results$gene_id),]
  t2 <- proc.time()
  # save time and markers
  time_file_name <- paste(data_name, 'time', 'monocle3', sep = '_')
  write.csv((t2 - t1)[3], stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{time_file_name}.csv"))
  marker_file_name <- paste(data_name,'markers', 'monocle3', sep = '_')
  write.csv(drop_dp_result, stringr::str_glue("/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData/{marker_file_name}.csv"))
  
}


source("/home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/load_sce.R")
setwd('/home/tdeng/SingleCell/scRNA-FeatureSelection/tempData')

args <- commandArgs(T)
data_name <- args[1]
method <- args[2]
n_genes <- as.numeric(args[3])
if("temp_X.csv" %in% list.files() & "temp_y.csv" %in% list.files()){
  sce <- load_sce("all")
} else if("temp_X_train.csv" %in% list.files() & "temp_y_train.csv" %in% list.files()){
  sce <- load_sce("train")
} else {
  stop("ERROR: There are no generating files.")
}


markergene.generate_marker_genes(method, sce, data_name, n_genes)



