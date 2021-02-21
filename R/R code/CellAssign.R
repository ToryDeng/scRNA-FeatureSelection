
cellassign <- function(data, n_genes){
  require(magrittr)
  require(limma)
  require(org.Hs.eg.db)
  require(edgeR)
  require(matrixStats)
  require(pheatmap)
  require(cellassign)
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
  return(marker_gene_mat)
}







