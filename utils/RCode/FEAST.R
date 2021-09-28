#' Calculate the gene-level F score and corresponding significance level.
#'
#' @param Y A gene expression matrix
#' @param classes The initial cluster labels NA values are allowed. This can directly from the \code{Consensus} function.
#' @return The score vector
#' @examples
#' data(Yan)
#' cal_F2(Y, classes = trueclass)
#' @export
cal_F2 = function(Y, classes){
  if (all(Y%%1 == 0)){
    L = colSums(Y) / median(colSums(Y))
    Ynorm = log2(sweep(Y, 2, L, FUN="/") + 1)
  }else{Ynorm = Y}
  cid = which(is.na(classes))
  if (length(cid) > 0){
    Ynorm = Ynorm[, -cid]
    classes = classes[-cid]
  }
  classes = as.factor(classes)
  unique_classes = unique(classes)
  k = length(unique(classes))
  row_class_mean = matrix(0, ncol = k, nrow = nrow(Y))
  
  row_mean = rowMeans(Ynorm)
  k = length(unique(classes))
  pb = txtProgressBar( min = 0, max = k, style = 3)
  for (i in seq_len(k)){
    setTxtProgressBar(pb, i)
    ix = which(classes == unique_classes[i])
    if (length(ix) > 1){
      tmp_mean =  rowMeans(Ynorm[, ix])
      row_class_mean[, i] = tmp_mean
    }else{
      row_class_mean[, i] = Ynorm[, ix]
    }
  }; close(pb)
  colnames(row_class_mean) = unique_classes
  ### make sure the classes are matched; otherwise, causing error ###
  table_class = table(classes)[unique_classes]
  BBS = table_class  %*% t((row_class_mean - row_mean)^2)
  BBS = BBS[1,]
  TSS = rowSums((Ynorm - row_mean) ^2)
  ESS = TSS - BBS
  df1 = k-1; df2 = ncol(Ynorm)-k
  F_scores = (BBS/df1) / (ESS/ df2)
  ps = pf(F_scores, df1, df2, lower.tail = FALSE)
  return(list(F_scores = as.numeric(F_scores), ps = ps))
}

#' function for convert a vector to a binary matrix
#' @param vec a vector.
#' @return a n by n binary matrix indicating the adjacency.
vector2matrix = function(vec){
  mat = matrix(0, nrow = length(vec), ncol = length(vec))
  diag(mat) = diag(mat) + 1
  classes = unique(vec)
  for (class in classes){
    tmp_ix = which(vec == class)
    # find all pair index of a class
    pair_ix = t(combn(tmp_ix, 2))
    pair_ix = rbind(pair_ix, pair_ix[,c(2,1)])
    mat[pair_ix] = 1
  }
  return(mat)
}

#' FEAST main function
#'
#' @param Y A expression matrix. Raw count matrix or normalized matrix.
#' @param k The number of input clusters (best guess).
#' @param num_pcs The number of top pcs that will be investigated through the consensus clustering.
#' @param dim_reduce dimension reduction methods chosen from pca, svd, or irlba.
#' @param split boolean. If T, using subsampling to calculate the gene-level significance.
#' @param batch_size when split is true, need to claim the batch size for spliting the cells.
#' @param BPPARAM parameters for BiocParallel. e.g. bpparam() and SnowParam.
#' @return the rankings of the gene-significance.
#' @examples
#' data(Yan)
#' k = length(unique(trueclass))
#' set.seed(123)
#' rixs = sample(nrow(Y), 500)
#' cixs = sample(ncol(Y), 40)
#' Y = Y[rixs, cixs]
#' ixs = FEAST(Y, k=k)
#' @export
FEAST = function (Y, k = 2, num_pcs = 10, dim_reduce = c("pca", "svd", "irlba"), split = FALSE, batch_size =1000, BPPARAM=bpparam()){
  require(mclust)
  genes = rownames(Y)
  Y = as.matrix(counts(Y))
  if (all(Y%%1 == 0)) {
    L = colSums(Y)/median(colSums(Y))
    Ynorm = log(sweep(Y, 2, L, FUN = "/") + 1)
  }else{
    Ynorm = Y
  }
  
  # dimention reduction part
  row_ms = rowMeans(Ynorm, na.rm = TRUE)
  row_sds = rowSds(Y, na.rm = TRUE)
  cv_scores = row_sds / row_ms
  gene_ranks = order(cv_scores, decreasing = TRUE, na.last = TRUE)
  # this top number of features for pca can be adjusted. Alternatively using top 1000 genes by rowMeans.
  top = round(nrow(Y)*0.33)
  gene_ixs = gene_ranks[seq_len(top)]
  YYnorm = Ynorm[gene_ixs, ]; ncells = ncol(Ynorm)
  YYnorm_scale = scale(YYnorm, scale=TRUE, center=TRUE)
  dim_reduce = match.arg(dim_reduce)
  
  # starting dimention reduction.
  pc_res = switch(dim_reduce,
                  pca = prcomp(t(YYnorm))$x,
                  svd = svd(t(YYnorm))$u,
                  irlba = prcomp_irlba(t(YYnorm_scale), n = num_pcs)$x)
  
  # setup for parallel computing. if it is SnowParam, it means on windows
  message("start consensus clustering ...")
  if (.Platform$OS.type=="windows"){
    con_mat = matrix(0, ncol = ncells, nrow = ncells)
    for (j in seq_len(num_pcs)){
      tmp_pca_mat = pc_res[, seq_len(j)]
      if (j == 1) {
        res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "V", verbose = FALSE))
      }
      else {
        res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "VVV", verbose = FALSE))
      }
      if (is.null(res)){
        res = suppressWarnings(Mclust(tmp_pca_mat, G = k, verbose = FALSE))
      }
      clusterid = apply(res$z, 1, which.max)
      con_mat = con_mat + vector2matrix(clusterid)
    }
    res = suppressWarnings(Mclust(con_mat, G = k, modelNames = "VII",  verbose = FALSE))
    if (is.null(res)) {
      res = suppressWarnings(Mclust(con_mat, G = k, verbose = FALSE))
    }
    cluster = apply(res$z, 1, function(x) {
      id = which(x > 0.95)
      if (length(id) == 0) {
        return(NA)
      }
      else {
        return(id)
      }
    })
    F_scores = cal_F2(Ynorm, classes = cluster)$F_scores
  }
  # consensus clustering for less cells (<5000)
  BPPARAM$progressbar = TRUE
  if (ncol(Ynorm) < 5000 & split == FALSE){
    mat_res = matrix(0, ncol = ncells, nrow = ncells)
    # write function for bplapply.
    bp_fun = function(i, pc_res, k){
      tmp_pca_mat = pc_res[,seq_len(i)]
      if (i == 1) {
        res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "V", verbose = FALSE))
      }
      else {
        res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "VVV", verbose = FALSE))
      }
      if (is.null(res)){
        res = suppressWarnings(Mclust(tmp_pca_mat, G = k, verbose = FALSE))
      }
      clusterid = apply(res$z, 1, which.max)
      return(clusterid)
    }
    pc_cluster = bplapply(seq_len(num_pcs), bp_fun, pc_res = pc_res, k=k, BPPARAM = BPPARAM)
    pc_mat = lapply(pc_cluster, vector2matrix)
    con_mat = Reduce("+", pc_mat)
    
    # final clustering
    res = suppressWarnings(Mclust(con_mat, G = k, modelNames = "VII",  verbose = FALSE))
    if (is.null(res)) {
      res = suppressWarnings(Mclust(con_mat, G = k, verbose = FALSE))
    }
    cluster = apply(res$z, 1, function(x) {
      id = which(x > 0.95)
      if (length(id) == 0) {
        return(NA)
      }
      else {
        return(id)
      }
    })
    F_scores = cal_F2(Ynorm, classes = cluster)$ps
  }
  else{
    split_k = round(ncells / batch_size)
    chunk_ixs = suppressWarnings(split(sample(ncol(Y)), seq_len(split_k)))
    # write function for bplapply.
    bp_fun = function(i){
      cell_ixs = chunk_ixs[[i]]
      tmp_pca_mat = pc_res[cell_ixs, seq_len(3)]
      res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "VVV", verbose = FALSE))
      if (is.null(res)){
        res = suppressWarnings(Mclust(tmp_pca_mat, G = k, verbose = FALSE))
      }
      clusterid = apply(res$z, 1, which.max)
      return(clusterid)
    }
    chunk_cluster = bplapply(seq_len(split_k), bp_fun, BPPARAM=BPPARAM)
    F_res_all = lapply(seq_len(split_k), function(j){
      cell_ixs = chunk_ixs[[j]]
      tmp_mat = Ynorm[, cell_ixs]
      tmp_cluster = chunk_cluster[[j]]
      F_scores = cal_F2(tmp_mat, classes = tmp_cluster)$ps
      return(F_scores)
    })
    F_mat = Reduce(cbind, F_res_all)
    F_scores = rowMeans(F_mat, na.rm = TRUE)
  }
  ixs = order(F_scores, decreasing = FALSE)
  
  result = data.frame(Gene=genes[ixs], p.value=F_scores[ixs])
  return(result)
}




