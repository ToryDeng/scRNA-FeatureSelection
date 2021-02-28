#compute deviance for each gene (row) in a matrix-like object...
#...or SummarizedExperiment

#' @importFrom methods as
sparseBinomialDeviance<-function(X,sz){
  #X has features in cols, observations in rows
  #assume X is a sparseMatrix object
  X<-as(X,"CsparseMatrix")
  LP<-L1P<-X/sz #recycling
  LP@x<-log(LP@x) #log transform nonzero elements only
  L1P@x<-log1p(-L1P@x) #rare case: -Inf if only a single gene nonzero in a cell
  ll_sat<-Matrix::colSums(X*(LP-L1P)+sz*L1P, na.rm=TRUE)
  sz_sum<-sum(sz)
  feature_sums<-Matrix::colSums(X)
  p<-feature_sums/sz_sum
  l1p<-log1p(-p)
  ll_null<-feature_sums*(log(p)-l1p)+sz_sum*l1p
  2*(ll_sat-ll_null)
}

denseBinomialDeviance<-function(X,sz){
  #X has features in cols, observations in rows
  P<-X/sz
  L1P<-log1p(-P)
  ll_sat<-DelayedArray::colSums(X*(log(P)-L1P)+sz*L1P, na.rm=TRUE)
  sz_sum<-sum(sz)
  feature_sums<-DelayedArray::colSums(X)
  p<-feature_sums/sz_sum
  l1p<-log1p(-p)
  ll_null<-feature_sums*(log(p)-l1p)+sz_sum*l1p
  2*(ll_sat-ll_null)
}

#' @importFrom methods as
sparsePoissonDeviance<-function(X,sz){
  #X has features in cols, observations in rows
  X<-as(X,"CsparseMatrix")
  LP<-X/sz #recycling
  LP@x<-log(LP@x) #log transform nonzero elements only
  ll_sat<-Matrix::colSums(X*LP, na.rm=TRUE)
  feature_sums<-Matrix::colSums(X)
  ll_null<-feature_sums*log(feature_sums/sum(sz))
  2*(ll_sat-ll_null)
}

densePoissonDeviance<-function(X,sz){
  #X has features in cols, observations in rows
  ll_sat<-DelayedArray::colSums(X*log(X/sz), na.rm=TRUE)
  feature_sums<-DelayedArray::colSums(X)
  ll_null<-feature_sums*log(feature_sums/sum(sz))
  2*(ll_sat-ll_null)
}

#' @importFrom Matrix t
#' @importFrom methods is
compute_deviance<-function(m,fam=c("binomial","poisson")){
  #m is either a Matrix or matrix object (later: support DelayedArrays)
  #m a data matrix with genes=rows
  fam <- match.arg(fam)
  sz <- colSums(m)
  if(fam=="poisson"){
    lsz<-log(sz)
    #make geometric mean of sz be 1 for poisson
    sz <- exp(lsz-mean(lsz))
  }
  m<-t(m) #column slicing faster than row slicing for matrix in R.
  #note: genes are now in the COLUMNS of m
  if(is(m,"sparseMatrix")){
    if(fam=="binomial"){
      return(sparseBinomialDeviance(m,sz))
    } else { #fam=="poisson"
      return(sparsePoissonDeviance(m,sz))
    }
  } else { #m is either 1) an ordinary dense array or matrix
    # 2) a non-sparse Matrix object like dgeMatrix
    # 3) a dense object like HDF5Array (on disk) or DelayedArray (in memory)
    if(fam=="binomial"){
      return(denseBinomialDeviance(m,sz))
    } else { #fam=="poisson"
      return(densePoissonDeviance(m,sz))
    }
  }
}