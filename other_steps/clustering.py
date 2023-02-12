import datetime

import anndata as ad
import anndata2ri
import numpy as np
import sc3s
import scanpy as sc
from rpy2.robjects import r, pandas2ri, globalenv
from rpy2.robjects.packages import importr
from sklearn.cluster import KMeans

from common_utils.utils import HiddenPrints, head
from config import cluster_cfg


def cluster_cells(adata: ad.AnnData):
    """
    cluster cells in adata. the clustering results are stored in adata.obs, with names of columns: {clustering_method}_{run}

    Parameters
    ----------
    adata
      the AnnData to be clustered
    Returns
    -------
    None
    """
    for method, n_runs in cluster_cfg.methods.items():
        print(f"{datetime.datetime.now()}: {method} clustering starts. {adata.n_obs} cells and {adata.n_vars} genes in data...")
        for run in range(n_runs):
            head(method, run)
            if method == 'SHARP':
                cluster_labels = SHARP_clustering(adata, random_seed=cluster_cfg.random_seed + run)
            elif method == 'Seurat_v4':
                cluster_labels = Seurat_v4_clustering(adata, random_seed=cluster_cfg.random_seed + run)
            elif method == 'TSCAN':
                cluster_labels = TSCAN_clustering(adata)
            elif method == 'CIDR':
                cluster_labels = CIDR_clustering(adata)
            elif method == 'KMeans':
                cluster_labels = KMeans_clustering(adata, random_seed=cluster_cfg.random_seed + run)
            elif method == 'SC3s':
                cluster_labels = SC3s_clustering(adata, random_seed=cluster_cfg.random_seed + run)
            else:
                raise NotImplementedError(f"{method} has not been implemented!")
            adata.obs[f'{method}_{run}'] = cluster_labels


def SHARP_clustering(adata: ad.AnnData, random_seed: int):
    """
    Clustering cells using SHARP.

    Parameters
    ----------
    adata : ad.AnnData
      AnnData object containing raw counts and cell types
    random_seed: int
      An integer for reproducibility

    Returns
    -------
    result : np.ndarray
      cluster labels
    """
    with HiddenPrints():
        anndata2ri.activate()
        importr('SHARP')
        globalenv['expr'] = pandas2ri.py2rpy(adata.to_df('log-normalized').T)
        globalenv['k'] = adata.obs['celltype'].unique().shape[0]
        globalenv['seed'] = random_seed
        r("res <- SHARP(expr, N.cluster = k, prep = FALSE, rN.seed = seed, n.cores = 1)")
        label_pred = np.squeeze(np.array(r('res$pred_clusters')))
        anndata2ri.deactivate()
    return label_pred


def Seurat_v4_clustering(adata: ad.AnnData, random_seed: int):
    """
    Clustering cells using Seurat v4.

    Parameters
    ----------
    adata : ad.AnnData
      AnnData object containing raw counts and cell types
    random_seed: int
      An integer for reproducibility

    Returns
    -------
    result : np.ndarray
      cluster labels
    """
    assert adata.raw is not None, "The raw counts in data must exist!"
    raw_adata = adata.raw.to_adata()
    with HiddenPrints():
        anndata2ri.activate()
        importr('Seurat')
        importr('future')
        importr('doParallel')
        globalenv['sce'] = anndata2ri.py2rpy(raw_adata)
        globalenv['seed'] = random_seed
        r("""
        options(future.globals.maxSize= Inf)
        plan("multicore", workers = 6)
        as(sce, 'SingleCellExperiment')
        seuset <- as.Seurat(sce, counts='X', data=NULL)
        seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", verbose = FALSE)
        #seuset <- FindVariableFeatures(seuset)
        seuset <- ScaleData(object = seuset, verbose = FALSE)
        seuset <- RunPCA(object = seuset, features = rownames(seuset), verbose = FALSE)
        seuset <- FindNeighbors(object = seuset, verbose = FALSE)
        seuset <- FindClusters(object = seuset, random.seed = seed, verbose = FALSE)
        """)
        label_pred = np.squeeze(np.array(r('as.integer(unname(seuset$seurat_clusters))')))
        anndata2ri.deactivate()
    return label_pred


def SC3s_clustering(adata: ad.AnnData, random_seed: int):
    """
    Clustering cells using SC3s. Need normalized data.

    Parameters
    ----------
    adata : ad.AnnData
      AnnData object containing normalized data and cell types
    random_seed: int
      An integer for reproducibility
    Returns
    -------
    result : np.ndarray
      cluster labels
    """
    with HiddenPrints():
        n_classes = adata.obs['celltype'].unique().shape[0]
        if 'X_pca' in adata.obsm:
            del adata.obsm['X_pca']  # do pca anyway
        sc.pp.pca(adata)
        cluster_adata = sc3s.tl.consensus(adata, n_clusters=[n_classes], random_state=random_seed)
    pred_labels = cluster_adata.obs[f"sc3s_{n_classes}"].to_numpy()
    del cluster_adata.obs[f"sc3s_{n_classes}"]
    return pred_labels


def TSCAN_clustering(adata: ad.AnnData):
    with HiddenPrints():
        anndata2ri.activate()
        importr('TSCAN')
        globalenv['expr'] = pandas2ri.py2rpy(adata.to_df('log-normalized').T)
        globalenv['k'] = adata.obs['celltype'].unique().shape[0]
        r("res <- exprmclust(expr, reduce = T, clusternum = k:k)")
        label_pred = np.squeeze(np.array(r('res$clusterid')))
        anndata2ri.deactivate()
    return label_pred


def KMeans_clustering(adata: ad.AnnData, random_seed: int):
    X_pca = sc.tl.pca(adata.X, return_info=False, random_state=random_seed)
    label_pred = KMeans(n_clusters=adata.obs['celltype'].unique().shape[0], random_state=random_seed).fit_predict(X_pca)
    return label_pred


def CIDR_clustering(adata: ad.AnnData):
    raw_adata = adata.raw.to_adata()
    with HiddenPrints():
        anndata2ri.activate()
        importr('cidr')
        globalenv['sce'], globalenv['k'] = anndata2ri.py2rpy(raw_adata), adata.obs['celltype'].unique().shape[0]
        r("""
        sData <- scDataConstructor(assay(sce, 'X'), tagType = "raw")
        sData <- determineDropoutCandidates(sData)
        sData <- wThreshold(sData)
        print('start to calculate dissimilarity...')
        sData <- scDissim(sData, threads = 0)
        sData <- scPCA(sData, plotPC = FALSE)
        sData <- nPC(sData)
        sDataC <- scCluster(object = sData, nCluster = k, nPC = sData@nPC, cMethod = "ward.D2")
        """)
        label_pred = np.squeeze(np.array(r("sDataC@clusters")))
        anndata2ri.deactivate()
    return label_pred
