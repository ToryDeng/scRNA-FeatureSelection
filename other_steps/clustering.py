import traceback

import anndata as ad
import anndata2ri
import numpy as np
import sc3s
import scanpy as sc
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.robjects import r, pandas2ri, globalenv
from rpy2.robjects.packages import importr

from common_utils.utils import HiddenPrints, head
from config import cluster_cfg


def cluster_cells(adata: ad.AnnData):
    """
    cluster cells in adata. the clustering results are stored in adata.obs, with names of columns: {clustering_method}_{run}

    Parameters
    ----------
    adata
      the anndata to be clustered
    Returns
    -------
    None
    """
    for method, n_runs in cluster_cfg.methods.items():
        try:
            print(f"{method} clustering starts. {adata.n_obs} cells and {adata.n_vars} genes in data...")
            for run in range(1, n_runs + 1):
                head(method, run)
                if method == 'SHARP':
                    cluster_labels = SHARP_clustering(adata, random_seed=cluster_cfg.random_seed + run)
                elif method == 'Seurat_v4':
                    if n_runs != 1:
                        raise RuntimeWarning("Seurat v4 clustering is not a random algorithm...")
                    cluster_labels = Seurat_v4_clustering(adata)
                elif method == 'SC3':
                    cluster_labels = SC3_clustering(adata, random_seed=cluster_cfg.random_seed + run)
                elif method == 'SC3s':
                    cluster_labels = SC3s_clustering(adata, random_seed=cluster_cfg.random_seed + run)
                else:
                    raise NotImplementedError(f"{method} has not been implemented!")
                adata.obs[f'{method}_{run}'] = cluster_labels
        except:
            print(f"{method} failed.")
            traceback.print_exc()


def SHARP_clustering(adata: ad.AnnData, random_seed: int = 0):
    """
    Need raw data

    Parameters
    ----------
    adata : ad.AnnData
      anndata object containing raw counts and cell types
    random_seed
      seed used to generate fixed numbers
    Returns
    -------
    result : np.ndarray
      cluster labels
    """
    assert adata.raw is not None, "The raw counts in data must exist!"
    raw_adata = adata.raw.to_adata()
    with HiddenPrints():
        anndata2ri.activate()
        importr('SHARP')
        exp_frame, n_classes = pandas2ri.py2rpy(raw_adata.to_df().T), raw_adata.obs['celltype'].unique().shape[0]
        globalenv['raw_counts'], globalenv['k'], globalenv['seed'] = exp_frame, n_classes, random_seed
        try:
            r("res <- SHARP(raw_counts, N.cluster = k, exp.type = 'count', rN.seed = seed)")
        except RRuntimeError:
            r("res <- SHARP(raw_counts, exp.type = 'count', rN.seed = seed)")
        labekl_pred = np.squeeze(np.array(r('res$pred_clusters')))
        anndata2ri.deactivate()
    return labekl_pred


def Seurat_v4_clustering(adata: ad.AnnData):
    """
    Need raw data

    Parameters
    ----------
    adata : ad.AnnData
      anndata object containing raw counts and cell types
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
        r("""
        options(future.globals.maxSize= 50 * 1024 ^ 3)
        plan("multicore", workers = 6)
        as(sce, 'SingleCellExperiment')
        seuset <- as.Seurat(sce, counts='X', data=NULL)
        seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", verbose = FALSE)
        #seuset <- FindVariableFeatures(seuset)
        seuset <- ScaleData(object = seuset, verbose = FALSE)
        seuset <- RunPCA(object = seuset, features = rownames(seuset), verbose = FALSE)
        seuset <- FindNeighbors(object = seuset, verbose = FALSE)
        seuset <- FindClusters(object = seuset, verbose = FALSE)
        """)
        label_pred = np.squeeze(np.array(r('as.integer(unname(seuset$seurat_clusters))')))
        anndata2ri.deactivate()
    return label_pred


def SC3s_clustering(adata: ad.AnnData, random_seed: int = 0):
    """
    Need norm data.

    Parameters
    ----------
    adata : ad.AnnData
      anndata object containing normalized data and cell types
    random_seed
      seed used to generate fixed numbers
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


def SC3_clustering(adata: ad.AnnData, random_seed: int = 0):
    """
    Need raw data. Very slow.

    Parameters
    ----------
    adata : ad.AnnData
      anndata object containing raw counts and cell types
    random_seed
      seed used to generate fixed numbers
    Returns
    -------
    result : np.ndarray
      cluster labels
    """
    assert adata.raw is not None, "The raw counts in data must exist!"
    raw_adata = adata.raw.to_adata()
    with HiddenPrints():
        anndata2ri.activate()
        importr('SC3')
        importr('stringr')
        importr('scater')
        globalenv['sce'], globalenv['seed'] = anndata2ri.py2rpy(raw_adata), random_seed
        r("""
        counts(sce) <- assay(sce, 'X')
        sce <- logNormCounts(sce)
        rowData(sce)$feature_symbol <- rownames(sce)
        n_classes <- dim(unique(colData(sce)['celltype']))[1]
        sink("/dev/null")
        sce <- sc3(sce, ks=n_classes, rand_seed = seed)
        sink()
        """)
        label_pred = np.squeeze(np.array(r("colData(sce)[str_glue('sc3_{n_classes}_clusters')]")))
        anndata2ri.deactivate()
    return label_pred
