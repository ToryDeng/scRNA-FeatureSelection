import numpy as np
import pandas as pd
import anndata as ad
import anndata2ri
import traceback
from scanpy.preprocessing._utils import _get_mean_var
from rpy2.robjects import r, globalenv
from rpy2.robjects.packages import importr
from typing import Literal
from config.experiments_config import batch_cfg


def correct_batch_effect(adata: ad.AnnData, method: Literal['Seurat_v4']) -> ad.AnnData:
    print(f"{method} correction starts.")
    try:
        if method == 'Seurat_v4':
            corrected = Seurat_v4_correct(adata)
        else:
            raise NotImplementedError(f"{method} has not been implemented!")
        return corrected
    except:
        traceback.print_exc()


def Seurat_v4_correct(adata: ad.AnnData) -> ad.AnnData:
    """
    Correct the anndata using Seurat v4

    Parameters
    ----------
    adata : ad.AnnData
      the combined data, containing norm and raw data

    Returns
    -------
    result : ad.AnnData
      corrected raw and norm data
    """
    raw_adata = adata.raw.to_adata()
    ubatches = raw_adata.obs['batch'].unique().categories.to_numpy()
    print('correcting data...')

    anndata2ri.activate()
    importr('Seurat')
    globalenv['sce'] = anndata2ri.py2rpy(raw_adata)
    r("""
    options(future.globals.maxSize= 50 * 1024 ^ 3)
    #plan("multiprocess", workers = 6)
    #saveRDS(sce, 'test.rds')

    batches_data <- as.Seurat(sce, counts='X', data=NULL)
    batch_list <- SplitObject(batches_data, split.by = 'batch')
    rm(sce)

    for(i in 1:length(batch_list)){
      batch_list[[i]] <- NormalizeData(object = batch_list[[i]], normalization.method = 'LogNormalize', verbose = FALSE)
    }
    cell_anchors <- FindIntegrationAnchors(object.list = batch_list, anchor.features = rownames(batches_data), nn.method = 'rann', verbose = FALSE)

    batches <- IntegrateData(anchorset = cell_anchors, verbose = FALSE)
    """)
    converted = r("to_converted <- as.SingleCellExperiment(batches, assay = 'integrated')")
    anndata2ri.deactivate()

    for key in converted.obs:
        if converted.obs[key].dtype == 'Int32':
            converted.obs[key] = converted.obs[key].astype(np.int)

    converted.X = converted.X.toarray()  # corrected norm data
    unnormalized = converted.copy()

    for ubatch in ubatches:
        batch_mask = (converted.obs['batch'] == ubatch).values
        batch_mean, batch_var = _get_mean_var(
            np.log1p(
                raw_adata.X[batch_mask, :] / np.expand_dims(converted.obs.loc[batch_mask, 'counts_per_cell'].values, 1) * batch_cfg.scale_factor
            ))
        batch_std = np.sqrt(batch_var)
        batch_std[batch_std == 0] = 1
        unnormalized.X[batch_mask, :] = unnormalized.X[batch_mask, :] * batch_std + batch_mean
    unnormalized.X = np.expm1(unnormalized.X)
    unnormalized.X = unnormalized.X / batch_cfg.scale_factor * np.expand_dims(converted.obs['counts_per_cell'].values, 1)
    unnormalized.X = np.where(unnormalized.X < 1e-4, 0, unnormalized.X)
    converted.raw = unnormalized  # corrected unnorm data
    converted.uns['data_name'] = adata.uns['data_name']
    converted.obs['batch'] = pd.Categorical(converted.obs['batch'])
    print('Correction complete!')
    return converted