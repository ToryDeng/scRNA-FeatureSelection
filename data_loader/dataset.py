import os
import re
from itertools import permutations

import anndata as ad
import numpy as np
import scanpy as sc
from scipy.sparse.csr import csr_matrix
from sklearn.model_selection import StratifiedKFold

from common_utils.utils import head, complexity
from config.datasets_config import data_cfg
from config.experiments_config import CellClassificationConfig
from data_loader.utils import standardize_adata, control_quality, sample_adata, store_markers, log_normalize,\
    set_rare_type


def _load_single_data(data_name: str) -> ad.AnnData:
    assert data_name in data_cfg.full_dataset_names.keys(), "Wrong argument 'data_name'!"
    adata = ad.read_h5ad(data_cfg.data_path + f'{data_cfg.full_dataset_names[data_name]}.h5ad')
    # convert csr_matrix to ndarray
    if isinstance(adata.X, csr_matrix):
        adata.X = adata.X.toarray()
    # check integers
    if not np.all(adata.X % 1 == 0):  # check raw counts
        raise ValueError(f"Dataset '{data_name}' may contain normalized data.")
    return adata


def show_data_info(adata: ad.AnnData):
    print(f"Dataset {adata.uns['data_name']} has {adata.n_obs} cells, {adata.n_vars} genes "
          f"and {adata.obs['celltype'].unique().shape[0]} classes after filtering.")
    if 'batch' not in adata.obs:
        print(f"Data complexity is {np.round(adata.uns['data_complexity'], 3)}.")
    else:
        ubatches = adata.obs['batch'].unique().categories.to_numpy()
        print(f"Dataset contains {ubatches.shape[0]} batches: {ubatches}")


def _load_data(data_name: str) -> ad.AnnData:
    """
    For specific data_name, get the corresponding normalized and raw data.

    Parameters
    ----------
    data_name
      the dataset name you want to get

    Returns
    -------
    anndata object with raw_data, norm_data and markers in data, or the concatenated batch data
    """
    if '+' not in data_name:
        dataset, number, sample_from = re.match(r"([a-zA-Z]*)([0-9]*)([a-zA-Z]*)?", data_name).groups()
        adata = _load_single_data(dataset)
        standardize_adata(adata)
        if dataset == 'Vento':  # at least 2 samples are in dataset, for computation time test
            adata = adata[adata.obs['celltype'].isin(adata.obs['celltype'].value_counts().index[:12].to_numpy()), :]
        # quality control
        adata = control_quality(adata)
        # sampling
        adata = sample_adata(adata, sample_from, number)
        # store marker genes and weights
        store_markers(adata)
        # store raw data
        adata.raw = adata
        # log-normalize anndata
        log_normalize(adata)
    else:
        batch_names = data_name.split('+')
        batch_data = [_load_data(batch_name) for batch_name in batch_names]
        adata = ad.concat(batch_data, join="inner", keys=batch_names, label='batch')
    set_rare_type(adata)
    return adata


def load_data(data_name: str) -> ad.AnnData:
    """
    load anndata object from cachedData/ directory, or the path to raw data. The adata.X is log-normalized, the raw data
    is stored in adata.raw, and the normalized data (without logarithm) is stored in adata.layers['normalized']

    Parameters
    ----------
    data_name
      the dataset name you want to get

    Returns
    -------
    adata
      anndata object with raw_data, norm_data and markers in data, or the concatenated batch data
    """
    file_name = data_name + '.h5ad'
    file_dir = os.path.join(data_cfg.cache_path, 'processedData')

    if os.path.exists(os.path.join(file_dir, file_name)):
        head(f"Loading processed {data_name} dataset", head_len=100)
        adata = sc.read_h5ad(os.path.join(file_dir, file_name))
    else:
        head(f"Loading original {data_name} dataset", head_len=100)
        adata = _load_data(data_name)
        if not os.path.exists(file_dir):
            os.makedirs(file_dir)
        adata.write_h5ad(os.path.join(file_dir, file_name))
    adata.uns['data_name'] = data_name
    adata.uns['data_complexity'] = complexity(adata, use_raw=False)  # use adata.X
    show_data_info(adata)
    return adata


def yield_train_test_data(adata: ad.AnnData, config: CellClassificationConfig):
    """
    For intra-dataset classification, this function yield (n_folds - 1) folds of data as training set, and the remaining
    fold as test set. For inter-dataset classification, this function yield a batch as training set, and another batch
    as test set.

    Parameters
    ----------
    adata
      the anndata object. For inter-dataset classification, it contains adata.obs['batch']
    config
      the configuration for cell classification
    Returns
    -------
    data_generator
      a generator that can generate train_data and test_data every time.
    """
    if config.is_intra:
        skf = StratifiedKFold(config.n_folds, random_state=config.random_seed, shuffle=True)
        for i, (train_idx, test_idx) in enumerate(skf.split(X=adata.X, y=adata.obs['celltype'].values)):
            train_data, test_data = adata[train_idx].copy(), adata[test_idx].copy()
            train_data.uns['fold'], test_data.uns['fold'] = i + 1, i + 1
            yield train_data, test_data
    else:
        indices = np.arange(adata.X.shape[0])
        for train_batch, test_batch in permutations(adata.obs['batch'].unique(), 2):
            train_mask, test_mask = adata.obs['batch'] == train_batch, adata.obs['batch'] == test_batch
            train_idx, test_idx = indices[train_mask], indices[test_mask]
            train_data, test_data = adata[train_idx].copy(), adata[test_idx].copy()
            name = ' to '.join([train_batch, test_batch])
            train_data.uns['data_name'], test_data.uns['data_name'] = name, name
            yield train_data, test_data
