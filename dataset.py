import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import besca as bc
import os
import re

from functools import reduce
from itertools import permutations
from typing import Literal, Tuple
from sklearn.model_selection import StratifiedKFold
from scipy.sparse.csr import csr_matrix
from config.datasets_config import data_cfg
from config.experiments_config import base_cfg, assign_cfg, CellClassificationConfig
from utils import head, complexity, HiddenPrints


def clean_var_names(gene_names: pd.Index) -> np.ndarray:
    regex = re.compile(pattern='[-_:+()|]')
    vreplace = np.vectorize(lambda x: regex.sub('.', x), otypes=[np.str])
    return vreplace(gene_names.to_numpy())


def standardize_adata(adata: ad.AnnData):
    """
    rename the columns in adata.obs and make obs and var names unique.
    Parameters
    ----------
    adata
        anndata object
    Returns
    -------
    inplace, return None
    """
    # rename genes
    adata.var_names = clean_var_names(adata.var_names)
    # rename celltype column
    adata.obs.rename(columns={adata.obs.columns[0]: 'celltype'}, inplace=True)
    # rename batch column
    if 'Batch' in adata.obs:
        adata.obs.rename(columns={'Batch': 'batch'}, inplace=True)
    # replace cell types
    for replaced, to_replace in data_cfg.replace_types.items():
        adata.obs['celltype'].replace(to_replace, replaced, inplace=True)
    # standardize cell types
    regex = re.compile(pattern="[- +']")
    vreplace = np.vectorize(lambda x: regex.sub('.', str(x)), otypes=[np.str])
    vladd = np.vectorize(lambda x: 'X' + x if x[0].isdigit() else x, otypes=[np.str])
    adata.obs['celltype'] = pd.Categorical(vladd(vreplace(adata.obs['celltype'].values)))
    # categorize batch column
    if 'batch' in adata.obs and not isinstance(adata.obs['batch'], pd.Categorical):
        adata.obs['batch'] = pd.Categorical(adata.obs['batch'])
    # make cells and genes unique
    adata.obs_names_make_unique(join='.')
    adata.var_names_make_unique(join='.')


def control_quality(adata: ad.AnnData) -> ad.AnnData:
    # filter unclear cells
    adata = adata[~adata.obs['celltype'].isin(data_cfg.remove_types)]

    adata.var['SYMBOL'] = adata.var_names
    with HiddenPrints():
        min_genes, min_cells, min_counts, n_genes, percent_mito, max_counts = bc.pp.valOutlier(adata)
        # filtering with thresholds of gene and cell counts
        adata = bc.st.filtering_cells_genes_min(adata, min_cells, min_genes, min_counts)
        # filtering with thresholds of proportion of mitochondrial genes and the upper limit of gene counts
        adata = bc.st.filtering_mito_genes_max(adata, percent_mito, n_genes, max_counts)
    adata.var.drop(columns=['SYMBOL'], inplace=True)
    return adata


def log_normalize(adata: ad.AnnData):
    """
    log-Normalize and scale data

    Parameters
    ----------
    adata
      anndata object

    Returns
    -------
    inplace, return None
    """
    #
    sc.pp.normalize_total(adata, target_sum=base_cfg.scale_factor, inplace=True, key_added='counts_per_cell')
    adata.layers['normalized'] = adata.X
    sc.pp.log1p(adata)
    if 'batch' in adata.obs:
        print('scaling data for each batch...')
        for batch in adata.obs['batch'].unique():
            batch_mask = adata.obs['batch'] == batch
            adata.X[batch_mask, :] = sc.pp.scale(adata.X[batch_mask, :], max_value=10, copy=True)
    else:
        sc.pp.scale(adata, max_value=10, copy=False)


def sample_adata(adata: ad.AnnData, sample_source: str, number: str) -> ad.AnnData:
    if sample_source != '' and number != '':
        if sample_source == 'cells':
            adata = sc.pp.subsample(adata, n_obs=int(number), random_state=0, copy=True)
        elif sample_source == 'genes':
            np.random.seed(base_cfg.random_seed)
            adata = adata[:, np.random.choice(adata.n_vars, size=int(number))]
        else:
            raise ValueError("You input an invalid  source to sample!")
    return adata


def load_markers(marker_type: Literal['PBMC', 'SimPBMCsmall', 'pancreas', 'brain']) -> Tuple[np.ndarray, np.ndarray]:
    marker_path = lambda x: os.path.join(data_cfg.marker_path, x)
    kwargs = {'skiprows': 1, 'dtype': np.str, 'delimiter': ','}
    if marker_type == 'PBMC':
        panglao = np.unique(np.loadtxt(marker_path('PBMC_panglaoDB.csv'), usecols=[0], **kwargs))
        cellmarker = np.unique(np.loadtxt(marker_path('PBMC_CellMarker.csv'), usecols=[1], **kwargs))
    elif marker_type == 'pancreas':
        panglao = np.loadtxt(marker_path('pancreas_panglaoDB.csv'), usecols=[0], **kwargs)
        cellmarker = np.loadtxt(marker_path('pancreas_CellMarker.csv'), usecols=[0], **kwargs)
    elif marker_type == 'brain':
        panglao = np.loadtxt(marker_path('MouseBrain_panglaoDB.csv'), usecols=[1], **kwargs)
        cellmarker = np.loadtxt(marker_path('MouseBrain_CellMarker.csv'), usecols=[0], **kwargs)
    else:
        raise ValueError(f"Wrong argument 'marker_type'!")
    return panglao, cellmarker


def store_markers(adata: ad.AnnData):
    if 'data_name' in adata.uns:
        panglao, cellmarker = None, None
        if adata.uns['data_name'] in ('PBMCsmall', 'PBMCbatchone', 'PBMCbatchtwo'):  # PBMC
            panglao, cellmarker = load_markers('PBMC')
        elif adata.uns['data_name'] in ('Segerstolpe', 'BaronHuman'):
            panglao, cellmarker = load_markers('pancreas')
        elif adata.uns['data_name'] == 'Zeisel':
            panglao, cellmarker = load_markers('brain')
        else:
            pass
        if panglao is not None and cellmarker is not None:  #
            all_markers = np.union1d(panglao, cellmarker)  # either in PanglaoDB or in CellMarker
            adata.var['is_marker'] = np.isin(adata.var_names, all_markers)  # only contains genes in data
            adata.var['marker_weight'] = np.where(np.isin(adata.var_names, np.intersect1d(panglao, cellmarker)), 2., 1.)
            adata.var['marker_weight'] = np.where(adata.var['is_marker'], adata.var['marker_weight'], 0)


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


def show_data_info(adata: ad.AnnData, is_batch: bool = False):
    if not is_batch:
        print(f"Dataset {adata.uns['data_name']} has {adata.n_obs} cells, {adata.n_vars} genes "
              f"and {adata.obs['celltype'].unique().shape[0]} classes after filtering.")
        if 'batch' not in adata.obs:

            print(f"Data complexity is {np.round(adata.uns['data_complexity'], 3)}.")
        else:
            ubatches = adata.obs['batch'].unique().categories.to_numpy()
            print(f"Dataset contains {ubatches.shape[0]} batches: {ubatches}")


def set_rare_type(adata: ad.AnnData):
    cell_type_counts = adata.obs['celltype'].value_counts(ascending=True)
    if 'batch' not in adata.obs:
        if cell_type_counts[0] >= assign_cfg.n_folds * 2:  # 10
            rare_type = cell_type_counts.index[0]
            print("Rare Cell Type:{:<15}     Rate:{:.2%}     Num:{}".format(
                rare_type, cell_type_counts[rare_type] / adata.n_obs, cell_type_counts[0])
            )
        else:
            rare_type = None
    else:
        u_batches = adata.obs['batch'].unique()
        bu_types = [adata[adata.obs['batch'] == b].obs['celltype'].unique() for b in u_batches]
        inter_types = np.intersect1d(*bu_types) if len(bu_types) == 2 else reduce(np.intersect1d, bu_types)
        try:
            rare_type = cell_type_counts.index[np.isin(cell_type_counts.index.to_numpy(), inter_types)][0]
            print("Rare Cell Type:{:<15}     Rate in the total dataset: {:.2%}     Num: {}".format(
                rare_type, cell_type_counts[rare_type] / adata.n_obs, cell_type_counts[rare_type]
            ))
            for b in u_batches:
                batch_mask = adata.obs['batch'] == b
                batch_counts = adata[batch_mask].obs['celltype'].value_counts()
                print("Rate in {}: {:.2%}     Num in {}: {}".format(
                    b, batch_counts[rare_type] / batch_mask.sum(), b, batch_counts[rare_type]
                ))
        except IndexError:
            rare_type = None
    adata.uns['rare_type'] = rare_type


def _load_data(data_name: str, is_batch=False) -> ad.AnnData:
    """
    For specific data_name, get the corresponding normalized and raw data.

    Parameters
    ----------
    data_name
      the dataset name you want to get
    is_batch
      whether it is a part of batch

    Returns
    -------
    anndata object with raw_data, norm_data and markers in data, or the concatenated batch data
    """
    if '+' not in data_name:
        dataset, number, sample_from = re.match(r"([a-zA-Z]*)([0-9]*)([a-zA-Z]*)?", data_name).groups()
        adata = _load_single_data(dataset)
        standardize_adata(adata)
        if not is_batch:
            head(dataset, head_len=100)
            adata.uns['data_name'] = dataset
            adata.uns['data_complexity'] = complexity(adata, use_raw=False)  # use adata.X
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
        head(data_name, head_len=100)
        batch_names = data_name.split('+')
        batch_data = [_load_data(batch_name, is_batch=True) for batch_name in batch_names]
        adata = ad.concat(batch_data, join="inner", keys=batch_names, label='batch')
        adata.uns['data_name'] = data_name
    set_rare_type(adata)
    show_data_info(adata, is_batch)
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
    anndata object with raw_data, norm_data and markers in data, or the concatenated batch data
    """
    file_name = data_name + '.h5ad'
    if os.path.exists(data_cfg.cache_path):
        if file_name in os.listdir(data_cfg.cache_path):
            adata = sc.read_h5ad(os.path.join(data_cfg.cache_path, file_name))
        else:
            adata = _load_data(data_name)
            adata.write_h5ad(os.path.join(data_cfg.cache_path, file_name))
    else:
        os.mkdir("cachedData/")
        adata = _load_data(data_name)
        adata.write_h5ad(os.path.join(data_cfg.cache_path, file_name))
    return adata


def yield_train_test_data(adata: ad.AnnData, config: CellClassificationConfig):
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
