import anndata as ad
import scanpy as sc
import numpy as np
import besca as bc
from typing import Literal, Tuple

import pandas as pd
from scipy.sparse.csr import csr_matrix
from config import data_cfg, exp_cfg
from utils.utils import head, complexity, HiddenPrints
import os
import re


def clean_var_names(gene_names: pd.Index) -> np.ndarray:
    regex = re.compile(pattern='[-_:+()|]')
    vreplace = np.vectorize(lambda x: regex.sub('.', x), otypes=[np.str])
    return vreplace(gene_names.to_numpy())


def standardize_adata(adata: ad.AnnData):
    adata.var_names = clean_var_names(adata.var_names)
    adata.obs.rename(columns={adata.obs.columns[0]: 'celltype'}, inplace=True)
    if 'Batch' in adata.obs:
        adata.obs.rename(columns={'Batch': 'batch'}, inplace=True)
    regex = re.compile(pattern="[- +']")
    vreplace = np.vectorize(lambda x: regex.sub('.', str(x)), otypes=[np.str])
    vladd = np.vectorize(lambda x: 'X' + x if x[0].isdigit() else x, otypes=[np.str])
    adata.obs['celltype'] = pd.Categorical(vladd(vreplace(adata.obs['celltype'].values)))
    if 'batch' in adata.obs and not isinstance(adata.obs['batch'], pd.Categorical):
        adata.obs['batch'] = pd.Categorical(adata.obs['batch'])
    adata.obs_names_make_unique(join='.')
    adata.var_names_make_unique(join='.')


def control_quality(adata: ad.AnnData) -> ad.AnnData:
    # filter unclear cells
    adata = adata[~adata.obs['celltype'].isin(data_cfg.remove_types)]
    # replace cell types
    for replaced, to_replace in data_cfg.replace_types.items():
        adata.obs['celltype'].replace(to_replace, replaced, inplace=True)
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
    sc.pp.normalize_total(adata, target_sum=exp_cfg.scale_factor, inplace=True, key_added='counts_per_cell')
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
            np.random.seed(exp_cfg.random_seed)
            adata = adata[:, np.random.choice(adata.n_vars, size=int(number))]
        else:
            raise ValueError("You input an invalid  source to sample!")
    return adata


def load_marekrs(marker_type: Literal['PBMC', 'SimPBMCsmall', 'pancreas', 'brain']) -> Tuple[np.ndarray, np.ndarray]:
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
            panglao, cellmarker = load_marekrs('PBMC')
        elif adata.uns['data_name'] in ('Segerstolpe', 'BaronHuman'):
            panglao, cellmarker = load_marekrs('pancreas')
        elif adata.uns['data_name'] == 'Zeisel':
            panglao, cellmarker = load_marekrs('brain')
        else:
            pass
        if panglao is not None and cellmarker is not None:  #
            markers = np.union1d(panglao, cellmarker)  # either in PanglaoDB or in CellMarker
            adata.uns['markers'] = np.intersect1d(markers, adata.var_names)  # only contains genes in data
            adata.uns['marker_weight'] = np.where(
                np.isin(adata.uns['markers'], np.intersect1d(panglao, cellmarker)), 2., 1.
            )


def load_single_data(data_name: str) -> ad.AnnData:
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
    print(f"Dataset {adata.uns['data_name']} has {adata.n_obs} cells, {adata.n_vars} genes "
          f"and {adata.obs['celltype'].unique().shape[0]} classes after filtering.")
    if not is_batch and 'batch' not in adata.obs:
        print(f"Data complexity is {np.round(adata.uns['data_complexity'], 3)}.")


def load_data(data_name: str, is_batch=False) -> ad.AnnData:
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
        adata = load_single_data(dataset)
        standardize_adata(adata)

        if dataset == 'Vento':  # at least 2 samples are in dataset, for computation time test
            adata = adata[adata.obs['celltype'].isin(adata.obs['celltype'].value_counts().index[:12].to_numpy()), :]

        if not is_batch:
            head(dataset, head_len=100)
            adata.uns['data_name'] = dataset
            adata.uns['data_complexity'] = complexity(adata, use_raw=False)  # use adata.X
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
        batch_data = [load_data(batch_name, is_batch=True) for batch_name in batch_names]
        adata = ad.concat(batch_data, join="inner", keys=batch_names, label='batch')
        adata.uns['data_name'] = data_name
    show_data_info(adata, is_batch)
    return adata
