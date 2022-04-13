import os
import re
from functools import reduce
from math import e
from typing import Literal, Tuple

import anndata as ad
import besca as bc
import numpy as np
import pandas as pd
import scanpy as sc

from common_utils.utils import HiddenPrints
from config import data_cfg, base_cfg


def complexity(adata: ad.AnnData, use_raw: bool = False):
    """
    Compute the complexity of the dataset.

    Parameters
    ----------
    adata
      anndata object
    use_raw
      whether to use raw data
    Returns
    -------
    complexity
      the complexity of the dataset
    """
    if use_raw:
        adata = adata.raw.to_adata()
    all_centroids = pd.concat(
        objs=[adata[adata.obs['celltype'] == t].to_df().mean(axis=0).rename(t) for t in adata.obs['celltype'].unique()],
        axis=1).T
    corr_mtx = np.corrcoef(all_centroids)
    return np.mean([corr_mtx[i, j] for i, j in enumerate(np.argsort(corr_mtx)[:, -2])])


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
    sc.pp.log1p(adata, base=e)
    if 'batch' in adata.obs:
        print('scaling data for each batch...')
        for batch in adata.obs['batch'].unique():
            batch_mask = adata.obs['batch'] == batch
            adata.X[batch_mask, :] = sc.pp.scale(adata.X[batch_mask, :], max_value=10, copy=True)
    else:
        sc.pp.scale(adata, max_value=10, copy=False)


def sample_adata(adata: ad.AnnData, sample_source: str, number: str) -> ad.AnnData:
    """
    sample the cells or genes from adata

    Parameters
    ----------
    adata
      the anndata object
    sample_source
      which dimension to sample from (cells or genes)
    number
      the number of cells or samples to sample
    Returns
    -------
    adata
     the sampled anndata
    """
    if sample_source != '' and number != '':
        if sample_source == 'cells':
            adata = sc.pp.subsample(adata, n_obs=int(number), random_state=0, copy=True)
        elif sample_source == 'genes':
            np.random.seed(base_cfg.random_seed)
            adata = adata[:, np.random.choice(adata.n_vars, size=int(number), replace=False)]
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
    """
    load and store marker genes (if the dataset has) in adata.

    Parameters
    ----------
    adata
      the anndata object
    Returns
    -------
    None
    """
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


def find_rare_cell_types(adata: ad.AnnData):
    cell_type_counts = adata.obs['celltype'].value_counts(ascending=True)
    cell_type_ratio = cell_type_counts / adata.n_obs
    rare_types = np.intersect1d(
        cell_type_counts[cell_type_counts >= data_cfg.rare_number].index,
        cell_type_ratio[cell_type_ratio <= data_cfg.rare_rate].index
    )
    return rare_types if rare_types.size != 0 else None


def set_rare_type(adata: ad.AnnData):
    if 'batch' not in adata.obs:
        rare_types = find_rare_cell_types(adata)
        rare_type = None if rare_types is None else rare_types[0]
    else:
        u_batches = adata.obs['batch'].unique()
        batch_rare_types = [find_rare_cell_types(adata[adata.obs['batch'] == b]) for b in u_batches]
        inter_types = np.intersect1d(*batch_rare_types) if len(batch_rare_types) == 2 else reduce(np.intersect1d, batch_rare_types)
        rare_type = None if inter_types is None else inter_types[0]
    if rare_type is not None:
        adata.uns['rare_type'] = rare_type


def show_data_info(adata: ad.AnnData):
    print(f"Dataset {adata.uns['data_name']} has {adata.n_obs} cells, {adata.n_vars} genes "
          f"and {adata.obs['celltype'].unique().shape[0]} classes after filtering.")
    if 'rare_type' in adata.uns:
        print(f"Rare cell type (> {data_cfg.rare_number} cells and < "
              f"{data_cfg.rare_rate * 100}% of all cells): {adata.uns['rare_type']}")
    if 'batch' not in adata.obs:
        print(f"Data complexity is {np.round(adata.uns['data_complexity'], 3)}.")
    else:
        ubatches = adata.obs['batch'].unique().categories.to_numpy()
        print(f"Dataset contains {ubatches.shape[0]} batches: {ubatches}")
