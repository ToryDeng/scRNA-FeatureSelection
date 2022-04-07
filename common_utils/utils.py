import os
import sys
import warnings
from typing import Optional, Literal

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('agg')
import scanpy as sc


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
        objs=[adata[adata.obs.celltype == t].to_df().mean(axis=0).rename(t) for t in adata.obs.celltype.unique()],
        axis=1).T
    corr_mtx = np.corrcoef(all_centroids)
    return np.mean([corr_mtx[i, j] for i, j in enumerate(np.argsort(corr_mtx)[:, -2])])


def head(name: str, n_times: int = None, head_len=65):
    """
    formatted printing

    Parameters
    ----------
    name
      head name
    n_times
      the number after the head name
    head_len
      length of head, including stars and spaces

    Returns
    -------
    None
    """
    title = ' ' + ' '.join([name, '-', str(n_times)]) + ' ' if n_times is not None else ' ' + name + ' '
    if head_len - len(title) <= 0:
        warnings.warn("title is too long.")
    print(title.center(head_len, '*'))


def plot_tsne(combined_adata: ad.AnnData, figure_name: str, **kwargs):
    sc.tl.tsne(combined_adata, use_rep='X_pca')
    sc.pl.tsne(combined_adata, save=figure_name, **kwargs)


def plot_umap(combined_adata: ad.AnnData, figure_name: str, **kwargs):
    sc.pp.neighbors(combined_adata, use_rep='X_pca')
    sc.tl.umap(combined_adata)
    sc.pl.umap(combined_adata, save=figure_name, **kwargs)


def plot_2D(combined_adata: ad.AnnData,
            fs_method: Optional[str] = None,
            bc_method: Optional[str] = None,
            mode: Literal['before', 'after'] = 'before',
            order: Optional[Literal['correction_first', 'selection_first']] = None
            ):
    """
    plot the umap and t-sne representations before and after feature selection and batch correction.

    Parameters
    ----------
    combined_adata
      the anndata object containing multiple batches
    fs_method
      the feature selection method
    bc_method
      the batch correction method
    mode
      'before' or 'after'
    order
      'correction_first' or 'selection_first'

    Returns
    -------
    None
    """
    assert mode in ('before', 'after'), "mode must be 'before' or 'after'."

    top_celltypes = combined_adata.obs['celltype'].value_counts(ascending=False).index.to_numpy()[:10]
    combined_adata.obs['celltype_to_plot'] = pd.Categorical(
        np.where(np.isin(combined_adata.obs['celltype'], top_celltypes), combined_adata.obs['celltype'], 'Others')
    )
    params = {'color': ['batch', 'celltype_to_plot'], 'hspace': 0.15, 'wspace': 0.15, 'frameon': False,
              'palette': sc.pl.palettes.vega_20_scanpy, 'show': False, 'legend_fontsize': 'xx-small'}
    sc.pp.pca(combined_adata)  # do PCA anyway
    if mode == 'before':
        fig_name = f" {combined_adata.uns['data_name']}-{mode}.png"
    else:
        fig_name = f" {combined_adata.uns['data_name']}-{mode}-{combined_adata.shape[1]}-{fs_method}-{bc_method}-{order}.png"
    plot_tsne(combined_adata, fig_name, **params)
    plot_umap(combined_adata, fig_name, **params)


class HiddenPrints:
    """
    Hide prints from terminal
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr
