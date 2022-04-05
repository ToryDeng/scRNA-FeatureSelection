import os
import sys
import warnings
import time
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
from functools import wraps
from typing import Optional, Literal
from itertools import permutations

from experiments.recorders import TimeRecorder

matplotlib.use('agg')
import scanpy as sc


def complexity(adata: ad.AnnData, use_raw: bool = False):
    """
    Compute the complexity of the dataset.
    """
    if use_raw:
        adata = adata.raw.to_adata()
    all_centroids = pd.concat(
        objs=[adata[adata.obs.celltype == t].to_df().mean(axis=0).rename(t) for t in adata.obs.celltype.unique()],
        axis=1).T
    corr_mtx = np.corrcoef(all_centroids)
    return np.mean([corr_mtx[i, j] for i, j in enumerate(np.argsort(corr_mtx)[:, -2])])


def head(name: str, fold='', head_len=65):
    """
    formatted printing

    :param name: head name
    :param fold: i-th fold in classification
    :param head_len: length of head, including stars
    :return: None
    """
    title = ' ' + ' '.join([name, '-', str(fold)]) + ' ' if fold != '' else ' ' + name + ' '
    if head_len - len(title) <= 0:
        warnings.warn("title is too long.")
    print(title.center(head_len, '*'))


def plot_tsne(combined_adata: ad.AnnData, figure_name:str, **kwargs):
    sc.tl.tsne(combined_adata, use_rep='X_pca')
    sc.pl.tsne(combined_adata, save=figure_name, **kwargs)


def plot_umap(combined_adata: ad.AnnData, figure_name:str, **kwargs):
    sc.pp.neighbors(combined_adata, use_rep='X_pca')
    sc.tl.umap(combined_adata)
    sc.pl.umap(combined_adata, save=figure_name, **kwargs)


def plot_2D(combined_adata: ad.AnnData,
            fs_method: Optional[str] = None,
            bc_method: Optional[str] = None,
            mode: Literal['before', 'after'] = 'before',
            order: Optional[Literal['correction_first', 'selection_first']] = None
            ):
    assert mode in ('before', 'after'), "mode must be 'before' or 'after'."

    top_celltypes = combined_adata.obs['celltype'].value_counts(ascending=False).index.to_numpy()[:10]
    combined_adata.obs['celltype_to_plot'] = pd.Categorical(
        np.where(np.isin(combined_adata.obs['celltype'], top_celltypes), combined_adata.obs['celltype'], 'Others')
    )
    params = {'color': ['batch', 'celltype_to_plot'], 'hspace': 0.15, 'wspace': 0.15, 'frameon': False,
              'palette': sc.pl.palettes.vega_20_scanpy, 'show': False, 'legend_fontsize': 'xx-small'}

    if mode == 'before':
        fig_name = f": {combined_adata.uns['data_name']}-{mode}.png"
    else:
        fig_name = f": {combined_adata.uns['data_name']}-{mode}-{combined_adata.shape[1]}-{fs_method}-{bc_method}-{order}.png "
    plot_tsne(combined_adata, fig_name, **params)
    plot_umap(combined_adata, fig_name, **params)





def save_genes(adata: ad.AnnData, method: str, selected_df: pd.DataFrame):
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    path = os.path.join('geneData', division_type, adata.uns['data_name'], method)
    if not os.path.exists(path):
        os.makedirs(path)
    selected_df.to_csv(os.path.join(path, f'{selected_df.shape[0]}-genes.csv'))
    print(f"{selected_df.shape[0]} Genes have been saved to {path}/")


def is_saved(adata: ad.AnnData, method: str, n_genes: int):
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    path = os.path.join('geneData', division_type, adata.uns['data_name'], method, f'{n_genes}-genes.csv')
    return os.path.exists(path)


def load_genes(adata: ad.AnnData, method: str, n_genes: int):
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    path = os.path.join('geneData', division_type, adata.uns['data_name'], method, f'{n_genes}-genes.csv')
    selection_df = pd.read_csv(path)
    return selection_df


def record_time(recorder: TimeRecorder):
    def timeit(func):
        @wraps(func)
        def func_timer(*args, **kwargs):
            t0 = time.perf_counter()
            result = func(*args, **kwargs)
            t1 = time.perf_counter()
            print(t1 - t0, 'seconds')
            for arg in args:
                if isinstance(arg, ad.AnnData):
                    print(arg.uns['data_name'])
                if isinstance(arg, str):
                    print(arg)
            return result
        return func_timer

    return timeit


class MyGroupSplit:
    def split(self, X: np.ndarray, y: np.ndarray, groups: np.ndarray):
        unique_groups = np.unique(groups)
        indices = np.arange(X.shape[0])
        for train_grp, test_grp in permutations(unique_groups, 2):
            train_mask, test_mask = groups == train_grp, groups == test_grp
            train_index, test_index = indices[train_mask], indices[test_mask]
            yield train_index, test_index


class HiddenPrints:
    """
    Hide prints

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
