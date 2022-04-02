import os
import sys
import warnings
import datetime
from itertools import permutations
import anndata as ad
import scib
import numpy as np
import pandas as pd
import matplotlib
from harmonypy import compute_lisi

matplotlib.use('agg')
import scanpy as sc
from config import exp_cfg


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


def filter_adata(adata: ad.AnnData, filter_gene: bool = False, filter_cell: bool = False):
    if filter_gene:
        gene_mask = sc.pp.filter_genes(adata, min_cells=exp_cfg.n_filter_cell, inplace=False)[0]
        adata = adata[:, gene_mask]
        if adata.raw is not None:
            adata.raw = adata.raw[:, gene_mask].to_adata()
        if adata.shape[1] < 10:
            warnings.warn(f"There are (is) only {adata.shape[1]} genes after filtering.")
    if filter_cell:
        cell_mask = sc.pp.filter_cells(adata, min_genes=exp_cfg.n_filter_gene, inplace=False)[0]
        adata = adata[cell_mask, :]
        if adata.shape[0] < 10:
            warnings.warn(f"There are (is) only {adata.shape[0]} cells after filtering.")
    if not filter_gene and not filter_cell:
        warnings.warn("function 'filter_adata' is actually not used.", RuntimeWarning)
    return adata


def plot_2D(combined_adata: ad.AnnData,
            dr_method: str = 'umap',
            fs_method: str = None,
            bc_method: str = None,
            mode: str = 'before',
            order: str = None,
            ):
    top_celltypes = combined_adata.obs['celltype'].value_counts(ascending=False).index.to_numpy()[:10]
    combined_adata.obs['celltype_to_plot'] = pd.Categorical(
        np.where(np.isin(combined_adata.obs['celltype'], top_celltypes), combined_adata.obs['celltype'], 'Others')
    )

    if mode == 'before':
        if dr_method == 'umap':
            sc.pp.neighbors(combined_adata, use_rep='X_pca')
            sc.tl.umap(combined_adata)
            sc.pl.umap(combined_adata,
                       color=['batch', 'celltype_to_plot'],
                       hspace=0.15,
                       wspace=0.15,
                       frameon=False,
                       palette=sc.pl.palettes.vega_20_scanpy,
                       save=f": {combined_adata.uns['data_name']}-{mode}.png",
                       show=False,
                       legend_fontsize='xx-small')
        elif dr_method == 'tsne':
            sc.tl.tsne(combined_adata, use_rep='X_pca')
            sc.pl.tsne(combined_adata,
                       color=['batch', 'celltype_to_plot'],
                       hspace=0.15,
                       wspace=0.15,
                       frameon=False,
                       palette=sc.pl.palettes.vega_20_scanpy,
                       save=f": {combined_adata.uns['data_name']}-{mode}.png",
                       show=False,
                       legend_fontsize='xx-small')
        else:
            raise NotImplementedError(f"{dr_method} have not been implemented yet.")

    elif mode == 'after':
        if dr_method == 'umap':
            sc.pp.neighbors(combined_adata, use_rep='X_pca')
            sc.tl.umap(combined_adata)
            sc.pl.umap(combined_adata,
                       color=['batch', 'celltype_to_plot'],
                       hspace=0.15,
                       wspace=0.15,
                       frameon=False,
                       palette=sc.pl.palettes.vega_20_scanpy,
                       save=f": {combined_adata.uns['data_name']}-{combined_adata.shape[1]}-{fs_method}-{bc_method}-{mode}-{order}.png",
                       show=False,
                       legend_fontsize='xx-small')
        elif dr_method == 'tsne':
            sc.tl.tsne(combined_adata, use_rep='X_pca')
            sc.pl.tsne(combined_adata,
                       color=['batch', 'celltype_to_plot'],
                       hspace=0.15,
                       wspace=0.15,
                       frameon=False,
                       palette=sc.pl.palettes.vega_20_scanpy,
                       save=f": {combined_adata.uns['data_name']}-{combined_adata.shape[1]}-{fs_method}-{bc_method}-{mode}-{order}.png",
                       show=False,
                       legend_fontsize='xx-small')
        else:
            raise NotImplementedError(f"{dr_method} have not been implemented yet. Must be umap or tsne.")
    else:
        raise ValueError(f"mode must be 'before' or 'after'.")


def now(lag_hours=0) -> str:
    return (datetime.datetime.now() + datetime.timedelta(hours=lag_hours)).strftime('%Y-%m-%d %H:%M:%S')


def compute_correction_metrics(adata: ad.AnnData) -> dict:
    correction_result = dict()

    # with HiddenPrints():
    kbet = scib.metrics.kBET(adata, batch_key='batch', label_key='celltype', embed=f'X_pca')
    cells_lisi = compute_lisi(adata.obsm[f'X_pca'], adata.obs, ['batch', 'celltype'])
    norm_lisi = (cells_lisi - cells_lisi.min(axis=0)) / (cells_lisi.max(axis=0) - cells_lisi.min(axis=0))
    median_lisi = np.median(norm_lisi, axis=0)

    correction_result[f'seuratV4_kBET'] = kbet
    correction_result[f'seuratV4_iLISI'] = median_lisi[0]
    correction_result[f'seuratV4_cLISI'] = median_lisi[1]
    correction_result[f'seuratV4_f1LISI'] = 2 * (1 - median_lisi[1]) * median_lisi[0] / (
            1 - median_lisi[1] + median_lisi[0])

    return correction_result


def save_genes(adata: ad.AnnData, method: str, selected_df: pd.DataFrame):
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'

    if 'geneData' not in os.listdir():
        os.mkdir('geneData')
        save_genes(adata, method, selected_df)
    else:
        if division_type not in os.listdir(os.path.join('geneData')):
            os.mkdir(os.path.join('geneData', division_type))
            save_genes(adata, method, selected_df)
        else:
            if adata.uns['data_name'] not in os.listdir(os.path.join('geneData', division_type)):
                os.mkdir(os.path.join('geneData', division_type, adata.uns['data_name']))
                save_genes(adata, method, selected_df)
            else:
                if method not in os.listdir(os.path.join('geneData', division_type, adata.uns['data_name'])):
                    os.mkdir(os.path.join('geneData', division_type, adata.uns['data_name'], method))
                    save_genes(adata, method, selected_df)
                else:
                    if selected_df.shape[0] not in exp_cfg.n_genes:
                        print(f"Saving {selected_df.shape[0]} genes! This result will be dismissed.")
                    else:
                        selected_df.to_csv(
                            os.path.join('geneData', division_type, adata.uns['data_name'], method, f'{selected_df.shape[0]}-genes.csv')
                        )
    print(f"Genes have been saved to {os.path.join('geneData', division_type, adata.uns['data_name'], method)}.")


def is_saved(adata: ad.AnnData, method: str, n_genes: int):
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    path = os.path.join('geneData', division_type, adata.uns['data_name'], method, f'{n_genes}-genes.csv')
    return os.path.exists(path)


def load_genes(adata: ad.AnnData, method: str, n_genes: int):
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    selection_df = pd.read_csv(
        os.path.join('geneData', division_type, adata.uns['data_name'], method, f'{n_genes}-genes.csv')
    )
    return selection_df


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
