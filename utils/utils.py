import os
import re
import sys
import warnings
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.model_selection import train_test_split

from config import data_cfg, exp_cfg


def get_gene_names(columns):
    """
    Get gene names array from features (dataframe).

    :param columns: dataframe.columns
    :return: an array contains gene names
    """
    if '__' in columns[0] and columns[0][0] != '"' and columns[0][-1] != '"':
        gene_names = pd.Series(columns).str.split('__', expand=True).iloc[:, 0].values
    elif '__' not in columns[0] and columns[0][0] == '"' and columns[0][-1] == '"':
        gene_names = np.char.strip(np.array(columns, dtype=np.str_), '"')
    elif '__' in columns[0] and columns[0][0] == '"' and columns[0][-1] == '"':
        gene_names = np.char.strip(np.array(columns, dtype=np.str_), '"')
        gene_names = pd.Series(gene_names).str.split('__', expand=True).iloc[:, 0].values
    elif '\t' in columns[0]:
        gene_names = pd.Series(columns).str.split('\t', expand=True).iloc[:, 1].values
    else:
        gene_names = np.array(columns, dtype=np.str_)
    regex = re.compile(pattern="[-_:+()|]")
    vreplace = np.vectorize(lambda x: regex.sub('.', x), otypes=[np.object])
    return vreplace(gene_names)


def load_data(data_name: str) -> ad.AnnData:
    """
    For specific data_name, get the corresponding raw data. Remove constant genes (=0) and normalize it
    using the method in Seurat.

    :param data_name: the dataset name you want to get
    :return: anndata object with raw_data, norm_data and markers in data, or the concatenated data.
    """
    sc.settings.verbosity = 1
    markers = None
    dataset, n_samples = re.match(r"([a-zA-Z]*)([0-9]*)", data_name).groups()

    project_path = os.getcwd()
    if '+' not in data_name:
        if dataset == 'PBMC':
            adata = ad.read_h5ad(data_cfg.data_path + 'ZhengPBMC92k.h5ad')
            os.chdir(data_cfg.PBMC_markers_path)
            panglao = np.loadtxt('PBMC_panglaoDB.csv', skiprows=1, usecols=[0], dtype=np.str, delimiter=',')
            cellmarker = np.loadtxt('PBMC_CellMarker.csv', skiprows=1, usecols=[1], dtype=np.str, delimiter=',')
            markers = np.union1d(panglao, cellmarker)
        elif dataset == 'PBMCsmall':
            adata = ad.read_h5ad(data_cfg.data_path + f'{dataset}.h5ad')
            os.chdir(data_cfg.PBMC_markers_path)
            panglao = np.loadtxt('PBMC_panglaoDB.csv', skiprows=1, usecols=[0], dtype=np.str, delimiter=',')
            cellmarker = np.loadtxt('PBMC_CellMarker.csv', skiprows=1, usecols=[1], dtype=np.str, delimiter=',')
            markers = np.union1d(panglao, cellmarker)
        elif dataset in ['muraro', 'segerstolpe', 'xin', 'baron']:
            adata = ad.read_h5ad(data_cfg.data_path + f"{dataset.capitalize() + 'HumanPancreas'}.h5ad")
            os.chdir(data_cfg.pancreas_markers_path)
            panglao = np.loadtxt('pancreas_panglaoDB.csv', skiprows=1, usecols=[0], dtype=np.str, delimiter=',')
            cellmarker = np.loadtxt('pancreas_CellMarker.csv', skiprows=1, usecols=[0], dtype=np.str_, delimiter=',')
            markers = np.union1d(panglao, cellmarker)
        elif dataset in ['ZilionisMouseLungCancer', 'AztekinTail', 'MarquesMouseBrain', 'ZeiselMouseBrain',
                         'BaronMousePancreas', 'LaMannoMouseAdult', 'LaMannoHumanEmbryonicMidbrain',
                         'VentoHumanPlacenta', 'LaMannoHumanEmbryonicStem', 'DengMouseEmbryoDevel',
                         'GoolamMouseEmbryoDevel', 'GuoHumanTestis', 'QuakeMouseHeart']:
            adata = ad.read_h5ad(data_cfg.data_path + f'{dataset}.h5ad')
        else:
            raise ValueError(f"data name {data_name} is wrong!")

        os.chdir(project_path)  # return to project dir

        adata.var_names = get_gene_names(adata.var_names)
        adata.obs.rename(columns={adata.obs.columns[0]: 'celltype'}, inplace=True)
        adata.obs_names_make_unique(join='.')
        adata.var_names_make_unique(join='.')

        if markers is not None:
            adata.uns['markers'] = np.intersect1d(markers, adata.var_names)  # only contains existing genes
        if dataset == 'VentoHumanPlacenta':  # at least 2 samples are in dataset, for computation time test
            adata = adata[adata.obs['celltype'].isin(adata.obs['celltype'].value_counts().index[:12].to_numpy()), :]
        head(data_name, head_len=100)
        adata.uns['data_name'] = data_name

        # filter almostly constant cells and genes
        sc.pp.filter_genes(data=adata, min_cells=exp_cfg.n_filter_cell, inplace=True)
        sc.pp.filter_cells(data=adata, min_genes=exp_cfg.n_filter_gene, inplace=True)
        # filter unclear cells
        adata = adata[~adata.obs['celltype'].isin(data_cfg.remove_types)]

        # filter mito, RP, ERCC
        adata.var['mt'] = adata.var_names.str.match(r'^MT-')
        adata.var['rp'] = adata.var_names.str.match(r'^RP[SL0-9]')
        adata.var['ERCC'] = adata.var_names.str.match(r'^ERCC-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'rp', 'ERCC'], percent_top=None, log1p=False, inplace=True)
        adata = adata[:, (~adata.var['mt']) & (~adata.var['rp']) & (~adata.var['ERCC'])]
        adata = adata[adata.obs.pct_counts_mt < 5, :]  # Cumulative percentage of counts for mito genes
        adata = adata[adata.obs.n_genes_by_counts < adata.obs.n_genes_by_counts.quantile(0.95), :]
        # sample data
        if n_samples != '':
            sc.pp.subsample(adata, n_obs=int(n_samples), random_state=exp_cfg.random_seed)
        # store raw data
        adata.raw = adata
        # normalize data
        sc.pp.normalize_total(adata, target_sum=exp_cfg.scale_factor, inplace=True)
        sc.pp.log1p(adata)
    else:
        batch_names = data_name.split('+')
        batch_data = [load_data(data_name) for data_name in batch_names]
        adata = ad.concat(batch_data, join="inner", keys=batch_names, label='batch')
        sc.pp.pca(adata)

        head(data_name, head_len=100)
        adata.uns['data_name'] = data_name

    return adata


def delete(path: str):
    """
    Recursively delete files in a temporary folder if 'path' is a dir, or delete the file when 'path'
    is a file path.

    :param path: folder or file path
    :return: None
    """
    if os.path.isdir(path):
        for file in os.listdir(path):
            p = os.path.join(path, file)
            delete(p) if os.path.isdir(p) else os.remove(p)  # recursion
    else:
        os.remove(path)

def save_data(adata_train: Optional[ad.AnnData] = None,
              adata_test: Optional[ad.AnnData] = None,
              use_rep: str = 'celltype',
              task: str = 'assign') -> None:
    """
    Save raw data in tempData dir for classification task and cluster task.

    :param adata_train: training cells when task == 'assign', or count matrix of all cells when task == 'cluster'.
    :param adata_test: test cells
    :param task: 'assign' or 'cluster'
    :return: None
    """
    assert task in ['assign', 'cluster', 'correct'], "Parameter 'task' is wrong."
    if task == 'assign':
        if adata_train is not None:
            adata_train.to_df().reset_index().to_feather('tempData/temp_X_train.feather')
            adata_train.obs[use_rep].to_frame().reset_index().to_feather('tempData/temp_y_train.feather')
        if adata_test is not None:
            adata_test.to_df().reset_index().to_feather('tempData/temp_X_test.feather')
            adata_test.obs[use_rep].to_frame().reset_index().to_feather('tempData/temp_y_test.feather')
        if adata_train is None and adata_test is None:
            raise ValueError("adata_train and adata_test are both None.")
    elif task == 'cluster':
        adata_train.to_df().reset_index().to_feather('tempData/temp_X.feather')
        adata_train.obs[use_rep].to_frame().reset_index().to_feather('tempData/temp_y.feather')
    else:  # batch correction
        if 'X_pca_harmony' in adata_train.obsm:
            latent_rep = pd.DataFrame(adata_train.obsm['X_pca_harmony'], index=adata_train.obs_names).reset_index()
        elif 'X_scanorama' in adata_train.obsm:
            latent_rep = pd.DataFrame(adata_train.obsm['X_scanorama'], index=adata_train.obs_names).reset_index()
        else:  # corrected data is in adata.X (scGen)
            latent_rep = pd.DataFrame(adata_train.X, index=adata_train.obs_names).reset_index()
        latent_rep.columns = latent_rep.columns.astype(np.str)
        latent_rep.to_feather('tempData/temp_X.feather')
        adata_train.obs['batch'].to_frame().reset_index().to_feather('tempData/temp_y.feather')


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
            data_name: str = None,
            mode: str = 'before'
            ):
    top_celltypes = combined_adata.obs['celltype'].value_counts(ascending=False).index.to_numpy()[:13]
    combined_adata.obs['celltype_to_plot'] = np.where(
        np.isin(combined_adata.obs['celltype'], top_celltypes),
        combined_adata.obs['celltype'],
        'Others'
    )

    if mode == 'before':
        if 'X_pca' not in combined_adata.obsm:
            sc.pp.pca(combined_adata)
        if dr_method == 'umap':
            sc.pp.neighbors(combined_adata, use_rep='X_pca')
            sc.tl.umap(combined_adata)
            sc.pl.umap(combined_adata,
                       color=['batch', 'celltype_to_plot'],
                       hspace=0.08,
                       wspace=0.01,
                       frameon=False,
                       palette=sc.pl.palettes.vega_20_scanpy,
                       save=f': {data_name}-{fs_method}-' + '+'.join(combined_adata.obs['batch'].unique().tolist()) + f'-{mode}.png', show=False)
        elif dr_method == 'tsne':
            sc.tl.tsne(combined_adata, use_rep='X_pca')
            sc.pl.tsne(combined_adata,
                       color=['batch', 'celltype_to_plot'],
                       hspace=0.08,
                       wspace=0.01,
                       frameon=False,
                       palette=sc.pl.palettes.vega_20_scanpy,
                       save=f': {data_name}-{fs_method}-' + '+'.join(combined_adata.obs['batch'].unique().tolist()) + f'-{mode}.png',
                       show=False)

    elif mode == 'after':
        if bc_method == 'harmony':
            representation = 'X_pca_harmony'
            print(bc_method, combined_adata.obsm[representation].shape)
        elif bc_method == 'scanorama':
            representation = 'X_scanorama'
            print(bc_method, combined_adata.obsm[representation].shape)
        elif bc_method == 'scgen':
            sc.pp.pca(combined_adata)  # execute pca anyway
            representation = 'X_pca'  # scGen
            print(bc_method, combined_adata.obsm[representation].shape)
        else:
            raise AttributeError(f"{bc_method} is wrong.")
        print(f"Use {representation} to plot.")

        if dr_method == 'umap':
            sc.pp.neighbors(combined_adata, use_rep=representation)
            sc.tl.umap(combined_adata)
            sc.pl.umap(combined_adata,
                       color=['batch', 'celltype_to_plot'],
                       hspace=0.08,
                       wspace=0.01,
                       frameon=False,
                       palette=sc.pl.palettes.vega_20_scanpy,
                       save=f': {data_name}-{fs_method}-{bc_method}-{combined_adata.shape[1]}-' +
                            '+'.join(combined_adata.obs['batch'].unique().tolist()) + f'-{mode}.png',
                       show=False)
        elif dr_method == 'tsne':
            sc.tl.tsne(combined_adata, use_rep=representation)
            sc.pl.tsne(combined_adata,
                       color=['batch', 'celltype_to_plot'],
                       hspace=0.08,
                       wspace=0.01,
                       frameon=False,
                       palette=sc.pl.palettes.vega_20_scanpy,
                       save=f': {data_name}-{fs_method}-{bc_method}-{combined_adata.shape[1]}-' +
                            '+'.join(combined_adata.obs['batch'].unique().tolist()) + f'-{mode}.png',
                       show=False)
        else:
            raise NotImplementedError(f"{dr_method} have not been implemented yet.")
    else:
        raise ValueError(f"mode must be 'before' or 'after'.")


def normalize_importances(importances: np.ndarray) -> np.ndarray:
    min_importance, max_importance = np.min(importances), np.max(importances)
    return (importances - min_importance) / (max_importance - min_importance)


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
