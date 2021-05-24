import os
import re
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
    regex = re.compile(pattern="[-_:+]")
    vreplace = np.vectorize(lambda x: regex.sub('.', x), otypes=[np.object])
    return vreplace(gene_names)


def load_data(data_name: str) -> ad.AnnData:
    """
    For specific data_name, get the corresponding raw data. Remove constant genes (=0) and normalize it
    using the method in Seurat.

    :param data_name: the dataset name you want to get
    :return:  anndata object with raw_data, norm_data and markers in data
    """
    sc.settings.verbosity = 1
    if data_name[:4] == 'PBMC':
        assert len(data_name) >= 4, ValueError("parameter 'data_name' is wrong!")
        data = pd.read_hdf(data_cfg.PBMC_path, key="AllCells")
        if len(data_name) > 4:
            if data_name[4:-1].isdigit() and data_name[-1] == '%':
                rate = float(data_name[4:-1])
                assert 0 < rate < 100, ValueError("The percentage is wrong!")
                X, y = data.iloc[:, :-1], data.iloc[:, -1]
                _, raw_features, __, labels = train_test_split(X, y, test_size=rate / 100,
                                                               random_state=exp_cfg.random_seed, stratify=y)
                data = pd.concat([raw_features, labels], axis=1)
            elif data_name[4:].isdigit():
                n_cells = int(data_name[4:])
                assert 0 < n_cells < data.shape[0], ValueError("The cell num is wrong!")
                X, y = data.iloc[:, :-1], data.iloc[:, -1]
                _, raw_features, __, labels = train_test_split(X, y, test_size=n_cells,
                                                               random_state=exp_cfg.random_seed, stratify=y)
                data = pd.concat([raw_features, labels], axis=1)
            else:
                raise ValueError("Wrong data name.")
        os.chdir(data_cfg.PBMC_markers_path)
        part1 = np.loadtxt('hsPBMC_markers_10x.txt', skiprows=1, usecols=[0], dtype=np.str_, delimiter=',')
        part2 = np.loadtxt('blood_norm_marker.txt', skiprows=1, usecols=[1], dtype=np.str_, delimiter=',')
        part3 = np.unique(np.loadtxt('CellMarker.csv', skiprows=1, usecols=[1], dtype=np.str_, delimiter=','))
        markers = np.union1d(np.union1d(part1, part2), part3)

    elif data_name in ['muraro', 'segerstolpe', 'xin']:
        data = pd.read_hdf(data_cfg.pancreas_path, key=data_name.capitalize() + '_pancreas')
        data = data.loc[~data.iloc[:, -1].isin(data_cfg.pancreas_remove_types).values]
        os.chdir(data_cfg.pancreas_markers_path)
        markers = np.loadtxt('pancreasMarkerGenes.csv', skiprows=1, usecols=[0], dtype=np.str_, delimiter=',')
    else:
        raise ValueError(f"data name {data_name} is wrong!")
    os.chdir("../../scRNA-FeatureSelection")  # return to scRNA-FeatureSelection dir
    count_matrix, labels = data.iloc[:, :-1].round(), data.iloc[:, -1]
    count_matrix.columns, labels.name = get_gene_names(count_matrix.columns), "type"

    adata = ad.AnnData(X=count_matrix, obs=labels.to_frame())
    adata.obs_names_make_unique(join='.')
    adata.var_names_make_unique(join='.')
    adata.uns['markers'] = np.intersect1d(markers, adata.var_names)  # only contains existing genes
    adata.uns['data_name'] = data_name
    # filter almostly constant cells and genes
    sc.pp.filter_genes(data=adata, min_cells=exp_cfg.n_filter_cell, inplace=True)
    sc.pp.filter_cells(data=adata, min_genes=exp_cfg.n_filter_gene, inplace=True)
    # store raw data
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=exp_cfg.scale_factor, inplace=True)
    sc.pp.log1p(adata)
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


def save_data(adata_train: Optional[ad.AnnData] = None, adata_test: Optional[ad.AnnData] = None,
              task: str = 'assign') -> None:
    """
    Save raw data in tempData dir for classification task and cluster task.

    :param adata_train: training cells when task == 'assign', or count matrix of all cells when task == 'cluster'.
    :param adata_test: test cells
    :param task: 'assign' or 'cluster'
    :return: None
    """
    assert task in ['assign', 'cluster'], "Parameter 'task' is wrong."
    if task == 'assign':
        if adata_train is not None:
            adata_train.to_df().to_csv('tempData/temp_X_train.csv')
            adata_train.obs['type'].to_csv('tempData/temp_y_train.csv')
        if adata_test is not None:
            adata_test.to_df().to_csv('tempData/temp_X_test.csv')
            adata_test.obs['type'].to_csv('tempData/temp_y_test.csv')
        if adata_train is None and adata_test is None:
            raise ValueError("adata_train and adata_test are both None.")
    else:
        adata_train.to_df().to_csv('tempData/temp_X.csv')
        adata_train.obs['type'].to_csv('tempData/temp_y.csv')


def head(name: str, fold='', head_len=65):
    """
    formatted printing

    :param name: head name
    :param fold: i-th fold
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
    if filter_cell:
        cell_mask = sc.pp.filter_cells(adata, min_genes=exp_cfg.n_filter_gene, inplace=False)[0]
        adata = adata[cell_mask, :]
    if not filter_gene and not filter_cell:
        warnings.warn("function 'filter_adata' is actually not used.", RuntimeWarning)
    return adata


def normalize_importances(importances:np.ndarray) -> np.ndarray:
    min_importance, max_importance = np.min(importances), np.max(importances)
    return (importances - min_importance) / (max_importance - min_importance)
