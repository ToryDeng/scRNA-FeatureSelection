import os
import numpy as np
import pandas as pd
import anndata as ad
from typing import Union
from config.datasets_config import data_cfg


def save_genes(adata: ad.AnnData, method: str, selected_df: pd.DataFrame):
    """
    save selected genes to disk.

    Parameters
    ----------
    adata
      anndata object
    method
      the feature selection method
    selected_df
     selected genes (and their importances)
    Returns
    -------
    None
    """
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    path = os.path.join(data_cfg.cache_path, 'geneData', division_type, adata.uns['data_name'], method)
    if not os.path.exists(path):
        os.makedirs(path)
    selected_df.to_csv(os.path.join(path, f'{selected_df.shape[0]}-genes.csv'))
    print(f"{selected_df.shape[0]} Genes have been saved to {path}/")


def is_saved(adata: ad.AnnData, method: str, n_genes: int):
    """
    check whether genes selected by certain method have been saved.

    Parameters
    ----------
    adata
      anndata object
    method
      the feature selection method
    n_genes
      number of genes to be selected
    Returns
    -------
    bool
      True represents saved. False represents unsaved.
    """
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    path = os.path.join(data_cfg.cache_path, 'geneData', division_type, adata.uns['data_name'], method, f'{n_genes}-genes.csv')
    return os.path.exists(path)


def load_genes(adata: ad.AnnData, method: str, n_genes: int):
    """
    load previously saved genes from disk

    Parameters
    ----------
    adata
      anndata object
    method
      the feature selection method
    n_genes
      number of genes to be selected
    Returns
    -------
    selection_df
      A dataframe. The first column contains gene names, and the second column contains gene importances.
    """
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    path = os.path.join(data_cfg.cache_path, 'geneData', division_type, adata.uns['data_name'], method, f'{n_genes}-genes.csv')
    selection_df = pd.read_csv(path)
    return selection_df


def subset_adata(adata: ad.AnnData, selected_genes: Union[np.ndarray, pd.Index], inplace=False):
    if isinstance(selected_genes, pd.Index):
        selected_genes = selected_genes.to_numpy()
    gene_mask = adata.var_names.isin(selected_genes)
    if inplace:
        if adata.raw is not None:
            adata.raw = adata.raw[:, gene_mask].to_adata()
        adata._inplace_subset_var(selected_genes)
        if adata.raw.shape[1] != selected_genes.shape[0] or adata.shape[1] != selected_genes.shape[0]:
            raise RuntimeError(
                f"{adata.raw.shape[1]} genes in raw data and {adata.shape[1]} in norm data were selected,"
                f" not {selected_genes.shape[0]} genes. Please check the gene names."
            )
    else:
        copied_adata = adata.copy()
        subset_adata(copied_adata, selected_genes, inplace=True)
        return copied_adata