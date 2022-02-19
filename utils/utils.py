import os
import re
import sys
import warnings
import datetime
from itertools import count, permutations
from typing import Optional

import anndata as ad
import scib
import numpy as np
import pandas as pd
import matplotlib
from harmonypy import compute_lisi
from scipy.sparse.csr import csr_matrix

matplotlib.use('agg')
import scanpy as sc
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


def standardize_cell_types(cell_types: np.ndarray) -> np.ndarray:
    regex = re.compile(pattern="[- +']")
    vreplace = np.vectorize(lambda x: regex.sub('.', str(x)), otypes=[np.object])
    vladd = np.vectorize(lambda x: 'X' + x if x[0].isdigit() else x, otypes=[np.object])
    return vladd(vreplace(cell_types))


def load_data(data_name: str,
              is_batch=False,
              min_cells: int = 3,
              min_genes: int = 200,
              ) -> ad.AnnData:
    """
    For specific data_name, get the corresponding raw data. Remove constant genes (=0) and normalize it
    using the method in Seurat.

    :param data_name: the dataset name you want to get
    :param is_batch: whether it is a part of batch
    :param min_cells: minimum cells in filtering genes
    :param min_genes: minimum cells in filtering cells
    :return: anndata object with raw_data, norm_data and markers in data, or the concatenated data.
    """
    sc.settings.verbosity = 0
    panglao, cellmarker, markers = None, None, None
    dataset, number, sample_from = re.match(r"([a-zA-Z]*)([0-9]*)([a-zA-Z]*)?", data_name).groups()

    project_path = os.getcwd()
    if '+' not in data_name:
        if dataset.startswith('PBMC'):
            os.chdir(data_cfg.marker_path)
            panglao = np.unique(np.loadtxt('PBMC_panglaoDB.csv', skiprows=1, usecols=[0], dtype=np.str, delimiter=','))
            cellmarker = np.unique(
                np.loadtxt('PBMC_CellMarker.csv', skiprows=1, usecols=[1], dtype=np.str, delimiter=','))
            markers = np.union1d(panglao, cellmarker)  # either in PanglaoDB or in CellMarker
            if dataset == 'PBMC':
                adata = ad.read_h5ad(data_cfg.data_path + 'ZhengPBMC92k.h5ad')
            elif dataset == 'PBMCsmall':
                adata = ad.read_h5ad(data_cfg.data_path + f'{dataset}.h5ad')
                celltype_convert_dict = {'Naive CD4 T': 'CD4.T.cell', 'Memory CD4 T': 'CD4.T.cell',
                                         'CD14+ Mono': 'Monocyte_CD14', 'B': 'B.cell', 'CD8 T': 'CD8.T.cell',
                                         'FCGR3A+ Mono': 'Monocyte_FCGR3A', 'NK': 'NK.cell'}
                adata.obs['seurat_annotations'] = adata.obs['seurat_annotations'].replace(celltype_convert_dict)
            elif dataset == 'PBMCbatches':
                adata = ad.read_h5ad(data_cfg.data_path + f'{dataset}.h5ad')
            elif dataset == 'PBMCbatchone':
                adata = ad.read_h5ad(data_cfg.data_path + f'PBMCbatches.h5ad')
                adata = adata[adata.obs['Batch'] == 'Batch1']
                del adata.obs['Batch']
            elif dataset == 'PBMCbatchtwo':
                adata = ad.read_h5ad(data_cfg.data_path + f'PBMCbatches.h5ad')
                adata = adata[adata.obs['Batch'] == 'Batch2']
                del adata.obs['Batch']
            else:
                raise AttributeError(f"data name '{dataset}' is wrong.")
        elif dataset == 'simulatingPBMCsmall':
            adata = ad.read_h5ad(data_cfg.data_path + f'{dataset}.h5ad')
            os.chdir(data_cfg.marker_path)
            markers = np.unique(
                np.loadtxt("pbmc3k_depr0.0045_DEGenes.csv", delimiter=',', skiprows=1, dtype=np.object).flatten())
            markers = markers[markers != 'NA']
        elif dataset in ['segerstolpe', 'baron']:
            adata = ad.read_h5ad(data_cfg.data_path + f"{dataset.capitalize() + 'HumanPancreas'}.h5ad")
            os.chdir(data_cfg.marker_path)
            panglao = np.loadtxt('pancreas_panglaoDB.csv', skiprows=1, usecols=[0], dtype=np.str, delimiter=',')
            cellmarker = np.loadtxt('pancreas_CellMarker.csv', skiprows=1, usecols=[0], dtype=np.str_, delimiter=',')
            markers = np.union1d(panglao, cellmarker)
        elif dataset == 'ZeiselMouseBrain':
            adata = ad.read_h5ad(data_cfg.data_path + f'{dataset}.h5ad')
            os.chdir(data_cfg.marker_path)
            panglao = np.loadtxt('MouseBrain_panglaoDB.csv', skiprows=1, usecols=[1], dtype=np.str_, delimiter=',')
            cellmarker = np.loadtxt('MouseBrain_CellMarker.csv', skiprows=1, usecols=[0], dtype=np.str_, delimiter=',')
            markers = np.union1d(panglao, cellmarker)
        elif dataset in ['ZilionisMouseLungCancer', 'AztekinTail', 'MarquesMouseBrain',
                         'BaronMousePancreas', 'LaMannoMouseAdult', 'LaMannoHumanEmbryonicMidbrain',
                         'VentoHumanPlacenta', 'LaMannoHumanEmbryonicStem', 'DengMouseEmbryoDevel',
                         'GuoHumanTestis', 'QuakeMouseHeart', 'QuakeMouseSpleen', 'QuakeMouseTongue',
                         'DarmanisBrain', 'Alles', 'Ariss', 'HaberIntestine', 'ToschesLizard',
                         'MouseAtlas', 'MouseRetina', 'CellLine', 'MouseHematopoieticStemProgenitor']:
            adata = ad.read_h5ad(f'{data_cfg.data_path}{dataset}.h5ad')
        else:
            raise ValueError(f"data name {data_name} is wrong!")

        os.chdir(project_path)  # return to project dir

        adata.var_names = get_gene_names(adata.var_names)
        adata.obs.rename(columns={adata.obs.columns[0]: 'celltype'}, inplace=True)
        if 'Batch' in adata.obs:
            adata.obs.rename(columns={'Batch': 'batch'}, inplace=True)
        adata.obs['celltype'] = pd.Categorical(standardize_cell_types(adata.obs['celltype'].values))
        adata.obs_names_make_unique(join='.')
        adata.var_names_make_unique(join='.')

        if isinstance(adata.X, csr_matrix):
            adata.X = adata.X.toarray()

        if dataset == 'VentoHumanPlacenta':  # at least 2 samples are in dataset, for computation time test
            adata = adata[adata.obs['celltype'].isin(adata.obs['celltype'].value_counts().index[:12].to_numpy()), :]
        if not is_batch:
            head(data_name, head_len=100)
            adata.uns['data_name'] = data_name

        # filter almostly constant cells and genes
        sc.pp.filter_genes(data=adata, min_cells=min_cells, inplace=True)
        sc.pp.filter_cells(data=adata, min_genes=min_genes, inplace=True)
        # filter unclear cells
        adata = adata[~adata.obs['celltype'].isin(data_cfg.remove_types)]
        for replaced, to_replace in data_cfg.replace_types.items():
            adata.obs['celltype'].replace(to_replace, replaced, inplace=True)

        # qc metrics
        adata.var['mito'] = adata.var_names.str.match(r'^MT.')
        adata.var['ERCC'] = adata.var_names.str.match(r'^ERCC.')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True)
        adata = adata[adata.obs.n_genes < adata.obs.n_genes.quantile(0.99), :]
        adata = adata[adata.obs.pct_counts_mito < 5, :]  # Cumulative percentage of counts for mito genes
        adata = adata[:, (~adata.var['ERCC'])]  # removing ERCC (external RNA controls consortium) spike-ins

        # sample data, only on large data sets which do not contain batches
        if number != '':
            if sample_from is not None:
                if sample_from == 'cells':
                    sc.pp.subsample(adata, n_obs=int(number), random_state=0)
                elif sample_from == 'genes':
                    np.random.seed(exp_cfg.random_seed)
                    adata = adata[:, np.random.choice(adata.n_vars, size=int(number))]
                    adata.var_names_make_unique()
                else:
                    raise AttributeError(
                        f"{sample_from} is a wrong source of sample. Please choose from ['cells', 'genes']")
            else:
                raise AttributeError(f"{data_name} did not specify source of sample.")

        # store marker genes and weights
        if markers is not None:  # baron, segerstolpe, PBMC3k, simulated PBMC3k
            adata.uns['markers'] = np.intersect1d(markers,
                                                  adata.var_names)  # only contains existing genes in data
            if panglao is not None and cellmarker is not None:  # only in this case the marker weight exists
                adata.uns['marker_weight'] = np.where(
                    np.isin(adata.uns['markers'], np.intersect1d(panglao, cellmarker)), 2., 1.)
            else:  # simulated PBMC3k
                adata.uns['marker_weight'] = np.ones_like(adata.uns['markers'])

        # store raw data
        adata.raw = adata
        # normalize and scale data
        sc.pp.normalize_total(adata, target_sum=exp_cfg.scale_factor, inplace=True, key_added='counts_per_cell')
        sc.pp.log1p(adata)
        if 'batch' in adata.obs:
            print('scaling data for each batch...')
            for batch in adata.obs['batch'].unique():
                batch_mask = adata.obs['batch'] == batch
                adata.X[batch_mask, :] = sc.pp.scale(adata.X[batch_mask, :], max_value=10, copy=True)
        else:
            sc.pp.scale(adata, max_value=10)
    else:
        head(data_name, head_len=100)

        batch_names = data_name.split('+')
        batch_data = [load_data(data_name, is_batch=True) for data_name in batch_names]
        adata = ad.concat(batch_data, join="inner", keys=batch_names, label='batch')

        adata.uns['data_name'] = data_name

    if 'batch' in adata.obs and not isinstance(adata.obs['batch'], pd.Categorical):
        adata.obs['batch'] = pd.Categorical(adata.obs['batch'])
    print(f"data set {data_name} has {adata.n_obs} cells, {adata.n_vars} genes "
          f"and {adata.obs['celltype'].unique().shape[0]} classes after filtering.")

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
        if os.path.exists(path):
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
    if task == 'assign':
        if adata_train is not None:
            print(f"There are {adata_train.n_obs} cells and {adata_train.n_vars} genes in train data.")
            adata_train.to_df().reset_index().to_feather('tempData/temp_X_train.feather', compression='uncompressed')
            adata_train.obs[use_rep].to_frame().reset_index().to_feather('tempData/temp_y_train.feather',
                                                                         compression='uncompressed')
        if adata_test is not None:
            print(f"There are {adata_test.n_obs} cells and {adata_test.n_vars} genes in test data.")
            adata_test.to_df().reset_index().to_feather('tempData/temp_X_test.feather', compression='uncompressed')
            adata_test.obs[use_rep].to_frame().reset_index().to_feather('tempData/temp_y_test.feather',
                                                                        compression='uncompressed')
        if adata_train is None and adata_test is None:
            raise ValueError("adata_train and adata_test are both None.")
    elif task == 'cluster':
        print(f"There are {adata_train.n_obs} cells and {adata_train.n_vars} genes in data.")
        adata_train.to_df().reset_index().to_feather('tempData/temp_X.feather', compression='uncompressed')
        adata_train.obs[use_rep].to_frame().reset_index().to_feather('tempData/temp_y.feather',
                                                                     compression='uncompressed')
    else:
        raise ValueError(f"Parameter 'task' is wrong. It should be 'assign' or 'cluster'.")


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


def normalize_importances(importances: np.ndarray) -> np.ndarray:
    min_importance, max_importance = np.min(importances), np.max(importances)
    return (importances - min_importance) / (max_importance - min_importance)


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


def save_genes(adata: ad.AnnData, method: str, feature_importance_list: list):
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'

    if 'geneData' not in os.listdir():
        os.mkdir('geneData')
        save_genes(adata, method, feature_importance_list)
    else:
        if division_type not in os.listdir(os.path.join('geneData')):
            os.mkdir(os.path.join('geneData', division_type))
            save_genes(adata, method, feature_importance_list)
        else:
            if adata.uns['data_name'] not in os.listdir(os.path.join('geneData', division_type)):
                os.mkdir(os.path.join('geneData', division_type, adata.uns['data_name']))
                save_genes(adata, method, feature_importance_list)
            else:
                if method not in os.listdir(os.path.join('geneData', division_type, adata.uns['data_name'])):
                    os.mkdir(os.path.join('geneData', division_type, adata.uns['data_name'], method))
                    save_genes(adata, method, feature_importance_list)
                else:
                    for genes, importances in feature_importance_list:
                        if genes.shape[0] not in exp_cfg.n_genes:
                            print(f"Saving genes with shape {genes.shape}! This result will not be saved.")
                        else:
                            np.save(
                                os.path.join('geneData', division_type, adata.uns['data_name'], method,
                                             f'{genes.shape[0]}-genes.npy'),
                                arr=genes
                            )
                            np.save(
                                os.path.join('geneData', division_type, adata.uns['data_name'], method,
                                             f'{genes.shape[0]}-impos.npy'),
                                arr=importances
                            )
    print(f"Genes have been saved to {os.path.join('geneData', division_type, adata.uns['data_name'], method)}.")


def is_saved(adata: ad.AnnData, method: str):
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    path = os.path.join('geneData', division_type, adata.uns['data_name'], method)
    is_fully_saved = os.path.exists(path) and len(os.listdir(path)) != 0 and np.char.endswith(os.listdir(path),
                                                                                              '.npy').sum() / 2 == len(
        exp_cfg.n_genes)
    if is_fully_saved:
        max_selected = np.load(os.path.join(path, f"{exp_cfg.n_genes[-1]}-genes.npy"), allow_pickle=True)
        if np.isin(max_selected, adata.var_names).sum() != exp_cfg.n_genes[-1]:
            is_fully_saved = False
    return is_fully_saved


def load_genes(adata: ad.AnnData, method: str):
    division_type = str(adata.uns['fold']) if 'fold' in adata.uns_keys() else 'all'
    feature_importance_list = []
    for n_gene in exp_cfg.n_genes:
        genes = np.load(os.path.join('geneData', division_type, adata.uns['data_name'], method, f'{n_gene}-genes.npy'),
                        allow_pickle=True)
        impos = np.load(os.path.join('geneData', division_type, adata.uns['data_name'], method, f'{n_gene}-impos.npy'))
        feature_importance_list.append((genes, impos))
    return feature_importance_list


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
