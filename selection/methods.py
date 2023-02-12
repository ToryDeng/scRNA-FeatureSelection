import datetime
import math
import traceback
from collections import defaultdict, Counter
from typing import Optional, List

import anndata as ad
import anndata2ri
import pandas as pd
import scanpy as sc
import triku as tk
from lightgbm import LGBMClassifier
from scGeneClust import scGeneClust
from rpy2.robjects import r, globalenv
from rpy2.robjects.packages import importr
from scGeneFit.functions import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import NearestCentroid
from sklearn.preprocessing import LabelEncoder
from xgboost import XGBClassifier

from config import base_cfg
from selection.fisher_score import fisher_score
from selection.giniclust3 import loessRegression, giniIndex
from selection.nearest_centroid import nearest_centroid_select
from selection.utils import is_saved, load_genes, save_genes, subset_adata, rpy2_wrapper


def single_select_by_batch(adata: ad.AnnData,
                           method: str,
                           n_selected_genes: int = base_cfg.n_genes[-1]) -> pd.DataFrame:
    """
    Calculate feature importances for each gene, and then rank and select genes.

    Parameters
    ----------
    adata: AnnData
      The AnnData object
    method: str
      The feature selection method
    n_selected_genes: int
      Number of genes to be selected
    Returns
    -------
    selected_genes_df
      A dataframe with 1 or 2 columns. The first column with name 'Gene' contains the selected genes. The second column
      (if exists) with name 'Importance' contains gene importances.
    """
    if method == 'rf':
        selected_genes_df = random_forest_compute_importance(adata)
    elif method == 'lgb':
        selected_genes_df = random_compute_importance(adata)
    elif method == 'xgb':
        selected_genes_df = xgboost_compute_importance(adata)
    elif method == 'seurat':
        selected_genes_df = seurat_compute_importance(adata)
    elif method == 'seurat_v3':
        selected_genes_df = seurat_v3_compute_importance(adata)
    elif method == 'cellranger':
        selected_genes_df = cellranger_compute_importance(adata)
    elif method == 'var':
        selected_genes_df = variance_compute_importance(adata)
    elif method == 'cv':
        selected_genes_df = cv_compute_importance(adata)
    elif method == 'mi':
        selected_genes_df = mutual_info_compute_importance(adata)
    elif method == 'scGeneFit':
        selected_genes_df = scGeneFit(adata, n_selected_genes)
    elif method == 'fisher_score':
        selected_genes_df = fisher_score_compute_importance(adata)
    elif method == 'nsc':
        selected_genes_df = nearest_shrunken_centroid_compute_importance(adata)
    elif method == 'feast':
        selected_genes_df = FEAST(adata)
    elif method == 'random':
        selected_genes_df = random_compute_importance(adata)
    elif method == 'm3drop':
        selected_genes_df = M3Drop_compute_importance(adata)
    elif method == 'scmap':
        selected_genes_df = scmap_compute_importance(adata)
    elif method == 'deviance':
        selected_genes_df = deviance_compute_importance(adata)
    elif method == 'scran':
        selected_genes_df = scran(adata)
    elif method == 'triku':
        selected_genes_df = triku_compute_importance(adata.copy())  # .copy() prevent modification on adata object
    elif method == 'sct':
        selected_genes_df = sctransform_compute_importance(adata)
    elif method == 'giniclust3':
        selected_genes_df = giniclust3_compute_importance(adata)
    elif method == 'sc3':
        selected_genes_df = sc3_compute_importance(adata)
    elif method == 'pagest1w':
        selected_genes = scGeneClust(adata.raw.to_adata(), 'fast', random_stat=base_cfg.random_seed, verbosity=2)
        selected_genes_df = pd.DataFrame({'Gene': selected_genes})
    elif method == 'pagest2w':
        selected_genes = scGeneClust(adata.raw.to_adata(), 'hc',
                                     n_cell_clusters=adata.obs.celltype.unique().shape[0],
                                     random_stat=base_cfg.random_seed, verbosity=2)
        selected_genes_df = pd.DataFrame({'Gene': selected_genes})
    else:
        raise NotImplementedError(f"No implementation of {method}!")

    # rank features if importance exist
    selected_genes_df.dropna(inplace=True)
    if selected_genes_df.shape[1] == 1:  # doesn't have importances, only genes
        selected_genes_df.rename(columns={selected_genes_df.columns[0]: 'Gene'}, inplace=True)
    elif selected_genes_df.shape[1] == 2:  # genes, importances
        col_names = selected_genes_df.columns
        selected_genes_df.rename(columns={col_names[0]: 'Gene', col_names[1]: 'Importance'}, inplace=True)
        selected_genes_df.sort_values(by='Importance', ascending=False, inplace=True)
    else:
        raise ValueError(f"There are more than 2 columns in results given by {method}: {selected_genes_df.columns}")

    # select features
    if selected_genes_df.shape[0] < n_selected_genes:
        raise RuntimeWarning("Genes are not enough to be selected!")
    else:
        selected_genes_df = selected_genes_df.iloc[:n_selected_genes, :]
    return selected_genes_df


def ensemble_select_by_batch(adata: ad.AnnData,
                             base_methods: List[str],
                             n_selected_genes: int = base_cfg.n_genes[-1],
                             ):
    # calculate ensemble feature importances
    gene_dict = defaultdict(float)
    if base_cfg.ensemble_mode == 'importance_sum':
        for base_method in base_methods:
            selected_df = single_select_by_batch(adata, base_method, n_selected_genes=n_selected_genes)
            if selected_df.shape[1] == 1:
                raise RuntimeError(f"{base_method} can't generate feature importances!")
            min_impo, max_impo = selected_df['Importance'].min(), selected_df['Importance'].max()  # min-max norm
            selected_df['Importance'] = (selected_df['Importance'] - min_impo) / (max_impo - min_impo)
            for gene, importance in zip(selected_df['Gene'], selected_df['Importance']):
                gene_dict[gene] += importance
    elif base_cfg.ensemble_mode == 'count_sum':
        for base_method in base_methods:
            selected_df = single_select_by_batch(adata, base_method, n_selected_genes=n_selected_genes)
            for gene in selected_df['Gene']:
                gene_dict[gene] += 1
    else:
        raise NotImplementedError(f"{base_cfg.ensemble_mode} has not been implemented yet.")
    # rank genes
    sorted_result = pd.Series(gene_dict).dropna().sort_values(ascending=False).reset_index()
    sorted_result.rename(columns={sorted_result.columns[0]: 'Gene', sorted_result.columns[1]: 'Importance'},
                         inplace=True)
    # select features
    if sorted_result.shape[0] < n_selected_genes:
        raise RuntimeWarning("Genes are not enough to be selected!")
    else:
        sorted_result = sorted_result.iloc[:n_selected_genes, :]
    return sorted_result


def select_genes_by_batch(adata: ad.AnnData,
                          method: str,
                          n_selected_genes: int = base_cfg.n_genes[-1],
                          use_saved: bool = True
                          ):
    if use_saved and is_saved(adata, method, n_selected_genes):
        final_result = load_genes(adata, method, n_selected_genes)
        print(f'{n_selected_genes} genes are selected by {method} using previously saved genes and importances...')
    else:  # do not use saved genes or genes have not been saved
        try:
            if '+' in method:
                final_result = ensemble_select_by_batch(adata, method.split('+'), n_selected_genes)
            else:
                final_result = single_select_by_batch(adata, method, n_selected_genes)
        except:
            traceback.print_exc()
            final_result = None
        if final_result is not None:
            save_genes(adata, method, final_result)
    return final_result


def select_genes(
        adata: ad.AnnData,
        method: str,
        n_selected_genes: int,
        inplace: bool = False,
        use_saved: bool = True,
        select_by_batch=True
) -> Optional[ad.AnnData]:
    print(f'{method} starts to select {n_selected_genes} genes...')
    if 'batch' in adata.obs and select_by_batch:
        # calculate feature importances
        ubatches = adata.obs['batch'].unique().categories.to_numpy()
        print(f"{adata.uns['data_name']} contains {ubatches.shape[0]} batches: {ubatches}")
        batch_features = []
        for ubatch in ubatches:
            batch_adata = adata[adata.obs['batch'] == ubatch].copy()
            # rename batch data set
            if '+' not in batch_adata.uns['data_name']:
                batch_adata.uns['data_name'] = '_'.join([batch_adata.uns['data_name'], ubatch])
            else:
                batch_adata.uns['data_name'] = ubatch

            print(
                f"{datetime.datetime.now()}: Start to apply {method} on batch {ubatch}, batch shape: {batch_adata.shape}.")
            batch_result = select_genes_by_batch(batch_adata, method, n_selected_genes, use_saved)
            if batch_result is None:
                raise RuntimeWarning(f"batch {ubatch} is not used.")
            else:
                if batch_result.shape[1] == 1:
                    raise RuntimeError(
                        f"{method} can't generate gene importances, so it can't select genes on batches!")
                batch_features.append(batch_result.set_index('Gene')['Importance'])
        batch_features_rank = [batch_feature.rank() for batch_feature in batch_features]
        cnt = Counter()
        for batch_feature in batch_features:
            cnt.update(batch_feature.index)
        ties = pd.Series(cnt).sort_values(ascending=False)
        untied_rank = []
        for tie in np.unique(ties.values):
            rank_dict = defaultdict(list)
            for gene in ties[ties == tie].index:
                for batch_rank in batch_features_rank:
                    if gene in batch_rank.index:
                        rank_dict[gene].append(batch_rank[gene])
            median_rank_dict = {key: np.median(value) for key, value in rank_dict.items()}
            untied = pd.Series(median_rank_dict).sort_values(ascending=True)
            untied_rank.append(untied)
        # rank features
        untied_features = pd.concat(untied_rank[::-1]).to_frame().reset_index()  # have ascending order
        # select features
        if untied_features.shape[0] < n_selected_genes:
            raise RuntimeWarning(f"Only {untied_features.shape[0]} genes! Genes are not enough to be selected!")
        else:
            untied_features = untied_features.iloc[:n_selected_genes, :]
            untied_features.rename(columns={untied_features.columns[0]: 'Gene', untied_features.columns[1]: 'Rank'},
                                   inplace=True)
        df = untied_features  # most important genes
    else:
        df = select_genes_by_batch(adata, method, n_selected_genes, use_saved)

    if df is not None:  # have enough genes
        if inplace:
            subset_adata(adata, df['Gene'].values, inplace=True)
        else:
            return subset_adata(adata, df['Gene'].values, inplace=False)


# python methods
def random_forest_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    y = adata.obs['celltype'].values
    forest = RandomForestClassifier(n_jobs=-1, random_state=0, verbose=0).fit(adata.raw.X, y)
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': forest.feature_importances_})


def lightgbm_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    y = adata.obs['celltype'].values
    lgb = LGBMClassifier(n_jobs=16, random_state=0).fit(adata.raw.X, y)
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': lgb.feature_importances_})


def xgboost_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    le = LabelEncoder()
    y = le.fit_transform(adata.obs['celltype'].values)
    xgb = XGBClassifier(eval_metric='mlogloss', use_label_encoder=False, nthread=-1).fit(adata.raw.X, y)
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': xgb.feature_importances_})


def seurat_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=adata.n_vars // 2)
    return pd.DataFrame({'Gene': adata.var.dispersions_norm.index, 'Importance': adata.var.dispersions_norm})


def seurat_v3_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    raw_adata = adata.raw.to_adata()
    sc.pp.highly_variable_genes(raw_adata, flavor='seurat_v3', n_top_genes=adata.n_vars // 2)
    return pd.DataFrame({'Gene': raw_adata.var['variances_norm'].index, 'Importance': raw_adata.var['variances_norm']})


def cellranger_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=adata.n_vars // 2)
    return pd.DataFrame({'Gene': adata.var['dispersions_norm'].index, 'Importance': adata.var.dispersions_norm})


def variance_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': np.var(adata.layers['log-normalized'], axis=0)})


def cv_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    std, mean = np.std(adata.layers['log-normalized'], axis=0), np.mean(adata.layers['log-normalized'], axis=0)
    print('Number of genes whose mean are less than 1e-4: {}'.format(np.sum(np.squeeze(mean) < 1e-4)))
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': std / mean})


def mutual_info_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    y = adata.obs['celltype'].values
    mutual_info = mutual_info_classif(adata.raw.X, y, discrete_features=False)
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': mutual_info})


def scGeneFit(adata: ad.AnnData, n_selected_genes: int) -> pd.DataFrame:
    markers = get_markers(adata.X, adata.obs['celltype'].values, n_selected_genes)
    return pd.DataFrame({'Gene': adata.var_names[markers]})


def fisher_score_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    y = adata.obs['celltype'].values
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': fisher_score(adata.raw.X, y)})


def nearest_shrunken_centroid_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    y = adata.obs['celltype'].values
    th_dict = {'shrink_threshold': np.arange(0.0, 1.01, 0.01)}
    gs = GridSearchCV(NearestCentroid(), param_grid=th_dict, cv=3, scoring='balanced_accuracy').fit(adata.X, y)
    print('best score:{}, best threshold:{}'.format(gs.best_score_, gs.best_params_['shrink_threshold']))
    importance = nearest_centroid_select(adata.X, y, shrink_threshold=gs.best_params_['shrink_threshold'])
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': importance})


def triku_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata)
    if 'distances' not in adata.obsp or 'connectivities' not in adata.obsp:
        sc.pp.neighbors(adata)
    tk.tl.triku(adata, verbose='error')  # the larger the distance, the more important the gene is
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': adata.var['triku_distance']})


def random_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    np.random.seed(base_cfg.random_seed)
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': np.random.rand(adata.X.shape[1])})


def giniclust3_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    geneGini, geneMax = giniIndex(adata.layers['normalized'].T)
    geneGini, geneLabel = np.array(geneGini), adata.var.index.tolist()
    logGeneMax = np.array(list(map(lambda x: math.log(x + 0.1, 2), geneMax)))
    _, sigGiniPvalue = loessRegression(geneGini, logGeneMax, geneLabel)
    result = pd.DataFrame.from_dict(sigGiniPvalue, orient='index', columns=['p.value']).reset_index()
    result['Importance'] = 1 - result['p.value']
    result.drop(columns=['p.value'], inplace=True)
    return result


def sc3_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    importance = np.zeros(shape=(adata.raw.shape[1],))  # the importances of filtered genes are set to zero
    dropouts = np.sum(adata.raw.X == 0, axis=0) / adata.raw.shape[1] * 100
    sc3_gene_filter = np.logical_and(10 < dropouts, dropouts < 90)
    print(f"Filtering {adata.raw.shape[1] - sc3_gene_filter.sum()} rarely and ubiquitously expressed genes...")
    importance[sc3_gene_filter] = adata.raw.X[:, sc3_gene_filter].mean(axis=0)  # expression levels
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': importance})


# R methods
@rpy2_wrapper
def FEAST(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    """
    Select features by FEAST. The input AnnData object contains norm and raw data. The raw data is normalized in R.

    Parameters
    ----------
    adata :
        anndata object.
    Returns
    -------
    A 1-column dataframe. The sole column contains genes sorted according to their F scores.
    The top genes are the most important.
    """
    importr('FEAST')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""
    Y <- process_Y(assay(sce, 'X'), thre = 2)
    Ynorm <- Norm_Y(Y)
    n_classes <- dim(unique(colData(sce)['celltype']))[1]
    rm(sce)
    idxs = FEAST_fast(Ynorm, k = n_classes)
    """)
    result = r("data.frame(Gene=rownames(Y)[idxs])")
    return result


@rpy2_wrapper
def M3Drop_compute_importance(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    """
    Select features by M3Drop. The input AnnData object contains norm and raw data. The raw data is normalized in R.

    Parameters
    ----------
    adata
      anndata object.
    Returns
    -------
    result
      A dataframe. The first column contains gene names, and the second column contains 1 - p.value.
    """
    importr('M3Drop')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""
    norm <- M3DropConvertData(assay(sce, 'X'), is.counts=TRUE)
    DE_genes <- M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=1, suppress.plot = TRUE)
    """)
    result = r("DE_genes").drop(columns=['effect.size', 'q.value'])
    result['Importance'] = 1 - result['p.value']
    result.drop(columns=['p.value'], inplace=True)
    return result


@rpy2_wrapper
def scmap_compute_importance(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    importr('scmap')
    importr('scater')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""
    counts(sce) <- as.matrix(assay(sce, 'X'))
    sce <- logNormCounts(sce)
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- selectFeatures(sce, dim(sce)[1], suppress_plot = TRUE)
    res <- na.omit(rowData(sce)['scmap_scores'])
    """)
    result = r("res").reset_index()  # index  scmap_scores
    return result


@rpy2_wrapper
def deviance_compute_importance(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    importr('scry')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())  # deviance need raw data
    result = r("rowData(devianceFeatureSelection(sce, assay='X'))['binomial_deviance']").reset_index()
    return result


@rpy2_wrapper
def scran(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    importr('scuttle')
    importr('scran')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""
    assay(sce, 'X') <- as.matrix(assay(sce, 'X'))
    sce <- logNormCounts(sce, assay.type = "X")
    HVGs <- getTopHVGs(sce)
    """)
    result = r("data.frame(Gene=HVGs)")
    return result


@rpy2_wrapper
def sctransform_compute_importance(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    importr('sctransform')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""res <- vst(as(assay(sce, 'X'), "sparseMatrix"), method = "glmGamPoi", verbosity = 0)""")
    result = r("res$gene_attr['residual_variance']").reset_index()
    return result
