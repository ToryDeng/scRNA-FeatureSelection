import anndata as ad
import bcubed
import numpy as np
import scanpy as sc
import scib
from harmonypy import compute_lisi
from sklearn.metrics import classification_report, cohen_kappa_score, adjusted_rand_score, v_measure_score, \
    silhouette_score, normalized_mutual_info_score

from common_utils.utils import HiddenPrints
from config import BatchCorrectionConfig, assign_cfg, cluster_cfg


def marker_discovery_rate(selected_adata: ad.AnnData, original_adata: ad.AnnData):
    assert 'is_marker' and 'marker_weight' in selected_adata.var, "'selected_adata' doesn't contain marker genes!"
    return selected_adata.var['marker_weight'].sum() / original_adata.var['marker_weight'].sum()


def BCubed_fbeta_score(adata: ad.AnnData, use_cluster_rep: str, only_rare: bool = False, beta: float = 1.0):
    if not only_rare:
        ldict = {cell: {cell_type} for cell, cell_type in adata.obs['celltype'].to_dict().items()}
        cdict = {cell: {cluster} for cell, cluster in adata.obs[use_cluster_rep].to_dict().items()}
    else:
        if adata.uns['rare_type'] is not None:
            if adata.uns['rare_type'] not in adata.obs['celltype'].unique():
                raise RuntimeError(f"{adata.uns['rare_type']} not in adata.obs['celltype']")
            rare_cell_clusters = adata.obs.loc[adata.obs['celltype'] == adata.uns['rare_type'], use_cluster_rep].to_dict()
            ldict = {cell: {adata.uns['rare_type']} for cell, _ in rare_cell_clusters.items()}
            cdict = {cell: {cluster} for cell, cluster in rare_cell_clusters.items()}
        else:
            raise RuntimeError("Only compute BCubed F1 score on rare cells but the rare cell type is None.")
    precision, recall = bcubed.precision(cdict, ldict), bcubed.recall(cdict, ldict)
    return bcubed.fscore(precision, recall, beta=beta)


def correction_metrics(adata: ad.AnnData, config: BatchCorrectionConfig) -> dict:
    sc.pp.pca(adata)  # do PCA anyway
    correction_result = dict()

    with HiddenPrints():
        if 'kBET' in config.metrics:
            correction_result[f'kBET'] = scib.metrics.kBET(adata, batch_key='batch', label_key='celltype', embed=f'X_pca')
        if 'cLISI' not in config.metrics and 'iLISI' not in config.metrics and 'f1LISI' not in config.metrics:  # only kBET
            cells_lisi = compute_lisi(adata.obsm[f'X_pca'], adata.obs, ['batch', 'celltype'])
            norm_lisi = (cells_lisi - cells_lisi.min(axis=0)) / (cells_lisi.max(axis=0) - cells_lisi.min(axis=0))
            median_lisi = np.median(norm_lisi, axis=0)

            correction_result['iLISI'] = median_lisi[0]
            correction_result['cLISI'] = median_lisi[1]
            correction_result['f1LISI'] = 2 * (1 - median_lisi[1]) * median_lisi[0] / (1 - median_lisi[1] + median_lisi[0])

    return correction_result


def classification_metrics(test_adata: ad.AnnData):
    assert test_adata.obs.columns.str.endswith('_label').sum() >= 1, "Run at least once for certain classification algorithm!"
    classification_result = dict()
    for method in assign_cfg.methods:
        per_method = dict()
        labels_true, labels_pred = test_adata.obs['celltype'].values, test_adata.obs[f'{method}_label']
        if 'f1' in assign_cfg.metrics:
            report = classification_report(labels_true, labels_pred, output_dict=True, zero_division=0)
            per_method['f1'] = report['macro avg']['f1-score']
        if 'ck' in assign_cfg.metrics:
            per_method['ck'] = cohen_kappa_score(labels_true, labels_pred)
        classification_result[method] = per_method
    return classification_result


def clustering_metrics(adata: ad.AnnData):
    assert adata.obs.columns.str.endswith('_0').sum() >= 1, "Run at least once for certain clustering algorithm!"
    clustering_result = dict()

    for method, n_runs in cluster_cfg.methods.items():
        per_method = dict()
        for run in range(n_runs):
            labels_true, labels_pred = adata.obs['celltype'].values, adata.obs[f'{method}_{run}'].values
            if 'ARI' in cluster_cfg.metrics:
                per_method[f'ARI_{run}'] = adjusted_rand_score(labels_true, labels_pred)
            if 'V' in cluster_cfg.metrics:
                per_method[f'V_{run}'] = v_measure_score(labels_true, labels_pred)
            if 'NMI' in cluster_cfg.metrics:
                per_method[f'NMI_{run}'] = normalized_mutual_info_score(labels_true, labels_pred)
            if 'SI' in cluster_cfg.metrics:
                if 'X_pca' in adata.obsm:
                    del adata.obsm['X_pca']  # do pca anyway
                sc.pp.pca(adata)
                per_method[f'SI_{run}'] = silhouette_score(adata.obsm['X_pca'], labels_pred)
            if 'bcubed' in cluster_cfg.metrics:
                if 'rare_type' in adata.uns:
                    per_method[f'bcubed_{run}'] = BCubed_fbeta_score(adata, f'{method}_{run}', True)
        clustering_result[method] = per_method
    return clustering_result


