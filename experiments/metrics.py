import anndata as ad
import bcubed
import numpy as np
import scib
from harmonypy import compute_lisi
from sklearn.metrics import classification_report, cohen_kappa_score, adjusted_rand_score, v_measure_score

from common_utils.utils import HiddenPrints
from config.experiments_config import BatchCorrectionConfig, CellClassificationConfig, CellClusteringConfig


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


def classification_metrics(test_adata: ad.AnnData, config: CellClassificationConfig):
    assert 'assign_label' in test_adata.obs, "'assign_label' not in adata.obs!"
    classification_result = dict()

    if 'f1_score' in config.metrics:
        report = classification_report(test_adata.obs['celltype'], test_adata.obs['assign_label'], output_dict=True, zero_division='warn')
        classification_result['f1_score'] = report['macro avg']['f1-score']
        if 'rare_type' in test_adata.uns:
            if test_adata.uns['rare_type'] not in report.keys():
                raise RuntimeWarning(f"rare cell type {test_adata.uns['rare_type']} not in {report.keys()}!")
            else:
                classification_result['f1_rare'] = report[test_adata.uns['rare_type']]['f1-score']
    if 'cohen_kappa' in config.metrics:
        classification_result['cohen_kappa'] = cohen_kappa_score(test_adata.obs['celltype'], test_adata.obs['assign_label'])
    return classification_result


def clustering_metrics(adata: ad.AnnData, clustering_method: str, config: CellClusteringConfig):
    assert adata.obs.columns.str.endswith('_1').sum() >= 1, "Run at least once for certain clustering algorithm!"
    clustering_result = dict()
    clustering_labels = adata.obs.loc[:, adata.obs.columns.str.startswith(clustering_method)]

    for run in range(1, config.methods[clustering_method] + 1):
        labels_true, labels_pred = adata.obs['celltype'].values, clustering_labels[f'{clustering_method}_{run}']
        if 'ARI' in config.metrics:
            clustering_result[f'ARI_{run}'] = adjusted_rand_score(labels_true, labels_pred)
        if 'V' in config.metrics:
            clustering_result[f'V_{run}'] = v_measure_score(labels_true, labels_pred)
        if 'bcubed' in config.metrics:
            if 'rare_type' in adata.uns:
                clustering_result[f'bcubed_{run}'] = BCubed_fbeta_score(adata, f'{clustering_method}_{run}', True)

    return clustering_result



