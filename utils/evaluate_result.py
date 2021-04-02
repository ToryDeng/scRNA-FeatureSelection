import numpy as np
from sklearn.metrics import adjusted_rand_score, f1_score
from utils.utils import delAll
import warnings


def evaluate_method(trusted_features, selected_result, mode='all', rank=True):
    """
    Calculate MRR (if rank=True), number of markers found, ARI and many other classification metrics.

    :param trusted_features: marker genes read from marker files
    :param selected_result: selected genes and importance
    :param rank: whether to calculate ARI metric
    :param mode: 'eval' for evaluation
    :return: a dict. If mode = 'all', return ARI and F1-score; If mode = 'eval', return MRR, marker genes found,
    ARI and F1-score
    """
    result = dict()
    if mode == 'eval':
        if rank:
            rank_list = np.argwhere(np.isin(selected_result[0], trusted_features))
            if len(rank_list) == 0:
                warnings.warn("MRR: Can not calculate MRR because no marker gene is selected!", RuntimeWarning)
            else:
                MRR = np.sum(1 / (rank_list + 1)) / rank_list.shape[0]
                result['MRR'] = MRR
        else:  # len(selected_result) == 1
            warnings.warn("The method can't obtain gene importance! MRR can't be calculated and is set to 0.", RuntimeWarning)
        marker_genes_found = np.intersect1d(trusted_features, selected_result[0]).shape[0]
        result['marker_genes_found'] = marker_genes_found

    temp_path = 'scRNA-FeatureSelection/tempData/'
    # evaluate using ARI
    label_true = np.loadtxt(temp_path + 'temp_y.csv', delimiter=',', skiprows=1, usecols=[1], dtype=np.str)
    for clustering_method in ['seurat', 'sc3']:
        label_pred = np.loadtxt(
            ''.join([temp_path, 'temp_', clustering_method, '.csv']),
            delimiter=',', skiprows=1, dtype=np.str)
        ari = adjusted_rand_score(np.squeeze(label_true), np.squeeze(label_pred))
        result[clustering_method + '_ARI'] = ari
    # evaluate using classification methods
    label_test = np.loadtxt(temp_path + 'temp_y_test.csv', delimiter=',', skiprows=1, dtype=np.object)[:, 1]
    for assign_method in ['scmap_cluster', 'scmap_cell', 'singlecellnet']:
        try:
            label_pred = np.loadtxt(''.join([temp_path, 'temp_', assign_method, '.csv']), delimiter=',', skiprows=1, dtype=np.str)
            f1 = f1_score(np.squeeze(label_test), np.char.strip(label_pred, '"'), average='weighted')
            result[assign_method + '_F1'] = f1
        except IOError:
            print('{} failed to execute successfully. Please check the R output in console.'.format(assign_method))
    # delAll(temp_path)
    return result
