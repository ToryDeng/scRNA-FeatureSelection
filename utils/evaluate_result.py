import numpy as np
from sklearn.metrics import adjusted_rand_score, classification_report
from utils.utils import delAll
import warnings


def evaluate_method(trusted_features, selected_features, rank=True):
    """
    Calculate MRR (if rank=True), number of markers found, ARI and many other classification metrics.

    :param trusted_features: marker genes read from marker files
    :param selected_features: selected genes
    :param rank: whether to calculate ARI metric
    :return: None
    """
    if rank:
        rank_list = np.argwhere(np.isin(selected_features, trusted_features))
        if len(rank_list) == 0:
            warnings.warn("MRR: Can not calculate MRR because no marker gene is selected!", RuntimeWarning)
        else:
            print('MRR:{:.5f}'.format(np.sum(1 / (rank_list + 1)) / rank_list.shape[0]))
    print('marker genes found:{}'.format(np.intersect1d(trusted_features, selected_features).shape[0]))

    temp_path = 'scRNA-FeatureSelection/tempData/'
    # evaluate using ARI
    label_true = np.loadtxt(temp_path + 'temp_y.csv', delimiter=',', skiprows=1, usecols=[1], dtype=np.str)
    for clustering_method in ['seurat', 'sc3']:
        label_pred = np.loadtxt(
            ''.join([temp_path, 'temp_', clustering_method, '.csv']),
            delimiter=',', skiprows=1, dtype=np.str)
        ari = adjusted_rand_score(np.squeeze(label_true), np.squeeze(label_pred))
        print("ARI of {}: {:.5f}".format(clustering_method, ari))
    # evaluate using classification methods
    label_test = np.loadtxt(temp_path + 'temp_y_test.csv', delimiter=',', skiprows=1, dtype=np.object)[:, 1]
    for assign_method in ['scmap_cluster', 'scmap_cell', 'singlecellnet']:
        try:
            label_pred = np.loadtxt(
                ''.join([temp_path, 'temp_', assign_method, '.csv']),
                delimiter=',', skiprows=1, dtype=np.str)
            print("classification report of {}:\n".format(assign_method))
            print(classification_report(np.squeeze(label_test), np.char.strip(label_pred, '"')))
        except IOError:
            print('{} failed to execute successfully. Please check your R output.'.format(assign_method))
    delAll(temp_path)
    print('\n')
