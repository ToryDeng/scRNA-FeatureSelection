import numpy as np
from sklearn.metrics import adjusted_rand_score, classification_report
from utils.utils import delAll
import os


def evaluate_method(trusted_features, selected_features, label_true, rank=True):
    """
    Calculate MRR (if rank=True), number of markers found, ARI and many other classification metrics.

    :param trusted_features: marker genes read from marker files
    :param selected_features: selected genes
    :param label_true: true cell types
    :param rank: whether to calculate ARI metric
    :return: None
    """
    if rank:
        rank = np.argwhere(np.isin(selected_features, trusted_features))
        print('MRR:{:.5f}'.format(np.sum(1 / (rank + 1)) / rank.shape[0]))
    print('marker genes found:{}'.format(np.intersect1d(trusted_features, selected_features).shape[0]))

    os.chdir(r'../TempData')
    # evaluate using ARI
    for clustering_method in ['seurat', 'sc3']:
        label_pred = np.loadtxt('temp_' + clustering_method + '.csv', delimiter=',', skiprows=1, dtype=np.object)
        ari = adjusted_rand_score(np.squeeze(label_true), np.squeeze(label_pred))
        print("ari of {}: {:.5f}".format(clustering_method, ari))
    # evaluate using classification methods
    label_test = np.loadtxt('temp_y_test.csv', delimiter=',', skiprows=1,dtype=np.object)[:, 1]
    for assign_method in ['scmap_cluster', 'scmap_cell', 'singlecellnet']:
        label_pred = np.loadtxt('temp_' + assign_method + '.csv', delimiter=',', skiprows=1, dtype=np.object)
        print("classification report of {}:\n".format(assign_method))
        print(classification_report(np.squeeze(label_test), np.squeeze(label_pred)))
#     delAll(r'/home/tdeng/SingleCell/FeatureSelection/R/TempData')
    print('\n')


