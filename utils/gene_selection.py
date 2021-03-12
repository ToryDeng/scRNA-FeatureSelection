from utils.utils import get_gene_names, save_filtered_data, load_data
from utils.evaluate_result import evaluate_method
from utils.importance import select_features
import numpy as np
import os


def evaluate_gene_selection_method(dataset=None, methods=None, datatype=None):
    """
    Evaluate certain gene selection methods on specific dataset.

    :param dataset: string, the name of dataset
    :param methods: string list, names of gene selection methods
    :param datatype: string, the method is used on raw or norm data
    :return: None
    """
    X_raw, X_norm, y, trusted_markers = load_data(dataset)
    gene_names = get_gene_names(X_raw.columns)
    for method in methods:
        if datatype == 'raw':
            result = select_features(dataset, 1000, method, gene_names, X=X_raw.values, y=np.squeeze(y.values))
        elif datatype == 'norm':
            result = select_features(dataset, 1000, method, gene_names, X=X_raw.values, y=np.squeeze(y.values))
        else:
            print("The parameter 'datatype' is wrong. Please check again.")
            return None
        save_filtered_data(X_raw, y, gene_names, result)
        # TODO same numbers of genes and clusters
        os.system('Rscript scRNA-FeatureSelection/utils/RCode/main.R')
        evaluate_method(trusted_markers, result, y)
