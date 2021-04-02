from utils.utils import get_gene_names, save_data, load_data
from utils.evaluate_result import evaluate_method
from utils.importance import select_features
from collections import defaultdict
import warnings
import pandas as pd
import numpy as np
import os


def evaluate_gene_selection_method(dataset=None, methods=None, data_type=None):
    """
    Evaluate certain gene selection methods on specific dataset.

    :param dataset: string, the name of dataset
    :param methods: string list, names of gene selection methods
    :param data_type: string, the method is used on raw or norm data
    :return: None
    """
    if dataset[:4] == 'PBMC' and 'scGeneFit' in methods:
        dataset = 'PBMC5'
        warnings.warn("Using 5% of PBMC cells because scGeneFit needs lots of system resource.", RuntimeWarning)
    X_raw, X_norm, y, trusted_markers = load_data(dataset)

    ix = ['MRR', 'marker_genes_found', 'seurat_ARI', 'sc3_ARI', 'scmap_cluster_F1', 'scmap_cell_F1', 'singlecellnet_F1']
    performance_record = pd.DataFrame(np.zeros((len(ix), len(methods))), index=ix, columns=methods)
    print("*************** Dataset Information ***************")
    print("Name:{}  Type:{}  Cell(s):{}  Gene(s):{}".format(dataset, data_type, X_raw.shape[0], X_raw.shape[1]))
    gene_names = get_gene_names(X_raw.columns)
    for method in methods:
        if data_type == 'raw':
            result = select_features(dataset, 1000, method, gene_names, X=X_raw.values, y=np.squeeze(y.values))
        elif data_type == 'norm':
            result = select_features(dataset, 1000, method, gene_names, X=X_norm.values, y=np.squeeze(y.values))
        else:
            print("The parameter 'data_type' is wrong. Please check again.")
            return None
        # calculate ARI and F1-score before gene selection
        save_data(X_raw, y, filtering=False)
        os.system('Rscript scRNA-FeatureSelection/utils/RCode/main.R')
        if len(result) == 2:
            before = evaluate_method(trusted_markers, result[0], mode='all')
        else:  # len(result) == 1, i. e. no feature importance
            before = evaluate_method(trusted_markers, result, mode='all', rank=False)
        # calculate ARI and F1-score after gene selection
        save_data(X_raw, y, gene_names, result)
        os.system('Rscript scRNA-FeatureSelection/utils/RCode/main.R')
        if len(result) == 2:
            after = evaluate_method(trusted_markers, result[0], mode='eval')
        else:
            after = evaluate_method(trusted_markers, result, mode='eval', rank=False)
        for key in after.keys():
            if key not in before.keys():
                performance_record.loc[key, method] = after[key]
            else:
                performance_record.loc[key, method] = after[key] - before[key]
    performance_record.to_csv(''.join(['scRNA-FeatureSelection/results/', dataset, '_', data_type, '_record.csv']))


def ensemble_gene_selection(dataset=None, base_methods=None):
    # loading data
    X_raw, X_norm, y, trusted_markers = load_data(dataset)
    gene_names = get_gene_names(X_raw.columns)

    gene_importance_dict = defaultdict(int)
    for method in base_methods:  # TODO: 1.data_type, also a dict? 2.No importance?
        genes, importances = select_features(dataset, 1000, method, gene_names, X=X_raw.values, y=np.squeeze(y.values))
        importances = (importances - np.mean(importances)) / np.std(importances)  # Normalization

        for gene, importance in zip(genes, importances):
            gene_importance_dict[gene] += importance
    sorted_result = pd.Series(gene_importance_dict).sort_values(ascending=False).iloc[:1000].reset_index().T.to_numpy()
    print(sorted_result)
    ix = ['MRR', 'marker_genes_found', 'seurat_ARI', 'sc3_ARI', 'scmap_cluster_F1', 'scmap_cell_F1', 'singlecellnet_F1']
    ensemble_method_name = '-'.join(base_methods)
    performance_record = pd.DataFrame(np.zeros((len(ix), 1)), index=ix, columns=[ensemble_method_name])

    # calculate ARI and F1-score before gene selection
    save_data(X_raw, y, filtering=False)
    os.system('Rscript scRNA-FeatureSelection/utils/RCode/main.R')
    before = evaluate_method(trusted_markers, sorted_result, mode='all')
    # calculate ARI and F1-score after gene selection
    save_data(X_raw, y, gene_names, sorted_result)
    os.system('Rscript scRNA-FeatureSelection/utils/RCode/main.R')
    after = evaluate_method(trusted_markers, sorted_result, mode='eval')
    for key in after.keys():
        if key not in before.keys():
            performance_record.loc[key, ensemble_method_name] = after[key]
        else:
            performance_record.loc[key, ensemble_method_name] = after[key] - before[key]
    performance_record.to_csv(''.join(['scRNA-FeatureSelection/results/', dataset, '_', ensemble_method_name, '_record.csv']))
