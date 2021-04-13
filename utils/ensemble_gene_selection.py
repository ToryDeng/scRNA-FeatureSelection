from utils.utils import get_gene_names, load_data, cal_marker_num_MRR, delete, filter_const_cells, PerformanceRecord,\
    save_raw_data, filter_const_genes, now, head
from utils.evaluation import evaluate_classification_result, evaluate_clustering_result
from utils.importance import select_features
from config import ClassificationConfig, ClusteringConfig
from collections import defaultdict
from sklearn.model_selection import StratifiedKFold
import pandas as pd
import numpy as np


def integrateResults(dataset: str, base_methods: list, gene_names, X_raw, X_norm, y, config):
    gene_importance_dict = defaultdict(int)
    for method in base_methods:
        # delete current files in tempData
        delete('scRNA-FeatureSelection/tempData/')
        if config.method_on[method] == 'raw':
            save_raw_data(X_train=X_raw, X_test=None, y_train=y, y_test=None, task='clustering')  # for consistency
            genes, importances = select_features(dataset, 1000, method, gene_names, X_raw.values,
                                                 np.squeeze(y.values))
        else:  # norm data
            save_raw_data(X_train=X_norm, X_test=None, y_train=y, y_test=None, task='clustering')  # for consistency
            genes, importances = select_features(dataset, 1000, method, gene_names, X_norm.values,
                                                 np.squeeze(y.values))
        importances = (importances - np.min(importances)) / (np.max(importances) - np.min(importances))  # Normalization

        for gene, importance in zip(genes, importances):
            gene_importance_dict[gene] += importance
    sorted_result = pd.Series(gene_importance_dict).sort_values(ascending=False).iloc[:1000].reset_index().T.to_numpy()
    return sorted_result


def ensemble_gene_selection(dataset: str, base_methods: list, task: str):
    # loading data
    X_raw, X_norm, y, trusted_markers = load_data(dataset)
    gene_names = get_gene_names(X_raw.columns)
    ensemble_method_name = '-'.join(base_methods)
    # TODO: No importance?
    performance_record = PerformanceRecord(methods=[ensemble_method_name], task=task)

    # output dataset information
    head(name='Dataset Information')
    print("Name:{}  Cell(s):{}  Gene(s):{}\nMarker Gene(s):{}".format(
        dataset, X_raw.shape[0], X_raw.shape[1], np.intersect1d(gene_names, trusted_markers).shape[0])
    )

    if task == 'assign':
        cfg = ClassificationConfig()
        skf = StratifiedKFold(n_splits=5, random_state=2020, shuffle=True)
        for train_idx, test_idx in skf.split(X_raw, y):
            X_raw_train, X_raw_test = X_raw.iloc[train_idx], X_raw.iloc[test_idx]
            X_norm_train, X_norm_test = X_norm.iloc[train_idx], X_norm.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
            delete('scRNA-FeatureSelection/tempData/')

            # calculate F1-score before feature selection
            save_raw_data(X_raw_train, X_raw_test, y_train, y_test, task='assign')
            before = evaluate_classification_result()
            print(before)

            # remove const genes before feature selection
            X_raw_train, X_norm_train = filter_const_genes(X_raw_train), filter_const_genes(X_norm_train)
            X_raw_test, X_norm_test = X_raw_test.loc[:, X_raw_train.columns], X_norm_test.loc[:, X_norm_train.columns]
            gene_names = get_gene_names(X_raw_train.columns)

            # integrate feature selection results from training set
            sorted_result = integrateResults(dataset, base_methods, gene_names, X_raw_train, X_norm_train, y_train, cfg)
            mask = np.isin(gene_names, sorted_result[0])

            # calculate and store markers_found & MRR
            markers_found, MRR = cal_marker_num_MRR(trusted_markers, sorted_result[0], rank=True)
            performance_record.loc['marker_genes_found', ensemble_method_name] += markers_found
            performance_record.loc['MRR', ensemble_method_name] += MRR

            # filter out non-markers
            X_train_selected, y_train_selected = filter_const_cells(X_raw_train.loc[:, mask], y_train)
            X_test_selected = X_raw_test.loc[:, mask]

            # save X_train and X_test after gene selection
            save_raw_data(X_train_selected, X_test_selected, y_train_selected, y_test, task='assign')

            # execute R script to generate classification result
            after = evaluate_classification_result()

            # update performance record
            print(after)
            for key in after.keys():
                performance_record.loc[key, ensemble_method_name] += after[key] - before[key]
        performance_record.divide(5).to_csv(
            ''.join(['scRNA-FeatureSelection/results/', dataset, '_', ensemble_method_name, '_assign.csv']))
    elif task == 'clustering':
        # save raw data and generate clustering result before feature selection
        save_raw_data(X_train=X_raw, X_test=None, y_train=y, y_test=None, task='clustering')
        before = evaluate_clustering_result()
        cfg = ClusteringConfig()

        # integrate feature selection results
        sorted_result = integrateResults(dataset, base_methods, gene_names, X_raw, X_norm, y, cfg)
        mask = np.isin(gene_names, sorted_result[0])

        # calculate and store markers_found & MRR
        markers_found, MRR = cal_marker_num_MRR(trusted_markers, sorted_result[0], rank=True)
        performance_record.loc['marker_genes_found', ensemble_method_name] = markers_found
        performance_record.loc['MRR', ensemble_method_name] = MRR

        # filter out non-markers
        X_selected, y_selected = filter_const_cells(X_raw.loc[:, mask], y)

        # save raw data and generate clustering result before feature selection
        save_raw_data(X_train=X_selected, X_test=None, y_train=y_selected, y_test=None, task='clustering')
        after = evaluate_clustering_result()

        # update performance record
        for key in after.keys():
            performance_record.loc[key, ensemble_method_name] = after[key] - before[key]

        # save performance record
        performance_record.to_csv(
            ''.join(['scRNA-FeatureSelection/results/', dataset, '_', ensemble_method_name,  '_clustering.csv']))
    print(now() + ": Evaluation is done!")
