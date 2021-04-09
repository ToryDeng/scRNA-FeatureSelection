from sklearn.metrics import adjusted_rand_score, f1_score
from utils.utils import get_gene_names, load_data, cal_marker_num_MRR, delete, filter_const_cells, PerformanceRecord, \
    save_raw_data, filter_const_genes, now, head
from utils.importance import select_features
from sklearn.model_selection import KFold
import numpy as np
import warnings
import os


def evaluate_classification_result():
    """
    Evaluate classification result using the F1-scores generating from three assign methods.

    :return: a dict containing F1-scores of three assign methods
    """
    os.system('Rscript scRNA-FeatureSelection/utils/RCode/classification.R')
    temp_path = 'scRNA-FeatureSelection/tempData/'
    classification_result = dict()
    label_test = np.loadtxt(temp_path + 'temp_y_test.csv', delimiter=',', skiprows=1, dtype=np.object)[:, 1]
    for assign_method in ['scmap_cluster', 'scmap_cell', 'singlecellnet']:
        try:
            label_pred = np.loadtxt(''.join([temp_path, 'temp_', assign_method, '.csv']), delimiter=',', skiprows=1,
                                    dtype=np.str)
            f1 = f1_score(np.squeeze(label_test), np.char.strip(label_pred, '"'), average='weighted')
            classification_result[assign_method + '_F1'] = f1
        except IOError:
            print('{} failed to execute. Please check the R output in console.'.format(assign_method))
    return classification_result


def evaluate_clustering_result():
    """
    Evaluate clustering result using the ARI generating from two clustering methods.

    :return: a dict containing ARI of two clustering methods
    """
    os.system('Rscript scRNA-FeatureSelection/utils/RCode/clustering.R')
    temp_path = 'scRNA-FeatureSelection/tempData/'
    clustering_result = dict()
    label_true = np.loadtxt(temp_path + 'temp_y.csv', delimiter=',', skiprows=1, usecols=[1], dtype=np.str)
    for clustering_method in ['seurat', 'sc3']:
        label_pred_file_name = ''.join([temp_path, 'temp_', clustering_method, '.csv'])
        label_pred = np.loadtxt(label_pred_file_name, delimiter=',', skiprows=1, dtype=np.str)
        ari = adjusted_rand_score(np.squeeze(label_true), np.squeeze(label_pred))
        clustering_result[clustering_method + '_ARI'] = ari
    return clustering_result


def evaluate_classification_methods(dataset: str, methods: list, data_type: str):
    """
    Evaluate classification methods on specific dataset (raw or norm) using 5-fold cross validation.
    The result is saved in results folder.

    :param dataset: string, data name
    :param methods: list of string, names of feature selection methods
    :param data_type: string, 'raw' or 'norm'
    :return: None
    """
    # load raw and norm data
    if dataset[:4] == 'PBMC' and 'scGeneFit' in methods:
        dataset = 'PBMC5'
        warnings.warn("Using 5% of PBMC cells because scGeneFit needs lots of system resource.", RuntimeWarning)
    X_raw, X_norm, y, trusted_markers = load_data(dataset)

    # prepare performance record
    performance_record = PerformanceRecord(methods=methods, task='assign')

    # output dataset information
    head(name='Dataset Information')
    print("Name:{}  Type:{}  Cell(s):{}  Gene(s):{}\nMarker Gene(s):{}".format(
        dataset, data_type, X_raw.shape[0], X_raw.shape[1],
        np.intersect1d(get_gene_names(X_raw.columns), trusted_markers).shape[0])
    )

    # 5-fold CV
    kf = KFold(n_splits=5, random_state=2020, shuffle=True)
    for train_idx, test_idx in kf.split(X_raw):
        # split and save raw data
        X_raw_train, X_raw_test = X_raw.iloc[train_idx], X_raw.iloc[test_idx]
        X_norm_train, X_norm_test = X_norm.iloc[train_idx], X_norm.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        delete('scRNA-FeatureSelection/tempData/')

        # calculate F1-score before feature selection
        save_raw_data(X_raw_train, X_raw_test, y_train, y_test, task='assign')
        before = evaluate_classification_result()

        # remove const genes before feature selection
        X_raw_train, X_norm_train = filter_const_genes(X_raw_train), filter_const_genes(X_norm_train)
        X_raw_test, X_norm_test = X_raw_test.loc[:, X_raw_train.columns], X_norm_test.loc[:, X_norm_train.columns]
        gene_names = get_gene_names(X_raw_train.columns)

        # feature selection
        for method in methods:
            # delete current files in tempData
            delete('scRNA-FeatureSelection/tempData/')
            if data_type == 'raw':
                result = select_features(dataset, 1000, method, gene_names, X_raw_train.values,
                                         np.squeeze(y_train.values))
            elif data_type == 'norm':
                result = select_features(dataset, 1000, method, gene_names, X_norm_train.values,
                                         np.squeeze(y_train.values))
            else:
                result = None
                warnings.warn("The parameter 'data_type' is wrong. Please check again.", RuntimeWarning)

            if len(result) == 2:  # method can generate feature importance
                markers_found, MRR = cal_marker_num_MRR(trusted_markers, result[0], rank=True)
                performance_record.loc['marker_genes_found', method] += markers_found
                performance_record.loc['MRR', method] += MRR
                mask = np.isin(gene_names, result[0])
            elif len(result) == 1:  # method can not generate feature importance
                markers_found = cal_marker_num_MRR(trusted_markers, result, rank=False)
                performance_record.loc['marker_genes_found', method] += markers_found
                warnings.warn("The method can not generate gene importance! MRR can not be calculated and is set to 0.",
                              RuntimeWarning)
                mask = np.isin(gene_names, result)
            else:
                warnings.warn("The length of feature selection result is 0 or more than 2.", RuntimeWarning)
                return None
            if mask.sum() == 0:
                warnings.warn("No gene is selected!", RuntimeWarning)
            # filter out non-markers
            X_train_selected, y_train_selected = filter_const_cells(X_raw_train.loc[:, mask], y_train)
            X_test_selected = X_raw_test.loc[:, mask]

            # save X_train and X_test after gene selection
            save_raw_data(X_train_selected, X_test_selected, y_train_selected, y_test, task='assign')

            # execute R script to generate classification result
            after = evaluate_classification_result()

            # update performance record
            print(before, after)
            for key in after.keys():
                performance_record.loc[key, method] += after[key] - before[key]

    # save performance record
    performance_record.divide(5).to_csv(
        ''.join(['scRNA-FeatureSelection/results/', dataset, '_', data_type, '_assign_record.csv']))
    print(now() + ": Evaluation is done!")


def evaluate_clustering_methods(dataset, methods, data_type):
    """
    Evaluate clustering methods on specific dataset (raw or norm). The result is saved in results folder.

    :param dataset: string, data name
    :param methods: list of string, names of feature selection methods
    :param data_type: string, 'raw' or 'norm'
    :return: None
    """
    # load raw and norm data
    if dataset[:4] == 'PBMC' and 'scGeneFit' in methods:
        dataset = 'PBMC5'
        warnings.warn("Using 5% of PBMC cells because scGeneFit needs lots of system resource.", RuntimeWarning)
    X_raw, X_norm, y, trusted_markers = load_data(dataset)
    gene_names = get_gene_names(X_raw.columns)

    # prepare performance record
    performance_record = PerformanceRecord(methods=methods, task='clustering')

    # output dataset information
    head(name='Dataset Information')
    print("Name:{}  Type:{}  Cell(s):{}  Gene(s):{}\nMarker Gene(s):{}".format(
        dataset, data_type, X_raw.shape[0], X_raw.shape[1], np.intersect1d(gene_names, trusted_markers).shape[0])
    )

    # save raw data and generate clustering result before feature selection
    delete('scRNA-FeatureSelection/tempData/')
    save_raw_data(X_train=X_raw, y_train=y, task='clustering')
    before = evaluate_clustering_result()

    for method in methods:
        delete('scRNA-FeatureSelection/tempData/')
        if data_type == 'raw':
            X = X_raw
        elif data_type == 'norm':
            X = X_norm
        else:
            warnings.warn("The parameter 'data_type' is wrong. Please check again.", RuntimeWarning)
            return None
        result = select_features(dataset, 1000, method, gene_names, X=X.values, y=np.squeeze(y.values))  # consistency
        if len(result) == 2:  # method can generate feature importance
            markers_found, MRR = cal_marker_num_MRR(trusted_markers, result[0], rank=True)
            performance_record.loc['marker_genes_found', method] = markers_found
            performance_record.loc['MRR', method] = MRR
            mask = np.isin(gene_names, result[0])
        elif len(result) == 1:  # method can not generate feature importance
            markers_found = cal_marker_num_MRR(trusted_markers, result, rank=False)
            performance_record.loc['marker_genes_found', method] = markers_found
            warnings.warn("The method can not generate gene importance! MRR can not be calculated and is set to 0.",
                          RuntimeWarning)
            mask = np.isin(gene_names, result)
        else:
            warnings.warn("The length of feature selection result is 0 or more than 2.", RuntimeWarning)
            return None
        if mask.sum() == 0:
            warnings.warn("No gene is selected!", RuntimeWarning)
        # filter out non-markers
        X_selected, y = filter_const_cells(X_raw.loc[:, mask], y)
        # save raw data and generate clustering result before feature selection
        save_raw_data(X_train=X_selected, y_train=y, task='clustering')
        after = evaluate_clustering_result()

        # update performance record
        print(before, after)
        for key in after.keys():
            performance_record.loc[key, method] = after[key] - before[key]

    # save performance record
    performance_record.to_csv(
        ''.join(['scRNA-FeatureSelection/results/', dataset, '_', data_type, '_clustering_record.csv']))
    print(now() + ": Evaluation is done!")
