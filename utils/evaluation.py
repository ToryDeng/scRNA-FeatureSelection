import os
import time

import numpy as np
import scanpy as sc
import anndata as ad
from sklearn.metrics import adjusted_rand_score, classification_report
from sklearn.model_selection import StratifiedKFold

from config import assign_cfg, cluster_cfg, exp_cfg
from utils.importance import select_genes
from utils.record import PerformanceRecorder
from utils.itclust import ItClust_predict
from utils.utils import load_data, delete, save_data, filter_adata


def evaluate_assign_result(recorder: PerformanceRecorder = None,
                           adata_train: ad.AnnData = None,
                           adata_test: ad.AnnData = None) -> dict:
    """
    Evaluate assign result using the F1-scores generating from three assign methods.

    :return: a dict containing F1-scores of three assign methods
    """
    os.system('Rscript utils/RCode/classification.R >& /dev/null')
    assign_result = dict()
    time_record = dict()
    # load test labels
    try:
        label_test = np.loadtxt('tempData/temp_y_test.csv', delimiter=',', skiprows=1, dtype=np.object_)[:, 1]
    except IndexError:  # if label_test is an empty array
        print(f"After selecting {recorder.n_genes} genes: y_test is empty.")
        label_test = None

    for assign_method in ['singlecellnet', 'singleR', 'itclust']:
        start_time = time.time()
        # generate pred labels
        if assign_method == 'itclust':
            try:
                label_pred = ItClust_predict(train_data=adata_train, test_data=adata_test)
            except (IndexError, ValueError, AttributeError) as e:
                print(f"After selecting {recorder.n_genes} genes: {e}")
                label_pred = np.empty(shape=(0,))
        else:
            try:
                label_pred = np.loadtxt(''.join(['tempData/', 'temp_', assign_method, '.csv']),
                                        delimiter=',', skiprows=1, dtype=np.str_)
            except OSError:
                label_pred = np.empty(shape=(0,))

        # calculate F1 score
        if label_pred.shape[0] == 0:  # label_pred is an empty array
            print(f"After selecting {recorder.n_genes} genes: {assign_method} failed. Set F1 (rare) score to 0.")
            f1_all, f1_rare = 0, 0
        else:
            if label_test is not None:
                report = classification_report(np.squeeze(label_test), np.char.strip(label_pred, '"'),
                                               output_dict=True, zero_division=0)
                f1_all = report['weighted avg']['f1-score']
                if hasattr(recorder, 'rare_type') and recorder.rare_type in report.keys():
                    if recorder.rare_type not in report.keys():
                        print(f"After selecting {recorder.n_genes} genes: "
                              f"{recorder.rare_type} is not in report. Set F1-rare to 0.")
                        f1_rare = 0
                    else:
                        f1_rare = report[recorder.rare_type]['f1-score']
                else:
                    f1_rare = 0
            else:  # test labels do not exist
                f1_all, f1_rare = 0, 0
        assign_result[assign_method + '_F1'] = f1_all
        if hasattr(recorder, 'rare_type'):
            assign_result[assign_method + '_F1_rare'] = f1_rare
        time_record[assign_method] = np.round(time.time() - start_time, decimals=2)
    print(f"{recorder.n_genes}: Evaluation of gene selection results costs: {time_record}")
    return assign_result


def evaluate_cluster_result(recorder: PerformanceRecorder = None) -> dict:
    """
    Evaluate clustering result using the ARI generating from two clustering methods.

    :return: a dict containing ARI of two clustering methods
    """
    os.system('Rscript utils/RCode/clustering.R >& /dev/null')
    cluster_result = dict()
    time_record = dict()
    label_true = np.squeeze(np.loadtxt('tempData/temp_y.csv', delimiter=',', skiprows=1, usecols=[1], dtype=np.str_))
    for cluster_method in ['seurat', 'sc3', 'cidr']:
        start_time = time.time()
        label_pred_file_name = ''.join(['tempData/', 'temp_', cluster_method, '.csv'])
        try:
            label_pred = np.squeeze(np.loadtxt(label_pred_file_name, delimiter=',', skiprows=1, dtype=np.str_))
            cluster_result[cluster_method + '_ARI'] = adjusted_rand_score(label_true, label_pred)
            if hasattr(recorder, 'rare_type'):
                if cluster_cfg.rare_type_detection_metric == 'f-measure':
                    cluster_result[cluster_method + '_F1_rare'] = recorder.fbeta_score(label_true, label_pred)
                elif cluster_cfg.rare_type_detection_metric == 'bcubed':
                    cluster_result[cluster_method + '_F1_rare'] = recorder.BCubed_score(label_true, label_pred)
                else:
                    raise NotImplementedError(f"{cluster_cfg.rare_type_detection_metric} has not been implemented yet.")
        except OSError:
            print(f"After selecting {recorder.n_genes} genes: {cluster_method} failed."
                  f" Set ARI and F1 of the rare cell type to 0.")
        time_record[cluster_method] = np.round(time.time() - start_time, decimals=2)
    print(f"{recorder.n_genes}: Evaluation of gene selection results costs: {time_record}")
    return cluster_result


def evaluate_assign_methods(dataset: str, methods: list) -> None:
    """
    Only this function can handle Timeout error because scGeneFit is a supervised method.

    :param dataset: str, name of data
    :param methods: list, methods to evaluate
    :return: None
    """
    # load raw and norm data
    adata = load_data(dataset)
    recorders = [PerformanceRecorder('assign', dataset, methods, adata, n_gene=n_gene) for n_gene in exp_cfg.n_genes]
    recorders[0].show_dataset_info()

    skf = StratifiedKFold(n_splits=assign_cfg.n_folds, random_state=assign_cfg.random_seed, shuffle=True)
    for i, (train_idx, test_idx) in enumerate(skf.split(adata.X, adata.obs['celltype'].values)):
        # clean directory
        delete('tempData/')
        # split train and test data
        adata_train, adata_test = adata[train_idx], adata[test_idx]
        # remove const genes and cells before feature selection
        gene_mask = sc.pp.filter_genes(adata_train, min_cells=exp_cfg.n_filter_cell, inplace=False)[0]
        adata_train, adata_test = adata_train[:, gene_mask], adata_test[:, gene_mask]
        adata_train.raw = adata_train.raw[:, gene_mask].to_adata()
        adata_test.raw = adata_test.raw[:, gene_mask].to_adata()
        adata_train, adata_test = filter_adata(adata_train, filter_cell=True), filter_adata(adata_test,
                                                                                            filter_cell=True)
        for method in methods:
            # init each gene-specific recorder
            for recorder in recorders:
                recorder.init_current_feature_selection_method(method, i)
            # select features using a kind of single or ensemble method
            result_list = select_genes(method, adata_train, recorders=recorders, config=assign_cfg)
            if result_list is not None:
                for recorder, result in zip(recorders, result_list):
                    # record markers_found_and_MRR
                    recorder.record_markers_found_and_MRR(result)
                    # get gene mask
                    marker_mask = recorder.get_mask(adata_train.var_names.values, result)
                    # filter out non-markers and save raw data
                    selected_train = filter_adata(adata_train.raw[:, marker_mask].to_adata(), filter_cell=True)
                    selected_test = filter_adata(adata_test.raw[:, marker_mask].to_adata(), filter_cell=True)
                    save_data(selected_train, selected_test)
                    # metric of downstream analysis
                    assign_metric = evaluate_assign_result(recorder, selected_train, selected_test)
                    recorder.record_metrics_from_dict(assign_metric)
    for recorder in recorders:
        recorder.summary(save=True, show=True)
        recorder.save()


def evaluate_cluster_methods(dataset: str, methods: list) -> None:
    # load raw and norm data
    adata = load_data(dataset)
    recorders = [PerformanceRecorder('cluster', dataset, methods, adata, n_gene=n_gene) for n_gene in exp_cfg.n_genes]
    recorders[0].show_dataset_info()

    for i in range(cluster_cfg.n_loops):
        skf = StratifiedKFold(n_splits=cluster_cfg.n_folds, random_state=cluster_cfg.random_seed + i, shuffle=True)
        left_index, right_index = list(skf.split(adata.X, adata.obs['celltype'].values))[0]

        # remove const genes before feature selection
        adata_left = filter_adata(adata[left_index], filter_gene=True, filter_cell=True)
        adata_right = filter_adata(adata[right_index], filter_gene=True, filter_cell=True)
        # clean directory
        delete('tempData/')
        for method in methods:
            # init each gene-specific recorder
            for recorder in recorders:
                recorder.init_current_feature_selection_method(method, i)
            # select features
            all_result_list = select_genes(method, adata, recorders=recorders, config=cluster_cfg)
            # record markers_found_and_MRR
            for recorder, all_result in zip(recorders, all_result_list):
                # markers found and MRR are generated from intact dataset
                recorder.record_markers_found_and_MRR(all_result)

            double_metrics = [[], []]  # left_metrics, right_metrics
            for j, data in enumerate([adata_left, adata_right]):
                half_result_list = select_genes(method, data, recorders=recorders, config=cluster_cfg)
                for recorder, half_result in zip(recorders, half_result_list):
                    half_mask = recorder.get_mask(data.var_names.values, half_result)
                    half_selected = filter_adata(data.raw[:, half_mask].to_adata(), filter_gene=True, filter_cell=True)
                    save_data(half_selected, task='cluster')
                    half_metric = evaluate_cluster_result(recorder)
                    double_metrics[j].append(half_metric)

            for recorder, left_metric, right_metric in zip(recorders, double_metrics[0], double_metrics[1]):
                recorder.record_metrics_from_dict([left_metric, right_metric])

    for recorder in recorders:
        recorder.summary(save=True, show=True)
        recorder.save()
