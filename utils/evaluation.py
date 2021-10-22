import os
import traceback
from datetime import datetime
from typing import List

import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import scgen
import warnings
from sklearn.metrics import adjusted_rand_score, classification_report, cohen_kappa_score, v_measure_score
from sklearn.model_selection import StratifiedKFold

from config import assign_cfg, cluster_cfg, exp_cfg
from utils.importance import select_genes
from utils.record import MasterRecorder, ClassificationRecorder, ClusteringRecorder
from utils.classification import scANVI_predict
from utils.utils import load_data, delete, save_data, filter_adata, plot_2D, HiddenPrints


def evaluate_assign_result(recorder: ClassificationRecorder = None,
                           n_gene: int = None,
                           adata_train: ad.AnnData = None,
                           adata_test: ad.AnnData = None) -> dict:
    """
    Evaluate assign result using the F1-scores generating from three assign methods.

    :return: a dict containing F1-scores of three assign methods
    """
    os.system('Rscript utils/RCode/classification.R >& /dev/null')  #
    assign_result = dict()
    time_record = dict()
    # load test labels
    try:
        label_test = np.char.replace(np.squeeze(
            pd.read_feather("tempData/temp_y_test.feather").set_index('index', drop=True).values
        ).astype(np.str), '-', '.')
        # label_test = np.loadtxt('tempData/temp_y_test.csv', delimiter=',', skiprows=1, dtype=np.object_)[:, 1]
    except IndexError:  # if label_test is an empty array
        print(f"After selecting {n_gene} genes: y_test is empty.")
        label_test = np.empty(shape=(0,))

    for assign_method in ['singlecellnet', 'singleR', 'scanvi']:
        start_time = datetime.now()
        # generate pred labels
        if assign_method == 'scanvi':
            # try:
            #     label_pred = scANVI_predict(train_data=adata_train.copy(), test_data=adata_test.copy())
            # except (IndexError, ValueError, AttributeError) as e:
            #     print(f"After selecting {n_gene} genes: {e}")
            label_pred = np.empty(shape=(0,))
        else:
            try:
                label_pred = np.char.strip(np.loadtxt(''.join(['tempData/', 'temp_', assign_method, '.csv']),
                                                      delimiter=',', skiprows=1, dtype=np.str), '"')
            except OSError:
                print(f"temp_{assign_method}.csv not found.")
                label_pred = np.empty(shape=(0,))

        # calculate F1 score
        if label_pred.shape[0] == 0:  # label_pred is an empty array
            print(f"After selecting {n_gene} genes: {assign_method} failed. Set F1 (rare) score to 0.")
            f1_all, ck_score, f1_rare = 0, 0, 0
        else:  # label_pred is not empty
            if label_test.shape[0] != 0:
                if label_test.shape[0] != label_pred.shape[0]:
                    raise ValueError(f"label_test is {label_test.shape[0]}, and label_pred is {label_pred.shape[0]}")
                report = classification_report(label_test, label_pred, output_dict=True, zero_division=0)
                f1_all = report['macro avg']['f1-score']  # scANVI is bad in rare cell detection
                ck_score = cohen_kappa_score(label_test, np.char.strip(label_pred, '"'))
                if hasattr(recorder, 'rare_type'):
                    if recorder.rare_type not in report.keys():
                        print(f"After selecting {n_gene} genes: "
                              f"{recorder.rare_type} is not in report. Set F1-rare to 0.")
                        f1_rare = 0
                    else:
                        f1_rare = report[recorder.rare_type]['f1-score']
                else:
                    f1_rare = 0
            else:  # test labels do not exist
                print('label_test does not exist!')
                f1_all, ck_score, f1_rare = 0, 0, 0
        assign_result[assign_method + '_f1'] = f1_all
        assign_result[assign_method + '_ck'] = ck_score
        if hasattr(recorder, 'rare_type'):
            assign_result[assign_method + '_f1_rare'] = f1_rare
        time_record[assign_method] = np.round((datetime.now() - start_time).total_seconds(), decimals=2)

    if n_gene == exp_cfg.n_genes[-1]:
        print(f"{n_gene}: Evaluation of gene selection results costs: {time_record}")
        print(f"{n_gene}: Evaluation results: {assign_result}")

    return assign_result


def evaluate_cluster_result(recorder: ClusteringRecorder = None,
                            n_gene: int = None
                            ) -> dict:
    """
    Evaluate clustering result using the ARI generating from three clustering methods.

    :return: a dict containing ARI of two clustering methods
    """
    os.system('Rscript utils/RCode/clustering.R >& /dev/null')
    cluster_result = dict()
    time_record = dict()
    label_true = np.squeeze(np.loadtxt('tempData/temp_y.csv', delimiter=',', skiprows=1, usecols=[1], dtype=np.str_))
    for cluster_method in ['seurat', 'sc3', 'cidr']:
        start_time = datetime.now()
        label_pred_file_name = ''.join(['tempData/', 'temp_', cluster_method, '.csv'])
        try:
            label_pred = np.squeeze(np.loadtxt(label_pred_file_name, delimiter=',', skiprows=1, dtype=np.str_))
            cluster_result[cluster_method + '_ARI'] = adjusted_rand_score(label_true, label_pred)
            cluster_result[cluster_method + '_BCubed'] = recorder.BCubed_fbeta_score(label_true, label_pred)
            cluster_result[cluster_method + '_V'] = v_measure_score(label_true, label_pred)
            if hasattr(recorder, 'rare_type'):
                cluster_result[cluster_method + '_BCubed_rare'] = recorder.BCubed_fbeta_score_rare(label_true,
                                                                                                   label_pred)
        except OSError:
            print(f"After selecting {n_gene} genes: {cluster_method} failed."
                  f" Set ARI and F1 of the rare cell type to 0.")
        time_record[cluster_method] = np.round((datetime.now() - start_time).total_seconds(), decimals=2)
    print(f"{n_gene}: Evaluation of gene selection results costs: {time_record}")
    return cluster_result


def evaluate_correction_result(combined_adata: ad.AnnData, fs_method: str = None, data_name: str = None):
    combined_adata = combined_adata[np.argsort(combined_adata.obs['batch']), :]
    correction_result = dict()
    with HiddenPrints():
        for correct_method in ['harmony', 'scanorama']:  # 'scgen',
            if correct_method == 'scgen':
                train = scgen.setup_anndata(combined_adata.raw.to_adata(), batch_key="batch", labels_key="celltype",
                                            copy=True)
                model = scgen.SCGEN(train)
                model.train(max_epochs=exp_cfg.epochs, batch_size=exp_cfg.batch_size,
                            early_stopping=exp_cfg.use_early_stopping,
                            early_stopping_patience=exp_cfg.patience, use_gpu=False)
                corrected = model.batch_removal()
            elif correct_method == 'harmony':
                sc.external.pp.harmony_integrate(combined_adata, key='batch', random_state=exp_cfg.random_seed,
                                                 verbose=0, nclust=50)  # basis: X_pca, so gene expression after batch
                corrected = combined_adata.copy()  # correction is similar, and thus iLISI values are similar
            elif correct_method == 'scanorama':
                sc.external.pp.scanorama_integrate(combined_adata, key='batch', verbose=0)  # basis: X_pca
                corrected = combined_adata.copy()
            else:
                raise NotImplementedError(f"{correct_method} have not been implemented yet.")
            plot_2D(corrected, dr_method='umap', fs_method=fs_method, bc_method=correct_method, data_name=data_name,
                    mode='after')
            plot_2D(corrected, dr_method='tsne', fs_method=fs_method, bc_method=correct_method, data_name=data_name,
                    mode='after')
            save_data(corrected, task='correct')
            os.system('Rscript utils/RCode/correction.R >& /dev/null')  #
            metrics = np.loadtxt('tempData/temp_correction.csv', skiprows=1, usecols=[1, 2], delimiter=',', dtype=np.float_)
            correction_result[correct_method + '_kBET'] = metrics[0]
            correction_result[correct_method + '_iLISI'] = metrics[1]
            print(correction_result)
    return correction_result


def evaluate_feature_selection_methods(measurements: List[str], methods: List[str]):
    print(f"measurements: {measurements}")
    master_recorder = MasterRecorder(measurements, methods)
    try:
        for measurement, datasets in exp_cfg.measurements.items():
            if measurement in measurements:
                if measurement == 'population_demixing':
                    for dataset in datasets:
                        delete('tempData/')
                        adata = load_data(dataset)
                        for method in methods:
                            master_recorder.mixture.set_current_method(method)
                            result_list = select_genes(method, adata, config=assign_cfg)
                            if result_list is not None:
                                for n_gene, result in zip(exp_cfg.n_genes, result_list):
                                    gene_mask = master_recorder.mixture.get_mask(adata.var_names.values, result, n_gene)
                                    selected_adata = filter_adata(adata.raw[:, gene_mask].to_adata(), filter_cell=True)
                                    master_recorder.mixture.record(dataset, n_gene, selected_adata)
                    master_recorder.mixture.save()
                elif measurement == 'marker_discovery':
                    for dataset in datasets:
                        delete('tempData/')
                        adata = load_data(dataset)
                        for method in methods:
                            master_recorder.marker.set_current_method(method)
                            result_list = select_genes(method, adata, config=assign_cfg)
                            if result_list is not None:
                                for n_gene, result in zip(exp_cfg.n_genes, result_list):
                                    master_recorder.marker.record(dataset, n_gene, adata.uns['markers'], adata.uns['marker_weight'], result)
                    master_recorder.marker.save()
                elif measurement == 'computation_time':
                    for dataset in datasets:
                        delete('tempData/')
                        adata = load_data(dataset)
                        for method in methods:
                            master_recorder.time.set_current_method(method)
                            select_genes(method, adata, recorder=master_recorder.time, config=assign_cfg)
                    master_recorder.time.save()
                elif measurement == 'batch_correction':
                    for comb_dataset in datasets:
                        comb_adata = load_data(comb_dataset)
                        for method in methods:
                            master_recorder.correct.set_current_method(method)
                            delete('tempData/')
                            # select features
                            all_result_list = select_genes(method, comb_adata, config=assign_cfg)
                            if all_result_list is not None:
                                for n_gene, result in zip(exp_cfg.n_genes, all_result_list):
                                    gene_mask = master_recorder.correct.get_mask(comb_adata.var_names.values, result,
                                                                                 n_gene)
                                    selected_comb = filter_adata(comb_adata[:, gene_mask], filter_cell=True)
                                    selected_comb.raw = comb_adata.raw.to_adata()[:, gene_mask]
                                    # print(n_gene, selected_comb.shape, selected_comb.raw.shape)
                                    if selected_comb.shape[0] < 1e3:
                                        warnings.warn("This iteration will be skipped.")
                                        continue
                                    plot_2D(selected_comb, dr_method='umap', fs_method=method, data_name=comb_dataset)
                                    plot_2D(selected_comb, dr_method='tsne', fs_method=method, data_name=comb_dataset)
                                    correction_result = evaluate_correction_result(selected_comb, fs_method=method,
                                                                                   data_name=comb_dataset)
                                    master_recorder.correct.record(comb_dataset, n_gene, correction_result)
                            master_recorder.correct.save()
                elif measurement == 'classification':
                    for dataset in datasets:
                        adata = load_data(dataset)
                        master_recorder.assign.set_rare_type(adata)
                        skf = StratifiedKFold(n_splits=assign_cfg.n_folds, random_state=assign_cfg.random_seed,
                                              shuffle=True)
                        for i, (train_idx, test_idx) in enumerate(skf.split(adata.X, adata.obs['celltype'].values)):

                            # split train and test data
                            adata_train, adata_test = adata[train_idx], adata[test_idx]
                            # remove const genes and cells before feature selection
                            gene_mask = sc.pp.filter_genes(adata_train, min_cells=exp_cfg.n_filter_cell, inplace=False)[0]
                            adata_train, adata_test = adata_train[:, gene_mask], adata_test[:, gene_mask]
                            adata_train.raw = adata_train.raw[:, gene_mask].to_adata()
                            adata_test.raw = adata_test.raw[:, gene_mask].to_adata()
                            adata_train, adata_test = filter_adata(adata_train, filter_cell=True), \
                                                      filter_adata(adata_test, filter_cell=True)
                            for method in methods:
                                master_recorder.assign.set_current_method(method, i + 1)
                                # select features
                                all_result_list = select_genes(method, adata_train, config=assign_cfg)
                                # store selected genes, 4 n_gene * 5 folds
                                if all_result_list is not None:
                                    master_recorder.assign.store_selected_genes(dataset, all_result_list)
                                    for n_gene, result in zip(exp_cfg.n_genes, all_result_list):
                                        # clean directory
                                        delete('tempData/')
                                        gene_mask = master_recorder.assign.get_mask(adata_train.var_names.values,
                                                                                    result,
                                                                                    n_gene)
                                        # filter out non-HVGs and save raw data
                                        selected_train = filter_adata(adata_train.raw[:, gene_mask].to_adata(),
                                                                      filter_cell=True)  # raw
                                        selected_test = filter_adata(adata_test.raw[:, gene_mask].to_adata(),
                                                                     filter_cell=True)  # raw
                                        print(selected_test.shape)
                                        save_data(selected_train, selected_test)
                                        # metric of downstream analysis
                                        assign_result = evaluate_assign_result(master_recorder.assign, n_gene,
                                                                               selected_train, selected_test)
                                        master_recorder.assign.record(dataset, n_gene, i, assign_result)
                        master_recorder.assign.cal_mean_rbo(dataset)
                        master_recorder.assign.save()
                elif measurement == 'clustering':
                    for dataset in datasets:
                        adata = load_data(dataset)
                        master_recorder.cluster.set_rare_type(adata)

                        for i in range(cluster_cfg.n_loops):
                            skf = StratifiedKFold(n_splits=cluster_cfg.n_folds,
                                                  random_state=cluster_cfg.random_seed + i,
                                                  shuffle=True)
                            left_index, right_index = list(skf.split(adata.X, adata.obs['celltype'].values))[0]

                            # remove const genes before feature selection
                            adata_left = filter_adata(adata[left_index], filter_gene=True, filter_cell=True)
                            adata_right = filter_adata(adata[right_index], filter_gene=True, filter_cell=True)
                            # clean directory
                            delete('tempData/')
                            for method in master_recorder.cluster.methods:
                                master_recorder.cluster.set_current_method(method, i + 1)
                                for j, data in enumerate([adata_left, adata_right]):
                                    half_result_list = select_genes(method, data, config=cluster_cfg)
                                    if half_result_list is not None:
                                        for n_gene, half_result in zip(exp_cfg.n_genes, half_result_list):
                                            half_mask = master_recorder.cluster.get_mask(
                                                all_genes=data.var_names.values,
                                                single_selected_result=half_result,
                                                n_gene=n_gene)
                                            half_selected = filter_adata(data.raw[:, half_mask].to_adata(),
                                                                         filter_gene=True, filter_cell=True)
                                            save_data(half_selected, task='cluster')
                                            # downstream analysis: clustering
                                            half_metric = evaluate_cluster_result(master_recorder.cluster, n_gene)
                                            # double_metrics[j].append(half_metric)
                                            master_recorder.cluster.record(dataset, n_gene, i, j, half_metric)
                        master_recorder.cluster.save()
    except:
        master_recorder.save_all()
        traceback.print_exc()
