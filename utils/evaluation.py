import os
import traceback
from datetime import datetime
from typing import List, Union
import warnings
import numpy as np
import scanpy as sc
import anndata as ad
from sklearn.metrics import adjusted_rand_score, classification_report, cohen_kappa_score, v_measure_score
from sklearn.model_selection import StratifiedKFold
from config import assign_cfg, cluster_cfg, exp_cfg
from utils.importance import select_genes
from utils.record import MasterRecorder, ClassificationRecorder, ClusteringRecorder, BatchCorrectionRecorder
from utils.utils import load_data, delete, filter_adata, plot_2D, compute_correction_metrics, head, MyGroupSplit
from otherSteps.clustering import seurat_v4_clustering, SC3s_clustering, SHARP_clustering, SC3_clustering2
from otherSteps.classification import SingleR
from otherSteps.correction import Seurat_v4_correct


def evaluate_assign_result(adata_train: ad.AnnData,
                           adata_test: ad.AnnData,
                           recorder: ClassificationRecorder = None,
                           n_gene: Union[int, str] = None,
                           ) -> dict:
    """
    Evaluate assign result using the F1-scores generating from three assign methods.

    :return: a dict containing F1-scores of three assign methods
    """
    if n_gene == 'all':
        print('Calculating baseline...')

    assign_result = dict()
    time_record = dict()
    # load test labels
    label_test = adata_test.obs['celltype'].to_numpy()

    for assign_method in assign_cfg.evaluation_method:
        start_time = datetime.now()
        # generate pred labels
        try:
            if assign_method == 'SingleR':
                label_pred = SingleR(adata_train=adata_train.raw.to_adata(), adata_test=adata_test.raw.to_adata())
            else:
                print(f"{assign_method} is not implemented.")
                label_pred = np.empty(shape=(0,))
        except:
            print(f"{assign_method} failed after selecting {n_gene} genes.")
            traceback.print_exc()
            label_pred = np.empty(shape=(0,))

        # calculate F1 score
        if label_pred.shape[0] == 0:  # label_pred is an empty array
            f1_all, ck_score, f1_rare = np.nan, np.nan, np.nan
        else:  # label_pred is not empty
            if label_test.shape[0] != 0:
                if label_test.shape[0] != label_pred.shape[0]:
                    raise ValueError(f"label_test is {label_test.shape[0]}, and label_pred is {label_pred.shape[0]}")
                report = classification_report(label_test, label_pred, output_dict=True, zero_division=0)
                f1_all = report['macro avg']['f1-score']
                ck_score = cohen_kappa_score(label_test, np.char.strip(label_pred, '"'))
                if recorder.rare_type is not None:
                    if recorder.rare_type not in report.keys():
                        print(f"After selecting {n_gene} genes: "
                              f"{recorder.rare_type} is not in report.")
                        f1_rare = np.nan
                    else:
                        f1_rare = report[recorder.rare_type]['f1-score']
                else:
                    f1_rare = np.nan
            else:  # test labels do not exist
                print('label_test does not exist!')
                f1_all, ck_score, f1_rare = np.nan, np.nan, np.nan
        assign_result[assign_method + '_f1'] = f1_all
        assign_result[assign_method + '_ck'] = ck_score
        if hasattr(recorder, 'rare_type'):
            assign_result[assign_method + '_f1_rare'] = f1_rare
        time_record[assign_method] = np.round((datetime.now() - start_time).total_seconds(), decimals=2)

    if n_gene == 'all':
        print(f"Baseline: {assign_result}")
    if n_gene == exp_cfg.n_genes[-1]:
        print(f"{n_gene} genes: Evaluation of gene selection results costs: {time_record}")
        print(f"{n_gene} genes: Evaluation results: {assign_result}")

    return assign_result


def evaluate_cluster_result(adata: ad.AnnData, recorder: ClusteringRecorder = None) -> dict:
    """
    Evaluate clustering result using the ARI generating from three clustering methods.

    :return: a dict containing ARI and v-measure of  clustering methods
    """
    cluster_result = dict()
    time_record = dict()
    label_true = adata.obs['celltype'].to_numpy()
    for cluster_method in cluster_cfg.evaluation_method:
        start_time = datetime.now()
        try:
            if cluster_method == 'seurat':
                label_pred = seurat_v4_clustering(adata.raw.to_adata())  # need raw data
            elif cluster_method == 'sc3s':
                label_pred = SC3s_clustering(adata)  # need norm data
            elif cluster_method == 'sc3':
                label_pred = SC3_clustering2(adata.raw.to_adata())
            elif cluster_method == 'sharp':
                label_pred = SHARP_clustering(adata.raw.to_adata())
            else:
                print(f"{cluster_method} is not implemented.")
                label_pred = np.empty(shape=(0,))
        except:
            print(f"After selecting {adata.n_vars} genes: {cluster_method} failed.")
            traceback.print_exc()
            label_pred = np.empty(shape=(0,))
        else:
            print(f"{cluster_method} has been finished.")
        if label_pred.shape[0] != 0:
            if all(label_pred == label_pred):
                cluster_result[cluster_method + '_ARI'] = adjusted_rand_score(label_true, label_pred)
                # cluster_result[cluster_method + '_BCubed'] = recorder.BCubed_fbeta_score(label_true, label_pred)
                cluster_result[cluster_method + '_V'] = v_measure_score(label_true, label_pred)
                if hasattr(recorder, 'rare_type'):
                    cluster_result[cluster_method + '_BCubed_rare'] = recorder.BCubed_fbeta_score_rare(label_true,
                                                                                                       label_pred)
            else:
                print("Prediction result contains abnormal values (np.nan).")
        time_record[cluster_method] = np.round((datetime.now() - start_time).total_seconds(), decimals=2)

    if adata.n_vars == exp_cfg.n_genes[-1]:
        print(f"{adata.n_vars} genes: Evaluation of selection results costs: {time_record}")
        print(f"{adata.n_vars} genes: Evaluation results: {cluster_result}\n")
    return cluster_result


def evaluate_batch_correction(combined_adata: ad.AnnData, fs_method: str, recorder: BatchCorrectionRecorder,
                              pre_corrected: ad.AnnData):
    warnings.filterwarnings('ignore')
    # do batch correction first
    # then select features
    delete('tempData/')
    print('selecting features...')
    all_result_list_correction_first = select_genes(fs_method, pre_corrected, config=assign_cfg, select_by_batch=False)
    # evaluate the results
    correction_first_result = list()
    for n_gene, result in zip(exp_cfg.n_genes, all_result_list_correction_first):
        gene_mask = recorder.get_mask(pre_corrected.var_names.values, result, n_gene)
        selected_correct = pre_corrected[:, gene_mask]
        selected_correct.raw = pre_corrected.raw.to_adata()[:, gene_mask]
        sc.pp.pca(selected_correct)
        correction_first_result.append(compute_correction_metrics(selected_correct))
        plot_2D(selected_correct, dr_method='umap', fs_method=fs_method, bc_method='seuratV4', mode='after',
                order='correction_first')
        plot_2D(selected_correct, dr_method='tsne', fs_method=fs_method, bc_method='seuratV4', mode='after',
                order='correction_first')
    del pre_corrected, selected_correct
    # do feature selection first
    delete('tempData/')
    print('selecting features...')
    all_result_list_selection_first = select_genes(fs_method, combined_adata, config=assign_cfg)
    selection_first_result = list()
    for n_gene, result in zip(exp_cfg.n_genes, all_result_list_selection_first):
        gene_mask = recorder.get_mask(combined_adata.var_names.values, result, n_gene)
        selected_comb = combined_adata[:, gene_mask]
        selected_comb.raw = combined_adata.raw.to_adata()[:, gene_mask]
        # then correct data
        print('correcting data...')
        corrected_adata = Seurat_v4_correct(selected_comb)
        # evaluate the results
        sc.pp.pca(corrected_adata)
        selection_first_result.append(compute_correction_metrics(corrected_adata))
        plot_2D(corrected_adata, dr_method='umap', fs_method=fs_method, bc_method='seuratV4', mode='after',
                order='selection_first')
        plot_2D(corrected_adata, dr_method='tsne', fs_method=fs_method, bc_method='seuratV4', mode='after',
                order='selection_first')

    return correction_first_result, selection_first_result  # [{}, {}],[{}, {}]


def evaluate_feature_selection_methods(measurements: List[str], methods: List[str]):
    sc.settings.verbosity = 0
    master_recorder = MasterRecorder(measurements, methods)
    print(f"Process pid: {os.getpid()}")
    try:
        for measurement, datasets in exp_cfg.measurements.items():
            if measurement in measurements:
                print(f"Current measurement: {measurement}")
                if measurement == 'population_demixing':
                    for dataset in datasets:
                        delete('tempData/')
                        adata = load_data(dataset)
                        master_recorder.mixture.record(dataset, filtered_adata=adata)
                        for method in methods:
                            master_recorder.mixture.set_current_method(method)
                            result_list = select_genes(method, adata, config=assign_cfg)
                            if result_list is not None:
                                for n_gene, result in zip(exp_cfg.n_genes, result_list):
                                    gene_mask = master_recorder.mixture.get_mask(adata.var_names.values, result, n_gene)
                                    selected_adata = adata.raw[:, gene_mask].to_adata()
                                    master_recorder.mixture.record(dataset, selected_adata, n_gene)
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
                                    master_recorder.marker.record(
                                        dataset, adata.uns['markers'], adata.uns['marker_weight'], result, n_gene
                                    )
                        master_recorder.marker.save()
                elif measurement == 'computation_time':
                    delete(os.path.join(exp_cfg.record_path, 'pkl', 'TimeRecorder.pkl'))
                    for dataset in datasets:
                        delete('tempData/')
                        adata = load_data(dataset)
                        for method in methods:
                            master_recorder.time.set_current_method(method)
                            select_genes(method, adata, master_recorder.time, config=assign_cfg, use_saved=False)
                elif measurement == 'batch_correction':
                    for comb_dataset in datasets:
                        comb_adata = load_data(comb_dataset)
                        # before batch correction and feature selection
                        sc.pp.pca(comb_adata)
                        plot_2D(comb_adata, dr_method='umap')
                        plot_2D(comb_adata, dr_method='tsne')
                        del comb_adata.obsm['X_pca']

                        # do batch correction first
                        print('correcting data first...')
                        if f"{comb_adata.uns['data_name']}.h5ad" not in os.listdir('correctedData/'):
                            precorrected = Seurat_v4_correct(comb_adata)  # 'batch' in obs
                            print('Saving pre-corrected data...')
                            sc.write(f"correctedData/{precorrected.uns['data_name']}.h5ad", precorrected)
                        else:
                            print('Using previously saved data...')
                            precorrected = sc.read_h5ad(f"correctedData/{comb_adata.uns['data_name']}.h5ad")

                        for method in methods:
                            master_recorder.correct.set_current_method(method)
                            correction_first, selection_first = evaluate_batch_correction(
                                comb_adata, method, master_recorder.correct, pre_corrected=precorrected)
                            for i, n_gene in enumerate(exp_cfg.n_genes):
                                if i < len(correction_first):
                                    master_recorder.correct.record(comb_dataset, n_gene, correction_first[i],
                                                                   'correction_first')
                                if i < len(selection_first):
                                    master_recorder.correct.record(comb_dataset, n_gene, selection_first[i],
                                                                   'selection_first')
                            master_recorder.correct.save()
                elif measurement.endswith('classification'):
                    assign_type = measurement.split('-')[0]
                    if assign_type == 'intra':
                        select_by_batch = True
                        split_obj = StratifiedKFold(assign_cfg.n_folds, random_state=assign_cfg.random_seed, shuffle=True)
                    else:  # inter
                        select_by_batch = False
                        split_obj = MyGroupSplit()

                    for dataset in datasets:
                        adata = load_data(dataset)
                        getattr(master_recorder, assign_type).set_rare_type_and_tables(adata)
                        for i, (train_idx, test_idx) in enumerate(split_obj.split(adata.X, adata.obs['celltype'].values, adata.obs['batch'].values)):
                            i = master_recorder.inter.perms[i] if assign_type == 'inter' else i
                            # split train and test data
                            adata_train, adata_test = adata[train_idx].copy(), adata[test_idx].copy()
                            adata_train.uns['fold'], adata_test.uns['fold'] = i, i
                            # remove const genes and cells before feature selection
                            gene_mask = sc.pp.filter_genes(adata_train, min_cells=exp_cfg.n_filter_cell, inplace=False)[0]
                            adata_train, adata_test = adata_train[:, gene_mask], adata_test[:, gene_mask]
                            adata_train.raw = adata_train.raw[:, gene_mask].to_adata()
                            adata_test.raw = adata_test.raw[:, gene_mask].to_adata()
                            adata_train, adata_test = filter_adata(adata_train, filter_cell=True), \
                                                      filter_adata(adata_test, filter_cell=True)
                            # baseline
                            baseline_result = evaluate_assign_result(adata_train, adata_test, getattr(master_recorder, assign_type), 'all')
                            getattr(master_recorder, assign_type).record(dataset, i, baseline_result)
                            # do feature selection
                            for method in methods:
                                getattr(master_recorder, assign_type).set_current_method(method, i)
                                # select features
                                all_result_list = select_genes(method, adata_train, config=assign_cfg, select_by_batch=select_by_batch)
                                # store selected genes, 4 n_gene * 5 folds
                                if all_result_list is not None:
                                    for n_gene, result in zip(exp_cfg.n_genes, all_result_list):
                                        # clean directory
                                        delete('tempData/')
                                        gene_mask = getattr(master_recorder, assign_type).get_mask(
                                            adata_train.var_names.values, result, n_gene)
                                        # filter out non-HVGs and save raw data
                                        selected_train, selected_test = adata_train[:, gene_mask], adata_test[:,
                                                                                                   gene_mask]
                                        selected_train.raw = adata_train.raw[:, gene_mask].to_adata()  # norm and raw
                                        selected_test.raw = adata_test.raw[:, gene_mask].to_adata()  # norm and raw
                                        # metric of downstream analysis
                                        assign_result = evaluate_assign_result(selected_train, selected_test,
                                                                               getattr(master_recorder, assign_type), n_gene)
                                        getattr(master_recorder, assign_type).record(dataset, i, assign_result, n_gene)
                        getattr(master_recorder, assign_type).save()
                elif measurement == 'clustering':
                    for dataset in datasets:
                        adata = load_data(dataset)
                        master_recorder.cluster.set_rare_type(adata)
                        # run baseline first
                        for i in range(cluster_cfg.n_loops):
                            print('Calculating baseline...', end=' ')
                            baseline_result = evaluate_cluster_result(adata.copy(), master_recorder.cluster)
                            master_recorder.cluster.record(dataset, baseline_result, i)
                        # evaluate each method
                        for method in master_recorder.cluster.methods:
                            master_recorder.cluster.set_current_method(method)
                            all_result_list = select_genes(method, adata, config=cluster_cfg)
                            if all_result_list is not None:
                                for i in range(cluster_cfg.n_loops):
                                    head(name=method, fold=i + 1)
                                    for n_gene, result in zip(exp_cfg.n_genes, all_result_list):
                                        # clean directory
                                        delete('tempData/')
                                        gene_mask = master_recorder.cluster.get_mask(adata.var_names.values, result,
                                                                                     n_gene)
                                        selected_adata = adata[:, gene_mask]
                                        selected_adata.raw = adata.raw[:, gene_mask].to_adata()  # norm and raw
                                        # downstream analysis: clustering
                                        cluster_result = evaluate_cluster_result(selected_adata, master_recorder.cluster)
                                        master_recorder.cluster.record(dataset, cluster_result, i, n_gene)
                        master_recorder.cluster.save()
    except:
        master_recorder.save_all()
        traceback.print_exc()
