import time
from typing import List

from common_utils.utils import plot_2D, head
from config.experiments_config import marker_cfg, assign_cfg, cluster_cfg, batch_cfg, time_cfg
from data_loader.dataset import yield_train_test_data, load_data
from experiments.metrics import marker_discovery_rate, correction_metrics, classification_metrics, clustering_metrics
from experiments.recorders import init_recorder
from other_steps.classification import classify_cells
from other_steps.clustering import cluster_cells
from other_steps.correction import correct_batch_effect
from selection.methods import select_genes
from selection.utils import subset_adata


def run_marker_discovery(fs_methods: List[str]):
    recorder = init_recorder(fs_methods, marker_cfg)
    for dataset_name in marker_cfg.datasets:
        adata = load_data(dataset_name)
        for fs_method in fs_methods:
            for n_genes in marker_cfg.n_genes:
                selected_adata = select_genes(adata, fs_method, n_genes)
                rate = marker_discovery_rate(selected_adata, adata)
                recorder.record(dataset_name, n_genes, fs_method, rate)
        recorder.sink()  # sink every dataset


def run_computation_time(fs_methods: List[str]):
    recorder = init_recorder(fs_methods, time_cfg)
    for dataset_name in time_cfg.datasets:
        adata = load_data(dataset_name)
        for fs_method in fs_methods:
            t0 = time.perf_counter()
            select_genes(adata, fs_method, time_cfg.n_genes[-1], use_saved=False)
            t1 = time.perf_counter()
            recorder.record(dataset_name, fs_method, t1 - t0)
        recorder.sink()  # sink every dataset


def run_batch_correction(fs_methods: List[str]):
    recorder = init_recorder(fs_methods, batch_cfg)
    for dataset_name in batch_cfg.datasets:
        adata = load_data(dataset_name)
        plot_2D(adata, mode='before')  # before the feature selection and batch correction
        for fs_method in fs_methods:
            for n_genes in batch_cfg.n_genes:
                for bc_method in batch_cfg.methods:
                    # feature selection + batch correction
                    selected_adata = select_genes(adata, fs_method, n_genes, select_by_batch=True)
                    corrected_adata = correct_batch_effect(selected_adata, bc_method)
                    results = correction_metrics(corrected_adata, batch_cfg)
                    recorder.record(dataset_name, 'selection_first', n_genes, bc_method, fs_method, results)
                    plot_2D(corrected_adata, fs_method, bc_method, 'after', 'selection_first')

                    # batch correction + feature selection
                    corrected_adata = correct_batch_effect(adata, bc_method)
                    selected_adata = select_genes(corrected_adata, fs_method, n_genes, select_by_batch=False)
                    results = correction_metrics(selected_adata, batch_cfg)
                    recorder.record(dataset_name, 'correction_first', n_genes, bc_method, fs_method, results)
                    plot_2D(corrected_adata, fs_method, bc_method, 'after', 'correction_first')
        recorder.sink()  # sink every dataset


def run_cell_classification(fs_methods: List[str]):
    recorder = init_recorder(fs_methods, assign_cfg)
    for dataset_name in assign_cfg.intra_datasets if assign_cfg.is_intra else assign_cfg.inter_datasets:
        adata = load_data(dataset_name)  # load the whole dataset
        recorder.extend_record_table(adata)  # add new part of table to the existing table, or create new table
        data_generator = yield_train_test_data(adata, assign_cfg)  # create a data generator
        for i, (train_adata, test_adata) in enumerate(data_generator):  # fold / batch:  split the data
            for assign_method in assign_cfg.methods:
                head(assign_method, i + 1)
                # baseline: all genes;
                classify_cells(train_adata, test_adata, assign_method)
                baseline_results = classification_metrics(test_adata, assign_cfg)
                recorder.record(test_adata.uns['data_name'], 'AllGenes', assign_method, baseline_results,
                                n_fold=i + 1 if assign_cfg.is_intra else None)
                for n_genes in assign_cfg.n_genes:
                    for fs_method in fs_methods:
                        # select genes on training set and do classification
                        selected_train_adata = select_genes(train_adata, fs_method, n_genes, select_by_batch=False)
                        selected_test_adata = subset_adata(test_adata, selected_train_adata.var_names)
                        classify_cells(selected_train_adata, selected_test_adata, assign_method)
                        results = classification_metrics(test_adata, assign_cfg)
                        recorder.record(test_adata.uns['data_name'], n_genes, assign_method, results, fs_method,
                                        n_fold=i + 1 if assign_cfg.is_intra else None)
        recorder.sink()


def run_cell_clustering(fs_methods: List[str]):
    recorder = init_recorder(fs_methods, cluster_cfg)
    for dataset_name in cluster_cfg.datasets:
        adata = load_data(dataset_name)
        for cluster_method in cluster_cfg.methods:
            # baseline: all genes;
            cluster_cells(adata, cluster_method, cluster_cfg.methods[cluster_method])
            baseline_results = clustering_metrics(adata, cluster_method, cluster_cfg)
            recorder.extend_record_table(dataset_name, 'AllGenes', cluster_method)
            recorder.record(dataset_name, 'AllGenes', cluster_method, baseline_results)
            for fs_method in fs_methods:
                for n_genes in cluster_cfg.n_genes:
                    recorder.extend_record_table(dataset_name, n_genes, cluster_method)
                    selected_adata = select_genes(adata, fs_method, n_genes)
                    cluster_cells(selected_adata, cluster_method, cluster_cfg.methods[cluster_method])
                    results = clustering_metrics(selected_adata, cluster_method, cluster_cfg)
                    recorder.record(dataset_name, n_genes, cluster_method, results, fs_method)
        recorder.sink()  # sink every dataset
