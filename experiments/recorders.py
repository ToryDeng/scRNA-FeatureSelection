from typing import List, Union
from abc import abstractmethod
from config.experiments_config import MarkerDiscoveryConfig, CellClassificationConfig, CellClusteringConfig, \
    BasicExperimentConfig, ComputationTimeConfig, BatchCorrectionConfig
from config.methods_config import method_cfg
from itertools import permutations
from collections.abc import Iterable
import anndata as ad
import datetime
import pickle
import os
import pandas as pd
import numpy as np


class BasicExperimentRecorder:
    def __init__(self, fs_methods: List[str], config: BasicExperimentConfig):
        self.config = config
        self.__check_fs_methods_exist(fs_methods)
        self.__create_sink_dir(config)

    @staticmethod
    def __check_fs_methods_exist(fs_methods: List[str]):
        for fs_method in fs_methods:
            if '+' not in fs_method and fs_method not in method_cfg.formal_names.keys():
                raise ValueError(f"Wrong name of feature selection method: {fs_method}")
            if '+' in fs_method:
                base_methods = fs_method.split('+')
                if np.isin(base_methods, method_cfg.formal_names.keys()).sum() < len(base_methods):
                    raise ValueError(f"Wrong name of ensemble method {fs_method}")

    @staticmethod
    def __create_sink_dir(config: BasicExperimentConfig):
        pkl_dir, xlsx_dir = os.path.join(config.sink_dir, 'pkl/'), os.path.join(config.sink_dir, 'xlsx/')
        if not os.path.exists(pkl_dir):
            os.makedirs(pkl_dir)
        if not os.path.exists(xlsx_dir):
            os.makedirs(xlsx_dir)

    @staticmethod
    def _create_record_table(*row_grps, fs_methods: List[str]):
        n_rows, row_idxs = 1, []
        for row_grp in row_grps:
            if isinstance(row_grp, str) or isinstance(row_grp, int):
                row_idxs.append([row_grp])
            else:  # non-string and non-integer
                if isinstance(row_grp, Iterable):
                    n_rows *= sum(1 for _ in row_grp)
                    row_idxs.append(row_grp)
                else:
                    raise ValueError(f"type {type(row_grp)} can not be used to create dataframe!")
        return pd.DataFrame(
            data=np.full(shape=(n_rows, len(fs_methods)), fill_value=np.nan),
            index=pd.MultiIndex.from_product(row_idxs),
            columns=fs_methods
        ).sort_index()

    @abstractmethod
    def record(self, *args, **kwargs):
        pass

    def sink(self):
        # sink pickle file
        with open(os.path.join(self.config.sink_dir, 'pkl', f'{self.__class__.__name__}.pkl'), 'wb') as f:
            pickle.dump(self, f)
        # sink tables
        for name, attr in self.__dict__.items():
            if isinstance(attr, pd.DataFrame):
                attr.to_excel(os.path.join(os.path.join(self.config.sink_dir, 'xlsx', f'{name}.xlsx')), encoding='utf-8')
        print(f"{datetime.datetime.now()}: Finished and saved record to {self.config.sink_dir}\n\n\n")


class MarkerRecorder(BasicExperimentRecorder):
    def __init__(self, fs_methods: List[str], config: MarkerDiscoveryConfig):
        super(MarkerRecorder, self).__init__(fs_methods, config)
        self.marker_discovery_rate = self._create_record_table(config.datasets, config.n_genes, fs_methods=fs_methods)

    def record(self, dataset: str, n_genes: Union[int, str], fs_method: str, rate: float):
        self.marker_discovery_rate.loc[(dataset, n_genes), fs_method] = rate


class TimeRecorder(BasicExperimentRecorder):
    def __init__(self, fs_methods: List[str], config: ComputationTimeConfig):
        super(TimeRecorder, self).__init__(fs_methods, config)
        self.computation_time = self._create_record_table(config.datasets, fs_methods=fs_methods)

    def record(self, dataset: str, fs_method: str, seconds: float):
        self.computation_time.loc[dataset, fs_method] = seconds


class BatchRecorder(BasicExperimentRecorder):
    def __init__(self, fs_methods: List[str], config: BatchCorrectionConfig):
        super(BatchRecorder, self).__init__(fs_methods, config)
        self.config = config
        self.correction = self._create_record_table(
            config.datasets, config.orders, config.n_genes + ['AllGenes'], config.methods, config.metrics,  # rows
            fs_methods=fs_methods
        )

    def record(self, dataset: str, order: str, n_genes: Union[int, str], method: str, fs_method: str, metrics: dict):
        for metric, value in metrics.items():
            if metric in self.config.metrics:
                self.correction.loc[(dataset, order, n_genes, method, metric), fs_method] = value


class AssignRecorder(BasicExperimentRecorder):
    def __init__(self, fs_methods: List[str], config: CellClassificationConfig):
        super(AssignRecorder, self).__init__(fs_methods, config)
        self.config = config
        self.fs_methods = fs_methods
        if self.config.is_intra:
            self.classification = self._create_record_table(
                config.intra_datasets, range(1, config.n_folds + 1), config.n_genes + ['AllGenes'], config.methods, config.metrics,
                fs_methods=fs_methods
            )

    def extend_record_table(self, adata: ad.AnnData):
        if not self.config.is_intra:
            assert 'batch' in adata.obs, "No batch in anndata!"
            datasets = [' to '.join([train_name, test_name]) for train_name, test_name in
                        permutations(adata.obs['batch'].unique(), 2)]
            part = self._create_record_table(
                datasets, self.config.n_genes + ['AllGenes'], self.config.methods, self.config.metrics,
                fs_methods=self.fs_methods
            )
            self.classification = self.classification.append(part).sort_index() if hasattr(self, 'classification') else part
        else:
            pass

    def record(self, dataset: str, n_genes: Union[int, str], method: str, metrics: dict, fs_method: str = None, n_fold: int = None):
        if n_fold is not None:  # intra-dataset classification
            for metric, value in metrics.items():
                if metric in self.config.metrics:
                    if isinstance(n_genes, str):  # all genes
                        self.classification.loc[(dataset, n_fold, n_genes, method, metric), :] = value
                    else:
                        self.classification.loc[(dataset, n_fold, n_genes, method, metric), fs_method] = value
        else:  # inter-dataset classification
            for metric, value in metrics.items():
                if metric in self.config.metrics:
                    if isinstance(n_genes, str):  # all genes
                        self.classification.loc[(dataset, n_genes, method, metric), :] = value
                    else:
                        self.classification.loc[(dataset, n_genes, method, metric), fs_method] = value


class ClusterRecorder(BasicExperimentRecorder):
    def __init__(self, fs_methods: List[str], config: CellClusteringConfig):
        super(ClusterRecorder, self).__init__(fs_methods, config)
        self.config = config
        self.fs_methods = fs_methods
        self.__check_fs_methods_all_unsupervised()

    def extend_record_table(self, dataset: str, n_genes: Union[int, str], clustering_method: str):
        part = self._create_record_table(
            dataset, n_genes, clustering_method, range(1, self.config.methods[clustering_method]), self.config.metrics,
            fs_methods=self.fs_methods
        )
        self.clustering = self.clustering.append(part).sort_index() if hasattr(self, 'clustering') else part

    def __check_fs_methods_all_unsupervised(self):
        for fs_method in self.fs_methods:
            if '+' not in fs_method and fs_method not in method_cfg.unsupervised:
                raise ValueError(f"{fs_method} not in the list of unsupervised methods in method config")
            if '+' in fs_method:
                base_methods = fs_method.split('+')
                if np.isin(base_methods, method_cfg.unsupervised.keys()).sum() < len(base_methods):
                    raise ValueError(f"One or more base methods in {fs_method}"
                                     f" are not in the list of unsupervised methods in method config")

    def record(self, dataset: str, n_genes: Union[int, str], method: str, metrics: dict, fs_method: str = None):
        for metric_run, value in metrics.items():
            metric, run = metric_run.split('_')
            if metric in self.config.metrics:
                if isinstance(n_genes, str):  # all genes
                    self.clustering.loc[(dataset, n_genes, method, int(run), metric), :] = value
                else:
                    self.clustering.loc[(dataset, n_genes, method, int(run), metric), fs_method] = value


def init_recorder(fs_methods: List[str], config: BasicExperimentConfig):
    if isinstance(config, MarkerDiscoveryConfig):
        recorder = MarkerRecorder(fs_methods, config)
    elif isinstance(config, CellClassificationConfig):
        recorder = AssignRecorder(fs_methods, config)
    elif isinstance(config, CellClusteringConfig):
        recorder = ClusterRecorder(fs_methods, config)
    elif isinstance(config, BatchCorrectionConfig):
        recorder = BatchRecorder(fs_methods, config)
    elif isinstance(config, ComputationTimeConfig):
        recorder = TimeRecorder(fs_methods, config)
    else:
        raise ValueError("Wrong config type!")
    return recorder