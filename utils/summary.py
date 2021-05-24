import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import seaborn as sns
import warnings

from config import data_cfg, assign_cfg, cluster_cfg, exp_cfg
from typing import Union, List
import pickle
import os


class PerformanceSummary:
    def __init__(self, record_path):
        self.record_path = record_path
        self.files = [file for file in os.listdir(record_path) if file.endswith('pkl')]
        # if data_cfg.n_datasets * 2 * 3 != len(self.files):  # 2 assign, cluster   3 n_genes
        #     raise RuntimeError("There are some missing records. Visualization stopped.")
        self.datasets = np.unique([name.split('-')[0] for name in self.files])

        self.assign_record_names = [name for name in self.files if name.split('-')[1] == 'assign']
        self.assign_records = self.load_record(self.assign_record_names)  # a dict
        self.assign_methods = self.assign_records[self.datasets[0]].scmap_cluster_F1.columns.tolist()

        self.cluster_record_names = [name for name in self.files if name.split('-')[1] == 'cluster']
        self.cluster_records = self.load_record(self.cluster_record_names)
        self.cluster_methods = self.cluster_records[self.datasets[0]].seurat_ARI.columns.tolist()

    def load_record(self, file_name: Union[str, List[str]]):
        if isinstance(file_name, str):
            with open(self.record_path + file_name, 'rb') as f:
                return pickle.load(f)
        else:
            if len(file_name) > 0:
                return {name.split('-')[0]: self.load_record(name) for name in file_name}
            else:
                warnings.warn("file_name doesn't have an element. Return null dict.")
                return dict()

    def save_specific_F1(self, method: str):
        """
        Save summary of specific F1 values to summary/

        :param method: can be 'scmap_cluster', 'scmap_cell' or 'singleR'.
        :return: None
        """
        specific_F1 = pd.DataFrame(index=self.datasets, columns=self.assign_methods)
        for dataset, record in self.assign_records.items():
            if method == 'scmap_cluster':
                f1_mean = record.scmap_cluster_F1.mean()
            elif method == 'scmap_cell':
                f1_mean = record.scmap_cell_F1.mean()
            elif method == 'singleR':
                f1_mean = record.singleR.mean()
            else:
                raise ValueError(f"{method} is an invalid argument.")
            f1_mean.name = dataset
            specific_F1.loc[dataset] = f1_mean
        specific_F1.name = method
        print(specific_F1)

    def save_specific_ARI(self, method: str):
        specific_ARI = pd.DataFrame(index=self.datasets, columns=self.cluster_methods)
        for dataset, record in self.cluster_records.items():
            if method == 'seurat':
                ari_mean = record.seurat_ARI.applymap(lambda x: np.mean(eval(x))).mean()
            elif method == 'sc3':
                ari_mean = record.sc3_ARI.mean()
            else:
                raise ValueError(f"{method} is an invalid argument.")
            ari_mean.name = dataset
            specific_ARI.loc[dataset] = ari_mean
        specific_ARI.name = method
        print(specific_ARI)


class PerformanceVisualizer:
    def __init__(self, task: str):
        assert task in ('assign', 'cluster'), ValueError(f"{task} is an invalid argument.")
        self.record_path = exp_cfg.record_path
        self.task = task
        self.record_names = [f for f in os.listdir(self.record_path) if f.endswith('pkl') and f.split('-')[1] == task]
        self.task_records = self.load_record(self.record_names)  # {dataset-n_gene: record}
        self.datasets = np.unique([file[:-4].split('-')[0] for file in self.record_names])
        self.n_genes = np.unique([file[:-4].split('-')[2] for file in self.record_names])

        if task == 'assign':
            self.task_methods = self.task_records[self.datasets[0] + '-' + self.n_genes[0]].singleR_F1.columns.tolist()
        else:
            self.task_methods = self.task_records[self.datasets[0] + '-' + self.n_genes[0]].seurat_ARI.columns.tolist()

    def load_record(self, file_name: Union[str, List[str]]):
        if isinstance(file_name, str):
            with open(self.record_path + file_name, 'rb') as f:
                return pickle.load(f)
        else:
            if len(file_name) > 0:
                return {name[:-4].split('-')[0] + '-' + name[:-4].split('-')[2]: self.load_record(name) for name in
                        file_name}
            else:
                warnings.warn("file_name doesn't have an element. Return null dict.")
                return dict()

    def plot_metric_heatmap(self, metric='MRR'):
        fig, axes = plt.subplots(2, 7, sharey='all')
        for i, method in enumerate(self.task_methods):
            spec_method_table = pd.DataFrame(np.zeros((len(self.datasets), len(self.n_genes))),
                                             index=self.datasets, columns=self.n_genes)
            ratio_factor = 1.0

            for dataset in self.datasets:
                for n_gene in self.n_genes:
                    record = self.task_records[dataset + '-' + n_gene]
                    if metric == 'markers_found':
                        ratio_factor = record.n_marker_contain
                    if hasattr(record, metric):
                        spec_method_table.loc[dataset, n_gene] = record.summary().loc[metric, method] / ratio_factor
                    else:
                        raise ValueError(f"record of {dataset + '-' + n_gene} doesn't have metric {metric}.")
            if i != len(self.task_methods) - 1:
                sns.heatmap(spec_method_table, annot=True, fmt='.3f' if metric != 'markers_found' else '.1%',
                            ax=axes[i // 7, i % 7], cmap="YlGnBu", cbar=False, vmin=0, vmax=1)
            else:
                sns.heatmap(spec_method_table, annot=True, fmt='.3f' if metric != 'markers_found' else '.1%',
                            ax=axes[i // 7, i % 7], cmap="YlGnBu", cbar=True, vmin=0, vmax=1,
                            cbar_kws={
                                'format': tick.FuncFormatter(lambda x, y: '{:.0%}'.format(x))
                            } if metric != 'markers_found' else None)
                axes[i // 7, i % 7 + 1].set_visible(False)
            axes[i // 7, i % 7].text(0.5, 4.5, method.replace('+', '+\n'), ha='center', va='center')
            axes[i // 7, i % 7].set_ylim([0, 5])
            axes[i // 7, i % 7].fill_between(np.arange(0.0, 1.0, 0.01), 4, 5, facecolor='grey', alpha=0.5)
        plt.show()

    def plot_task_bar(self):
        if self.task == 'assign':
            fig, axes = plt.subplots(3, 1)
            metrics = ['scmap_cluster_F1', 'scmap_cell_F1', 'singleR_F1']
        else:
            fig, axes = plt.subplots(2, 1)
            metrics = ['seurat_ARI', 'sc3_ARI']

        for i, metric in enumerate(metrics):
            spec_metric_table = pd.DataFrame(np.zeros((len(self.task_methods), len(self.n_genes))),
                                             index=self.task_methods, columns=self.n_genes)  # row method, # col n_gene
            for dataset in self.datasets:
                for n_gene in self.n_genes:
                    record = self.task_records[dataset + '-' + n_gene]
                    spec_metric_table.loc[:, n_gene] += record.summary().loc[metric, :]
            spec_metric_table.plot.bar(ax=axes[i])
        plt.show()
