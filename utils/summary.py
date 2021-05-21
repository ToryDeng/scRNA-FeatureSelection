import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings

from config import data_cfg, assign_cfg, cluster_cfg
from typing import Union, List
import pickle
import os


class PerformanceSummary:
    def __init__(self, record_path):
        self.record_path = record_path
        self.files = [file for file in os.listdir(record_path) if file.endswith('pkl')]
        # if data_cfg.n_datasets * 2 * 3 != len(self.files):
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



