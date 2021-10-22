import numpy as np
import anndata as ad
import warnings
import datetime
import pickle
from typing import List, Union
from abc import abstractmethod
from itertools import combinations, product
from collections import defaultdict
from sklearn.metrics import silhouette_score
from rbo import rbo
from utils.utils import head
from config import assign_cfg, cluster_cfg, exp_cfg
import pandas as pd


class PerformanceRecorder:
    def __init__(self, datasets: List[str], methods: List[str]):
        self.datasets = datasets
        self.methods = methods
        self.current_method = None

    def set_current_method(self, method, fold=''):
        self.current_method = method
        if '+' in method:
            self.base_methods = method.split('+')
        head(name=method, fold=fold)

    def get_mask(self, all_genes: np.ndarray, single_selected_result: tuple, n_gene: Union[float, int]) -> np.ndarray:
        """
        Get the mask that indicates which genes are selected.

        :param all_genes: all gene names
        :param single_selected_result: the result of a selection
        :param n_gene: number of selected genes
        :return: a bool array
        """
        selected_genes, all_genes = np.char.strip(single_selected_result[0].astype(np.str), "'"), np.char.strip(all_genes.astype(np.str), "'")
        missing_genes = self._check_selection_result(selected_genes, all_genes, n_gene)
        while missing_genes is not None:
            if np.char.startswith(missing_genes, 'X').sum() == missing_genes.shape[0]:  # gene names start with X
                print("Removing 'X' at start...")
                prefixed = np.char.strip(np.char.lstrip(missing_genes, 'X'), '.')
                selected_genes = np.union1d(np.intersect1d(selected_genes, all_genes), prefixed)
            else:
                raise RuntimeError(f"There are genes not starting with X.")
                # break
            missing_genes = self._check_selection_result(selected_genes, all_genes, n_gene)
        mask = np.isin(all_genes, selected_genes)
        return mask

    @staticmethod
    def _check_selection_result(selected: np.ndarray, total: np.ndarray, need: Union[float, int]):
        msg, missing_genes = '', None
        if selected.shape[0] < need:
            expected_num = selected.shape[0]
            msg += f"The number of selected genes is {selected.shape[0]}, not {need}."
        else:
            expected_num = need
        mask = np.isin(total, selected)
        selected_markers_num = mask.sum()
        if selected_markers_num == 0:
            raise RuntimeError("No gene is selected!")
        elif 0 < selected_markers_num < expected_num:
            missing_genes = np.setdiff1d(selected, total)
            msg += f"{expected_num - selected_markers_num} selected genes are not used: {missing_genes}"
        else:  # selected_markers_num == expected_num
            print(f"Selected genes ({selected_markers_num}) are all successfully masked.")
        if msg != '':
            warnings.warn(msg, RuntimeWarning)
        return missing_genes

    @abstractmethod
    def save(self):
        pass


class MixtureRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(MixtureRecorder, self).__init__(datasets, methods)
        self.silhouette_coef = pd.DataFrame(
            data=np.zeros((len(exp_cfg.n_genes) * len(datasets), len(methods))),
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes]), columns=methods
        ).sort_index()

    def record(self, dataset: str, n_gene: int, filtered_adata: ad.AnnData):
        try:
            self.silhouette_coef.loc[(dataset, n_gene), self.current_method] = silhouette_score(
                X=filtered_adata.X, labels=filtered_adata.obs['celltype'].values
            )
        except ValueError as e:
            print(e)
            self.silhouette_coef.loc[(dataset, n_gene), self.current_method] = 0.

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        file_path = exp_cfg.record_path + 'MixtureRecorder.pkl'
        with open(file_path, 'wb') as f:
            pickle.dump(self, f)
        self.silhouette_coef.to_excel(exp_cfg.record_path + 'silhouette_coefficient.xlsx')
        print(f"{(datetime.datetime.now() + datetime.timedelta(hours=8)).strftime('%Y-%m-%d %H:%M:%S')}: "
              f"Finished and saved record to {file_path}\n\n\n")  # UTC + 8 hours = Beijing time


class MarkerRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(MarkerRecorder, self).__init__(datasets, methods)
        self.markers_found = pd.DataFrame(
            data=np.zeros((len(datasets) * len(exp_cfg.n_genes), len(methods))),
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes]), columns=methods
        ).sort_index()
        self.MRR = pd.DataFrame(
            data=np.zeros((len(datasets) * len(exp_cfg.n_genes), len(methods))),
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes]), columns=methods
        ).sort_index()

    @staticmethod
    def cal_markers_found_and_MRR(total_markers: np.ndarray,
                                  total_marker_weight: np.ndarray,
                                  single_selected_result: tuple):
        """
        Calculate proportion of selected marker genes and MRR if importances exist, otherwise only calculate the proportion.

        :param total_markers: all marker genes the dataset contains
        :param single_selected_result: result from function 'select_features'
        :return: proportion of marker genes found, and MRR
        """
        selected_genes = single_selected_result[0]
        rank = False if len(single_selected_result) == 1 else True
        markers_found_rate = total_marker_weight[np.isin(total_markers, selected_genes)].sum() / total_marker_weight.sum()
        if rank:
            rank_list = np.argwhere(np.isin(selected_genes, total_markers))
            if len(rank_list) == 0:
                warnings.warn("No marker gene is selected! MRR is set to 0.", RuntimeWarning)
                MRR = 0
            else:
                MRR = np.sum(1 / (rank_list + 1)) / rank_list.shape[0]
        else:  # len(selected_result) == 1
            warnings.warn("The method can't obtain gene importance! MRR is set to 0.", RuntimeWarning)
            MRR = 0
        return markers_found_rate, MRR

    def record(self, dataset: str,
               n_gene: int,
               all_markers: np.ndarray,
               all_marker_weight: np.ndarray,
               single_selected_result: tuple):
        markers_found_rate, MRR = self.cal_markers_found_and_MRR(all_markers, all_marker_weight, single_selected_result)
        self.markers_found.loc[(dataset, n_gene), self.current_method] = markers_found_rate
        self.MRR.loc[(dataset, n_gene), self.current_method] = MRR

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        file_path = exp_cfg.record_path + 'MarkerRecorder.pkl'
        with open(file_path, 'wb') as f:
            pickle.dump(self, f)
        self.markers_found.to_excel(exp_cfg.record_path + 'markers_found_rate.xlsx')
        self.MRR.to_excel(exp_cfg.record_path + 'MRR.xlsx')
        print(f"{(datetime.datetime.now() + datetime.timedelta(hours=8)).strftime('%Y-%m-%d %H:%M:%S')}: "
              f"Finished and saved record to {file_path}\n\n\n")  # UTC + 8 hours = Beijing time


class TimeRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(TimeRecorder, self).__init__(datasets, methods)
        self.computation_time = pd.DataFrame(
            data=np.zeros((len(datasets), len(methods))), index=datasets, columns=methods
        ).sort_index()

    def py_method_start(self):
        self.start_time = datetime.datetime.now()

    def py_method_end(self):
        self.end_time = datetime.datetime.now()

    @staticmethod
    def R_method_time(data_name, method):
        time_path = f'tempData/{data_name}_time_{method}.csv'
        seconds = np.round(np.loadtxt(time_path, skiprows=1, usecols=[1], delimiter=',').tolist(), decimals=5)
        return seconds

    def record(self, dataset, method):
        if assign_cfg.method_lan[method] == 'python':
            self.computation_time.loc[dataset, self.current_method] += (
                    self.end_time - self.start_time).total_seconds()
        else:
            self.computation_time.loc[dataset, self.current_method] += self.R_method_time(dataset, method)

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        file_path = exp_cfg.record_path + 'TimeRecorder.pkl'
        with open(file_path, 'wb') as f:
            pickle.dump(self, f)
        self.computation_time.to_excel(exp_cfg.record_path + 'computation_time.xlsx')
        print(f"{(datetime.datetime.now() + datetime.timedelta(hours=8)).strftime('%Y-%m-%d %H:%M:%S')}: "
              f"Finished and saved record to {file_path}\n\n\n")  # UTC + 8 hours = Beijing time


class BatchCorrectionRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(BatchCorrectionRecorder, self).__init__(datasets, methods)
        correct_methods = ['harmony', 'scanorama']  # 'scgen',
        correct_metrics = ['kBET', 'iLISI']
        self.correction = pd.DataFrame(
            np.zeros((len(datasets) * len(exp_cfg.n_genes) * len(correct_methods) * len(correct_metrics), len(methods))),
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes, correct_methods, correct_metrics]),
            columns=methods
        ).sort_index()

    def record(self, dataset: str, n_gene: int, correction_result: dict):
        for key, value in correction_result.items():
            correct_method, metric = key.split('_')
            self.correction.loc[(dataset, n_gene, correct_method, metric), self.current_method] = value

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        file_path = exp_cfg.record_path + 'BatchCorrectionRecorder.pkl'
        with open(file_path, 'wb') as f:
            pickle.dump(self, f)
        self.correction.to_excel(exp_cfg.record_path + 'batch_correction_metrics.xlsx')
        print(f"{(datetime.datetime.now() + datetime.timedelta(hours=8)).strftime('%Y-%m-%d %H:%M:%S')}: "
              f"Finished and saved record to {file_path}\n\n\n")  # UTC + 8 hours = Beijing time


class ClassificationRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(ClassificationRecorder, self).__init__(datasets, methods)
        assign_methods = ['singlecellnet', 'singleR']  # , 'scanvi'
        assign_metrics = ['f1_score', 'cohen_kappa_score']
        self.overall_metrics = pd.DataFrame(
            np.zeros((len(datasets) * len(exp_cfg.n_genes) * assign_cfg.n_folds * len(assign_methods) * len(
                assign_metrics), len(methods))),
            index=pd.MultiIndex.from_product(
                [datasets, exp_cfg.n_genes, range(assign_cfg.n_folds), assign_methods, assign_metrics]),
            columns=methods
        ).sort_index()
        self.rare_metric = pd.DataFrame(
            np.zeros((len(datasets) * len(exp_cfg.n_genes) * assign_cfg.n_folds * len(assign_methods), len(methods))),
            # rare F1-score
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes, range(assign_cfg.n_folds), assign_methods]),
            columns=methods
        ).sort_index()
        self.selected_genes = defaultdict(list)  # for stability
        self.stability = pd.DataFrame(
            np.zeros((len(datasets) * len(exp_cfg.n_genes), len(methods))),
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes]),
            columns=methods
        ).sort_index()

    def set_rare_type(self, adata: ad.AnnData):
        cell_type_counts = adata.obs['celltype'].value_counts(ascending=True)
        if cell_type_counts[0] >= assign_cfg.n_folds * 2:
            self.rare_type = cell_type_counts.index[0]
            print("Rare Cell Type:{:<15}  Rate:{:.2%}     Num:{}".format(
                self.rare_type, cell_type_counts[0] / adata.n_obs, cell_type_counts[0])
            )

    def store_selected_genes(self, dataset: str, result_list: list):
        for n_gene, result in zip(exp_cfg.n_genes, result_list):
            self.selected_genes[(dataset, self.current_method, n_gene)].append(result[0])

    def cal_mean_rbo(self, dataset: str):
        for method, n_gene in product(self.methods, exp_cfg.n_genes):
            self.stability.loc[(dataset, n_gene), method] = np.median(
                [rbo.RankingSimilarity(a, b).rbo(p=0.999) for a, b in combinations(self.selected_genes[(dataset, method, n_gene)], 2)]
            )  # when p = 0.999: the top 500, 1000, 1500, 2000 genes have 67%, 85%, 93%, 96%  of the weight
            # print(f"RBO (dataset: {dataset}, n_gene: {n_gene}): {self.stability.loc[(dataset, n_gene), method]}")

    def record(self, dataset: str, n_gene: int, n_fold: int, assign_result: dict):
        for key, value in assign_result.items():
            if key.endswith('rare'):
                assign_method, *_ = key.split('_')
                self.rare_metric.loc[(dataset, n_gene, n_fold, assign_method), self.current_method] = value
            else:
                assign_method, assign_metric = key.split('_')
                if assign_metric == 'f1':
                    self.overall_metrics.loc[(dataset, n_gene, n_fold, assign_method, 'f1_score'), self.current_method] = value
                else:  # cohen_kappa_score
                    self.overall_metrics.loc[(dataset, n_gene, n_fold, assign_method, 'cohen_kappa_score'), self.current_method] = value

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        file_path = exp_cfg.record_path + 'ClassificationRecorder.pkl'
        with open(file_path, 'wb') as f:
            pickle.dump(self, f)
        self.overall_metrics.to_excel(exp_cfg.record_path + 'classification_overall_metrics.xlsx')
        self.stability.to_excel(exp_cfg.record_path + 'selection_stability.xlsx')
        self.rare_metric.to_excel(exp_cfg.record_path + 'classification_rare_metrics.xlsx')
        print(f"{(datetime.datetime.now() + datetime.timedelta(hours=8)).strftime('%Y-%m-%d %H:%M:%S')}: "
              f"Finished and saved record to {file_path}\n\n\n")  # UTC + 8 hours = Beijing time


class ClusteringRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(ClusteringRecorder, self).__init__(datasets, methods)
        cluster_methods = ['seurat', 'sc3', 'cidr']
        cluster_metrics = ['ARI', 'BCubed', 'V']
        # get unsupervised methods
        self.methods = []
        for method in methods:
            if '+' in method:
                isUnsupervised = True
                base_methods = method.split('+')
                for base_method in base_methods:
                    if base_method not in cluster_cfg.method_lan.keys():
                        isUnsupervised = False
                if isUnsupervised:
                    self.methods.append(method)
            else:
                if method in cluster_cfg.method_lan.keys():
                    self.methods.append(method)

        self.overall_metrics = pd.DataFrame(
            np.zeros((len(datasets) * len(exp_cfg.n_genes) * cluster_cfg.n_loops * cluster_cfg.n_folds *
                      len(cluster_methods) * len(cluster_metrics), len(self.methods))),
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes, range(cluster_cfg.n_loops), range(cluster_cfg.n_folds), cluster_methods, cluster_metrics]),
            columns=self.methods
        ).sort_index()

        self.rare_metric = pd.DataFrame(
            np.zeros((len(datasets) * len(exp_cfg.n_genes) * cluster_cfg.n_loops * cluster_cfg.n_folds *
                      len(cluster_methods), len(self.methods))),
            # rare bcubed f-score
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes, range(cluster_cfg.n_loops), range(cluster_cfg.n_folds), cluster_methods]),
            columns=self.methods
        ).sort_index()

    def set_rare_type(self, adata: ad.AnnData):
        cell_type_counts = adata.obs['celltype'].value_counts(ascending=True)
        if cell_type_counts[0] >= assign_cfg.n_folds * 2:
            self.rare_type = cell_type_counts.index[0]
            print("Rare Cell Type:{:<15}  Rate:{:.2%}     Num:{}".format(
                self.rare_type, cell_type_counts[0] / adata.n_obs, cell_type_counts[0])
            )

    @staticmethod
    def BCubed_fbeta_score(labels_true: np.ndarray, labels_pred: np.ndarray, beta=1.):
        assert len(labels_true.shape) == len(labels_pred.shape) == 1, \
            ValueError("Dimensions are inconsistent or not equal to 1.")

        item_precision_list, item_recall_list = [], []
        for item_type, item_cluster in zip(labels_true, labels_pred):
            type_mask = labels_true == item_type
            cluster_mask = labels_pred == item_cluster

            cluster_item_type = labels_true[cluster_mask]
            item_precision_list.append((cluster_item_type == item_type).sum() / cluster_mask.sum())  # item precision
            item_recall_list.append((cluster_item_type == item_type).sum() / type_mask.sum())  # item recall

        bcubed_precision = np.array(item_precision_list).sum() / len(item_precision_list)
        bcubed_recall = np.array(item_recall_list).sum() / len(item_recall_list)
        return (1 + beta ** 2) * bcubed_precision * bcubed_recall / (beta ** 2 * bcubed_precision + bcubed_recall)

    def BCubed_fbeta_score_rare(self, labels_true: np.ndarray, labels_pred: np.ndarray, beta=1.):
        """
        BCubed_fbeta_score for rare cell type detection.

        :param labels_true: category annotation
        :param labels_pred: cluster annotation
        :param beta: the weighting hyperparameter
        :return: BCubed score for rare type
        """
        rare_type_mask = labels_true == self.rare_type
        rare_cell_num = np.sum(rare_type_mask)
        precision, recall = [], []
        for cluster in np.unique(labels_pred[rare_type_mask]):  # clusters that rare cells are divided to
            cluster_mask = labels_pred == cluster
            cluster_true = labels_true[cluster_mask]

            n_true_positive = np.sum(cluster_true == self.rare_type)
            precision.append(n_true_positive ** 2 / np.sum(cluster_mask))  # sum of precision in this cluster
            recall.append(n_true_positive ** 2 / rare_cell_num)  # sum of recall in this cluster

        ave_precision, ave_recall = np.sum(precision) / rare_cell_num, np.sum(recall) / rare_cell_num
        fbeta_BCubed_rare = (beta ** 2 + 1) * ave_precision * ave_recall / (beta ** 2 * ave_precision + ave_recall)
        return fbeta_BCubed_rare

    def record(self, dataset: str, n_gene: int, n_loop: int, n_fold: int, cluster_result: dict):
        for key, value in cluster_result.items():
            if key.endswith('rare'):
                cluster_method, *_ = key.split('_')
                self.rare_metric.loc[(dataset, n_gene, n_loop, n_fold, cluster_method), self.current_method] = value
            else:
                cluster_method, cluster_metric = key.split('_')
                self.overall_metrics.loc[(dataset, n_gene, n_loop, n_fold, cluster_method, cluster_metric), self.current_method] = value

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        file_path = exp_cfg.record_path + 'ClusteringRecorder.pkl'
        with open(file_path, 'wb') as f:
            pickle.dump(self, f)
        self.overall_metrics.to_excel(exp_cfg.record_path + 'clustering_overall_metrics.xlsx')
        self.rare_metric.to_excel(exp_cfg.record_path + 'clustering_rare_metrics.xlsx')
        print(f"{(datetime.datetime.now() + datetime.timedelta(hours=8)).strftime('%Y-%m-%d %H:%M:%S')}: "
              f"Finished and saved record to {file_path}\n\n\n")  # UTC + 8 hours = Beijing time


class MasterRecorder:
    """
    Contains all base recorders.
    """
    def __init__(self, measurements: List[str], methods: List[str]):
        print(f"{(datetime.datetime.now() + datetime.timedelta(hours=8)).strftime('%Y-%m-%d %H:%M:%S')}: "
              f"evaluation starts.")  # UTC + 8 hours = Beijing time
        self.measurements = []
        for measurement, datasets in exp_cfg.measurements.items():
            if measurement in measurements:
                self.measurements.append(measurement)
                if measurement == 'population_demixing':
                    self.mixture = MixtureRecorder(datasets=datasets, methods=methods)
                elif measurement == 'computation_time':
                    self.time = TimeRecorder(datasets=datasets, methods=methods)
                elif measurement == 'classification':
                    self.assign = ClassificationRecorder(datasets=datasets, methods=methods)
                elif measurement == 'marker_discovery':
                    self.marker = MarkerRecorder(datasets=datasets, methods=methods)
                elif measurement == 'batch_correction':
                    self.correct = BatchCorrectionRecorder(datasets=datasets, methods=methods)
                elif measurement == 'clustering':
                    self.cluster = ClusteringRecorder(datasets=datasets, methods=methods)
                else:
                    raise ValueError(f"{measurement} is not appropriate.")

    def save_all(self):
        for measurement in self.measurements:
            if measurement == 'population_demixing':
                self.mixture.save()
            elif measurement == 'computation_time':
                self.time.save()
            elif measurement == 'classification':
                self.assign.save()
            elif measurement == 'clustering':
                self.cluster.save()
            elif measurement == 'marker_discovery':
                self.marker.save()
            elif measurement == 'batch_correction':
                self.correct.save()
            else:
                raise ValueError(f"{measurement} is not appropriate.")
