import numpy as np
import anndata as ad
import scanpy as sc
import warnings
import datetime
import pickle
import os
from typing import List, Union, Literal
from functools import reduce
from abc import abstractmethod
from sklearn.metrics import silhouette_score
from itertools import permutations
from utils.utils import head, now
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
        selected_genes, all_genes = np.char.strip(single_selected_result[0].astype(np.str), "'"), np.char.strip(
            all_genes.astype(np.str), "'")
        missing_genes = self._check_selection_result(selected_genes, all_genes, n_gene)
        while missing_genes is not None:
            if np.char.startswith(missing_genes, 'X').sum() == missing_genes.shape[0]:  # gene names start with X
                print("Removing 'X' at start...")
                prefixed = np.char.strip(np.char.lstrip(missing_genes, 'X'), '.')
                selected_genes = np.union1d(np.intersect1d(selected_genes, all_genes), prefixed)
            else:
                raise RuntimeError(f"There are genes not starting with 'X'.")
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
            print(f"Selected genes ({selected_markers_num}) are all successfully masked.", end=' ')
        if msg != '':
            warnings.warn(msg, RuntimeWarning)
        return missing_genes

    @staticmethod
    def _generate_table(*args):
        n_rows, idx, col = 1, [], args[-1]
        for arg in args[:-1]:
            if isinstance(arg, int):
                n_rows *= arg
                idx.append(range(arg))
            else:
                n_rows *= len(arg)
                idx.append(arg)
        return pd.DataFrame(np.full((n_rows, len(col)), np.nan), index=pd.MultiIndex.from_product(idx), columns=col)

    @abstractmethod
    def save(self):
        pass


class MixtureRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(MixtureRecorder, self).__init__(datasets, methods)
        self.silhouette_coef = pd.DataFrame(
            data=np.full(shape=(len(exp_cfg.n_genes) * len(datasets), len(methods)), fill_value=np.nan),
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes]), columns=methods
        ).sort_index()
        self.baseline = pd.DataFrame(
            data=np.full(shape=(len(datasets), 1), fill_value=np.nan),
            index=datasets, columns=['AllGenes']
        ).sort_index()

    def record(self, dataset: str, filtered_adata: ad.AnnData, n_gene: int = None):
        try:
            if n_gene is None:
                self.baseline.loc[dataset, 'AllGenes'] = silhouette_score(
                    X=sc.pp.pca(filtered_adata.X), labels=filtered_adata.obs['celltype'].values
                )
            else:
                self.silhouette_coef.loc[(dataset, n_gene), self.current_method] = silhouette_score(
                    X=sc.pp.pca(filtered_adata.X), labels=filtered_adata.obs['celltype'].values
                )
        except ValueError as e:
            print(e)

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        with open(os.path.join(exp_cfg.record_path, 'pkl', 'MarkerRecorder.pkl'), 'wb') as f:
            pickle.dump(self, f)

        self.silhouette_coef.to_excel(os.path.join(exp_cfg.record_path, 'xlsx', 'silhouette_coefficient.xlsx'))
        # baseline results
        self.baseline.to_excel(os.path.join(exp_cfg.record_path, 'xlsx', 'silhouette_coefficient_baseline.xlsx'))
        print(f"{now()}: Finished and saved record to {exp_cfg.record_path}\n\n\n")  # UTC + 8 hours = Beijing time


class MarkerRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(MarkerRecorder, self).__init__(datasets, methods)
        self.markers_found = pd.DataFrame(
            data=np.full(shape=(len(datasets) * len(exp_cfg.n_genes), len(methods)), fill_value=np.nan),
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes]), columns=methods
        ).sort_index()
        self.MRR = pd.DataFrame(
            data=np.full(shape=(len(datasets) * len(exp_cfg.n_genes), len(methods)), fill_value=np.nan),
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes]), columns=methods
        ).sort_index()

        self.markers_found_baseline = pd.DataFrame(
            data=np.full(shape=(len(datasets), 1), fill_value=np.nan),
            index=datasets, columns=['AllGenes']
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
        markers_found_rate = total_marker_weight[
                                 np.isin(total_markers, selected_genes)].sum() / total_marker_weight.sum()
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
               all_markers: np.ndarray,
               all_marker_weight: np.ndarray,
               single_selected_result: tuple,
               n_gene: int = None
               ):
        if n_gene is None:
            markers_found_rate = self.cal_markers_found_and_MRR(all_markers, all_marker_weight, single_selected_result)
            self.markers_found_baseline.loc[dataset, 'AllGenes'] = markers_found_rate[0]
        else:
            markers_found_rate, MRR = self.cal_markers_found_and_MRR(all_markers, all_marker_weight,
                                                                     single_selected_result)
            self.markers_found.loc[(dataset, n_gene), self.current_method] = markers_found_rate
            self.MRR.loc[(dataset, n_gene), self.current_method] = MRR

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        with open(os.path.join(exp_cfg.record_path, 'pkl', 'MarkerRecorder.pkl'), 'wb') as f:
            pickle.dump(self, f)

        self.markers_found.to_excel(os.path.join(exp_cfg.record_path, 'xlsx', 'markers_found_rate.xlsx'))
        self.MRR.to_excel(os.path.join(exp_cfg.record_path, 'xlsx', 'MRR.xlsx'))
        # baseline results
        self.markers_found_baseline.to_excel(
            os.path.join(exp_cfg.record_path, 'xlsx', 'markers_found_rate_baseline.xlsx'))

        print(f"{now()}: Finished and saved record to {exp_cfg.record_path}\n\n\n")  # UTC + 8 hours = Beijing time


class TimeRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(TimeRecorder, self).__init__(datasets, methods)
        self.computation_time = pd.DataFrame(
            data=np.zeros((len(datasets), len(methods))), index=datasets, columns=methods
        )

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
        if os.path.exists(os.path.join(exp_cfg.record_path, 'pkl', 'TimeRecorder.pkl')):
            with open(os.path.join(exp_cfg.record_path, 'pkl', 'TimeRecorder.pkl'), 'rb') as f:
                self.computation_time = pickle.load(f).computation_time

        if assign_cfg.method_lan[method] == 'python':
            self.computation_time.loc[dataset, self.current_method] = (
                    self.end_time - self.start_time).total_seconds()
            self.start_time, self.end_time = None, None
        else:
            self.computation_time.loc[dataset, self.current_method] = self.R_method_time(dataset, method)
        print(f"computation time: {self.computation_time.loc[dataset, self.current_method]} seconds.")
        self.save()

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        with open(os.path.join(exp_cfg.record_path, 'pkl', 'TimeRecorder.pkl'), 'wb') as f:
            pickle.dump(self, f)

        self.computation_time.to_excel(os.path.join(exp_cfg.record_path, 'xlsx', 'computation_time.xlsx'))

        print(f"{now()}: Finished and saved record to {exp_cfg.record_path}\n\n\n")  # UTC + 8 hours = Beijing time


class BatchCorrectionRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(BatchCorrectionRecorder, self).__init__(datasets, methods)
        correct_methods = ['seuratV4']
        correct_metrics = ['kBET', 'iLISI', 'cLISI', 'f1LISI']
        correct_orders = ['correction_first', 'selection_first']
        self.correction = pd.DataFrame(
            np.full(
                (len(datasets) * len(exp_cfg.n_genes) * len(correct_methods) * len(correct_metrics) * len(
                    correct_orders), len(methods)),
                fill_value=np.nan),
            index=pd.MultiIndex.from_product(
                [datasets, exp_cfg.n_genes, correct_methods, correct_metrics, correct_orders]),
            columns=methods
        ).sort_index()

    def record(self, dataset: str, n_gene: int, result: dict, result_type: str):
        for key, value in result.items():
            correct_method, metric = key.split('_')
            self.correction.loc[(dataset, n_gene, correct_method, metric, result_type), self.current_method] = value

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        with open(os.path.join(exp_cfg.record_path, 'pkl', 'BatchCorrectionRecorder.pkl'), 'wb') as f:
            pickle.dump(self, f)

        self.correction.to_excel(os.path.join(exp_cfg.record_path, 'xlsx', 'batch_correction_metrics.xlsx'))

        print(f"{now()}: Finished and saved record to {exp_cfg.record_path}\n\n\n")  # UTC + 8 hours = Beijing time


class ClassificationRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str], assign_type: Literal['intra', 'inter']):
        super(ClassificationRecorder, self).__init__(datasets, methods)
        self.assign_methods = assign_cfg.evaluation_method
        self.assign_metrics = ['f1_score', 'cohen_kappa_score']
        self.methods = methods
        self.assign_type = assign_type
        if self.assign_type == 'intra':
            self.classification_type = ''
            self.overall_metrics = self._generate_table(datasets, exp_cfg.n_genes, assign_cfg.n_folds, self.assign_methods, self.assign_metrics, methods).sort_index()
            self.rare_metric = self._generate_table(datasets, exp_cfg.n_genes, assign_cfg.n_folds, self.assign_methods, methods).sort_index()
            self.overall_metrics_baseline = self._generate_table(datasets, assign_cfg.n_folds, self.assign_methods, self.assign_metrics, ['AllGenes']).sort_index()
            self.rare_metric_baseline = self._generate_table(datasets, assign_cfg.n_folds, self.assign_methods, ['AllGenes']).sort_index()
        else:  # inter
            self.overall_metrics, self.rare_metric = None, None
            self.overall_metrics_baseline, self.rare_metric_baseline = None, None

        self.rare_type = None

    def set_record_tables(self, adata: ad.AnnData):
        if self.assign_type == 'inter':
            ubatches = adata.obs['batch'].unique().categories.to_numpy()
            self.perms = [' to '.join(perm) for perm in permutations(ubatches, 2)]
            dataset = [adata.uns['data_name']]
            if self.overall_metrics is None:  # inter-dataset
                self.overall_metrics = self._generate_table(dataset, exp_cfg.n_genes, self.perms, self.assign_methods, self.assign_metrics, self.methods)
                self.rare_metric = self._generate_table(dataset, exp_cfg.n_genes, self.perms, self.assign_methods, self.methods)
                self.overall_metrics_baseline = self._generate_table(dataset, self.perms, self.assign_methods, self.assign_metrics, ['AllGenes'])
                self.rare_metric_baseline = self._generate_table(dataset, self.perms, self.assign_methods, ['AllGenes'])
            else:
                overall_metrics = self._generate_table(dataset, exp_cfg.n_genes, self.perms, self.assign_methods, self.assign_metrics, self.methods)
                rare_metric = self._generate_table(dataset, exp_cfg.n_genes, self.perms, self.assign_methods, self.methods)
                overall_metrics_baseline = self._generate_table(dataset, self.perms, self.assign_methods, self.assign_metrics, ['AllGenes'])
                rare_metric_baseline = self._generate_table(dataset, self.perms, self.assign_methods, ['AllGenes'])

                self.overall_metrics = pd.concat([self.overall_metrics, overall_metrics])
                self.rare_metric = pd.concat([self.rare_metric, rare_metric])
                self.overall_metrics_baseline = pd.concat([self.overall_metrics_baseline, overall_metrics_baseline])
                self.rare_metric_baseline = pd.concat([self.rare_metric_baseline, rare_metric_baseline])

    def set_rare_type(self, adata: ad.AnnData):
        cell_type_counts = adata.obs['celltype'].value_counts(ascending=True)
        if 'batch' not in adata.obs:
            if cell_type_counts[0] >= assign_cfg.n_folds * 2:  # 10
                self.rare_type = cell_type_counts.index[0]
                print("Rare Cell Type:{:<15}     Rate:{:.2%}     Num:{}".format(
                    self.rare_type, cell_type_counts[self.rare_type] / adata.n_obs, cell_type_counts[0])
                )
            else:
                self.rare_type = None
        else:
            u_batches = adata.obs['batch'].unique()
            bu_types = [adata[adata.obs['batch'] == b].obs['celltype'].unique() for b in u_batches]
            inter_types = np.intersect1d(*bu_types) if len(bu_types) == 2 else reduce(np.intersect1d, bu_types)
            try:
                self.rare_type = cell_type_counts.index[np.isin(cell_type_counts.index.to_numpy(), inter_types)][0]
                print("Rare Cell Type:{:<15}     Rate in the total dataset: {:.2%}     Num: {}".format(
                    self.rare_type, cell_type_counts[self.rare_type] / adata.n_obs, cell_type_counts[self.rare_type]
                ))
                for b in u_batches:
                    batch_mask = adata.obs['batch'] == b
                    batch_counts = adata[batch_mask].obs['celltype'].value_counts()
                    print("Rate in {}: {:.2%}     Num in {}: {}".format(
                        b, batch_counts[self.rare_type] / batch_mask.sum(), b, batch_counts[self.rare_type]
                    ))
            except IndexError:
                self.rare_type = None

    def set_rare_type_and_tables(self, adata: ad.AnnData):
        self.set_rare_type(adata)
        self.set_record_tables(adata)

    def record(self, dataset: str, n_split: Union[int, str], assign_result: dict, n_gene: int = None):
        for key, value in assign_result.items():
            if key.endswith('rare'):
                assign_method, *_ = key.split('_')
                if n_gene is None:  # baseline
                    self.rare_metric_baseline.loc[(dataset, n_split, assign_method), 'AllGenes'] = value
                else:
                    self.rare_metric.loc[(dataset, n_gene, n_split, assign_method), self.current_method] = value
            else:
                assign_method, assign_metric = key.split('_')
                if assign_metric == 'f1':
                    if n_gene is None:
                        self.overall_metrics_baseline.loc[(dataset, n_split, assign_method, 'f1_score')] = value
                    else:
                        self.overall_metrics.loc[
                            (dataset, n_gene, n_split, assign_method, 'f1_score'), self.current_method] = value
                else:  # cohen_kappa_score
                    if n_gene is None:
                        self.overall_metrics_baseline.loc[(dataset, n_split, assign_method, 'cohen_kappa_score')] = value
                    else:
                        self.overall_metrics.loc[
                            (dataset, n_gene, n_split, assign_method, 'cohen_kappa_score'), self.current_method] = value

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """
        with open(os.path.join(exp_cfg.record_path, 'pkl', 'ClassificationRecorder.pkl'), 'wb') as f:
            pickle.dump(self, f)

        self.overall_metrics.to_excel(
            os.path.join(exp_cfg.record_path, 'xlsx', f'{self.assign_type}-classification_overall_metrics.xlsx'))
        self.rare_metric.to_excel(
            os.path.join(exp_cfg.record_path, 'xlsx', f'{self.assign_type}-classification_rare_metrics.xlsx'))
        # baseline results
        self.overall_metrics_baseline.to_excel(
            os.path.join(exp_cfg.record_path, 'xlsx', f'{self.assign_type}-classification_overall_metrics_baseline.xlsx'))
        self.rare_metric_baseline.to_excel(
            os.path.join(exp_cfg.record_path, 'xlsx', f'{self.assign_type}-classification_rare_metrics_baseline.xlsx'))

        print(f"{now()}: Finished and saved record to {exp_cfg.record_path}\n\n\n")  # UTC + 8 hours = Beijing time


class ClusteringRecorder(PerformanceRecorder):
    def __init__(self, datasets: List[str], methods: List[str]):
        super(ClusteringRecorder, self).__init__(datasets, methods)
        cluster_methods = cluster_cfg.evaluation_method
        cluster_metrics = ['ARI', 'V']
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
                    print(f"method {method} is not an unsupervised method!")
            else:
                if method in cluster_cfg.method_lan.keys():
                    self.methods.append(method)
                else:
                    print(f"method {method} is not an unsupervised method!")

        self.overall_metrics = pd.DataFrame(
            np.full(shape=(len(datasets) * len(exp_cfg.n_genes) * cluster_cfg.n_loops *
                           len(cluster_methods) * len(cluster_metrics), len(self.methods)),
                    fill_value=np.nan),
            index=pd.MultiIndex.from_product(
                [datasets, exp_cfg.n_genes, range(cluster_cfg.n_loops), cluster_methods, cluster_metrics]),
            columns=self.methods
        ).sort_index()

        self.rare_metric = pd.DataFrame(
            np.full(shape=(
                len(datasets) * len(exp_cfg.n_genes) * cluster_cfg.n_loops * len(cluster_methods), len(self.methods)),
                fill_value=np.nan),  # rare bcubed f-score
            index=pd.MultiIndex.from_product([datasets, exp_cfg.n_genes, range(cluster_cfg.n_loops), cluster_methods]),
            columns=self.methods
        ).sort_index()

        self.overall_metrics_baseline = pd.DataFrame(
            np.full(shape=(len(datasets) * cluster_cfg.n_loops * len(cluster_methods) * len(cluster_metrics), 1),
                    fill_value=np.nan),
            index=pd.MultiIndex.from_product([datasets, range(cluster_cfg.n_loops), cluster_methods, cluster_metrics]),
            columns=['AllGenes']
        ).sort_index()

        self.rare_metric_baseline = pd.DataFrame(
            np.full(shape=(len(datasets) * cluster_cfg.n_loops * len(cluster_methods), 1), fill_value=np.nan),
            index=pd.MultiIndex.from_product([datasets, range(cluster_cfg.n_loops), cluster_methods]),
            columns=['AllGenes']
        ).sort_index()

        self.rare_type = None

    def set_rare_type(self, adata: ad.AnnData):
        cell_type_counts = adata.obs['celltype'].value_counts(ascending=True)
        self.rare_type = cell_type_counts.index[0]
        print("Rare Cell Type:{:<15}  Rate:{:.2%}     Num:{}".format(
            self.rare_type, cell_type_counts[0] / adata.n_obs, cell_type_counts[0])
        )

    @staticmethod
    def BCubed_fbeta_score(labels_true: np.ndarray, labels_pred: np.ndarray, beta=1.):
        correctness = lambda x, y: 1 if labels_true[x] == labels_true[y] and labels_pred[x] == labels_pred[y] else 0

        BCubed_precision = np.mean([np.mean(
            [correctness(i, j) for j in range(labels_pred.shape[0]) if labels_pred[i] == labels_pred[j]]
        ) for i in range(labels_pred.shape[0])])

        BCubed_recall = np.mean([np.mean(
            [correctness(i, j) for j in range(labels_true.shape[0]) if labels_true[i] == labels_true[j]]
        ) for i in range(labels_true.shape[0])])
        return (1 + beta ** 2) * BCubed_precision * BCubed_recall / (beta ** 2 * BCubed_precision + BCubed_recall)

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
        if rare_cell_num == 0:
            print(
                f"There is no rare cell ({self.rare_type}) in the data! All cell types in data: f{np.unique(labels_true)}")
            fbeta_BCubed_rare = np.nan
        else:
            precision, recall = [], []
            for cluster in np.unique(labels_pred[rare_type_mask]):  # clusters that rare cells are divided to
                cluster_mask = labels_pred == cluster
                cluster_true = labels_true[cluster_mask]

                n_true_positive = np.sum(cluster_true == self.rare_type)
                precision.append(n_true_positive ** 2 / np.sum(cluster_mask))  # sum of precision in this cluster
                recall.append(n_true_positive ** 2 / rare_cell_num)  # sum of recall in this cluster

            ave_precision, ave_recall = np.sum(precision) / rare_cell_num, np.sum(recall) / rare_cell_num
            if ave_precision == ave_precision == 0:
                print("precision and recall are all equal to zero!")
                fbeta_BCubed_rare = np.nan
            else:
                fbeta_BCubed_rare = (beta ** 2 + 1) * ave_precision * ave_recall / (
                        beta ** 2 * ave_precision + ave_recall)
        return fbeta_BCubed_rare

    def record(self, dataset: str, cluster_result: dict, n_loop: int = None, n_gene: int = None):
        for key, value in cluster_result.items():
            if key.endswith('rare'):
                cluster_method, *_ = key.split('_')
                if n_gene is None:
                    self.rare_metric_baseline.loc[(dataset, n_loop, cluster_method), 'AllGenes'] = value
                else:
                    self.rare_metric.loc[(dataset, n_gene, n_loop, cluster_method), self.current_method] = value
            else:
                cluster_method, cluster_metric = key.split('_')
                if n_gene is None:
                    self.overall_metrics_baseline.loc[
                        (dataset, n_loop, cluster_method, cluster_metric), 'AllGenes'] = value
                else:
                    self.overall_metrics.loc[
                        (dataset, n_gene, n_loop, cluster_method, cluster_metric), self.current_method] = value

    def save(self):
        """
        Save this object to file_path.

        :return: None
        """

        with open(os.path.join(exp_cfg.record_path, 'pkl', 'ClusteringRecorder.pkl'), 'wb') as f:
            pickle.dump(self, f)

        self.overall_metrics.to_excel(os.path.join(exp_cfg.record_path, 'xlsx', 'clustering_overall_metrics.xlsx'))
        self.rare_metric.to_excel(os.path.join(exp_cfg.record_path, 'xlsx', 'clustering_rare_metrics.xlsx'))
        # baseline results
        self.overall_metrics_baseline.to_excel(
            os.path.join(exp_cfg.record_path, 'xlsx', 'clustering_overall_metrics_baseline.xlsx'))
        self.rare_metric_baseline.to_excel(
            os.path.join(exp_cfg.record_path, 'xlsx', 'clustering_rare_metrics_baseline.xlsx'))

        print(f"{now()}: Finished and saved record to {exp_cfg.record_path}\n\n\n")  # UTC + 8 hours = Beijing time


class MasterRecorder:
    """
    Contains all base recorders.
    """

    def __init__(self, measurements: List[str], methods: List[str]):
        print(f"{(datetime.datetime.now()).strftime('%Y-%m-%d %H:%M:%S')}: "
              f"evaluation starts.")  # UTC + 8 hours = Beijing time
        self._check_record_path()
        self.measurements = []
        for measurement, datasets in exp_cfg.measurements.items():
            if measurement in measurements:
                self.measurements.append(measurement)
                if measurement == 'population_demixing':
                    self.mixture = MixtureRecorder(datasets=datasets, methods=methods)
                elif measurement == 'computation_time':
                    self.time = TimeRecorder(datasets=datasets, methods=methods)
                elif measurement == 'intra-classification':
                    self.intra = ClassificationRecorder(datasets=datasets, methods=methods, assign_type='intra')
                elif measurement == 'inter-classification':
                    self.inter = ClassificationRecorder(datasets=datasets, methods=methods, assign_type='inter')
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
            elif measurement == 'intra-classification':
                self.intra.save()
            elif measurement == 'inter-classification':
                self.inter.save()
            elif measurement == 'clustering':
                self.cluster.save()
            elif measurement == 'marker_discovery':
                self.marker.save()
            elif measurement == 'batch_correction':
                self.correct.save()
            else:
                raise ValueError(f"{measurement} is not appropriate.")

    @staticmethod
    def _check_record_path():
        if 'xlsx' not in os.listdir(exp_cfg.record_path):
            os.mkdir(os.path.join(exp_cfg.record_path, 'xlsx'))

        if 'pkl' not in os.listdir(exp_cfg.record_path):
            os.mkdir(os.path.join(exp_cfg.record_path, 'pkl'))
