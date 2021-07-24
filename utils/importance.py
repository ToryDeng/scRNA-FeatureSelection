from typing import Tuple, Optional, Union, List
from collections import defaultdict
import anndata as ad
import pandas as pd
import timeout_decorator
from utils.utils import get_gene_names, save_data, delete, normalize_importances
from utils.record import TimeRecorder
from config import AssignConfig, ClusterConfig, exp_cfg, assign_cfg
# tree models
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier
# var, cv2
import numpy as np
# seurat
import scanpy as sc
# scGeneFit
from selection.scgenefit import get_importance
# fisher score
from selection.fisher_score import fisher_score
# nearest shrunken centroid
from sklearn.neighbors import NearestCentroid
from selection.nearest_centroid import nearest_centroid_select
from sklearn.model_selection import GridSearchCV
# execute R methods
import os


def most_important_genes(importances: Optional[np.ndarray],
                         all_features: Union[np.ndarray, List[Tuple]]
                         ) -> List[Tuple]:
    """
    Select max feature_num important genes according to importance array.

    :param importances: importance of each gene
    :param all_features: names of all genes
    :return: list of the most important genes and importance of which number equals to
    the corresponding element in experiment configuration (descending, not normalized)
    """
    if importances is None:
        return all_features
    else:
        if importances.shape[0] != all_features.shape[0]:
            raise RuntimeError(f"The length of importance({importances.shape[0]}) and gene names"
                               f"({all_features.shape[0]}) are not the same. Please check again!")
        all_idx = np.argsort(importances)[::-1]
        fea_impo_list = []
        for n_gene in exp_cfg.n_genes:
            if importances.shape[0] >= n_gene:
                fea_impo_list.append((all_features[all_idx[:n_gene]], importances[all_idx[:n_gene]]))
            else:
                fea_impo_list.append((all_features[all_idx], importances[all_idx]))
        return fea_impo_list


def execute_Rscript(data_name: str,
                    method: str,
                    feature_num=exp_cfg.n_genes[-1],
                    show_detail: bool = False
                    ) -> str:
    """
    Execute R script for methods in R.

    :param data_name: the dataset name
    :param method: feature selection method
    :param feature_num: number of features
    :param show_detail: whether to show details during execution
    :return: path of selected markers
    """
    cmd = ['Rscript utils/RCode/cal_feature_importance.R', data_name, method, str(feature_num)]
    if not show_detail:
        cmd.append('>& /dev/null')
    command = ' '.join(cmd)
    os.system(command) if method != 'cellassign' else os.system('export MKL_THREADING_LAYER=GNU && ' + command)
    gene_path = f'tempData/{data_name}_markers_{method}.csv'
    return gene_path


@timeout_decorator.timeout(seconds=exp_cfg.max_timeout)
def cal_feature_importance(method: str,
                           adata: ad.AnnData,
                           recorder: Optional[TimeRecorder],
                           config: Union[AssignConfig, ClusterConfig]
                           ):
    """
    Calculate importance of each feature using specified method.
    If running time is over the given seconds, this function will raise TimeoutError.

    :param method: feature selection method
    :param adata: the anndata instance containing raw and norm data
    :param recorder: performance recorder
    :param config: assign configuration or cluster configuration
    :return: importances and feature names, respectively.
    """
    assert config is not None, "config must exist."
    if config.method_lan[method] == 'python':  # python
        y, all_features = adata.obs['celltype'].values, adata.var_names.values
        if config.method_on[method] == 'raw':
            X = adata.raw.X
            adata = adata.raw.to_adata()
        else:
            X = adata.X
        if recorder is not None:
            recorder.py_method_start()
        if method == 'rf':  # random forest
            forest = RandomForestClassifier(n_jobs=15, random_state=0, verbose=0).fit(X, y)
            importance = forest.feature_importances_
        elif method == 'lgb':
            lgb = LGBMClassifier(n_jobs=15, random_state=0).fit(X, y)
            importance = lgb.feature_importances_
        elif method == 'xgb':
            le = LabelEncoder()
            y = le.fit_transform(y)
            xgb = XGBClassifier(eval_metric='mlogloss', use_label_encoder=False, nthread=15).fit(X, y)
            importance = xgb.feature_importances_
        elif method == 'seurat':
            sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=exp_cfg.n_genes[-1])
            importance, all_features = adata.var.variances_norm.values, adata.var.variances_norm.index.values
        elif method == 'cellranger':
            sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=exp_cfg.n_genes[-1])
            importance, all_features = adata.var.dispersions_norm.values, adata.var.dispersions_norm.index.values
        elif method == 'var':
            importance = np.var(X, axis=0)
        elif method == 'cv2':
            print('Number of genes whose mean are less than 1e-4: {}'.format(
                np.sum(np.squeeze(np.mean(X, axis=0)) < 1e-4)))
            importance = np.square(np.std(X, axis=0) / np.mean(X, axis=0))
        elif method == 'scGeneFit':
            importance = get_importance(X, y, exp_cfg.n_genes[-1], verbose=False)
        elif method == 'fisher_score':
            importance = fisher_score(X, y)
        elif method == 'nsc':
            th_dict = {'shrink_threshold': np.arange(0.0, 1.01, 0.01)}
            gs = GridSearchCV(NearestCentroid(), param_grid=th_dict, cv=3, scoring='balanced_accuracy').fit(X, y)
            print('best score:{}, best threshold:{}'.format(gs.best_score_, gs.best_params_['shrink_threshold']))
            importance = nearest_centroid_select(X, y, shrink_threshold=gs.best_params_['shrink_threshold'])
        else:
            raise NotImplementedError(f"{method} has not been implemented.")
        if recorder is not None:
            recorder.py_method_end()
            recorder.record(adata.uns['data_name'], method)
    else:  # R
        # save data first
        if config.method_on[method] == 'raw':
            save_data(adata.raw.to_adata(), task='assign' if isinstance(config, AssignConfig) else 'cluster')
        else:
            save_data(adata, task='assign' if isinstance(config, AssignConfig) else 'cluster')

        if method == 'm3drop':
            genes = pd.read_csv(execute_Rscript(adata.uns['data_name'], method),
                                usecols=[1, 3]
                                ).sort_values(by='p.value', ascending=True)
            all_features, importance = get_gene_names(genes['Gene'].values), 1 - genes['p.value'].values
        elif method == 'deviance':
            gene_importance = np.loadtxt(execute_Rscript(adata.uns['data_name'], method),
                                         dtype=np.object_, delimiter=',', usecols=[0, 1], skiprows=1)
            all_features, importance = get_gene_names(gene_importance[:, 0]), gene_importance[:, 1].astype(np.float_)
        elif method == 'scmap':
            gene_importance = pd.read_csv(execute_Rscript(adata.uns['data_name'], method), usecols=[0, 2]).dropna()
            all_features, importance = gene_importance.feature_symbol.values, gene_importance.scmap_scores.values
        else:
            raise NotImplementedError(f"{method} has not been implemented yet.")
        if recorder is not None:
            recorder.record(adata.uns['data_name'], method)
    return importance, all_features


def select_genes(method: str,
                 adata: ad.AnnData,
                 recorder: Optional[TimeRecorder] = None,
                 config: Union[AssignConfig, ClusterConfig] = assign_cfg
                 ):
    """
        Run ensemble selection method or single selection method according to the method name.

        :param method: feature selection method
        :param adata: the anndata instance containing raw and norm data
        :param recorder: performance recorder
        :param config: assign configuration or cluster configuration
        :return: importances and feature names, respectively.
    """
    try:
        if '+' in method:  # ensemble learning
            base_methods = method.split('+')
            gene_dict = defaultdict(int)
            if exp_cfg.ensemble_mode == 'importance_sum':
                for base_method in base_methods:
                    # delete current files in tempData
                    delete('tempData/')
                    result = cal_feature_importance(base_method, adata, recorder=recorder, config=config)
                    if isinstance(result, list):  # base method is cellassign
                        raise RuntimeError(f"{base_method} can't generate feature importances!")
                    else:
                        importances, genes = result
                        # Normalization
                        importances_norm = normalize_importances(importances)
                        for gene, importance in zip(genes, importances_norm):
                            gene_dict[gene] += importance
            elif exp_cfg.ensemble_mode == 'count_sum':
                for base_method in base_methods:
                    # delete current files in tempData
                    delete('tempData/')
                    result = cal_feature_importance(base_method, adata, recorder=recorder, config=config)
                    selected_genes = result[-1][0] if isinstance(result, list) else result[0]
                    for gene in selected_genes:
                        gene_dict[gene] += 1
            else:
                raise NotImplementedError(f"{exp_cfg.ensemble_mode} has not been implemented yet.")
            sorted_result = pd.Series(gene_dict).sort_values(ascending=False).reset_index().T.to_numpy()
            final_result = most_important_genes(sorted_result[1], sorted_result[0])  # importance, all features
        else:  # single method
            final_result = most_important_genes(
                *cal_feature_importance(method, adata, recorder=recorder, config=config))
    except timeout_decorator.timeout_decorator.TimeoutError:
        final_result = None
    return final_result
