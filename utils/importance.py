from typing import Tuple, Optional, Union
from collections import defaultdict
import anndata as ad
import pandas as pd
from utils.utils import get_gene_names, save_data, delete
from utils.record import PerformanceRecorder
from config import AssignConfig, ClusterConfig, exp_cfg
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
from models.scGeneFit import scGeneFit
# fisher score
from models.fisher_score import fisher_score
# nearest shrunken centroid
from sklearn.neighbors import NearestCentroid
from models.nearest_centroid import nearest_centroid_select
from sklearn.model_selection import GridSearchCV
# execute R methods
import os


def most_important_genes(importances: np.ndarray, feature_num: int, all_features: np.ndarray) \
        -> Tuple[np.ndarray, np.ndarray]:
    """
    Select max feature_num important genes according to importance array.

    :param importances: importance of each gene
    :param feature_num: the number of selected genes
    :param all_features: names of all genes
    :return: the max feature_num important genes
    """
    if importances.shape[0] != all_features.shape[0]:
        raise RuntimeError(f"The length of importance({importances.shape[0]}) and gene names({all_features.shape[0]}) "
                           f"are not the same. Please check again!")
    else:
        idx = np.argsort(importances)[::-1][:feature_num]
        return all_features[idx], importances[idx]


def execute_Rscript(data_name: str, method: str, show_detail: bool = False) -> str:
    cmds = ['Rscript utils/RCode/select_markers.R', data_name, method, str(exp_cfg.n_genes)]
    if not show_detail:
        cmds.append('>& /dev/null')
    command = ' '.join(cmds)
    os.system(command) if method != 'cellassign' else os.system('export MKL_THREADING_LAYER=GNU && ' + command)
    gene_path = f'tempData/{data_name}_markers_{method}.csv'
    return gene_path


def select_features(feature_num: int, method: str, adata: ad.AnnData = None,
                    recorder: Optional[PerformanceRecorder] = None,
                    config: Union[AssignConfig, ClusterConfig] = None) \
        -> Union[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, ], None]:
    """
    For python method, use adata to select features; For R method, save data as csv and call Rscript.

    :param feature_num:
    :param method:
    :param adata:
    :param recorder: record performance (time mainly)
    :param config: config for assign task or clustering task
    :return: selection result
    """
    assert recorder is not None and config is not None, "recorder and config must exist."
    if config.method_lan[method] == 'python':  # python
        y, all_features = adata.obs['type'].values, adata.var_names.values
        if config.method_on[method] == 'raw':
            X = adata.raw.X
            adata = adata.raw.to_adata()
        else:
            X = adata.X
        recorder.cmpt_time_start()
        if method == 'rf':  # random forest
            forest = RandomForestClassifier(n_estimators=100, n_jobs=15, random_state=0, verbose=0).fit(X, y)
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
            sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=feature_num)
            importance, all_features = adata.var.variances_norm.values, adata.var.variances_norm.index.values
        elif method == 'var':
            importance = np.var(X, axis=0)
        elif method == 'cv2':
            print('Number of genes whose mean are less than 1e-4: {}'.format(
                np.sum(np.squeeze(np.mean(X, axis=0)) < 1e-4)))
            importance = np.square(np.std(X, axis=0) / np.mean(X, axis=0))
        elif method == 'scGeneFit':
            le = LabelEncoder()
            y_encoded = le.fit_transform(y)
            importance = scGeneFit(X, y_encoded, 0)
        elif method == 'fisher_score':
            importance = fisher_score(X, y)
        elif method == 'nsc':
            th_dict = {'shrink_threshold': np.arange(0.0, 1.01, 0.01)}
            gs = GridSearchCV(NearestCentroid(), param_grid=th_dict, cv=3, scoring='balanced_accuracy').fit(X, y)
            print('best score:{}, best threshold:{}'.format(gs.best_score_, gs.best_params_['shrink_threshold']))
            importance = nearest_centroid_select(X, y, shrink_threshold=gs.best_params_['shrink_threshold'])
        else:
            raise NotImplementedError(f"{method} has not been implemented.")
        recorder.cmpt_time_end()
        return most_important_genes(importance, feature_num, all_features)
    else:  # R
        # save data first
        if config.method_on[method] == 'raw':
            save_data(adata.raw.to_adata(), task='assign' if isinstance(config, AssignConfig) else 'cluster')
        else:
            save_data(adata, task='assign' if isinstance(config, AssignConfig) else 'cluster')

        if method == 'm3drop':
            genes = pd.read_csv(execute_Rscript(adata.uns['data_name'], method), usecols=[1, 3]).sort_values(
                by='p.value', ascending=True)
            recorder.record_cmpt_time_from_csv()
            return get_gene_names(genes['Gene'].values[:feature_num]), 1 - genes['p.value'].values[:feature_num]
        elif method == 'cellassign':
            gene_path = execute_Rscript(adata.uns['data_name'], method)
            recorder.record_cmpt_time_from_csv()
            return (get_gene_names(np.loadtxt(gene_path, dtype=np.object, delimiter=',', usecols=[0], skiprows=1)),)
        elif method == 'deviance':
            gene_and_importance = np.loadtxt(execute_Rscript(adata.uns['data_name'], method), dtype=np.object_,
                                             delimiter=',', usecols=[0, 1], skiprows=1)
            recorder.record_cmpt_time_from_csv()
            return get_gene_names(gene_and_importance[:, 0]), gene_and_importance[:, 1].astype(np.float_)
        elif method == 'monocle3':
            gene_and_importance = pd.read_csv(execute_Rscript(adata.uns['data_name'], method),
                                              usecols=['gene_id', 'pseudo_R2']).values
            selected_genes_num = gene_and_importance.shape[0]
            recorder.record_cmpt_time_from_csv()
            if selected_genes_num > feature_num:
                return get_gene_names(gene_and_importance[:feature_num, 0]), \
                       gene_and_importance[:feature_num, 1].astype(np.float)
            else:
                return get_gene_names(gene_and_importance[:, 0]), gene_and_importance[:, 1]
        else:
            raise NotImplementedError(f"{method} has not been implemented.")


def select_genes(feature_num: int, method: str, adata: ad.AnnData = None,
                 recorder: Optional[PerformanceRecorder] = None,
                 config: Union[AssignConfig, ClusterConfig] = None)\
        -> Union[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, ], None]:
    if '+' in method:  # ensemble learning
        base_methods = method.split('+')
        gene_dict = defaultdict(int)
        if exp_cfg.ensemble_mode == 'importance_sum':
            for base_method in base_methods:
                # delete current files in tempData
                delete('tempData/')
                result = select_features(feature_num, base_method, adata, recorder=recorder, config=config)
                if len(result) != 2:
                    raise RuntimeError(f"{base_method} can't generate feature importances!")
                else:
                    genes, importances = result
                    importances_norm = (importances - importances.min()) / (importances.max() - importances.min())  # Normalization
                    for gene, importance in zip(genes, importances_norm):
                        gene_dict[gene] += importance
        elif exp_cfg.ensemble_mode == 'count_sum':
            for base_method in base_methods:
                # delete current files in tempData
                delete('tempData/')
                result = select_features(exp_cfg.n_genes, base_method, adata, recorder=recorder, config=config)
                for gene in result[0]:
                    gene_dict[gene] += 1
        else:
            raise NotImplementedError(f"{exp_cfg.ensemble_mode} has not been implemented yet.")
        sorted_result = pd.Series(gene_dict).sort_values(ascending=False).iloc[:1000].reset_index().T.to_numpy()
        return sorted_result[0], sorted_result[1]
    else:  # single method
        return select_features(feature_num, method, adata, recorder=recorder, config=config)

