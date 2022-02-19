from typing import Tuple, Optional, Union, List
from collections import defaultdict
import anndata as ad
import pandas as pd
import timeout_decorator
from collections import Counter
from utils.utils import get_gene_names, save_data, delete, normalize_importances, now, save_genes, is_saved, load_genes
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
from selection import scgenefit
# fisher score
from selection.fisher_score import fisher_score
# nearest shrunken centroid
from sklearn.neighbors import NearestCentroid
from selection.nearest_centroid import nearest_centroid_select
from sklearn.model_selection import GridSearchCV
# execute R methods
from selection.R_methods import FEAST_select
import os


def most_important_genes(importances: Optional[np.ndarray],
                         all_features: Union[np.ndarray, List[Tuple]],
                         ) -> Union[List[Tuple], Tuple]:
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
        na_mask = np.isnan(importances.astype(np.float))
        importances, all_features = importances[~na_mask].astype(np.float), all_features[~na_mask]
        importances = np.where(importances == np.inf, 10e11, importances)
        all_idx = np.argsort(importances)[::-1]

        fea_impo_list = []  # [(feature, importance), (feature, importance), (feature, importance), ()]
        for n_gene in exp_cfg.n_genes:
            if importances.shape[0] >= n_gene:
                fea_impo_list.append((all_features[all_idx[:n_gene]], importances[all_idx[:n_gene]]))
            else:
                print(f"No enough genes ({importances.shape[0]}) compared with {n_gene} to be selected.")
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
    os.system(command)
    gene_path = f'tempData/{data_name}_genes_{method}.csv'
    return gene_path


# @timeout_decorator.timeout(seconds=exp_cfg.max_timeout, use_signals=False)
def cal_feature_importance(method: str,
                           adata: ad.AnnData,
                           recorder: Optional[TimeRecorder],
                           config: Union[AssignConfig, ClusterConfig],
                           use_rep: str = 'celltype',
                           ):
    """
        Calculate importance of each feature using specified method.
        If running time is over the given seconds, this function will raise TimeoutError.

        :param method: feature selection method
        :param adata: the anndata instance containing raw and norm data
        :param recorder: performance recorder
        :param config: assign configuration or cluster configuration
        :param use_rep:  string that represents cell type
        :return: unsorted importances and feature names, respectively.
    """
    assert config is not None, "config must exist."
    if config.method_lan[method] == 'python':  # python
        y, all_features = adata.obs[use_rep].values, adata.var_names.values
        if config.method_on[method] == 'raw':
            X = adata.raw.X
            adata = adata.raw.to_adata()
        else:
            X = adata.X
        if recorder is not None:
            recorder.py_method_start()
        if method == 'rf':  # random forest
            forest = RandomForestClassifier(n_jobs=-1, random_state=0, verbose=0).fit(X, y)
            importance = forest.feature_importances_
        elif method == 'lgb':
            lgb = LGBMClassifier(n_jobs=16, random_state=0).fit(X, y)
            importance = lgb.feature_importances_
        elif method == 'xgb':
            le = LabelEncoder()
            y = le.fit_transform(y)
            xgb = XGBClassifier(eval_metric='mlogloss', use_label_encoder=False, nthread=-1).fit(X, y)
            importance = xgb.feature_importances_
        elif method == 'seurat':
            sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=exp_cfg.n_genes[-1])
            importance, all_features = adata.var.dispersions_norm.values, adata.var.dispersions_norm.index.values
        elif method == 'seurat_v3':
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
            importance = scgenefit.get_importance(X, y, exp_cfg.n_genes[-1], verbose=False)
        elif method == 'fisher_score':
            importance = fisher_score(X, y)
        elif method == 'nsc':
            th_dict = {'shrink_threshold': np.arange(0.0, 1.01, 0.01)}
            gs = GridSearchCV(NearestCentroid(), param_grid=th_dict, cv=3, scoring='balanced_accuracy').fit(X, y)
            print('best score:{}, best threshold:{}'.format(gs.best_score_, gs.best_params_['shrink_threshold']))
            importance = nearest_centroid_select(X, y, shrink_threshold=gs.best_params_['shrink_threshold'])
        elif method == 'random':
            np.random.seed(exp_cfg.random_seed)
            importance = np.random.rand(X.shape[1])
        else:
            raise NotImplementedError(f"{method} has not been implemented in python.")
        if recorder is not None:
            recorder.py_method_end()
            recorder.record(adata.uns['data_name'], method)
    else:  # R
        # save data first
        if config.method_on[method] == 'raw':
            save_data(adata.raw.to_adata(), task='assign' if isinstance(config, AssignConfig) else 'cluster',
                      use_rep=use_rep)
        else:
            save_data(adata, task='assign' if isinstance(config, AssignConfig) else 'cluster', use_rep=use_rep)

        if method == 'm3drop':
            genes = pd.read_csv(execute_Rscript(adata.uns['data_name'], method),
                                usecols=[1, 3]).sort_values(by='p.value', ascending=True)
            all_features, importance = get_gene_names(genes['Gene'].values), 1 - genes['p.value'].values
        elif method == 'feast':
            # genes = FEAST_select(adata.raw.to_adata())
            genes = pd.read_csv(execute_Rscript(adata.uns['data_name'], method), usecols=[1, 2])
            all_features, importance = get_gene_names(genes['Gene'].values), genes['F_scores'].values
        elif method == 'deviance':
            gene_importance = np.loadtxt(execute_Rscript(adata.uns['data_name'], method),
                                         dtype=np.object_, delimiter=',', usecols=[0, 1], skiprows=1)
            all_features, importance = get_gene_names(gene_importance[:, 0]), gene_importance[:, 1].astype(np.float_)
        elif method == 'scmap':
            gene_importance = pd.read_csv(execute_Rscript(adata.uns['data_name'], method), usecols=[0, 2]).dropna()
            all_features, importance = gene_importance.feature_symbol.values, gene_importance.scmap_scores.values
        else:
            raise NotImplementedError(f"{method} has not been implemented yet in R.")
        if recorder is not None:
            recorder.record(adata.uns['data_name'], method)
    return importance, all_features


def _select_genes(method: str,
                  adata: ad.AnnData,
                  recorder: Optional[TimeRecorder] = None,
                  config: Union[AssignConfig, ClusterConfig] = assign_cfg,
                  use_rep: str = 'celltype',
                  use_saved: bool = True
                  ):
    if use_saved and is_saved(adata, method):
        final_result = load_genes(adata, method)
        print('Using previously saved genes and importances...')
    else:  # do not use saved genes or genes have not been saved
        try:
            if '+' in method:  # ensemble learning
                base_methods = method.split('+')
                gene_dict = defaultdict(int)
                if exp_cfg.ensemble_mode == 'importance_sum':
                    for base_method in base_methods:
                        # delete current files in tempData
                        delete('tempData/')
                        result = cal_feature_importance(base_method, adata, recorder=recorder, config=config)
                        if isinstance(result, list):
                            raise RuntimeError(f"{base_method} can't generate feature importances!")
                        else:
                            to_dropna = pd.Series(result[0], index=result[1]).dropna()
                            importances, genes = to_dropna.values, to_dropna.index.to_numpy()
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
                sorted_result = pd.Series(gene_dict).dropna().sort_values(ascending=False).reset_index().T.to_numpy()
                final_result = most_important_genes(sorted_result[1], sorted_result[0])  # importance, all features
            else:  # single method
                final_result = most_important_genes(
                    *cal_feature_importance(method, adata, recorder=recorder, config=config, use_rep=use_rep),
                    )
        except timeout_decorator.timeout_decorator.TimeoutError:
            final_result = None
            print(f"{method} is running out of time.")
        except FileNotFoundError as e:
            final_result = None
            print(f"The R method {method} failed on dataset {adata.uns['data_name']}.")
            print(e)
        except MemoryError:
            final_result = None
            print(f"{method} is running out of memory.")
        except ValueError as e:
            final_result = None
            print(f"{method} failed: {e}")
        if final_result is not None:
            save_genes(adata, method, final_result)
    return final_result


def select_genes(method: str,
                 adata: ad.AnnData,
                 recorder: Optional[TimeRecorder] = None,
                 config: Union[AssignConfig, ClusterConfig] = assign_cfg,
                 use_rep: str = 'celltype',
                 use_saved: bool = True,
                 select_by_batch = True
                 ):
    """
        Run ensemble selection method or single selection method according to the method name.

        :param use_rep: string that represents cell type
        :param method: feature selection method
        :param adata: the anndata instance containing raw and norm data
        :param recorder: performance recorder
        :param config: assign configuration or cluster configuration
        :return: importances and feature names, respectively.
    """
    if 'batch' in adata.obs and select_by_batch:
        ubatches = adata.obs['batch'].unique().categories.to_numpy()
        print(f"{adata.uns['data_name']} contains {ubatches.shape[0]} batches: {ubatches}")
        batch_features = []
        for ubatch in ubatches:
            batch_adata = adata[adata.obs['batch'] == ubatch].copy()
            # rename batch data set
            if '+' not in batch_adata.uns['data_name']:
                batch_adata.uns['data_name'] = '_'.join([batch_adata.uns['data_name'], ubatch])
            else:
                batch_adata.uns['data_name'] = ubatch

            print(f"{now()}: Start to apply {method} on batch {ubatch}, batch shape: {batch_adata.shape}.")
            selected = _select_genes(method, batch_adata, recorder, config, use_rep, use_saved)
            if selected is None:
                raise RuntimeWarning(f"batch {ubatch} is not used.")
            else:
                batch_features.append(pd.Series(data=selected[-1][1], index=selected[-1][0]))
        batch_features_rank = [batch_feature.rank() for batch_feature in batch_features]
        cnt = Counter()
        for batch_feature in batch_features:
            cnt.update(batch_feature.index)
        ties = pd.Series(cnt).sort_values(ascending=False)
        untied_rank = []
        for tie in np.unique(ties.values):
            rank_dict = defaultdict(list)
            for gene in ties[ties == tie].index:
                for batch_rank in batch_features_rank:
                    if gene in batch_rank.index:
                        rank_dict[gene].append(batch_rank[gene])
            median_rank_dict = {key: np.median(value) for key, value in rank_dict.items()}
            untied = pd.Series(median_rank_dict).sort_values(ascending=True)
            untied_rank.append(untied)
        untied_features = pd.concat(untied_rank[::-1])
        return most_important_genes(untied_features.values, untied_features.index.to_numpy())
    else:
        return _select_genes(method, adata, recorder, config, use_rep, use_saved)


