import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from utils.utils import get_gene_names, head
import warnings
import os
# tree models
from sklearn.ensemble import RandomForestClassifier
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier
# nearest shrunken centroid
from sklearn.neighbors import NearestCentroid
from utils.nearest_centroid import nearest_centroid_select
from sklearn.model_selection import GridSearchCV
# seurat
from loess.loess_1d import loess_1d
# scGeneFit
from utils.scGeneFit import scGeneFit
# fisher score
from utils.fisher_score import fisher_score


def most_important_genes(importances, feature_num, all_features):
    """
    Select max feature_num important genes according to importance array.

    :param importances: importance of each gene
    :param feature_num: the number of selected genes
    :param all_features: names of all genes
    :return: the max feature_num important genes
    """
    if importances.shape[0] != all_features.shape[0]:
        warnings.warn('The length of importance and gene names are not the same. Please check again!', RuntimeWarning)
        return None
    else:
        idx = np.argsort(importances)[::-1][:feature_num]
        return all_features[idx], importances[idx]


def select_features(data_name, feature_num, method, all_features, X, y):
    """
    Select max feature_num important genes from data.

    :param data_name: dataset name
    :param feature_num: the number of selected genes
    :param method: name of feature selection method
    :param all_features: names of all genes
    :param X: count matrix
    :param y: cell types
    :return: the max feature_num important genes and the corresponding importance (except cellassign)
    """
    head(name=method)
    if method == 'rf':  # random forest
        forest = RandomForestClassifier(n_estimators=100, n_jobs=15, random_state=0, verbose=0).fit(X, y)
        return most_important_genes(forest.feature_importances_, feature_num, all_features)
    elif method == 'lgb':
        lgb = LGBMClassifier(n_jobs=15, random_state=0).fit(X, y)
        return most_important_genes(lgb.feature_importances_, feature_num, all_features)
    elif method == 'xgb':
        le = LabelEncoder()
        y = le.fit_transform(y)
        xgb = XGBClassifier(eval_metric='mlogloss', use_label_encoder=False, nthread=15).fit(X, y)
        return most_important_genes(xgb.feature_importances_, feature_num, all_features)
    elif method == 'seurat':
        mean, var = X.mean(axis=0), X.var(axis=0)
        mean_fit, var_fit, weights = loess_1d(np.log10(mean), np.log10(var), frac=0.3, degree=2)
        z = (X - mean) / (10 ** (var_fit / 2))
        z[z > X.shape[0] ** 0.5] = X.shape[0] ** 0.5
        z = np.var(np.squeeze(z), axis=0)
        return most_important_genes(z, feature_num, all_features)
    elif method == 'var':
        return most_important_genes(np.var(X, axis=0), feature_num, all_features)
    elif method == 'cv2':
        cv2 = np.square(np.std(X, axis=0) / np.mean(X, axis=0))
        print('Number of genes whose mean are less than 1e-4: {}'.format(np.sum(np.squeeze(np.mean(X, axis=0)) < 1e-4)))
        return most_important_genes(cv2, feature_num, all_features)
    elif method == 'nsc':
        th_dict = {'shrink_threshold': np.arange(0.0, 1.01, 0.01)}
        gs = GridSearchCV(NearestCentroid(), param_grid=th_dict, cv=3, scoring='balanced_accuracy').fit(X, y)
        print('best score:{}, best threshold:{}'.format(gs.best_score_, gs.best_params_['shrink_threshold']))
        var = nearest_centroid_select(X, y, shrink_threshold=gs.best_params_['shrink_threshold'])
        return most_important_genes(var, feature_num, all_features)
    elif method == 'm3drop':
        os.system("Rscript scRNA-FeatureSelection/utils/RCode/M3Drop.R " + data_name)
        gene_path = 'scRNA-FeatureSelection/tempData/' + data_name + '_markers_m3drop.csv'
        genes = pd.read_csv(gene_path, usecols=[1, 3]).sort_values(by='p.value', ascending=True)
        return get_gene_names(genes['Gene'].values[:feature_num]), 1 - genes['p.value'].values[:feature_num]
    elif method == 'cellassign':
        os.system("export MKL_THREADING_LAYER=GNU && Rscript scRNA-FeatureSelection/utils/RCode/CellAssign.R " + data_name)
        gene_path = 'scRNA-FeatureSelection/tempData/' + data_name + '_markers_cellassign.csv'
        return get_gene_names(np.loadtxt(gene_path, dtype=np.object, delimiter=',', usecols=[0], skiprows=1))
    elif method == 'deviance':
        os.system("Rscript scRNA-FeatureSelection/utils/RCode/Deviance.R " + data_name)
        gene_path = 'scRNA-FeatureSelection/tempData/' + data_name + '_markers_deviance.csv'
        gene_and_importance = np.loadtxt(gene_path, dtype=np.object, delimiter=',', usecols=[0, 1], skiprows=1)
        return get_gene_names(gene_and_importance[:, 0]), gene_and_importance[:, 1].astype(np.float)
    elif method == 'monocle3':
        os.system("Rscript scRNA-FeatureSelection/utils/RCode/monocle3.R " + data_name)
        gene_path = 'scRNA-FeatureSelection/tempData/' + data_name + '_markers_monocle3.csv'
        gene_and_importance = pd.read_csv(gene_path, usecols=['gene_id', 'pseudo_R2']).values
        selected_genes_num = gene_and_importance.shape[0]
        print(f'Monocle3 selected {selected_genes_num} marker genes.')
        if selected_genes_num > feature_num:
            return get_gene_names(gene_and_importance[:feature_num, 0]), gene_and_importance[:feature_num, 1].astype(np.float)
        else:
            return get_gene_names(gene_and_importance[:, 0]), gene_and_importance[:, 1]
    elif method == 'scGeneFit':
        le = LabelEncoder()
        y_encoded = le.fit_transform(y)
        importances = scGeneFit(X, y_encoded, 0)
        return most_important_genes(importances, feature_num, all_features)
    elif method == 'fisher_score':
        importances = fisher_score(X, y)
        return most_important_genes(importances, feature_num, all_features)
    else:
        warnings.warn("Parameter 'method' is wrong. Please check again.")
        return None
