import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_fscore_support
from sklearn.neighbors import NearestCentroid, NearestNeighbors
from sklearn.metrics import adjusted_rand_score, classification_report
import rpy2.robjects as robjects
from scipy.optimize import linprog
# LightGBM
from lightgbm import LGBMClassifier
# xgboost
from xgboost import XGBClassifier
from loess.loess_1d import loess_1d
import os


def get_data(data_name, with_marker=False, norm=False, scale_factor=1e4):
    """
    For specific data_name, get the corresponding data.

    ----------
    :param data_name: the dataset you want to get
    :param with_marker: whether to return marker genes
    :param norm: whether to normalize features(using the normalization in Seurat)
    :param scale_factor: size factor, default 1e4
    :return:if with_marker=False, features(row:cell, col:gene, dataframe) and labels(dataframe) of raw data.
    if with_marker=True, features(row:cell, col:gene, dataframe), labels(dataframe) and marker genes(array) of raw data
    """
    if data_name[:4] == 'PBMC':
        os.chdir('/home/tdeng/SingleCell/data/PBMC/integrated data')
        if data_name == 'PBMC':
            data = pd.read_hdf('PBMC_AllCells_withLables.h5', key='AllCells')
            features, labels = data.iloc[:, :-1], data.iloc[:, -1]
        elif 0 < int(data_name[4:]) < 100:
            features = pd.read_csv('raw_features_sample' + data_name[4:] + '.csv', index_col=0)
            labels = pd.read_csv('raw_labels_sample' + data_name[4:] + '.csv', usecols=[1])
        else:
            print("parameter 'data_name' is wrong!")
            return None

        if norm:
            features = np.log1p(features / features.sum(1).values.reshape(features.shape[0], 1) * scale_factor)
        if with_marker:
            part1 = np.squeeze(
                pd.read_csv('/home/tdeng/SingleCell/data/PBMC/hsPBMC_markers_10x.txt', usecols=[0]).values)
            part2 = np.squeeze(
                pd.read_csv('/home/tdeng/SingleCell/data/PBMC/blood_norm_marker.txt', usecols=[1]).values)
            markers = np.union1d(part1, part2)
            return features, labels, markers
        return features, labels
    elif data_name in ['muraro', 'segerstolpe', 'xin']:
        os.chdir('/home/tdeng/SingleCell/data/pancreas/separated data')
        features = pd.read_csv(data_name.capitalize() + '_pancreas_filtered.csv', index_col=0)
        labels = pd.read_csv(data_name.capitalize() + '_trueID_filtered.csv', usecols=[1])
        if norm:
            features = np.log1p(features / features.sum(1).values.reshape(features.shape[0], 1) * scale_factor)
        if with_marker:
            markers = np.squeeze(pd.read_csv('/home/tdeng/SingleCell/data/pancreas/pancreasMarkerGenes.csv',
                                             usecols=[0]).values)
            return features, labels, markers
        return features, labels
    elif data_name == 'all_pancreas':
        os.chdir('/home/tdeng/SingleCell/data/pancreas/integrated data')
        features = pd.read_csv('features.csv', index_col=0)
        labels = pd.read_csv('labels.csv', usecols=[1])
        if norm:
            features = np.log1p(features / features.sum(1).values.reshape(features.shape[0], 1) * scale_factor)
        if with_marker:
            markers = np.squeeze(pd.read_csv('/home/tdeng/SingleCell/data/pancreas/pancreasMarkerGenes.csv',
                                             usecols=[0]).values)
            return features, labels, markers
        return features, labels
    else:
        print("parameter 'data_name' is wrong!")
        return None


def filter_const_genes(X):
    """
    Remove constant genes.

    ----------
    :param X: Count matrix in  dataframe format (row:cell, col:gene)
    :return: Filtered count matrix
    """
    return X.loc[:, X.mean(axis=0) != 0]


def pancreas_filter_and_save_data():
    """
    Read separated data, remove cells which do not have a clear cell type and save filtered data.

    ----------
    :return: None
    """
    os.chdir('/home/tdeng/SingleCell/data/pancreas/separated data')
    not_cell_type = ['unclear', 'not applicable', 'unclassified', 'co-expression', 'beta.contaminated',
                     'alpha.contaminated', 'delta.contaminated', 'gamma.contaminated']

    for pan_type in ['muraro', 'segerstolpe', 'xin']:
        # get raw data
        features = pd.read_csv(pan_type.capitalize() + '_pancreas.csv', index_col=0).T
        labels = pd.read_csv(pan_type.capitalize() + '_trueID.csv', usecols=[1])
        # which cells should be kept
        is_clear = ~labels['x'].isin(not_cell_type).values
        filtered_features, filtered_labels = features.loc[is_clear, :], labels.loc[is_clear]
        # save filtered data
        print('The raw data has {} cells. The filtered data has {} cells.'.format(labels.shape[0],
                                                                                  filtered_labels.shape[0]))
        filtered_features.to_csv(pan_type.capitalize() + '_pancreas_filtered.csv')
        filtered_labels.reset_index(drop=True).to_csv(pan_type.capitalize() + '_trueID_filtered.csv')


def get_gene_names(columns):
    """
    Get gene names array from features (dataframe).
    :param columns: dataframe.columns
    :return: an array contains gene names
    """
    if '__' in columns[0]:
        gene_names = pd.Series(columns).str.split('__', expand=True).iloc[:, 0].values
    elif '\t' in columns[0]:
        gene_names = pd.Series(columns).str.split('\t', expand=True).iloc[:, 1].values
    else:
        gene_names = pd.Series(columns).values
    return gene_names


X, y = get_data('muraro')
print(X.shape)
X = filter_const_genes(X)
print(X.shape)