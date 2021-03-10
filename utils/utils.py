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
    description:
    for specific data_name, get the corresponding data.

    parameters:
    data_name: the dataset you want to get
    with_marker: whether to return marker genes
    norm: whether to normalize features(using the normalization in Seurat)
    scale_factor: size factor, default 1e4

    return:
    if with_marker=False, features(row:cell, col:gene, dataframe) and labels(dataframe) of raw data
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
