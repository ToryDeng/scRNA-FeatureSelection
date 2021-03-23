import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import warnings
import os


def filter_const_genes(X):
    """
    Remove constant genes.

    :param X: Count matrix in  dataframe format (row:cell, col:gene)
    :return: Filtered count matrix
    """
    mask = X.sum(axis=0) != 0
    print("Removed {} gene(s).".format(X.shape[1] - mask.sum()))
    return X.loc[:, mask]


def filter_const_cells(X, y):
    mask = X.sum(axis=1) != 0
    print("Removed {} cell(s).".format(X.shape[0] - mask.sum()))
    return X.loc[mask, :], y[mask.values]


def load_data(data_name, scale_factor=1e4):
    """
    For specific data_name, get the corresponding raw data. Remove constant genes (=0) and normalize it
    using the method in Seurat.

    :param data_name: the dataset name you want to get
    :param scale_factor: size factor, default 1e4
    :return: raw_features, norm_features, labels, markers
    """
    if data_name[:4] == 'PBMC':
        os.chdir('/home/tdeng/SingleCell/data/PBMC/integrated data')
        if len(data_name) == 4:
            data = pd.read_hdf('PBMC_AllCells_withLables.h5', key='AllCells')
            raw_features, labels = data.iloc[:, :-1], data.iloc[:, -1]
        elif len(data_name) > 4:
            raw_features = pd.read_csv('raw_features_sample' + data_name[4:] + '.csv', index_col=0)
            labels = pd.read_csv('raw_labels_sample' + data_name[4:] + '.csv', usecols=[1])
        else:
            print("parameter 'data_name' is wrong!")
            return None
        os.chdir('../')
        part1 = np.loadtxt('hsPBMC_markers_10x.txt', skiprows=1, usecols=[0], dtype=np.object, delimiter=',')
        part2 = np.loadtxt('blood_norm_marker.txt', skiprows=1, usecols=[1], dtype=np.object, delimiter=',')
        markers = np.union1d(part1, part2)
    elif data_name in ['muraro', 'segerstolpe', 'xin']:
        os.chdir('/home/tdeng/SingleCell/data/pancreas/separated data')
        raw_features = pd.read_csv(data_name.capitalize() + '_pancreas_filtered.csv', index_col=0)
        labels = pd.read_csv(data_name.capitalize() + '_trueID_filtered.csv', usecols=[1])
        os.chdir('../')
        markers = np.loadtxt('pancreasMarkerGenes.csv', skiprows=1, usecols=[0], dtype=np.object, delimiter=',')
    elif data_name == 'all_pancreas':
        os.chdir('/home/tdeng/SingleCell/data/pancreas/integrated data')
        raw_features = pd.read_csv('features.csv', index_col=0)
        labels = pd.read_csv('labels.csv', usecols=[1])
        os.chdir('../')
        markers = np.loadtxt('pancreasMarkerGenes.csv', skiprows=1, usecols=[0], dtype=np.object, delimiter=',')
    else:
        print("parameter 'data_name' is wrong!")
        return None
    raw_features = filter_const_genes(raw_features)
    raw_features, labels = filter_const_cells(X=raw_features, y=labels)
    norm_features = np.log1p(raw_features / raw_features.sum(1).values.reshape(raw_features.shape[0], 1) * scale_factor)
    os.chdir("../..")
    return raw_features, norm_features, labels, markers


def pancreas_filter_and_save_data():
    """
    Read separated data, remove cells which do not have a clear cell type and save filtered data.

    :return: None
    """
    path = '/home/tdeng/SingleCell/data/pancreas/separated data/'
    not_cell_type = ['unclear', 'not applicable', 'unclassified', 'co-expression', 'beta.contaminated',
                     'alpha.contaminated', 'delta.contaminated', 'gamma.contaminated', 'unassigned']

    for pan_type in ['muraro', 'segerstolpe', 'xin']:
        # get raw data
        features = pd.read_csv(path + pan_type.capitalize() + '_pancreas.csv', index_col=0).T
        labels = pd.read_csv(path + pan_type.capitalize() + '_trueID.csv', usecols=[1])
        # which cells should be kept
        is_clear = ~labels['x'].isin(not_cell_type).values
        filtered_features, filtered_labels = features.loc[is_clear, :], labels.loc[is_clear]
        # save filtered data
        filtered_features.to_csv(path + pan_type.capitalize() + '_pancreas_filtered.csv')
        filtered_labels.reset_index(drop=True).to_csv(path + pan_type.capitalize() + '_trueID_filtered.csv')
        print('{}: The raw data has {} cells. The filtered data has {} cells.'.format(
            pan_type.capitalize(), labels.shape[0], filtered_labels.shape[0]))


def get_gene_names(columns):
    """
    Get gene names array from features (dataframe).

    :param columns: dataframe.columns
    :return: an array contains gene names
    """
    if '__' in columns[0] and columns[0][0] != '"' and columns[0][-1] != '"':
        gene_names = pd.Series(columns).str.split('__', expand=True).iloc[:, 0].values
    elif '__' not in columns[0] and columns[0][0] == '"' and columns[0][-1] == '"':
        gene_names = np.char.strip(np.array(columns, dtype=np.str), '"')
    elif '__' in columns[0] and columns[0][0] == '"' and columns[0][-1] == '"':
        gene_names = np.char.strip(np.array(columns, dtype=np.str), '"')
        gene_names = pd.Series(gene_names).str.split('__', expand=True).iloc[:, 0].values
    elif '\t' in columns[0]:
        gene_names = pd.Series(columns).str.split('\t', expand=True).iloc[:, 1].values
    else:
        gene_names = np.array(columns, dtype=np.str)
    return gene_names


def save_filtered_data(X, Y, all_genes, selected_genes):
    """
    Save raw counts with only marker genes, labels and split data to training set and test set
    for evaluation using clustering and classification methods.

    :param X: raw counts
    :param Y: labels
    :param all_genes: all genes
    :param selected_genes: selected genes
    :return: None
    """
    # clustering
    mask = np.isin(all_genes, selected_genes)
    if mask.sum() == 0:
        warnings.warn("No gene is selected!", RuntimeWarning)
    X_selected, y = filter_const_cells(X.loc[:, mask], Y)  # some cells may have zero counts
    print('After gene selection, the dataset now has {} cells and {} genes.'.format(
        X_selected.shape[0], X_selected.shape[1]))
    X_selected.to_csv('scRNA-FeatureSelection/tempData/temp_X.csv')
    y.index = X_selected.index
    y.to_csv('scRNA-FeatureSelection/tempData/temp_y.csv')
    # classification
    X_train, X_test, y_train, y_test = train_test_split(X_selected, y, shuffle=True, test_size=0.3, random_state=2021)
    X_train.to_csv('scRNA-FeatureSelection/tempData/temp_X_train.csv')
    X_test.to_csv('scRNA-FeatureSelection/tempData/temp_X_test.csv')
    y_train.to_csv('scRNA-FeatureSelection/tempData/temp_y_train.csv')
    y_test.to_csv('scRNA-FeatureSelection/tempData/temp_y_test.csv')


def delAll(path):
    """
    Recursively delete files in a temporary folder.

    :param path: folder or file path
    :return: None
    """
    if os.path.isdir(path):
        files = os.listdir(path)
        for file in files:
            p = os.path.join(path, file)
            if os.path.isdir(p):
                # recursion
                delAll(p)
            else:
                os.remove(p)
    else:
        os.remove(path)


def PBMC_sample_to_csv(fraction):
    """
    Randomly sample some cells and save as csv file from PBMC dataset.

    :param fraction: double, sampling rate
    :return: None
    """
    path = '/home/tdeng/SingleCell/data/PBMC/integrated data/'
    data = pd.read_hdf(path + 'PBMC_AllCells_withLables.h5', key='AllCells')
    samples = data.sample(frac=fraction, random_state=2020).reset_index(drop=True)
    samples.columns.name = ''
    samples.columns = pd.Series(samples.columns).str.split('\t', expand=True).iloc[:, 1].values
    samples.iloc[:, :-1].to_csv(path + 'raw_features_sample' + str(int(fraction * 100)) + '.csv')
    samples.iloc[:, -1].to_csv(path + 'raw_labels_sample' + str(int(fraction * 100)) + '.csv')
    print('{} PBMC cells have been saved.'.format(samples.shape[0]))


def random_choice(X, y, size):
    """
    Randomly choose cells and labels.

    :param X: ndarray, count matrix
    :param y: ndarray, labels
    :param size: number of cells that are randomly chosen
    :return: random cells and corresponding labels
    """
    rand_idx = np.random.choice(a=X.shape[0], replace=False, size=size)
    return X[rand_idx], y[rand_idx]
