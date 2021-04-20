import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.model_selection import train_test_split
from config import data_cfg
import warnings
import datetime
import os


def filter_const_genes(X) -> pd.DataFrame:
    """
    Remove zero count genes.

    :param X: Count matrix in  dataframe format (row:cell, col:gene)
    :return: Filtered count matrix
    """
    mask = X.sum(axis=0) != 0
    print("Removed {} gene(s).".format(X.shape[1] - mask.sum()))
    return X.loc[:, mask]


def filter_const_cells(X, y):
    """
    Remove zero count cells.

    :param X: count matrix
    :param y: cell types
    :return: filtered count matrix and cell types
    """
    assert X.shape[0] == y.shape[0], "Inconsistency between X and y in dim 0."
    mask = X.sum(axis=1) != 0
    print("Removed {} cell(s).".format(X.shape[0] - mask.sum()))
    return X.loc[mask, :], y[mask.values]


def polish(df, dataset: str):
    """
    For PBMC dataset, simplify the gene names. For pancreas dataset, filter cells.

    :param df: the original dataframe
    :param dataset: dataset name
    :return: the polished dataframe
    """
    if dataset == 'PBMC':
        df.columns.name = ''
        names = get_gene_names(df.columns).tolist()  # may have duplicated names
        df.columns = [v + '-' + str(names[:i].count(v) + 1) if names.count(v) > 1 else v for i, v in enumerate(names)]
        return df
    elif dataset == 'pancreas':
        # which cells should be kept
        return df.loc[~df['label'].isin(data_cfg.remove_types).values]


def load_data(data_name: str, scale_factor=1e4):
    """
    For specific data_name, get the corresponding raw data. Remove constant genes (=0) and normalize it
    using the method in Seurat.

    :param data_name: the dataset name you want to get
    :param scale_factor: size factor, default 1e4
    :return: raw_features, norm_features, labels, markers
    """
    if data_name[:4] == 'PBMC':
        assert len(data_name) >= 4, "parameter 'data_name' is wrong!"
        data = pd.read_hdf(data_cfg.PBMC_path, key='AllCells').reset_index(drop=True)
        if len(data_name) == 4:
            data = polish(data, dataset='PBMC')
            raw_features, labels = data.iloc[:, :-1], data.iloc[:, -1]
        else:
            rate = float(data_name[4:])
            assert 0 < rate < 100, "The percentage is wrong!"
            data = polish(data, dataset='PBMC')
            X, y = data.iloc[:, :-1], data.iloc[:, -1]
            _, raw_features, __, labels = train_test_split(X, y, test_size=rate / 100, random_state=2020, stratify=y)
        os.chdir('/home/tdeng/SingleCell/data/PBMC/')
        part1 = np.loadtxt('hsPBMC_markers_10x.txt', skiprows=1, usecols=[0], dtype=np.object_, delimiter=',')
        part2 = np.loadtxt('blood_norm_marker.txt', skiprows=1, usecols=[1], dtype=np.object_, delimiter=',')
        markers = np.union1d(part1, part2)
    elif data_name in ['muraro', 'segerstolpe', 'xin']:
        data = polish(pd.read_hdf(data_cfg.pancreas_path, key=data_name.capitalize() + '_pancreas'), dataset='pancreas')
        raw_features, labels = data.iloc[:, :-1], data.iloc[:, -1]
        os.chdir('/home/tdeng/SingleCell/data/pancreas/')
        markers = np.loadtxt('pancreasMarkerGenes.csv', skiprows=1, usecols=[0], dtype=np.object_, delimiter=',')
    elif data_name == 'all_pancreas':
        os.chdir('/home/tdeng/SingleCell/data/pancreas/integrated data')
        raw_features = pd.read_csv('features.csv', index_col=0)
        labels = pd.read_csv('labels.csv', usecols=[1])
        os.chdir('../')
        markers = np.loadtxt('pancreasMarkerGenes.csv', skiprows=1, usecols=[0], dtype=np.object_, delimiter=',')
    else:
        print("parameter 'data_name' is wrong!")
        return None
    raw_features = filter_const_genes(raw_features).round()
    raw_features, labels = filter_const_cells(X=raw_features, y=labels)
    norm_features = np.log1p(raw_features / raw_features.sum(1).values.reshape(raw_features.shape[0], 1) * scale_factor)
    os.chdir("../..")  # return to SingleCell dir
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
        gene_names = np.char.strip(np.array(columns, dtype=np.str_), '"')
    elif '__' in columns[0] and columns[0][0] == '"' and columns[0][-1] == '"':
        gene_names = np.char.strip(np.array(columns, dtype=np.str_), '"')
        gene_names = pd.Series(gene_names).str.split('__', expand=True).iloc[:, 0].values
    elif '\t' in columns[0]:
        gene_names = pd.Series(columns).str.split('\t', expand=True).iloc[:, 1].values
    else:
        gene_names = np.array(columns, dtype=np.str_)
    return gene_names


def delete(path: str):
    """
    Recursively delete files in a temporary folder if 'path' is a dir, or delete the file when 'path'
    is a file path.

    :param path: folder or file path
    :return: None
    """
    if os.path.isdir(path):
        files = os.listdir(path)
        for file in files:
            p = os.path.join(path, file)
            if os.path.isdir(p):
                # recursion
                delete(p)
            else:
                os.remove(p)
    else:
        os.remove(path)


def PBMC_sample_to_csv(fraction: float):
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


def cal_marker_num_MRR(trusted_features, selected_features, rank=True):
    """
    Calculate number of selected marker genes and MRR if rank == True.

    :param trusted_features: all marker genes
    :param selected_features: selected marker genes
    :param rank: bool. Whether to calculate MRR
    :return: if rank == True, return marker num and MRR, Otherwise return marker num
    """
    marker_genes_found = np.intersect1d(trusted_features, selected_features).shape[0]
    if rank:
        rank_list = np.argwhere(np.isin(selected_features, trusted_features))
        if len(rank_list) == 0:
            warnings.warn("Can not calculate MRR because no marker gene is selected! MRR is set to 0.", RuntimeWarning)
            MRR = 0
        else:
            MRR = np.sum(1 / (rank_list + 1)) / rank_list.shape[0]
        return marker_genes_found, MRR
    else:  # len(selected_result) == 1
        warnings.warn("The method can't obtain gene importance! MRR can't be calculated and is set to 0.",
                      RuntimeWarning)
        return marker_genes_found


def PerformanceRecord(methods: list, task: str):
    """
    Generate a performance record.

    :param methods: the list of gene selection methods
    :param task: 'assign' or 'clustering'
    :return: If task == 'assign', return a dataframe with MRR, marker_genes_found and 3 F1-scores.
    If task == 'clustering', return a dataframe with MRR, marker_genes_found and 2 ARI.
    """
    assert len(methods) > 0, "methods list is null."
    assert task in ['assign', 'clustering'], "Parameter 'task' is wrong."
    if task == 'assign':
        index = ['MRR', 'marker_genes_found', 'scmap_cluster_F1', 'scmap_cell_F1', 'singlecellnet_F1']
    elif task == 'clustering':
        index = ['MRR', 'marker_genes_found', 'seurat_ARI', 'sc3_ARI']
    else:
        index = None
    return pd.DataFrame(data=np.zeros((len(index), len(methods))), index=index, columns=methods)


def save_raw_data(X_train, X_test, y_train, y_test, task: str):
    """
    Save raw data for classification task and clustering task in tempData dir.

    :param X_train: count matrix of training cells when task == 'assign',
     or count matrix of all cells when task == 'clustering'.
    :param X_test: count matrix of test cells
    :param y_train: types of training cells when task == 'assign' or types of  all cells when task == 'clustering'.
    :param y_test: types of test cells
    :param task:
    :return: None
    """
    assert (task in ['assign', 'clustering']), "Parameter 'task' is wrong."
    if task == 'assign':
        X_train.to_csv('scRNA-FeatureSelection/tempData/temp_X_train.csv')
        X_test.to_csv('scRNA-FeatureSelection/tempData/temp_X_test.csv')
        y_train.to_csv('scRNA-FeatureSelection/tempData/temp_y_train.csv')
        y_test.to_csv('scRNA-FeatureSelection/tempData/temp_y_test.csv')
    elif task == 'clustering' and X_train is not None and y_train is not None:
        X_train.to_csv('scRNA-FeatureSelection/tempData/temp_X.csv')
        y_train.index = X_train.index
        y_train.to_csv('scRNA-FeatureSelection/tempData/temp_y.csv')


def now():
    """
    UTC + 8 hours = Beijing time

    :return: Beijing time in string
    """
    return (datetime.datetime.now() + datetime.timedelta(hours=8)).strftime("%Y-%m-%d %H:%M:%S")


def head(name: str, head_len=50):
    """
    formatted printing

    :param name: head name
    :param head_len: length of head, including stars
    :return: None
    """
    stars_num = head_len - len(name)
    assert stars_num > 0, "name is too long."
    if stars_num % 2 == 0:
        print(''.join(['*' * (stars_num // 2), ' ', name, ' ', '*' * (stars_num // 2)]))
    else:
        print(''.join(['*' * (stars_num // 2), ' ', name, ' ', '*' * (stars_num // 2 + 1)]))


def plot_result(dataset: str, data_type: str, task: str, save=False):
    """
    Visualize results.

    :param dataset: dataset name
    :param data_type: 'norm' or 'raw'
    :param task: 'assign' or 'clustering'
    :param save: whether to save figure
    :return: None
    """
    assert (data_type in ['norm', 'raw']), "Parameter 'data_type' is wrong."
    assert (task in ['assign', 'clustering']), "Parameter 'task' is wrong."
    result_path = 'scRNA-FeatureSelection/results/plot'
    result = pd.read_csv(result_path + '_'.join([dataset, data_type, task, 'record.csv']), index_col=0)

    result[result.index.str.contains('F1') | result.index.str.contains('ARI')] *= 100

    fig, axes = plt.subplots(result.shape[0], 1, figsize=(6, result.shape[0] * 3))
    for i in range(result.shape[0]):
        if 'F1' in result.iloc[i, :].name or 'ARI' in result.iloc[i, :].name:
            axes[i].set_yscale('symlog')
            axes[i].axhline(0)
            axes[i].set_title(result.index[i] + '(%)')
        else:
            axes[i].set_title(result.index[i])
        cmap = mpl.cm.get_cmap('Set1', result.shape[1])
        colors = cmap(np.linspace(0, 1, result.shape[1]))
        result.iloc[i, :].plot(kind='bar', ax=axes[i], rot=15, color=colors)
    plt.tight_layout()

    if save:
        plt.savefig(result_path + '_'.join([dataset, data_type, task, 'graph.jpg']), dpi=150, bbox_inches='tight')
    plt.show()
