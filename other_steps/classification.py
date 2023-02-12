import datetime

import anndata as ad
import anndata2ri
import numpy as np
from rpy2.robjects import r, globalenv, pandas2ri
from rpy2.robjects.packages import importr

from common_utils.utils import HiddenPrints
from config import assign_cfg
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV


def classify_cells(train_adata: ad.AnnData, test_adata: ad.AnnData):
    """
    train a classifier on training set and predict labels on test set.

    Parameters
    ----------
    train_adata
      the training set
    test_adata
      the test set

    Returns
    -------
    None
    """
    for method in assign_cfg.methods:
        print(f"{datetime.datetime.now()}: {method} starts. {train_adata.n_obs} cells and {train_adata.n_vars} genes in train data; "
              f"{test_adata.n_obs} cells and {test_adata.n_vars} genes in test data...")
        if method == 'SingleR':
            assign_labels = SingleR(train_adata, test_adata)
        elif method == 'SVM':
            assign_labels = SVM_rejection(train_adata, test_adata)
        elif method == 'scmap':
            assign_labels = scmap_cluster(train_adata, test_adata)
        else:
            raise NotImplementedError(f"{method} has not been implemented!")
        test_adata.obs[f'{method}_label'] = assign_labels


def SingleR(adata_train: ad.AnnData, adata_test: ad.AnnData) -> np.ndarray:
    """
    Using raw data to train a classifier on training dataset and predict on test dataset

    Parameters
    ----------
    adata_train
      training dataset
    adata_test
      test dataset
    Returns
    -------
    label_pred
      the predicted labels of test samples
    """
    assert adata_train.raw is not None and adata_test.raw is not None, "datasets do not contain raw counts!"
    with HiddenPrints():
        anndata2ri.activate()
        importr('SingleR')
        importr('scater')
        importr('BiocParallel')
        importr('doParallel')
        globalenv['train_data'] = anndata2ri.py2rpy(adata_train.raw.to_adata())
        globalenv['test_data'] = anndata2ri.py2rpy(adata_test.raw.to_adata())
        r("""
        BPPARAM <- MulticoreParam(workers = detectCores(logical = F) - 1)
        labels <- train_data$celltype
    
        train_data <- scater::logNormCounts(train_data, assay.type = 'X') 
        test_data <- scater::logNormCounts(test_data, assay.type = 'X')
    
        # disable feature selection
        preds <- SingleR(test=test_data, ref=train_data, labels=labels, de.n = ncol(train_data), BPPARAM = BPPARAM)
        print(preds$labels)
        """)  # print() is necessary, but don't know why
        label_pred = np.array(r('preds$labels'))
        anndata2ri.deactivate()
    return label_pred


def SVM_rejection(adata_train: ad.AnnData, adata_test: ad.AnnData) -> np.ndarray:
    classifier = LinearSVC(random_state=assign_cfg.random_seed)
    clf = CalibratedClassifierCV(classifier, n_jobs=-1)
    clf.fit(adata_train.layers['log-normalized'], adata_train.obs['celltype'])
    label_pred = clf.predict(adata_test.X)
    prob = np.max(clf.predict_proba(adata_test.layers['log-normalized']), axis=1)
    unlabeled = np.where(prob < 0.7)
    label_pred[unlabeled] = 'Unknown'
    return label_pred


def scmap_cluster(adata_train: ad.AnnData, adata_test: ad.AnnData) -> np.ndarray:
    anndata2ri.activate()
    importr('scmap')
    importr('SingleCellExperiment')
    globalenv['train_data'] = pandas2ri.py2rpy(adata_train.to_df(layer='log-normalized').T)
    globalenv['test_data'] = pandas2ri.py2rpy(adata_test.to_df(layer='log-normalized').T)
    globalenv['ann'] = pandas2ri.py2rpy(adata_train.obs['celltype'].to_frame())
    r("""
    train_sce <- SingleCellExperiment(assays = list(logcounts = train_data), colData = ann)
    test_sce <- SingleCellExperiment(assays = list(logcounts = test_data))

    rowData(train_sce)$feature_symbol <- rownames(train_sce)
    rowData(test_sce)$feature_symbol <- rownames(test_sce)
    train_sce <- setFeatures(train_sce, rownames(train_sce))
    test_sce <- setFeatures(test_sce, rownames(test_sce))

    train_sce <- indexCluster(train_sce, cluster_col = 'celltype')
    res <- scmapCluster(test_sce, list(metadata(train_sce)$scmap_cluster_index), threshold = 0.5)
    """)
    label_pred = np.array(r('res[["combined_labs"]]'))
    anndata2ri.deactivate()
    return label_pred
