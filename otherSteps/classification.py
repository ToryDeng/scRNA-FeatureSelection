import numpy as np
import anndata as ad
import anndata2ri
import traceback
from rpy2.robjects import r, globalenv
from rpy2.robjects.packages import importr
from utils import HiddenPrints


def classify_cells(train_adata: ad.AnnData, test_adata: ad.AnnData, method: str, col_name: str = 'assign_label'):
    print(f"{method} clustering starts. {train_adata.n_obs} cells and {train_adata.n_vars} genes in train data..."
          f"{test_adata.n_obs} cells and {test_adata.n_vars} genes in test data...")
    try:
        if method == 'SingleR':
            assign_labels = SingleR(train_adata, test_adata)
        else:
            raise NotImplementedError(f"{method} has not been implemented!")
        test_adata.obs[col_name] = assign_labels
    except:
        traceback.print_exc()



def SingleR(adata_train: ad.AnnData, adata_test: ad.AnnData):
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
