import ItClust as ic
from utils.utils import HiddenPrints
import os
import warnings
import tensorflow as tf
import numpy as np
import anndata as ad
import scvi


# these two algorithms are not used
def ItClust_predict(train_data: ad.AnnData, test_data: ad.AnnData) -> np.ndarray:

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
    warnings.filterwarnings("ignore")

    del train_data.uns
    del test_data.uns

    with HiddenPrints():
        clf = ic.transfer_learning_clf()
        clf.fit(train_data, test_data, maxiter=100, pretrain_epochs=20, verbose=True)  # not print under HiddenPrints
        pred, prob, cell_type_pred = clf.predict(write=False)

    return pred.cluster.map({key: value[0] for key, value in cell_type_pred.items()}).values.astype(np.str_)


def scANVI_predict(train_data: ad.AnnData, test_data: ad.AnnData) -> np.ndarray:
    with HiddenPrints():
        test_data.obs['celltype'] = 'unknown'
        concat_data = ad.concat([train_data, test_data])
        is_test = concat_data.obs['celltype'] == 'unknown'
        scvi.model.SCANVI.setup_anndata(concat_data, batch_key=None, labels_key='celltype')

        scvi_model = scvi.model.SCANVI(concat_data, unlabeled_category='unknown')
        scvi_model.train(max_epochs=200, early_stopping=True, batch_size=512)
        y_pred = scvi_model.predict()
    return y_pred[is_test]

