import ItClust as ic
from utils.utils import HiddenPrints
import os
import warnings
import tensorflow as tf
import numpy as np
import anndata as ad


def ItClust_predict(train_data: ad.AnnData, test_data: ad.AnnData) -> np.ndarray:

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
    warnings.filterwarnings("ignore")

    del train_data.uns
    del test_data.uns

    with HiddenPrints():
        clf = ic.transfer_learning_clf()
        clf.fit(train_data, test_data, maxiter=200, pretrain_epochs=50, verbose=True)  # not print under HiddenPrints
        pred, prob, cell_type_pred = clf.predict(write=False)

    return pred.cluster.map({key: value[0] for key, value in cell_type_pred.items()}).values.astype(np.str_)

