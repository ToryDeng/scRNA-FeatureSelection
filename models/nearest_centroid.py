import numpy as np
from sklearn.preprocessing import LabelEncoder


def nearest_centroid_select(X, y, shrink_threshold):
    """
    Select genes using nearest shrunken centroid.

    :param X: count matrix
    :param y: labels
    :param shrink_threshold: hyperparameter shrink threshold
    :return: deviation of each gene
    """
    n_samples, n_features = X.shape
    le = LabelEncoder()
    y_ind = le.fit_transform(y)
    classes = le.classes_
    n_classes = classes.size
    centroids_ = np.empty((n_classes, n_features), dtype=np.float64)
    nk = np.zeros(n_classes)
    for cur_class in range(n_classes):
        center_mask = y_ind == cur_class
        nk[cur_class] = np.sum(center_mask)
    if n_classes < 2:
        raise ValueError('The number of classes has to be greater than'
                             ' one; got %d class' % (n_classes))
    if shrink_threshold is not None:
        if np.all(np.ptp(X, axis=0) == 0):
            raise ValueError("All features have zero variance. "
                                    "Division by zero.")
        dataset_centroid_ = np.mean(X, axis=0)

        # m parameter for determining deviation
        m = np.sqrt((1. / nk) - (1. / n_samples))
        # Calculate deviation using the standard deviation of centroids.
        variance = (X - centroids_[y_ind]) ** 2
        variance = variance.sum(axis=0)
        s = np.sqrt(variance / (n_samples - n_classes))
        s += np.median(s)  # To deter outliers from affecting the results.
        mm = m.reshape(len(m), 1)  # Reshape to allow broadcasting.
        ms = mm * s
        deviation = ((centroids_ - dataset_centroid_) / ms)
        # Soft thresholding: if the deviation crosses 0 during shrinking,
        # it becomes zero.
        signs = np.sign(deviation)
        deviation = (np.abs(deviation) - shrink_threshold)
        np.clip(deviation, 0, None, out=deviation)
        deviation *= signs
        return np.var(deviation, axis=0)