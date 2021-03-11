import numpy as np
from sklearn.preprocessing import LabelEncoder
# tree models
from sklearn.ensemble import RandomForestClassifier
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier

from loess.loess_1d import loess_1d


def most_important_genes(importances, feature_num, all_features):
    if importances.shape[0] != all_features.shape[0]:
        print('The length of importance and gene names are not the same! Please check again!')
        return None
    else:
        return all_features[np.argsort(importances)[::-1][:feature_num]]


def select_features(data_name, feature_num, method, all_features, X, y):
    print('the method you are using is {}.'.format(method))
    if method == 'rf':  # random forest
        forest = RandomForestClassifier(n_estimators=100, n_jobs=15, random_state=0, verbose=0).fit(X, y)
        importances = forest.feature_importances_
        result = most_important_genes(importances, feature_num, all_features)
        return result
    elif method == 'lgb':
        lgb = LGBMClassifier(n_jobs=15, random_state=0).fit(X, y)
        importances = lgb.feature_importances_
        result = most_important_genes(importances, feature_num, all_features)
        return result
    elif method == 'xgb':
        le = LabelEncoder()
        y = le.fit_transform(y)
        xgb = XGBClassifier(eval_metric='mlogloss', use_label_encoder=False, nthread=15).fit(X, y)
        importances = xgb.feature_importances_
        result = most_important_genes(importances, feature_num, all_features)
        return result
    elif method == 'seurat':
        mean, var = X.mean(axis=0), X.var(axis=0)
        mean_fit, var_fit, weigts = loess_1d(np.log10(mean), np.log10(var), frac=0.3, degree=2)

        z = (X - mean) / (10 ** (var_fit / 2))
        z[z > X.shape[0] ** 0.5] = X.shape[0] ** 0.5
        z = np.var(np.squeeze(z), axis=0)

        result = most_important_genes(z, feature_num, all_features)
        return result
    elif method == 'var':
        var = np.var(X, axis=0)
        result = most_important_genes(var, feature_num, all_features)
        return result
    elif method == 'cv2':
        cv = np.squeeze(np.std(X, axis=0) / np.squeeze(np.mean(X, axis=0)))
        print('num of genes whose mean values are less than 1e-4:{}'.format(
            np.sum(np.squeeze(np.mean(X, axis=0)) < 1e-4)))
        result = most_important_genes(cv ** 2, feature_num, all_features)
        return result

