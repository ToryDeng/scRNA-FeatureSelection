from utils.utils import get_gene_names, save_filtered_data, load_data
from utils.evaluate_result import evaluate_method
from utils.importance import select_features
import numpy as np
import os


X_raw, X_norm, y, trusted_markers = load_data('muraro')
gene_names = get_gene_names(X_raw.columns)

for method in ['rf']:
    result = select_features('muraro', 1000, method, all_features=gene_names, X=X_norm.values, y=np.squeeze(y.values))
    save_filtered_data(X_raw, y, gene_names, result)
    os.system('Rscript scRNA-FeatureSelection/utils/RCode/main.R')
    evaluate_method(trusted_markers, result, y)

