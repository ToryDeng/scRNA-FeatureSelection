from utils.utils import get_gene_names, filter_const_genes, save_filtered_data, load_data
from utils.evaluate_result import evaluate_method
from utils.importance import select_features
import numpy as np
import os


X_raw, X_norm, y, trusted_markers = load_data('PBMC10')
gene_names = get_gene_names(X_raw.columns)

for method in ['scGeneFit']:
    result = select_features('PBMC10', 1000, method, all_features=gene_names, X=X_norm.values, y=np.squeeze(y.values))
    save_filtered_data(X_raw, y, gene_names, result)
    print(result.shape)


# os.system("Rscript /home/tdeng/SingleCell/scRNA-FeatureSelection/utils/RCode/test.R 111")