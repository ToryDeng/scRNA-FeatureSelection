from utils.utils import get_gene_names, filter_const_genes, save_filtered_data, load_data
from utils.evaluate_result import evaluate_method
from utils.importance import select_features
import numpy as np



X_raw, X_norm, y, trusted_markers = load_data('muraro')

gene_names = get_gene_names(X_raw.columns)

result = select_features('muraro', 1000, 'cv2', all_features=gene_names, X=X_norm.values, y=np.squeeze(y.values))
save_filtered_data(X_raw, y, gene_names, result)
print(result.shape)