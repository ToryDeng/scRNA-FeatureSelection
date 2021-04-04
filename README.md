# scRNA-FeatureSelection
1. Evaluation of several gene selection methods
2. Ensemble gene selection

## Program Structure
    utils-           # functions about data processing and several feature selection methods
      RCode-               # R scripts
      evaluation.py        # evaluat specific feature selection methods 
      importance.py        # calculate importance of every gene and select the most important ones   
      nearest_centroid.py  # nearest shrunken centroid 
      scGeneFit.py         # scGeneFit method for feature selection
      utils.py             # data processing functions 
      ensemble_gene_selection.py    # the ensemble gene selection function
    tempData-        # store temporary data
    results-         # store results of evaluation
    examples.ipynb-  # some examples about usage

## Included Methods
| Method | Classification  | Clustering |
| :----: | :-------------: | :--------: |
| Random Forest | ✔ | ❌ |
| XGBoost    | ✔ | ❌ |
| LigthtGBM   | ✔ | ❌ |
| Variance    | ✔ | ✔ |
| CV2         | ✔ | ✔ |
| Nearest Shrunken Centroid | ✔ | ❌ |
| Seurat      | ✔ | ✔ |
| Deviance | ✔ | ✔ |
| M3Drop          | ✔ | ✔ |
| CellAssign | ✔ | ❌ |
| scGeneFit | ✔ | ❌ |



## Example

```python
from utils.evaluation import evaluate_classification_methods, evaluate_clustering_methods

evaluate_classification_methods(dataset='xin', methods=['rf', 'lgb', 'xgb', 'nsc', 'cv2', 'var'], data_type='raw')
evaluate_clustering_methods(dataset='muraro', methods=['cellassign', 'deviance', 'm3drop'], data_type='norm')
```
## Datasets
### Pancreas Datasets
- *'muraro'*
- *'segerstolpe'*
- *'xin'*
### PBMC Datasets
- *'PBMC50'*: 50% of the complete PBMC dataset
- *'PBMC20'*: 20% of the complete PBMC dataset
- *'PBMC10'*: 10% of the complete PBMC dataset
- *'PBMC5'*: 5% of the complete PBMC dataset

### Notes
1. You can choose which type of data the methods will be implemented on by specifying the parameter *'data_type'*: 'raw' means raw data, 'norm' means normalized data.
2. Parameter *'methods'* must be a list. If *'methods'* contains 'scGeneFit' and dataset is PBMC, then all the methods will only run on PBMC5. This is because the scGeneFit written by Python is very slow and occupies much system resource when the number of cells is more than 3000.
3. CV2 may select 0 marker gene and cause IO error.
