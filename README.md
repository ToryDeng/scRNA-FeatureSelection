# scRNA-FeatureSelection
Several feature selection methods for scRNA-seq analysis using python and R

## Program Structure
    utils-           # functions about data processing and several feature selection methods  
      evaluate_result.py   # evaluat the result of specific feature selection method   
      importance.py        # calculate importance of every gene and select the most important ones   
      nearest_centroid.py  # nearest shrunken centroid 
      scGeneFit.py         # scGeneFit method for feature selection
      utils.py             # data processing functions 
      gene_selection.py    # the main function 'evaluate_gene_selection_method'
    examples.ipynb-  # some examples

## Included Methods
| Method | Classification  | Clustering |
| ------------- | :-------------: | :--------: |
| Random Forest | ✅ | 🔶
- Random Forest  
- XGBoost
- LigthtGBM      
- Seurat
- Variance       
- Deviance
- M3Drop         
- CellAssign
- CV2            
- scGeneFit
- Nearest Shrunken Centroid


## Example

```python
from utils.evaluation import evaluate_classification_methods, evaluate_clustering_methods

evaluate_classification_methods(dataset='xin', methods=['rf', 'lgb', 'xgb', 'nsc', 'cv2', 'var'], data_type='raw')
evaluate_clustering_methods(dataset='muraro', methods=['cellassign', 'deviance', 'm3drop'], data_type='norm')
```
### Datasets
#### Pancreas Datasets
- *'muraro'*
- *'segerstolpe'*
- *'xin'*
#### PBMC Datasets
- *'PBMC50'*: 50% of the complete PBMC dataset
- *'PBMC20'*: 20% of the complete PBMC dataset
- *'PBMC10'*: 10% of the complete PBMC dataset
- *'PBMC5'*: 5% of the complete PBMC dataset

### Notes
1. You can choose which type of data the methods will be implemented on by specifying the parameter *'data_type'*: 'raw' means raw data, 'norm' means normalized data.
2. Parameter *'methods'* must be a list. If *'methods'* contains 'scGeneFit' and dataset is PBMC, then all the methods will only run on PBMC5. This is because the scGeneFit written by Python is very slow and occupies much system resource when the number of cells is more than 3000.
3. CV2 may select 0 marker gene and cause IO error.
