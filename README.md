# scRNA-FeatureSelection
Several feature selection methods for scRNA-seq analysis using python and R

## Program Structure
    utils-  # functions about data processing and several feature selection methods  
      evaluate_result.py   # evaluat the result of specific feature selection method   
      importance.py        # calculate importance of every gene and select the most important ones   
      nearest_centroid.py  # nearest shrunken centroid 
      scGeneFit.py         # scGeneFit method for feature selection
      utils.py             # data processing functions 
      gene_selection.py    # the main function 'evaluate_gene_selection_method'

## Included Methods
- Random Forest  
- XGBOOST
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
from utils.gene_selection import evaluate_gene_selection_method


evaluate_gene_selection_method(dataset='xin', methods=['rf', 'lgb', 'xgb', 'nsc', 'cv2', 'var'], data_type='raw')
evaluate_gene_selection_method(dataset='muraro', methods=['cellassign', 'deviance', 'm3drop'], data_type='norm')
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
