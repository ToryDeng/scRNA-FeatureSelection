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
- Square of Coefficient of variation
- Nearest Shrunken Centroid
- Deviance
- M3Drop
- CellAssign
- scGeneFit

## Example
```python
from utils.gene_selection import evaluate_gene_selection_method


evaluate_gene_selection_method(dataset='xin', methods=['rf', 'lgb', 'xgb', 'nsc', 'cv2', 'var'], data_type='raw')
evaluate_gene_selection_method(dataset='muraro', methods=['cellassign', 'deviance', 'm3drop'], data_type='norm')
```
