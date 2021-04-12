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
| Method | Classification  | Clustering |  Language  |  Reference |
| :----: | :-------------: | :--------: | :--------: | :--------: |
| Random Forest | ✔ | ❌ | Python | - |
| XGBoost    | ✔ | ❌ | Python | - |
| LigthtGBM   | ✔ | ❌ | Python | - |
| Variance    | ✔ | ✔ | Python | - |
| CV2         | ✔ | ✔ | Python | - |
| Nearest Shrunken Centroid | ✔ | ❌ | Python | [1] |
| Seurat       | ✔ | ✔ | Python | [2] |
| Deviance     | ✔ | ✔ | R | [3] |
| M3Drop       | ✔ | ✔ | R | [4] |
| CellAssign   | ✔ | ❌ | R | [5] |
| scGeneFit    | ✔ | ❌ | Python | [6] |
| Monocle3     | ✔ | ✔ |  R     | [7] |
| Fisher Score | ✔ | ❌ | Python | [8] |

>1. Klassen M, Kim N. Nearest Shrunken Centroid as Feature Selection of Microarray Data[C]//CATA. 2009: 227-232.
>2. Stuart T, Butler A, Hoffman P, et al. Comprehensive integration of single-cell data[J]. Cell, 2019, 177(7): 1888-1902. e21.  
>3. Townes F W, Hicks S C, Aryee M J, et al. Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model[J]. Genome biology, 2019, 20(1): 1-16.  
>4. Andrews T S, Hemberg M. M3Drop: dropout-based feature selection for scRNASeq[J]. Bioinformatics, 2019, 35(16): 2865-2867.  
>5. Zhang A W, O’Flanagan C, Chavez E A, et al. Probabilistic cell-type assignment of single-cell RNA-seq for tumor microenvironment profiling[J]. Nature methods, 2019, 16(10): 1007-1015.  
>6. Dumitrascu B, Villar S, Mixon D G, et al. Optimal marker gene selection for cell type discrimination in single cell analyses[J]. Nature communications, 2021, 12(1): 1-8.  
>7. Cao J, Spielmann M, Qiu X, et al. The single-cell transcriptional landscape of mammalian organogenesis[J]. Nature, 2019, 566(7745): 496-502.
>8. Li J, Cheng K, Wang S, et al. Feature selection: A data perspective[J]. ACM Computing Surveys (CSUR), 2017, 50(6): 1-45.

## Datasets
- **PBMC**  
  A certain proportion of samples are extracted from the data because it is large. 
  Different sample rates(5%, 10%, 20% and 50%) are used to ensure robustness of the evaluation.
- **pancreas**  
  Three datasets which contain some cells with unclear cell types (unclear, not applicable, unclassified,
  co-expression) or contamination(alpha.contaminated, delta.contaminated, gamma.contaminated, beta.contaminated).

## Marker Genes
| Dataset | All Marker Genes | Marker Genes in Data |
| :-----: | :-----------: | :-----------: | 
|'PBMC50', 'PBMC20', 'PBMC10', 'PBMC5'|101|42
| 'muraro'    | 365 | 313 |
|'segerstolpe'| 365 | 320 |
|'xin'        | 365 | 321 |

## Normalization
The normalization method in **Seurat**.

## Metrics
- **Number of Marker Genes Found**  
  The number of marker genes which are found in data.
- **Mean Reciprocal Rank (MRR)**  
  ![1](https://latex.codecogs.com/gif.latex?MRR=\frac{1}{\vert&space;Q&space;\vert}\sum_{i=1}^{\vert&space;Q&space;\vert}\frac{1}{rank_{i}})
- **Adjusted Rand Index (ARI)**  
  ARI increment of 2 clustering methods: Seurat and SC3.
- **F1-score**  
  F1-score increment of 2 classification methods: scmap and singlecellnet (using 5-fold CV).
## Examples
### Evaluation
```python
from utils.evaluation import evaluate_classification_methods, evaluate_clustering_methods

evaluate_classification_methods(dataset='xin', methods=['rf', 'xgb', 'nsc', 'var'], data_type='raw')
evaluate_clustering_methods(dataset='muraro', methods=['cellassign', 'deviance'], data_type='norm')
```
### Ensemble Gene Selection
```python
from utils.ensemble_gene_selection import ensemble_gene_selection


ensemble_gene_selection(dataset='xin', base_methods=['rf', 'fisher_score'], task='assign')
```


## Notes
1. The ***data*** directory is at the same level as this program.
2. You can choose which type of data the methods will be implemented on by specifying the parameter ***'data_type'***: 
   'raw' means raw data, 'norm' means normalized data.
3. Parameter ***'methods'*** must be a list. If ***'methods'*** contains 'scGeneFit' and dataset is PBMC, then all the 
   methods will only run on PBMC5. This is because the scGeneFit written by Python is very slow and occupies much system 
   resource when the number of cells is more than 3000.
4. Results will be stored in the directory ***results***.
5. CV2 may select 0 marker gene and cause IO error.
6. The parameter ***'n_features'*** is fixed at 1000.
