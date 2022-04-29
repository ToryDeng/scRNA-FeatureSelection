![GitHub](https://img.shields.io/github/license/ToryDeng/scRNA-FeatureSelection)
![GitHub Repo stars](https://img.shields.io/github/stars/ToryDeng/scRNA-FeatureSelection)
![GitHub repo size](https://img.shields.io/github/repo-size/ToryDeng/scRNA-FeatureSelection)
# scRNA-FeatureSelection
Evaluation of several gene selection methods (including ensemble gene selection methods).

## Program structure
```
│  main.py
│          
├─cache
│  │  
│  ├─geneData       # store selected genes
│  └─preprocessedData  # store preprocessed datasets
│                   
├─common_utils
│      utils.py           #  common utils
│ 
├─config                  
│      datasets_config.py   
│      experiments_config.py 
│      methods_config.py
│      
├─data_loader
│      dataset.py         # load and preprocess datasets
│      utils.py           # utils used in loading and preprocessing data
│      
├─experiments
│      metrics.py         # metrics used in batch correction, cell classification and cell clustering
│      recorders.py       # record the evaluation results and sink them to disk
│      run_experiments.py # run each experiment by calling the corresponding function
│      
├─figures                 # store the umap and t-sne figures
│      
├─other_steps
│      classification.py  # cell classification algorithms
│      clustering.py      # cell clustering algorithms
│      correction.py      # batch correction algorithms
│      
├─records                 # store the evaluation results and recorders
└─selection
        fisher_score.py
        methods.py        # all feature selection algorithms
        nearest_centroid.py
        scgenefit.py
        utils.py          # utils used in feature selection
```

## Included methods
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
| scmap        | ✔ | ✔ | R | [5] |
| scGeneFit    | ✔ | ❌ | Python | [6] |
| CellRanger   | ✔ | ✔ |  Python | [7] |
| Fisher Score | ✔ | ❌ | Python | [8] |
| FEAST        | ✔ | ✔ |  R     |  [9] |



## Quality control
The function that detects ouliers in [Besca](https://bedapub.github.io/besca/preprocessing/besca.pp.valOutlier.html).

## Normalization
The normalization method in Seurat and the implementation in [Scanpy](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.recipe_seurat.html).


## Examples
### Specify the paths to data and marker genes
Before the evaluation you should specify the paths to data and marker genes in `config/datasets_config.py`:
```python
self.data_path = "/your/path/to/datasets/"
self.marker_path = "/your/path/to/marker/genes/"
```

### Evaluation of single gene selection method
```python
from experiments.run_experiments import run_cell_clustering

run_cell_clustering(fs_methods=['var', 'feast'])
```
### Evaluation of ensemble gene selection method
```python
from experiments.run_experiments import run_cell_classification

run_cell_classification(fs_methods=['lgb+rf'])
```
All the records will be stored in the directory `records/`. The recorders are in `records/pkl/`, and the tables are in `records/xlsx/`.

## Evaluating new feature selection methods step by step
Here we present an easy way to evaluate new feature selection methods on all datasets we used. if you just
want to test on a few datasets, please check the [notebook](https://github.com/ToryDeng/scRNA-FeatureSelection/blob/main/feature_selection.ipynb) for examples.
1. Add new methods to the function `single_select_by_batch()` in `selection/methods.py`:
   ```python
   elif method == 'deviance':
       selected_genes_df = deviance_compute_importance(adata)
   elif method == 'abbreviation_1':
       selected_genes_df = your_new_fucntion_1(adata)
   elif method == 'abbreviation_2':
       selected_genes_df = your_new_fucntion_2(adata)
   else:
       raise NotImplementedError(f"No implementation of {method}!")
   ```
    - ***input*** of your new functions: an `anndata` object, in whcih the `anndata.X` is the log-normalized data, 
    the `anndata.raw` is the  data after quality control but before normalization, and the normalized data is in `adata.layers['normalized']`.
    - ***output*** of your new functions: a dataframe. The first column with name `Gene` contains gene names. The second column
    is not necessary. It contains scores of each genes (if they exist). The higher the score is, the more important the gene.
    
2. Modify the method configuration `config/methods_config.py`:
    - in `self.formal_names`
    ```python
    'feast': 'FEAST',
    'abbreviation_1': 'formal_name_1',
    'abbreviation_2': 'formal_name_2',
    'rf+fisher_score': 'RF+\nFisher Score',
    ```
    - unsupervised methods should be added in `self.unsupervised`, and supervised methods should be added in `self.supervised`
    ```python
    self.unsupervised = ['abbreviation_1', 'var', 'cv2', ...]
    self.supervised = ['abbreviation_2', 'rf', 'lgb', 'xgb', ...]
    ```
3. Then you can run the function as shown in examples!
    ```python
    from experiments.run_experiments import run_cell_clustering

    run_cell_clustering(fs_methods=['abbreviation_1', 'abbreviation_2'])
    ```
   

