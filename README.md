![GitHub](https://img.shields.io/github/license/ToryDeng/scRNA-FeatureSelection)
![GitHub Repo stars](https://img.shields.io/github/stars/ToryDeng/scRNA-FeatureSelection)
![GitHub repo size](https://img.shields.io/github/repo-size/ToryDeng/scRNA-FeatureSelection)
# scRNA-FeatureSelection
Evaluation of several gene selection methods (including ensemble gene selection methods).

## Program Structure
```
│  main.py
│          
├─cache
│  │  
│  ├─geneData       # store selected genes
│  └─preprocessedData  # store preprocessed datasets
│                   
├─common_utils
│      __init__.py
│      utils.py           #  common utils
│ 
├─config
│      __init__.py                  
│      datasets_config.py   
│      experiments_config.py 
│      methods_config.py
│      
├─data_loader
│      __init__.py
│      dataset.py         # load and preprocess datasets
│      utils.py           # utils used in loading and preprocessing data
│      
├─experiments
│      __init__.py
│      metrics.py         # metrics used in batch correction, cell classification and cell clustering
│      recorders.py       # record the evaluation results and sink them to disk
│      run_experiments.py # run each experiment by calling the corresponding function
│      
├─figures                 # store the umap and t-sne figures
│      
├─other_steps
│      __init__.py
│      classification.py  # cell classification algorithms
│      clustering.py      # cell clustering algorithms
│      correction.py      # batch correction algorithms
│      
├─records                 # store the evaluation results and recorders
└─selection
        __init__.py
        fisher_score.py
        methods.py        # all feature selection algorithms
        nearest_centroid.py
        utils.py          # utils used in feature selection
```

## Included Methods
| Method | Classification  | Clustering |  Language  |  Reference |
| :----: | :-------------: | :--------: | :--------: | :--------: |
| Random Forest | ✔ | ❌ | Python | [[1]](https://doi.org/10.1186/1471-2105-7-3) |
| XGBoost     | ✔ | ❌ | Python | [[2]](https://doi.org/10.1145/2939672.2939785) |
| LightGBM    | ✔ | ❌ | Python | [[3]](https://papers.nips.cc/paper/2017/hash/6449f44a102fde848669bdd9eb6b76fa-Abstract.html) |
| Variance    | ✔ | ✔ | Python | [[4]](https://doi.org/10.1038/s41586-020-2649-2) |
| CV2         | ✔ | ✔ | Python | [[4]](https://doi.org/10.1038/s41586-020-2649-2) |
| Nearest Shrunken Centroid | ✔ | ❌ | Python | [[5]](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.565.4073&rep=rep1&type=pdf) |
| Seurat       | ✔ | ✔ | Python | [[6]](https://doi.org/10.1016/j.cell.2019.05.031) |
| Deviance     | ✔ | ✔ | R | [[7]](https://doi.org/10.1186/s13059-019-1861-6) |
| M3Drop       | ✔ | ✔ | R | [[8]](https://doi.org/10.1093/bioinformatics/bty1044) |
| scmap        | ✔ | ✔ | R | [[9]](https://doi.org/10.1038/nmeth.4644) |
| scGeneFit    | ✔ | ❌ | Python | [[10]](https://doi.org/10.1038/s41467-021-21453-4) |
| CellRanger   | ✔ | ✔ |  Python | [[11]](https://doi.org/10.1038/ncomms14049) |
| Fisher Score | ✔ | ❌ | Python | [[12]](http://dx.doi.org/10.4238/gmr.15028798) |
| FEAST        | ✔ | ✔ |  R     |  [[13]](https://doi.org/10.1093/bib/bbab034) |
| Mutual Information| ✔ | ❌ | Python  | [[14]](https://doi.org/10.1016/j.neucom.2008.04.005) |
| scran        | ✔ | ✔ |  R     | [[15]](https://doi.org/10.1186/s13059-016-0947-7) |
| triku        | ✔ | ✔ | Python | [[16]](https://doi.org/10.1093/gigascience/giac017) |
| sctransform  | ✔ | ✔ |    R   | [[17]](https://doi.org/10.1186/s13059-021-02584-9) |
| GiniClust3   | ✔ | ✔ | Python |  [[18]](https://doi.org/10.1186/s12859-020-3482-1) |
| pagest       | ✔ | ✔ | Python |     |

## Quality Control
The function that detects ouliers in [Besca](https://bedapub.github.io/besca/preprocessing/besca.pp.valOutlier.html).

## Normalization
The normalization method in Seurat and the implementation in [Scanpy](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.recipe_seurat.html).

## Reproduce Our Results
Before the evaluation you should specify the paths to data (and marker genes if you want to run the marker discovery experiment) in `config/datasets_config.py`:
```python
class DatasetConfig:
    def __init__(self):
        self.data_path = "/path/to/datasets/"
        self.marker_path = "/path/to/marker/genes/"  # optional
```
Then you can run certain experiment with one line of code:
```python
from experiments.run_experiments import run_cell_clustering, run_cell_classification

run_cell_clustering(fs_methods=['var', 'feast'])  # single FS methods
run_cell_classification(fs_methods=['lgb+rf'])  # ensemble FS method
```
All the records will be stored in the directory `records/`. The recorders in `.pkl` format are in `records/pkl/`, and the tables are in `records/xlsx/`.

## Evaluating new feature selection methods step by step
Here we present an easy way to evaluate new feature selection methods on all datasets we used. if you just
want to test on only a few datasets, please check the [notebook](https://github.com/ToryDeng/scRNA-FeatureSelection/blob/main/feature_selection.ipynb) for examples.
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
    - ***input*** of your new functions: an `AnnData` object, in whcih the `AnnData.X` is the log-normalized data, 
    the `AnnData.raw` is the  data after quality control but before normalization, and the normalized data is in `adata.layers['normalized']`.
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
   

