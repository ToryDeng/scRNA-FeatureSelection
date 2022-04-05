# scRNA-FeatureSelection
Evaluation of several gene selection methods (including ensemble gene selection methods).

## Program Structure
    selection-       # some feature selection methods in python
      fisher_score.py      # calculate fisher score
      nearest_centroid.py  # nearest shrunken centroid 
      scgenefit.py         # scGeneFit method for feature selection
    utils-           # functions about data processing and several feature selection methods
      RCode-               # R scripts
      classification.py    # classification algorithms in Python
      evaluation.py        # evaluat specific feature selection methods 
      importance.py        # calculate importance of every gene and select the most important ones   
      record.py            # record performance during evaluation
      utils.py             # data processing functions
    tempData-        # store temporary data
    records-         # store records of evaluation
    main.py          # the main function
    examples.ipynb-  # some examples about usage
    config.py        # experiment and data configuration

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
| scmap        | ✔ | ✔ | R | [5] |
| scGeneFit    | ✔ | ❌ | Python | [6] |
| CellRanger   | ✔ | ✔ |  Python | [7] |
| Fisher Score | ✔ | ❌ | Python | [8] |
| FEAST        | ✔ | ✔ |  R     |  [9] |

>1. Klassen M, Kim N. Nearest Shrunken Centroid as Feature Selection of Microarray Data[C]//CATA. 2009: 227-232.
>2. Stuart T, Butler A, Hoffman P, et al. Comprehensive integration of single-cell data[J]. Cell, 2019, 177(7): 1888-1902. e21.  
>3. Townes F W, Hicks S C, Aryee M J, et al. Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model[J]. Genome biology, 2019, 20(1): 1-16.  
>4. Andrews T S, Hemberg M. M3Drop: dropout-based feature selection for scRNASeq[J]. Bioinformatics, 2019, 35(16): 2865-2867.  
>5. Kiselev V Y, Yiu A, Hemberg M. scmap: projection of single-cell RNA-seq data across data sets[J]. Nature methods, 2018, 15(5): 359-362.
>6. Dumitrascu B, Villar S, Mixon D G, et al. Optimal marker gene selection for cell type discrimination in single cell analyses[J]. Nature communications, 2021, 12(1): 1-8.  
>7. Zheng G X Y, Terry J M, Belgrader P, et al. Massively parallel digital transcriptional profiling of single cells[J]. Nature communications, 2017, 8(1): 1-12.
>8. Li J, Cheng K, Wang S, et al. Feature selection: A data perspective[J]. ACM Computing Surveys (CSUR), 2017, 50(6): 1-45.
>9. Su K, Yu T, Wu H. Accurate feature selection improves single-cell RNA-seq cell clustering[J]. Briefings in Bioinformatics, 2021.


## Normalization
The normalization method in **Seurat**.


## Examples
### Evaluation of Single Gene Selection Method
```python

```
### Evaluation of Ensemble Gene Selection Method
```python

```
All the records will be stored in the directory ***records/***.


## Notes

By specifying data type in **config.py**, you can choose which type of data the methods will use.


