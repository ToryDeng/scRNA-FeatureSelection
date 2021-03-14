# scRNA-FeatureSelection
Several feature selection methods for scRNA-seq analysis using python and R

## Program Structure
    utils-  # functions about data processing and several feature selection methods  
      evaluate_result.py   # evaluat the result of specific feature selection method   
      importance.py        # calculate importance of every gene and select the most important ones   
      nearest_centroid.py  # nearest shrunken centroid 
      scGeneFit.py         # scGeneFit method for feature selection <br>  
      utils.py             # data processing functions 
      gene_selection.py    # the main function 'evaluate_gene_selection_method'
