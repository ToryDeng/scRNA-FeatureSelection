from utils.evaluation import evaluate_assign_methods, evaluate_cluster_methods

#
assign_methods = ['rf', 'lgb', 'xgb', 'seurat', 'var', 'cv2', 'nsc', 'fisher_score', 'scGeneFit',
                  'cellassign', 'm3drop', 'deviance', 'rf+fisher_score', 'monocle3']
cluster_methods = ['cv2', 'var',  'm3drop',  'seurat', 'deviance', 'seurat+deviance', 'monocle3']

#
for dataset in ['xin', 'muraro', 'segerstolpe', 'PBMC3000']:
    evaluate_assign_methods(dataset=dataset, methods=assign_methods)
    evaluate_cluster_methods(dataset=dataset, methods=cluster_methods)
