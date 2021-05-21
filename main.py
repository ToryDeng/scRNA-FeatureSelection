from utils.evaluation import evaluate_assign_methods, evaluate_cluster_methods

#
assign_methods = ['rf', 'lgb', 'xgb', 'seurat', 'var', 'cv2', 'nsc', 'fisher_score',
                  'cellassign', 'm3drop', 'deviance', 'monocle3', 'rf+fisher_score']
cluster_methods = ['cv2', 'var',  'm3drop',  'seurat', 'deviance', 'monocle3', 'seurat+deviance']
#
for dataset in ['segerstolpe', 'xin', 'muraro', 'PBMC10000']:
    evaluate_assign_methods(dataset=dataset, methods=assign_methods)
    evaluate_cluster_methods(dataset=dataset, methods=cluster_methods)
