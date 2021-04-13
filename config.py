class ClassificationConfig:
    def __init__(self):
        self.method_on = {
            'rf': 'raw',
            'lgb': 'raw',
            'xgb': 'raw',
            'var': 'raw',
            'cv2': 'raw',
            'nsc': 'raw',
            'seurat': 'raw',
            'deviance': 'raw',
            'm3drop': 'raw',
            'cellassign': 'raw',
            'monocle3': 'raw',
            'fisher_score': 'raw',
            'scGeneFit': 'raw'
        }
        self.method_lan = {
            'rf': 'python',
            'lgb': 'python',
            'xgb': 'python',
            'var': 'python',
            'cv2': 'python',
            'nsc': 'python',
            'seurat': 'r',
            'deviance': 'r',
            'm3drop': 'r',
            'cellassign': 'r',
            'monocle3': 'r',
            'fisher_score': 'python',
            'scGeneFit': 'python'
        }


class ClusteringConfig:
    def __init__(self):
        self.method_on = {
            'var': 'raw',
            'cv2': 'raw',
            'seurat': 'raw',
            'deviance': 'raw',
            'm3drop': 'raw',
            'monocle3': 'raw'
        }
        self.method_lan = {
            'var': 'python',
            'cv2': 'python',
            'seurat': 'r',
            'deviance': 'r',
            'm3drop': 'r',
            'monocle3': 'r'
        }


class DataConfig:
    def __init__(self):
        # PBMC
        self.PBMC_path = "/volume/scRNA/python_data/PBMC_AllCells_withLabels.h5"

        # pancreas
        self.pancreas_path = "/volume/scRNA/python_data/pancreas.h5"
        self.remove_types = ['unclear', 'not applicable', 'unclassified', 'co-expression', 'beta.contaminated',
                             'alpha.contaminated', 'delta.contaminated', 'gamma.contaminated', 'unassigned']


classification_cfg, clustering_cfg, data_cfg = ClassificationConfig(), ClusteringConfig(), DataConfig()
