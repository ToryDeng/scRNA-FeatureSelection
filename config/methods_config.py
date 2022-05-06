class MethodConfig:
    def __init__(self):
        self.formal_names = {
            'var': 'Variance',
            'cv2': 'CV2',
            'seurat': 'Seurat',
            'seurat_v3': 'Seurat v3',
            'deviance': 'Deviance',
            'm3drop': 'M3Drop',
            'scmap': 'scmap',
            'rf': 'Random Forest',
            'lgb': 'LightGBM',
            'xgb': 'XGBoost',
            'nsc': 'NSC',
            'fisher_score': 'Fisher Score',
            'scGeneFit': 'scGeneFit',
            'cellranger': 'CellRanger',
            'feast': 'FEAST',
            'rf+fisher_score': 'Random Forest+\nFisher Score',
            'seurat_v3+deviance': 'Seurat v3+\nDeviance',
            'scran': 'scran',
            'geneclust': 'GeneClust',
            'mi': 'Mutual Information',
            'gestect': 'gestect'
        }
        self.unsupervised = ['var', 'cv2', 'seurat', 'seurat_v3', 'deviance', 'm3drop', 'feast', 'scmap', 'cellranger',
                             'scran', 'geneclust', 'gestect']
        self.supervised = ['rf', 'lgb', 'xgb', 'nsc', 'mi', 'fisher_score', 'scGeneFit']


method_cfg = MethodConfig()
