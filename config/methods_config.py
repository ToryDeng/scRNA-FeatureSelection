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
            'rf': 'RF',
            'lgb': 'LightGBM',
            'xgb': 'XGBoost',
            'nsc': 'NSC',
            'fisher_score': 'Fisher Score',
            'scGeneFit': 'scGeneFit',
            'cellranger': 'CellRanger',
            'rf+fisher_score': 'RF+\nFisher Score',
            'seurat_v3+deviance': 'Seurat v3+\nDeviance'
        }
        self.unsupervised = ['var', 'cv2', 'seurat', 'seurat_v3', 'deviance', 'm3drop', 'scmap', 'cellranger']
        self.supervised = ['rf', 'lgb', 'xgb', 'nsc', 'fisher_score', 'scGeneFit']


method_cfg = MethodConfig()
