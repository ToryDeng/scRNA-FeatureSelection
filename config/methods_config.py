class MethodConfig:
    def __init__(self):
        self.formal_names = {
            'var': 'Variance',
            'cv': 'CV',
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
            'mi': 'Mutual Information',
            'triku': 'triku',
            'sct': 'SCTransform',
            'giniclust3': 'GiniClust3',
            'sc3': 'SC3',
            'pagest1w': 'pagest_oneway',
            'pagest2w': 'pagest_twoway'
        }
        self.unsupervised = ['var', 'cv', 'seurat', 'seurat_v3', 'deviance', 'm3drop', 'feast', 'scmap', 'cellranger',
                             'scran', 'triku', 'sct', 'giniclust3', 'sc3', 'pagest1w', 'pagest2w']
        self.supervised = ['rf', 'lgb', 'xgb', 'nsc', 'mi', 'fisher_score', 'scGeneFit']


method_cfg = MethodConfig()
