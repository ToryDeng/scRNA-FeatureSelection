class ExperimentConfig:
    def __init__(self):
        self.method_on = {
            'var': 'raw',
            'cv2': 'raw',
            'seurat': 'raw',
            'deviance': 'raw',
            'm3drop': 'raw',
            'monocle3': 'raw',
            'DE': 'norm'
        }
        self.method_lan = {
            'var': 'python',
            'cv2': 'python',
            'seurat': 'python',
            'deviance': 'r',
            'm3drop': 'r',
            'monocle3': 'r',
            'DE': 'python'
        }
        self.random_seed = 2020
        self.scale_factor = 1e4
        self.n_filter_cell = 3
        self.n_filter_gene = 3
        self.n_genes = 2000  # TODO: sim 1000 genes
        self.max_timeout = 60 * 30  # 60 seconds per minute * 30 minutes
        self.ensemble_mode = 'importance_sum'

        self.record_path = 'records/'
        self.figure_path = 'figures/'


class AssignConfig(ExperimentConfig):
    def __init__(self):
        super().__init__()
        self.method_on.update({
            'rf': 'raw',
            'lgb': 'raw',
            'xgb': 'raw',
            'nsc': 'norm',
            'fisher_score': 'raw',
            'scGeneFit': 'raw',
            'cellassign': 'raw'
        })
        self.method_lan.update({
            'rf': 'python',
            'lgb': 'python',
            'xgb': 'python',
            'nsc': 'python',
            'fisher_score': 'python',
            'scGeneFit': 'python',
            'cellassign': 'r'
        })
        self.n_folds = 5  # 5 folds


class ClusterConfig(ExperimentConfig):
    def __init__(self):
        super().__init__()
        self.n_folds = 2
        self.n_loops = 1  # 10 loops


class DataConfig:
    def __init__(self):
        # how many datasets are in use
        self.n_datasets = 4

        # PBMC
        self.PBMC_path = "/volume/scRNA/python_data/PBMC_AllCells_withLabels.h5"
        self.PBMC_markers_path = "/home/tdeng/SingleCell/data/PBMC/"

        # pancreas
        self.pancreas_path = "/volume/scRNA/python_data/pancreas.h5"
        self.pancreas_markers_path = "/home/tdeng/SingleCell/data/pancreas/"

        # sim data
        self.sim_path = "/home/tdeng/SingleCell/data/sim/sim_raw.h5"
        self.sim_markers_path = "/home/tdeng/SingleCell/data/sim/"

        # remove_types
        self.pancreas_remove_types = ['unclear', 'not applicable', 'unclassified', 'co-expression', 'beta.contaminated',
                                      'alpha.contaminated', 'delta.contaminated', 'gamma.contaminated', 'unassigned',
                                      'MHC class II', 'unclassified endocrine']


formal_method_names = {
    'var': 'Variance',
    'cv2': 'CV2',
    'seurat': 'Seurat',
    'deviance': 'Deviance',
    'm3drop': 'M3Drop',
    'monocle3': 'Monocle3',
    'DE': 'DE',
    'rf': 'RF',
    'lgb': 'LightGBM',
    'xgb': 'XGBoost',
    'nsc': 'NSC',
    'fisher_score': 'Fisher Score',
    'scGeneFit': 'scGeneFit',
    'cellassign': 'CellAssign',
    'rf+fisher_score': 'RF+\nFisher Score',
    'seurat+deviance': 'seurat+\ndeviance'
}

exp_cfg = ExperimentConfig()
assign_cfg, cluster_cfg, data_cfg = AssignConfig(), ClusterConfig(), DataConfig()
