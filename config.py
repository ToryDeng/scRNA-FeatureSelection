class ExperimentConfig:
    def __init__(self):
        # unsupervised methods
        self.method_on = {
            'var': 'raw',
            'cv2': 'raw',
            'seurat': 'raw',
            'deviance': 'raw',
            'm3drop': 'raw',
            'scmap': 'raw',
            'cellranger': 'raw'
        }
        self.method_lan = {
            'var': 'python',
            'cv2': 'python',
            'seurat': 'python',
            'deviance': 'r',
            'm3drop': 'r',
            'scmap': 'r',
            'cellranger': 'python'
        }
        # measurements
        self.measurements = {
            'population_demixing': ['xin', 'muraro', 'segerstolpe', 'PBMC3000'],  # all datasets
            'marker_discovery': ['xin', 'muraro', 'segerstolpe', 'PBMC3000'],  # 4 datasets
            'selection_stability': ['xin', 'muraro', 'segerstolpe', 'PBMC3000'],  #

            'classification': ['xin', 'muraro', 'segerstolpe', 'PBMC3000'],  # singlecellnet singleR itclust, ck f1
            'clustering': ['xin', 'muraro', 'segerstolpe', 'PBMC3000'],  # seurat sc3 cidr, bcubed_f1 ARI v-measure
            'rare_cell_detection': ['xin', 'muraro', 'segerstolpe', 'PBMC3000'],  # f1_rare, bcubed_f1_rare
            'batch_correction': ['xin+muraro'],
            # differential expression analysis?

            'computation_time': ['xin', 'muraro']  # large-scale datasets

            #, 'muraro', 'segerstolpe', 'PBMC3000'
        }

        self.random_seed = 2020
        self.scale_factor = 1e4
        self.n_filter_cell = 5
        self.n_filter_gene = 5
        self.n_genes = [500, 1000, 1500, 2000]
        self.max_timeout = 60 * 30  # 60 seconds per minute * 30 minutes
        self.ensemble_mode = 'importance_sum'

        self.record_path = 'records/'
        self.figure_path = 'figures/'


class AssignConfig(ExperimentConfig):
    def __init__(self):
        super().__init__()
        # supervised methods
        self.method_on.update({
            'rf': 'raw',
            'lgb': 'raw',
            'xgb': 'raw',
            'nsc': 'norm',
            'fisher_score': 'raw',
            'scGeneFit': 'raw',
        })
        self.method_lan.update({
            'rf': 'python',
            'lgb': 'python',
            'xgb': 'python',
            'nsc': 'python',
            'fisher_score': 'python',
            'scGeneFit': 'python',
        })
        self.n_folds = 5  # 5 folds


class ClusterConfig(ExperimentConfig):
    def __init__(self):
        super().__init__()
        self.n_folds = 2
        self.n_loops = 10  # 10 loops
        self.rare_type_detection_metric = 'bcubed'  # f-measure


class DataConfig:
    def __init__(self):
        # how many datasets are in use
        self.n_datasets = 4

        # PBMC
        self.PBMC_path = "/volume/scRNA/python_data/PBMC_AllCells_withLabels.h5"
        self.PBMC3K_path = "/home/tdeng/SingleCell/data/PBMC/pbmc3k_raw.h5ad"
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
    'scmap': 'scmap',
    'rf': 'RF',
    'lgb': 'LightGBM',
    'xgb': 'XGBoost',
    'nsc': 'NSC',
    'fisher_score': 'Fisher Score',
    'scGeneFit': 'scGeneFit',
    'cellranger': 'CellRanger',
    'rf+fisher_score': 'RF+\nFisher Score',
    'seurat+deviance': 'Seurat+\nDeviance'
}

exp_cfg = ExperimentConfig()
assign_cfg, cluster_cfg, data_cfg = AssignConfig(), ClusterConfig(), DataConfig()
