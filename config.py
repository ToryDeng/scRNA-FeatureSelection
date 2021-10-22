import numpy as np


class ExperimentConfig:
    def __init__(self):
        # unsupervised methods
        self.method_on = {
            'var': 'raw',
            'cv2': 'raw',
            'seurat': 'raw',
            'feast': 'norm',
            'deviance': 'raw',
            'm3drop': 'raw',
            'scmap': 'raw',
            'cellranger': 'raw'
        }
        self.method_lan = {
            'var': 'python',
            'cv2': 'python',
            'seurat': 'python',
            'feast': 'r',
            'deviance': 'r',
            'm3drop': 'r',
            'scmap': 'r',
            'cellranger': 'python'
        }
        # measurements
        self.measurements = {
            'population_demixing': ['baron', 'segerstolpe', 'ZilionisMouseLungCancer', 'MarquesMouseBrain',
                                    'GuoHumanTestis', 'QuakeMouseHeart', 'ZeiselMouseBrain', 'BaronMousePancreas',
                                    'LaMannoMouseAdult', 'LaMannoHumanEmbryonicStem', 'LaMannoHumanEmbryonicMidbrain',
                                    'DengMouseEmbryoDevel', 'GoolamMouseEmbryoDevel', 'PBMCsmall'],  # 13 datasets
            #

            'marker_discovery': ['baron', 'simulatingPBMCsmall'],  # 4 datasets  , 'segerstolpe', 'PBMCsmall'
            'selection_stability': ['baron', 'segerstolpe', 'ZilionisMouseLungCancer', 'MarquesMouseBrain',
                                    'GuoHumanTestis', 'QuakeMouseHeart', 'ZeiselMouseBrain', 'BaronMousePancreas',
                                    'LaMannoMouseAdult', 'LaMannoHumanEmbryonicStem', 'LaMannoHumanEmbryonicMidbrain',
                                    'DengMouseEmbryoDevel', 'GoolamMouseEmbryoDevel', 'PBMCsmall'],  #

            'classification': ['BaronMousePancreas'],
            # singlecellnet singleR itclust, ck f1
            #
            #
            #
            #
            #
            'clustering': ['baron', 'segerstolpe'],  # seurat sc3 cidr, bcubed_f1 ARI v-measure
            'rare_cell_detection': ['baron', 'segerstolpe'],  # f1_rare, bcubed_f1_rare
            'batch_correction': ['baron+segerstolpe', 'AztekinTail'],  #
            # differential expression analysis?

            'computation_time': ['PBMC50000',  # large-scale datasets, only use 50000 samples
                                 'VentoHumanPlacenta200', 'VentoHumanPlacenta500', 'VentoHumanPlacenta1000',
                                 'VentoHumanPlacenta2000', 'VentoHumanPlacenta5000', 'VentoHumanPlacenta10000',
                                 'VentoHumanPlacenta20000', 'VentoHumanPlacenta50000']

        }

        self.random_seed = 2021
        # normalization
        self.scale_factor = 1e4
        # filtering
        self.n_filter_cell = 10
        self.n_filter_gene = 10
        self.n_genes = [500, 1000, 1500, 2000]  #
        self.max_timeout = 60 * 60 * 3  # 60 seconds per minute * 300 minutes
        self.ensemble_mode = 'importance_sum'
        # batch-effect removal: scGen
        self.epochs = 500
        self.batch_size = 128
        self.use_early_stopping = True
        self.patience = 20
        # saving directory
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

        # data path
        self.data_path = "/volume/scRNA/python_data/"

        # PBMC
        self.PBMC_markers_path = "/home/tdeng/SingleCell/data/PBMC/"

        # pancreas
        self.pancreas_markers_path = "/home/tdeng/SingleCell/data/pancreas/"

        # mouse brain
        self.mouse_brain_markers_path = "/home/tdeng/SingleCell/data/MouseBrain"

        # sim data
        self.sim_markers_path = "/home/tdeng/SingleCell/data/sim/"

        # remove_types
        self.remove_types = ['unclear', 'not applicable', 'unclassified', 'co-expression', 'beta.contaminated',
                             'alpha.contaminated', 'delta.contaminated', 'gamma.contaminated', 'unassigned',
                             'MHC class II', 'unclassified endocrine', 'PSC', 'mast', np.NaN]


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
