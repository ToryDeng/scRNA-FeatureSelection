import numpy as np


class ExperimentConfig:
    def __init__(self):
        # unsupervised methods
        self.method_on = {
            'var': 'raw',
            'cv2': 'norm',
            'random': 'raw',
            'seurat': 'norm',
            'seurat_v3': 'raw',
            'feast': 'raw',  # normalize in R
            'deviance': 'raw',  # normalize in R
            'm3drop': 'raw',  # normalize in R
            'scmap': 'raw',  # normalize in R
            'cellranger': 'raw'
        }
        self.method_lan = {
            'var': 'python',
            'cv2': 'python',
            'random': 'python',
            'seurat': 'python',
            'seurat_v3': 'python',
            'feast': 'r',
            'deviance': 'r',
            'm3drop': 'r',
            'scmap': 'r',
            'cellranger': 'python'
        }
        # measurements
        self.measurements = {
            'population_demixing': ['baron', 'segerstolpe', 'ZilionisMouseLungCancer', 'MarquesMouseBrain', 'PBMCsmall',
                                    'DarmanisBrain', 'GuoHumanTestis', 'QuakeMouseHeart', 'ZeiselMouseBrain',
                                    'BaronMousePancreas', 'LaMannoHumanEmbryonicStem', 'LaMannoHumanEmbryonicMidbrain',
                                    'QuakeMouseSpleen', 'QuakeMouseTongue', 'Alles', 'Ariss', 'ToschesLizard'],  #

            'marker_discovery': ['baron', 'segerstolpe', 'PBMCsmall', 'simulatingPBMCsmall', 'ZeiselMouseBrain'],

            'intra-classification': ['baron', 'segerstolpe', 'ZilionisMouseLungCancer', 'MarquesMouseBrain',
                                     'PBMCsmall', 'DarmanisBrain', 'GuoHumanTestis', 'QuakeMouseHeart',
                                     'ZeiselMouseBrain', 'BaronMousePancreas', 'LaMannoHumanEmbryonicStem',
                                     'LaMannoHumanEmbryonicMidbrain', 'QuakeMouseSpleen', 'QuakeMouseTongue',
                                     'Alles', 'Ariss', 'ToschesLizard'],
            'inter-classification': ['PBMCbatchone+PBMCsmall', 'PBMCbatchtwo+PBMCsmall'],

            'clustering': ['PBMCbatchone', 'PBMCbatchtwo'],
            # seurat sc3s, bcubed_f1 ARI v-measure
            'batch_correction': ['AztekinTail', 'MouseAtlas', 'MouseHematopoieticStemProgenitor', 'MouseRetina',
                                 'PBMCbatches', 'baron+segerstolpe'],  #

            'computation_time': ['VentoHumanPlacenta500cells', 'VentoHumanPlacenta1000cells',
                                 'VentoHumanPlacenta2000cells', 'VentoHumanPlacenta5000cells',
                                 'VentoHumanPlacenta10000cells', 'VentoHumanPlacenta20000cells',
                                 'VentoHumanPlacenta50000cells',
                                 'GuoHumanTestis5000genes', 'GuoHumanTestis10000genes',
                                 'GuoHumanTestis15000genes', 'GuoHumanTestis20000genes',
                                 'GuoHumanTestis25000genes',
                                 'PBMC50000cells']  # large-scale datasets, only use 50000 samples
            # 'VentoHumanPlacenta5000cells'
        }

        self.random_seed = 2021
        # normalization
        self.scale_factor = 1e4
        # filtering
        self.n_filter_cell = 10
        self.n_filter_gene = 10
        self.n_genes = [500, 1000, 1500, 2000]  # 500, 1000, 1500, 2000
        self.max_timeout = 60 * 60 * 24  # 60 seconds * 60 minutes * 24 hours
        self.ensemble_mode = 'importance_sum'
        # batch-effect removal
        self.plot = True
        self.metric = True
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
            'scGeneFit': 'norm',
        })
        self.method_lan.update({
            'rf': 'python',
            'lgb': 'python',
            'xgb': 'python',
            'nsc': 'python',
            'fisher_score': 'python',
            'scGeneFit': 'python',
        })
        self.n_folds = 5  # 5 folds, intra-dataset
        self.evaluation_method = ['SingleR']


class ClusterConfig(ExperimentConfig):
    def __init__(self):
        super().__init__()
        self.n_loops = 20  # 1 loops
        self.rare_type_detection_metric = 'bcubed'  # f-measure
        self.evaluation_method = ['sc3s']


class DataConfig:
    def __init__(self):
        # data path
        self.data_path = "/volume2/bioinfo/scRNA/python_data/"

        # marker path
        self.marker_path = "/volume1/home/tdeng/SingleCell/Data/MarkerGene/"

        # remove_types in pancreas data
        self.remove_types = ['unclear', 'not.applicable', 'unclassified', 'co.expression', 'beta.contaminated',
                             'alpha.contaminated', 'delta.contaminated', 'gamma.contaminated', 'unassigned',
                             'MHC.class.II', 'unclassified.endocrine', np.NaN, 'nan']
        self.replace_types = {'stellate': ['activated_stellate', 'quiescent_stellate', 'PSC']}


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
