import numpy as np


class ExperimentConfig:
    def __init__(self):
        # measurements
        self.measurements = {
            'marker_discovery': ['PBMCbatchone', 'PBMCbatchtwo', 'PBMCsmall', 'SimPBMCsmall',
                                 'BaronHuman', 'Segerstolpe', 'Zeisel'],

            'intra-classification': ['PBMCbatchone', 'PBMCbatchtwo'],
            'inter-classification': ['PBMCbatches'],

            'clustering': ['baron', 'segerstolpe', 'ZilionisMouseLungCancer', 'MarquesMouseBrain',
                           'DarmanisBrain', 'Guo', 'QuakeMouseHeart',
                           'ZeiselMouseBrain', 'BaronMousePancreas', 'LaMannoStem',
                           'LaMannoMidbrain', 'QuakeMouseSpleen', 'QuakeMouseTongue',
                           'Alles', 'Ariss', 'ToschesLizard', 'PBMCbatchone', 'PBMCbatchtwo'],
            # seurat sc3s, bcubed_f1 ARI v-measure
            'batch_correction': ['AztekinTail', 'MouseAtlas', 'MouseHematopoieticStemProgenitor', 'MouseRetina',
                                 'PBMCbatches', 'baron+segerstolpe'],  #

            'computation_time': ['Vento500cells', 'Vento1000cells', 'Vento2000cells', 'Vento5000cells',
                                 'Vento10000cells', 'Vento20000cells', 'Vento50000cells',
                                 'Guo5000genes', 'Guo10000genes', 'Guo15000genes', 'Guo20000genes', 'Guo25000genes',
                                 'PBMC50000cells']  # large-scale datasets
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

        # used datasets
        self.full_dataset_names = {'Aztekin': 'AztekinTail',
                                   'BaronHuman': 'BaronHumanPancreas',
                                   'BaronMouse': 'BaronMousePancreas',
                                   'Segerstolpe': 'SegerstolpeHumanPancreas',
                                   'Vento': 'VentoHumanPlacenta',
                                   'Zilionis': 'ZilionisMouseLungCancer',
                                   'LaMannoMidbrain': 'LaMannoHumanEmbryonicMidbrain',
                                   'LaMannoStem': 'LaMannoHumanEmbryonicStem',
                                   'Zeisel': 'ZeiselMouseBrain',
                                   'Marques': 'MarquesMouseBrain',
                                   'Ariss': 'Ariss',
                                   'ToschesLizard': 'ToschesLizard',
                                   'PBMCsmall': 'PBMCsmall',
                                   'QuakeHeart': 'QuakeMouseHeart',
                                   'SimPBMCsmall': 'simulatingPBMCsmall',
                                   'QuakeTongue': 'QuakeMouseTongue',
                                   'Guo': 'GuoHumanTestis',
                                   'Alles': 'Alles',
                                   'QuakeSpleen': 'QuakeMouseSpleen',
                                   'Darmanis': 'DarmanisBrain',
                                   'MouseCellAtlas': 'MouseAtlas',
                                   'PBMCbatchone': 'PBMCbatch1',
                                   'PBMCbatchtwo': 'PBMCbatch2',
                                   'MouseRetina': 'MouseRetina',
                                   'MouseHSP': 'MouseHematopoieticStemProgenitor'}

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
