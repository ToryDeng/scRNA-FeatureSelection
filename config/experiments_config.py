class BasicExperimentConfig:
    def __init__(self):
        self.random_seed = 2022
        self.scale_factor = 1e4  # normalization
        self.ensemble_mode = 'count_sum'  # ensemble selection
        self.n_genes = [500, 1000, 1500, 2000]  # selected genes  
        self.sink_dir = 'records/'


class MarkerDiscoveryConfig(BasicExperimentConfig):
    def __init__(self):
        super(MarkerDiscoveryConfig, self).__init__()
        self.datasets = ['PBMCbatchone', 'PBMCbatchtwo', 'PBMCsmall', 'BaronHuman', 'Segerstolpe', 'Zeisel']


class CellClassificationConfig(BasicExperimentConfig):
    def __init__(self):
        super(CellClassificationConfig, self).__init__()
        self.methods = ['SingleR']
        self.is_intra = False
        self.intra_datasets = ['PBMCbatchone', 'PBMCbatchtwo']
        self.n_folds = 5  # 5 folds, in intra-dataset
        self.inter_datasets = ['PBMCbatchone+PBMCbatchtwo']
        self.metrics = ['f1', 'ck']


class CellClusteringConfig(BasicExperimentConfig):
    def __init__(self):
        super(CellClusteringConfig, self).__init__()
        self.datasets = ['PBMCbatchone', 'PBMCbatchtwo', 'Alles', 'Ariss', 'ToschesLizard',
                         'BaronHuman', 'Segerstolpe', 'Zilionis', 'Marques', 'Darmanis', 'Guo', 'QuakeHeart',
                         'Zeisel', 'BaronMouse', 'LaMannoStem', 'LaMannoMidbrain', 'QuakeSpleen', 'QuakeTongue'
                         ]

        self.methods = {'Seurat_v4': 1}  # clustering_method: number of runs  # , 'SC3s': 1
        self.metrics = ['ARI', 'V', 'bcubed']


class BatchCorrectionConfig(BasicExperimentConfig):
    def __init__(self):
        super(BatchCorrectionConfig, self).__init__()
        self.datasets = ['Aztekin', 'MouseAtlas', 'MouseHSP', 'MouseRetina',
                         'PBMCbatchone+PBMCbatchtwo', 'BaronHuman+Segerstolpe']
        self.methods = ['Seurat_v4']
        self.orders = ['correction_first', 'selection_first']
        self.plot = True
        self.metrics = ['kBET', 'iLISI', 'cLISI', 'f1LISI']
        self.figure_dir = 'figures/'


class ComputationTimeConfig(BasicExperimentConfig):
    def __init__(self):
        super(ComputationTimeConfig, self).__init__()
        self.datasets = ['Vento500cells', 'Vento1000cells', 'Vento2000cells', 'Vento5000cells',
                         'Vento10000cells', 'Vento20000cells', 'Vento50000cells',
                         'Guo5000genes', 'Guo10000genes', 'Guo15000genes', 'Guo20000genes', 'Guo25000genes']


base_cfg = BasicExperimentConfig()
marker_cfg = MarkerDiscoveryConfig()
assign_cfg = CellClassificationConfig()
cluster_cfg = CellClusteringConfig()
batch_cfg = BatchCorrectionConfig()
time_cfg = ComputationTimeConfig()
