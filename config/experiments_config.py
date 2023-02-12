class BasicExperimentConfig:
    def __init__(self):
        self.random_seed = 0
        self.scale_factor = 1e4  # normalization
        self.ensemble_mode = 'count_sum'  # ensemble selection
        self.n_genes = [2000, 3000, 4000, 5000]  # selected genes  , 1000, 1500, 2000
        self.sink_dir = 'records/'


class MarkerDiscoveryConfig(BasicExperimentConfig):
    def __init__(self):
        super(MarkerDiscoveryConfig, self).__init__()
        self.PBMC_markers = ['PBMCSLEA', 'PBMCSLEB', 'PBMCSLEC', 'PBMCSLEctrl', 'PBMCSLEstim']
        self.pancreas_markers = ['BaronHumanDonorTwo']
        self.mouse_brain_marekrs = []


class CellClassificationConfig(BasicExperimentConfig):
    def __init__(self):
        super(CellClassificationConfig, self).__init__()
        self.methods = ['scmap']  # 'SingleR', 'SVM', 'scmap'
        self.is_intra = False
        self.intra_datasets = ['PBMCSLEA', 'PBMCSLEB', 'PBMCSLEC', 'PBMCSLEctrl', 'PBMCSLEstim']
        self.n_folds = 5  # 5 folds, in intra-dataset
        self.inter_datasets = ['PBMCSLEctrl+PBMCSLEstim']
        self.metrics = ['f1', 'ck']


class CellClusteringConfig(BasicExperimentConfig):
    def __init__(self):
        super(CellClusteringConfig, self).__init__()
        self.datasets = [
            'QuakeTrachea']

        # 'PBMCSLEA', 'PBMCSLEB', 'PBMCSLEC', 'PBMCSLEctrl', 'PBMCSLEstim',
        # 'Adam', 'Chen', 'Guo', 'Plasschaert', 'QuakeTrachea', 'ToschesLizard',
        # 'PBMCeightkilo', 'PBMCsevenkilo', 'ZeiselBrain', 'ZilionisLung'

        # 'Seurat_v4': 1, 'KMeans': 100, 'TSCAN': 1, 'SC3s': 100s
        self.methods = {'Seurat_v4': 1}  # clustering_method: number of runs 'Seurat_v4', 'KMeans', 'TSCAN', 'SC3s'
        self.metrics = ['ARI', 'V']  #, 'SI', 'NMI',


class BatchCorrectionConfig(BasicExperimentConfig):
    def __init__(self):
        super(BatchCorrectionConfig, self).__init__()
        self.datasets = ['Aztekin', 'Campbell', 'Shekhar', 'PBMCbatchone+PBMCbatchtwo', 'BaronHuman+Segerstolpe']
        self.methods = ['Seurat_v4']
        self.orders = ['correction_first', 'selection_first']
        self.plot = True
        self.metrics = ['kBET', 'iLISI', 'cLISI', 'f1LISI']
        self.figure_dir = 'figures/'


class ComputationTimeConfig(BasicExperimentConfig):
    def __init__(self):
        super(ComputationTimeConfig, self).__init__()
        self.datasets = ['Zhao500cells', 'Zhao1000cells', 'Zhao5000cells', 'Zhao10000cells', 'Zhao50000cells',
                         'Guo5000genes', 'Guo10000genes', 'Guo15000genes', 'Guo20000genes', 'Guo25000genes']


base_cfg = BasicExperimentConfig()
marker_cfg = MarkerDiscoveryConfig()
assign_cfg = CellClassificationConfig()
cluster_cfg = CellClusteringConfig()
batch_cfg = BatchCorrectionConfig()
time_cfg = ComputationTimeConfig()
