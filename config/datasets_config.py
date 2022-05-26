import numpy as np


class DatasetConfig:
    def __init__(self):
        # data path
        self.data_path = "/volume2/bioinfo/scRNA/python_data/"
        self.cache_path = 'cache/'

        # marker path
        self.marker_path = "/volume1/home/tdeng/SingleCell/Data/MarkerGene/"

        # used datasets
        self.full_dataset_names = {
            'Ariss': 'Ariss',
            'Aztekin': 'AztekinTail',
            'Adam': 'Adam',
            'Bach': 'Bach',
            'BachDonorTwo': 'Bach_donor_2',
            'BachDonorThree': 'Bach_donor_3',
            'BachDonorFour': 'Bach_donor_4',
            'BaronHuman': 'BaronHuman',
            'BaronHumanDonorTwo': 'BaronHuman_GSM2230759',
            'Campbell': 'Campbell',
            'Chen': 'Chen',
            'ChenHungry': 'Chen_Hungry',
            'ChenHungryTwo': 'Chen_Hungry2',
            'ChenNormal': 'Chen_Normal',
            'Darmanis': 'DarmanisBrain',
            'Guo': 'Guo',
            'Hochane': 'Hochane',
            'HochaneDonorTwo': 'Hochane_donor2',
            'HochaneDonorFour': 'Hochane_donor4',
            'PBMCSLEA': 'PBMCSLEA',
            'PBMCSLEB': 'PBMCSLEB',
            'PBMCSLEC': 'PBMCSLEC',
            'PBMCSLEctrl': 'PBMCSLEctrl',
            'PBMCSLEstim': 'PBMCSLEstim',
            'Plasschaert': 'Plasschaert',
            'QuakeTrachea': 'QuakeTrachea',
            'QuakeTracheaDonorFS': 'QuakeTrachea_3-M-5-6',
            'QuakeTracheaDonorSE': 'QuakeTrachea_3-M-7-8',
            'QuakeSpleen': 'QuakeSpleen',
            'QuakeSpleenDonorFS': 'QuakeSpleen_3-F-56',
            'QuakeSpleenDonorE': 'QuakeSpleen_3-M-8',
            'Shekhar': 'Shekhar',
            'ToschesLizard': 'ToschesLizard',
            'Zhao': 'ZhaoImmuneLiver'
        }

        # remove_types in pancreas data
        self.remove_types = [
            'unclear', 'not.applicable', 'unclassified.cell', 'co.expression.cell', 'beta.contaminated',
            'alpha.contaminated', 'delta.contaminated', 'gamma.contaminated', 'unassigned',
            'MHC.class.II', 'unclassified.endocrine.cell', np.NaN, 'nan'
        ]
        self.replace_types = {
            'stellate': ['activated_stellate', 'quiescent_stellate', 'PSC', 'PSC cell'],
            'MHC class II cell': 'MHC.class.II',
            'acinar': 'acinar cell',
            'alpha': 'alpha cell',
            'beta': 'beta cell',
            'delta': 'delta cell',
            'ductal': 'ductal cell',
            'endothelial': 'endothelial cell',
            'epsilon': 'epsilon cell',
            'gamma': 'gamma cell',
            'mast': 'mast cell',
            'CD4.T.cell': ['Naive CD4 T', 'Memory CD4 T'],
            'Monocyte_CD14': 'CD14+ Mono',
            'Monocyte_FCGR3A': 'FCGR3A+ Mono',
            'NK.cell': 'NK',
            'B.cell': 'B',
            'CD8.T.cell': 'CD8 T'
        }
        self.force_min_cells = 3
        self.rare_rate = 0.05  # at most 5% of all cells
        self.rare_number = 30  # at least 30 cells


data_cfg = DatasetConfig()
