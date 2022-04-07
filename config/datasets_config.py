import numpy as np


class DatasetConfig:
    def __init__(self):
        # data path
        self.data_path = "/volume2/bioinfo/scRNA/python_data/"
        self.cache_path = 'cache/'

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
                                   'MouseAtlas': 'MouseCellAtlas',
                                   'PBMCbatchone': 'PBMCbatch1',
                                   'PBMCbatchtwo': 'PBMCbatch2',
                                   'MouseRetina': 'MouseRetina',
                                   'MouseHSP': 'MouseHematopoieticStemProgenitor'}

        # remove_types in pancreas data
        self.remove_types = ['unclear', 'not.applicable', 'unclassified', 'co.expression', 'beta.contaminated',
                             'alpha.contaminated', 'delta.contaminated', 'gamma.contaminated', 'unassigned',
                             'MHC.class.II', 'unclassified.endocrine', np.NaN, 'nan']
        self.replace_types = {'stellate': ['activated_stellate', 'quiescent_stellate', 'PSC'],
                              'CD4.T.cell': ['Naive CD4 T', 'Memory CD4 T'],
                              'Monocyte_CD14': 'CD14+ Mono',
                              'Monocyte_FCGR3A': 'FCGR3A+ Mono',
                              'NK.cell': 'NK',
                              'B.cell': 'B',
                              'CD8.T.cell': 'CD8 T'
                              }


data_cfg = DatasetConfig()