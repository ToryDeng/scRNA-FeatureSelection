import warnings
warnings.filterwarnings('ignore')

from utils.evaluation import evaluate_feature_selection_methods
# evaluate_feature_selection_methods(measurements=['population_demixing', 'marker_discovery',
#                                                  'classification', 'clustering',
#                                                  'batch_correction', 'computation_time'],
#                                    methods=['var', 'cv2', 'seurat_v3', 'deviance', 'm3drop', 'scmap',
#                                             'rf', 'lgb', 'xgb', 'nsc', 'fisher_score', 'scGeneFit',
#                                             'cellranger', 'feast', 'rf+fisher_score', 'feast+deviance'])

evaluate_feature_selection_methods(measurements=['clustering'],
                                   methods=['var', 'cv2', 'seurat_v3', 'deviance', 'm3drop', 'scmap',
                                            'rf', 'lgb', 'xgb', 'nsc', 'fisher_score', 'scGeneFit',
                                            'cellranger', 'feast', 'rf+fisher_score', 'feast+deviance'])
